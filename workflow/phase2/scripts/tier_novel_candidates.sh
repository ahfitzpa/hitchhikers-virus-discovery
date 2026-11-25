#!/usr/bin/env bash
# Tier the phase2 viral positives into novel-candidate confidence classes.
#
# Reads every data/05.final_selection/<family>/<sample>/taxonomy_merged.tsv
# under --root and writes one row per contig that passes the gate, with a tier
# and the explicit reasons behind it. Pure bash + awk: no conda env, no network,
# safe to run on a login node while the pipeline is still going.
#
# The tiers encode the signals the pipeline already computes:
#
#   NOVEL signal   kraken2 assigned no family while the tool consensus did ->
#                  viral by consensus, but matching no known species. This is
#                  the divergence signature; a contig kraken2 names is a known
#                  virus, not a candidate.
#   A  strong      NOVEL + tools agree + >=1 geNomad hallmark gene
#                  + CheckV High-quality/Complete + plausible genome length
#   B  review      NOVEL + any two of {agreement, hallmark, CheckV good}
#   C  weak        NOVEL but fragmentary, single-tool, or hallmark-free
#   -  known       no NOVEL signal (kraken2 placed it) -- reported only with
#                  --include-known
#
# Flags are additive and independent of tier: CHIMERA_RISK, LOW_COMPLETENESS,
# NO_HALLMARK, NO_CHECKV, SINGLE_TOOL, CONFLICT.
#
# NOT a dereplication step. Identical viruses recovered from several libraries
# appear as several rows; on the Hepeviridae set twelve ERR libraries carried
# RdRp blocks that were identical at the amino acid level. Cluster the
# representatives (mmseqs2/cd-hit) before counting anything, or every total is
# inflated. --exact-dup-check catches only byte-identical contigs, which is a
# strict subset of that problem.
#
# Usage:
#   tier_novel_candidates.sh --root $PWD/workflow/phase2 \
#       [--family Hepeviridae] [--out candidates.tsv] [--min-len 1000]
#       [--genome-min 3000] [--genome-max 20000] [--include-known]
#       [--exact-dup-check] [--summary-only]
set -euo pipefail

ROOT="."
FAMILY=""
OUT="-"
MIN_LEN=1000
GENOME_MIN=3000
GENOME_MAX=20000
INCLUDE_KNOWN=0
DUP_CHECK=0
SUMMARY_ONLY=0

die() { printf 'error: %s\n' "$*" >&2; exit 2; }

usage() { sed -n '2,40p' "$0" | sed 's/^# \{0,1\}//'; exit 0; }

while [ $# -gt 0 ]; do
    case "$1" in
        --root)          ROOT="${2:?--root needs a path}"; shift 2 ;;
        --family)        FAMILY="${2:?--family needs a name}"; shift 2 ;;
        --out)           OUT="${2:?--out needs a path}"; shift 2 ;;
        --min-len)       MIN_LEN="${2:?}"; shift 2 ;;
        --genome-min)    GENOME_MIN="${2:?}"; shift 2 ;;
        --genome-max)    GENOME_MAX="${2:?}"; shift 2 ;;
        --include-known) INCLUDE_KNOWN=1; shift ;;
        --exact-dup-check) DUP_CHECK=1; shift ;;
        --summary-only)  SUMMARY_ONLY=1; shift ;;
        -h|--help)       usage ;;
        *)               die "unknown argument: $1 (try --help)" ;;
    esac
done

for v in MIN_LEN GENOME_MIN GENOME_MAX; do
    case "${!v}" in (*[!0-9]*|"") die "--${v,,} must be a non-negative integer, got '${!v}'";; esac
done
[ -d "$ROOT" ] || die "--root is not a directory: $ROOT"
SEL="$ROOT/data/05.final_selection"
[ -d "$SEL" ] || die "no data/05.final_selection under $ROOT"

# Collect the per-sample merges. A family filter restricts the directory level,
# so an unmatched name yields a clear "no input" rather than an empty table.
mapfile -d '' -t TABLES < <(
    if [ -n "$FAMILY" ]; then
        find "$SEL/$FAMILY" -name taxonomy_merged.tsv -print0 2>/dev/null
    else
        find "$SEL" -name taxonomy_merged.tsv -print0 2>/dev/null
    fi
)
[ "${#TABLES[@]}" -gt 0 ] || die "no taxonomy_merged.tsv found under $SEL${FAMILY:+/$FAMILY}"
printf 'scanning %d sample table(s)\n' "${#TABLES[@]}" >&2

TMP="$(mktemp -d "${TMPDIR:-/tmp}/tier_novel.XXXXXX")"
trap 'rm -rf "$TMP"' EXIT

# One awk pass over every table. Columns are resolved BY NAME from each file's
# own header, never by fixed index: taxonomy_merged.tsv has gained columns
# before (CheckV, then VirSorter2) and positional parsing would silently
# misread older or newer runs.
awk -v FS='\t' -v OFS='\t' \
    -v min_len="$MIN_LEN" -v gmin="$GENOME_MIN" -v gmax="$GENOME_MAX" \
    -v want_family="$FAMILY" -v include_known="$INCLUDE_KNOWN" '
function col(name,   i) { return (name in H) ? $H[name] : "" }
function num(v) { return (v == "" || v ~ /^[[:space:]]*$/) ? -1 : v + 0 }

FNR == 1 {
    delete H
    for (i = 1; i <= NF; i++) { h = $i; gsub(/^[ \t]+|[ \t]+$/, "", h); H[h] = i }
    # .../<family>/<sample>/taxonomy_merged.tsv
    n = split(FILENAME, P, "/")
    sample = P[n-1]; family = P[n-2]
    next
}

{
    contig = col("contig_id")
    if (contig == "") next

    cons = col("consensus_family")
    if (cons == "") next
    if (want_family != "" && cons != want_family) next

    # SPAdes encodes length in the name: NODE_5_length_4981_cov_59.4
    len = 0
    if (match(contig, /length_[0-9]+/)) len = substr(contig, RSTART+7, RLENGTH-7) + 0
    if (len < min_len) next

    kraken_fam = col("kraken2_family")
    agreement  = col("family_agreement")
    hallmarks  = num(col("genomad_n_hallmarks"))
    quality    = col("checkv_quality")
    cv_compl      = num(col("checkv_completeness"))
    ntools     = num(col("n_tools"))

    novel     = (kraken_fam == "")
    agreed    = (agreement == "agree")
    hallmark  = (hallmarks >= 1)
    checkv_ok = (quality == "High-quality" || quality == "Complete")
    len_ok    = (len >= gmin && len <= gmax)

    if (!novel && !include_known) next

    reasons = ""; flags = ""
    if (!hallmark)                       flags = flags "NO_HALLMARK;"
    if (quality == "" || quality == "Not-determined") flags = flags "NO_CHECKV;"
    if (cv_compl >= 0 && cv_compl < 50)        flags = flags "LOW_COMPLETENESS;"
    if (ntools == 1 || agreement == "single_tool") flags = flags "SINGLE_TOOL;"
    if (agreement == "conflict")         flags = flags "CONFLICT;"
    # Well over a plausible genome for the family: usually a misassembly or a
    # host-fused contig rather than a large novel genome.
    if (gmax > 0 && len > gmax)          flags = flags "CHIMERA_RISK;"

    if (!novel) {
        tier = "-"; reasons = "kraken2 assigned " kraken_fam
    } else {
        support = agreed + hallmark + checkv_ok
        if (agreed && hallmark && checkv_ok && len_ok) {
            tier = "A"
            reasons = "no kraken2 family;tools agree;hallmark;CheckV " quality ";length in range"
        } else if (support >= 2) {
            tier = "B"
            reasons = "no kraken2 family"
            if (agreed)    reasons = reasons ";tools agree"
            if (hallmark)  reasons = reasons ";hallmark"
            if (checkv_ok) reasons = reasons ";CheckV " quality
            if (!len_ok)   reasons = reasons ";length outside " gmin "-" gmax
        } else {
            tier = "C"
            reasons = "no kraken2 family;weak support"
            if (!agreed)    reasons = reasons ";agreement=" (agreement == "" ? "NA" : agreement)
            if (!hallmark)  reasons = reasons ";no hallmark"
            if (!checkv_ok) reasons = reasons ";CheckV " (quality == "" ? "NA" : quality)
        }
    }
    sub(/;$/, "", flags)

    print tier, family, sample, contig, len, cons, agreement, \
          (hallmarks < 0 ? "" : hallmarks), quality, \
          (cv_compl < 0 ? "" : cv_compl), kraken_fam, col("genomad_family"), \
          col("rdrpcatch_family"), col("detected_by"), flags, reasons
}
' "${TABLES[@]}" | LC_ALL=C sort -t"$(printf '\t')" -k1,1 -k5,5nr > "$TMP/rows.tsv"

# Byte-identical contigs across libraries: the cheap, strict subset of
# dereplication. Equal-length rows are hashed from each sample's
# viral_positives.fasta; anything subtler needs real clustering.
DUPCOL=""
if [ "$DUP_CHECK" -eq 1 ]; then
    : > "$TMP/hashes.tsv"
    while IFS=$'\t' read -r _ family sample contig _; do
        fa="$SEL/$family/$sample/viral_positives.fasta"
        [ -s "$fa" ] || continue
        printf '%s\t%s\t%s\n' "$family" "$sample" "$contig" >> "$TMP/want.tsv"
    done < "$TMP/rows.tsv"
    if [ -s "${TMP}/want.tsv" ]; then
        awk -v FS='\t' -v sel="$SEL" '
        { want[$1 FS $2 FS $3] = 1; fams[$1 FS $2] = 1 }
        END {
            for (fs in fams) {
                split(fs, a, FS)
                fa = sel "/" a[1] "/" a[2] "/viral_positives.fasta"
                id = ""; seq = ""
                while ((getline line < fa) > 0) {
                    if (substr(line, 1, 1) == ">") {
                        if (id != "" && ((a[1] FS a[2] FS id) in want)) print a[1], a[2], id, seq
                        id = substr(line, 2); sub(/[ \t].*$/, "", id); seq = ""
                    } else seq = seq line
                }
                if (id != "" && ((a[1] FS a[2] FS id) in want)) print a[1], a[2], id, seq
                close(fa)
            }
        }' OFS='\t' "$TMP/want.tsv" \
        | while IFS=$'\t' read -r family sample contig seq; do
              printf '%s\t%s\t%s\t%s\n' "$family" "$sample" "$contig" \
                     "$(printf '%s' "$seq" | md5sum | cut -d' ' -f1)"
          done > "$TMP/hashes.tsv"
    fi
    DUPCOL="exact_dup_group"
fi

{
    printf 'tier\tfamily\tsample\tcontig_id\tlength\tconsensus_family\tfamily_agreement'
    printf '\tgenomad_hallmarks\tcheckv_quality\tcheckv_completeness\tkraken2_family'
    printf '\tgenomad_family\trdrpcatch_family\tdetected_by\tflags\treasons'
    [ -n "$DUPCOL" ] && printf '\t%s\t%s' "$DUPCOL" "exact_dup_n"
    printf '\n'
    if [ -n "$DUPCOL" ] && [ -s "$TMP/hashes.tsv" ]; then
        awk -v FS='\t' -v OFS='\t' '
        NR == FNR { h[$1 FS $2 FS $3] = $4; n[$4]++; next }
        { k = $2 FS $3 FS $4
          if (k in h) { g = substr(h[k], 1, 8); print $0, g, n[h[k]] }
          else        { print $0, "", "" } }
        ' "$TMP/hashes.tsv" "$TMP/rows.tsv"
    else
        cat "$TMP/rows.tsv"
    fi
} > "$TMP/out.tsv"

{
    printf '\n== novel-candidate tiers ==\n'
    awk -v FS='\t' 'NR>1 { c[$1]++ } END {
        order = "A B C -"; n = split(order, o, " ")
        for (i = 1; i <= n; i++) if (o[i] in c) printf "  tier %-2s %6d\n", o[i], c[o[i]]
    }' "$TMP/out.tsv"
    printf '  %-9s %6d\n' "total" "$(( $(wc -l < "$TMP/out.tsv") - 1 ))"
    printf '\n== flags (independent of tier) ==\n'
    awk -v FS='\t' 'NR>1 && $15 != "" { n = split($15, f, ";"); for (i=1;i<=n;i++) if (f[i] != "") c[f[i]]++ }
         END { for (k in c) printf "  %-18s %6d\n", k, c[k] }' "$TMP/out.tsv" | LC_ALL=C sort -k2,2nr
    if [ -n "$DUPCOL" ]; then
        printf '\n== byte-identical contig groups ==\n'
        awk -v FS='\t' 'NR>1 && $18 != "" && $18+0 > 1 { g[$17] = $18 }
             END { n=0; for (k in g) n++; printf "  %d group(s) with >1 member\n", n }' "$TMP/out.tsv"
        printf '  NOTE: exact duplicates only; near-identical viruses need mmseqs2/cd-hit.\n'
    fi
    printf '\n'
} >&2

if [ "$SUMMARY_ONLY" -eq 1 ]; then
    exit 0
elif [ "$OUT" = "-" ]; then
    cat "$TMP/out.tsv"
else
    mkdir -p "$(dirname "$OUT")"
    cp "$TMP/out.tsv" "$OUT"
    printf 'wrote %s\n' "$OUT" >&2
fi
