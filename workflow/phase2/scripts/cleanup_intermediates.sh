#!/usr/bin/env bash
# One-off purge of intermediates left by runs that predate the per-rule cleanup
# added on 2026-11-16 (rules/03.assembly.smk, rules/04.taxonomy.smk). Those
# rules now clean up after themselves; this catches what is already on disk.
#
# Only removes files no rule declares as output and nothing downstream reads,
# and only for a sample whose stage result already exists. Safe to re-run.
#
#   bash scripts/cleanup_intermediates.sh            # dry run: list + sizes
#   bash scripts/cleanup_intermediates.sh --delete   # actually delete
#
# Run from the phase2 workflow dir, with no phase2 snakemake running.
set -euo pipefail
cd "$(dirname "$0")/.."

DELETE=false
[ "${1:-}" = "--delete" ] && DELETE=true

if squeue -u "$USER" -h -o "%j" 2>/dev/null | grep -qE '^(snakemake_p2|p2_stage3|serratus_[a-z0-9_]+)$|^[0-9a-f-]{36}$'; then
    echo "phase2 jobs are queued/running (squeue) -- refusing to run." >&2
    exit 1
fi

LIST=$(mktemp)
trap 'rm -f "$LIST"' EXIT
add() { for p in "$@"; do [ -e "$p" ] && printf '%s\0' "$p" >> "$LIST"; done; return 0; }

S=data/03.rnaviralSPAdes
T=data/04.taxonomy

for d in "$S"/*/*/; do
    d=${d%/}
    # bwa-mem3 index: only ever needed to build sorted.bam
    add "$d"/contigs.{0123,amb,ann,bwt.2bit.64,pac}
    # samtools sort spill files from killed sorts
    add "$d"/sorted.bam.tmp.*.bam
    if [ -s "$d/filtered_contigs.fasta" ]; then
        # temp() outputs + the undeclared .bai, left by pre-temp() runs.
        # The BAM only once coverage has been taken from it: this script
        # deleting sorted.bam before contig_coverage_stats existed is why 96
        # runs had to be re-downloaded (2026-11).
        add "$d"/contigs.fasta
        if [ -s "$d/coverage_stats.tsv" ]; then
            add "$d"/sorted.bam "$d"/sorted.bam.bai
        fi
        # SPAdes working files from before -o moved to node-local /tmp
        add "$d"/assembly_graph*.gfa "$d"/scaffolds.fasta "$d"/scaffolds.paths \
            "$d"/contigs.paths "$d"/before_rr.fasta "$d"/first_pe_contigs.fasta \
            "$d"/misc "$d"/.bin_reads "$d"/tmp "$d"/corrected "$d"/split_input \
            "$d"/pipeline_state "$d"/K[0-9]*
    fi
    if [ -f "$d/virsorter2/final-viral-score.tsv" ]; then
        while IFS= read -r -d '' p; do add "$p"; done < <(find "$d/virsorter2" -mindepth 1 -maxdepth 1 \
            ! -name 'final-viral-*' ! -name 'config.yaml' ! -name 'log' ! -name '.snakemake_timestamp' -print0)
    fi
done

for d in "$T"/*/*/; do
    d=${d%/}
    if [ -f "$d/genomad/filtered_contigs_summary/filtered_contigs_virus_summary.tsv" ]; then
        while IFS= read -r -d '' p; do add "$p"; done < <(find "$d/genomad" -mindepth 1 -maxdepth 1 -type d \
            ! -name 'filtered_contigs_summary' -print0)
    fi
    # rdrpcatch scratch: run_rdrpcatch recreates both on its next run
    add "$d"/rdrpcatch/scan "$d"/rdrpcatch/contigs_le_*.fasta
    add "$d"/checkv/tmp
done

# Cleaned reads of samples dropped from the manifest (configs/removed_samples.tsv)
if [ -f configs/removed_samples.tsv ]; then
    while IFS=$'\t' read -r sid _; do
        case "$sid" in ''|'#'*|sample_id) continue ;; esac
        add data/02.fastp/*/"$sid"/"$sid"_cleaned_{1,2}.fastq.gz
    done < configs/removed_samples.tsv
fi

# Never delete what later analysis reads (configs/retained_outputs.tsv): drop
# any protected path, or directory containing one, from the list.
SAFE=$(mktemp)
trap 'rm -f "$LIST" "$SAFE"' EXIT
python3 -I scripts/retention_guard.py filter-deletions < "$LIST" > "$SAFE" || {
    echo "retention guard failed -- refusing to delete anything." >&2; exit 1; }
mv "$SAFE" "$LIST"

N=$(tr -cd '\0' < "$LIST" | wc -c)
if [ "$N" -eq 0 ]; then echo "Nothing to clean."; exit 0; fi

echo "== by category (GB, count)"
xargs -0 du -sb < "$LIST" | awk -F'\t' '{
    n=$2; sub(/.*\//,"",n)
    if (n ~ /^contigs\.(0123|amb|ann|bwt\.2bit\.64|pac)$/) k="bwa-mem3 index"
    else if (n ~ /^sorted\.bam/ || n=="contigs.fasta") k="temp leftovers (bam/bai/contigs.fasta)"
    else if ($2 ~ /virsorter2/) k="virsorter2 working files"
    else if ($2 ~ /genomad/) k="genomad module dirs"
    else if ($2 ~ /rdrpcatch/) k="rdrpcatch scratch"
    else if ($2 ~ /checkv/) k="checkv tmp"
    else if ($2 ~ /02\.fastp/) k="fastp reads of removed samples"
    else k="old SPAdes working files"
    a[k]+=$1; c[k]++; t+=$1
} END { for (k in a) printf "  %8.1f  %6d  %s\n", a[k]/1e9, c[k], k; printf "  %8.1f  total\n", t/1e9 }' | sort -rn

if $DELETE; then
    xargs -0 rm -rf < "$LIST"
    echo "Deleted $N paths."
    python3 -I scripts/retention_guard.py audit
else
    echo "Dry run: $N paths. Re-run with --delete to remove them."
fi
