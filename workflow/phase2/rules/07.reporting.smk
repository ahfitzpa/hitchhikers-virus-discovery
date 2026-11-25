# Machine-readable capture for the conifer_pest_virome run (spec items 1-6 and
# the candidate contig table). Nothing here changes what the pipeline calls
# viral -- these rules only read what the taxonomy path already produced.
#
# Items 1, 2, 3 and the candidate table aggregate across every library, so they
# take get_all_final_files() (defined in rules/05.merge_tax.smk) as input and
# therefore run once, after the last sample's merge.

REPORT_DIR = "results/reporting"


# -------------------------------------------------------------------------
# Item 6a: RdRp domains as protein FASTA, per library.
#
# Reads RdRpCATCH's six-frame translation from its scan working directory if
# it is still there. In Serratus it never is: run_rdrpcatch's EXIT trap
# removes scan/ on success (rules/04.taxonomy.smk), unlike conifer, where
# taxonomy_cleanup does the purge and this rule is ordered before it.
#
# So when no translation is found, it is rebuilt: the hit contigs (Contig_name
# column) are pulled from filtered_contigs.fasta -- RdRpCATCH's own input,
# minus contigs over 299,999 bp, which never carry a hit -- and translated
# with RdRpCATCH's exact call (rdrpcatch/run_seqkit.py, defaults):
#   seqkit translate --clean --append-frame -f 6 -T 1
# Each contig is translated independently, so the frames, names
# (..._frame=N) and amino-acid coordinates match the hit table exactly.
#
# A header-only hit table (RdRpCATCH found nothing) is passed as /dev/null,
# which the script treats as "no hits"; its column check would otherwise
# reject the bare "Contig_name" header.
# -------------------------------------------------------------------------
rule rdrp_domains:
    input:
        rdrp = "data/04.taxonomy/{family}/{sample}/rdrpcatch/rdrp_annotated.tsv",
        contigs = "data/03.rnaviralSPAdes/{family}/{sample}/filtered_contigs.fasta"
    output:
        faa = "data/06.rdrp/{family}/{sample}/{sample}_rdrp_domains.faa",
        tsv = "data/06.rdrp/{family}/{sample}/{sample}_rdrp_domains.tsv"
    benchmark:
        "benchmarks/rdrp_domains/{family}_{sample}.tsv"
    threads: 1
    resources: mem_gb=8, runtime=60
    conda: get_env("TAX_CORE_ENV")
    shell:
        """
        SCAN_DIR="data/04.taxonomy/{wildcards.family}/{wildcards.sample}/rdrpcatch/scan"
        OUT_DIR=$(dirname {output.faa})
        mkdir -p "$OUT_DIR"
        IDS="$OUT_DIR/.hit_contigs.txt"
        REBUILT="$OUT_DIR/.hit_contigs_translated.fasta"
        trap 'rm -f "$IDS" "$REBUILT"' EXIT

        RDRP_TSV={input.rdrp}
        if [ "$(wc -l < {input.rdrp})" -le 1 ]; then
            echo "[rdrp_domains] {wildcards.sample}: no RdRpCATCH hits" >&2
            RDRP_TSV=/dev/null
        fi

        # run_rdrpcatch deletes scan/ on success, so for almost every sample
        # this find exits 1. Under `set -euo pipefail` that aborted the whole
        # rule at this line, before any work: 197 of 199 samples failed in
        # ~5 s. Guard the directory and swallow the status.
        TRANS=""
        if [ -d "$SCAN_DIR" ]; then
            TRANS=$(find "$SCAN_DIR" -path "*seqkit_translate_output*" -name "*.fasta" 2>/dev/null | head -1 || true)
        fi
        if [ -z "$TRANS" ] && [ "$RDRP_TSV" != /dev/null ]; then
            awk -F'\\t' 'NR==1{{for(i=1;i<=NF;i++) if($i=="Contig_name") c=i; next}} c&&$c!=""{{print $c}}' \
                {input.rdrp} | sort -u > "$IDS"
            echo "[rdrp_domains] {wildcards.sample}: scan translation gone; re-translating $(wc -l < "$IDS") hit contig(s)" >&2
            seqkit grep -n -f "$IDS" {input.contigs} 2>/dev/null \
                | seqkit translate --clean --append-frame -f 6 -T 1 -o "$REBUILT"
            if [ ! -s "$REBUILT" ]; then
                # -n matches whole headers; SPAdes headers have no description,
                # but fall back to ID matching rather than lose the sample.
                seqkit grep -f "$IDS" {input.contigs} \
                    | seqkit translate --clean --append-frame -f 6 -T 1 -o "$REBUILT"
            fi
            TRANS="$REBUILT"
        fi
        if [ -z "$TRANS" ]; then
            TRANS=/dev/null
        fi
        "${{CONDA_PREFIX:+$CONDA_PREFIX/bin/}}python" scripts/extract_rdrp_domains.py \
            --rdrp-tsv "$RDRP_TSV" \
            --translated "$TRANS" \
            --family {wildcards.family} \
            --sample {wildcards.sample} \
            --out-faa {output.faa} \
            --out-tsv {output.tsv}
        """


def get_all_rdrp_domains(wildcards):
    return [
        f"data/06.rdrp/{row.family}/{row.sample_id}/{row.sample_id}_rdrp_domains.faa"
        for _, row in MANIFEST.iterrows()
    ]


def get_all_coverage_stats(wildcards):
    return [
        f"data/03.rnaviralSPAdes/{row.family}/{row.sample_id}/coverage_stats.tsv"
        for _, row in MANIFEST.iterrows()
    ]


# -------------------------------------------------------------------------
# Item 6b: MetaCerberus ORF calls, lifted out of its output directory.
#
# MetaCerberus is run with --fraggenescan; its internal layout is not a stable
# contract, so the ORF FASTA is located by extension rather than a fixed path.
# -------------------------------------------------------------------------
rule export_metacerberus_orfs:
    input:
        mc = "results/annotation/metacerberus"
    output:
        faa = REPORT_DIR + "/metacerberus_orfs.faa",
        index = REPORT_DIR + "/metacerberus_orfs_sources.tsv"
    benchmark:
        "benchmarks/export_metacerberus_orfs/export.tsv"
    threads: 1
    resources: mem_gb=4, runtime=60
    shell:
        """
        mkdir -p $(dirname {output.faa})
        : > {output.faa}
        printf 'source_file\\tn_sequences\\n' > {output.index}
        FOUND=0
        while IFS= read -r f; do
            [ -s "$f" ] || continue
            N=$(grep -c '^>' "$f" || true)
            [ "$N" -gt 0 ] || continue
            cat "$f" >> {output.faa}
            printf '%s\\t%s\\n' "$f" "$N" >> {output.index}
            FOUND=$((FOUND+1))
        done < <(find {input.mc} -type f \\( -name '*.faa' -o -name '*.ffn' -o -name '*protein*.fa' \\) 2>/dev/null | sort)
        if [ "$FOUND" -eq 0 ]; then
            echo "WARNING: no ORF FASTA found under {input.mc}; {output.faa} is empty" >&2
        fi
        """


# -------------------------------------------------------------------------
# Item 1: per-library, per-stage counts (tidy long + wide pivot).
# -------------------------------------------------------------------------
rule stage_counts:
    input:
        finals = get_all_final_files,
        coverage = get_all_coverage_stats
    output:
        long = REPORT_DIR + "/stage_counts_long.tsv",
        wide = REPORT_DIR + "/stage_counts_wide.tsv"
    params:
        manifest = config.get("manifest", "data/intake/manifest.csv"),
        min_len = config["annotation"]["min_contig_len"]
    benchmark:
        "benchmarks/stage_counts/stage_counts.tsv"
    threads: 1
    resources: mem_gb=4, runtime=60
    conda: get_env("TAX_CORE_ENV")
    shell:
        """
        "${{CONDA_PREFIX:+$CONDA_PREFIX/bin/}}python" scripts/collect_stage_counts.py \
            --manifest {params.manifest} \
            --root . \
            --min-len {params.min_len} \
            --out-long {output.long} \
            --out-wide {output.wide}
        """


# -------------------------------------------------------------------------
# Item 3: per-library escalation fractions.
# -------------------------------------------------------------------------
rule escalation_fractions:
    input:
        finals = get_all_final_files,
        coverage = get_all_coverage_stats
    output:
        tsv = REPORT_DIR + "/escalation_fractions.tsv"
    params:
        manifest = config.get("manifest", "data/intake/manifest.csv")
    benchmark:
        "benchmarks/escalation_fractions/escalation.tsv"
    threads: 1
    resources: mem_gb=4, runtime=60
    conda: get_env("TAX_CORE_ENV")
    shell:
        """
        "${{CONDA_PREFIX:+$CONDA_PREFIX/bin/}}python" scripts/escalation_fractions.py \
            --manifest {params.manifest} \
            --root . \
            --out {output.tsv}
        """


# -------------------------------------------------------------------------
# OUTPUT: candidate contig table + the metric dictionary that names and
# defines every reported quantity (spec: the previous run's distance was
# unnamed and uninterpretable).
# -------------------------------------------------------------------------
rule candidate_contig_table:
    input:
        finals = get_all_final_files,
        coverage = get_all_coverage_stats
    output:
        tsv = REPORT_DIR + "/candidate_contigs.tsv",
        metrics = REPORT_DIR + "/metrics_definitions.tsv"
    params:
        manifest = config.get("manifest", "data/intake/manifest.csv")
    benchmark:
        "benchmarks/candidate_contig_table/candidates.tsv"
    threads: 1
    resources: mem_gb=8, runtime=120
    conda: get_env("TAX_CORE_ENV")
    shell:
        """
        "${{CONDA_PREFIX:+$CONDA_PREFIX/bin/}}python" scripts/build_candidate_table.py \
            --manifest {params.manifest} \
            --root . \
            --out {output.tsv} \
            --out-metrics {output.metrics}
        """


# -------------------------------------------------------------------------
# Item 2: all Snakemake benchmarks combined.
#
# Last rule in the run: it reads the benchmark tree, so everything whose timing
# should appear must already have finished. Its own benchmark, and those of the
# reporting rules it waits on, are necessarily absent from its output.
# -------------------------------------------------------------------------
rule combined_benchmarks:
    input:
        finals = get_all_final_files,
        coverage = get_all_coverage_stats,
        rdrp = get_all_rdrp_domains,
        counts = REPORT_DIR + "/stage_counts_long.tsv",
        escalation = REPORT_DIR + "/escalation_fractions.tsv",
        candidates = REPORT_DIR + "/candidate_contigs.tsv"
    output:
        tsv = REPORT_DIR + "/benchmarks_combined.tsv"
    params:
        manifest = config.get("manifest", "data/intake/manifest.csv")
    threads: 1
    resources: mem_gb=4, runtime=60
    conda: get_env("TAX_CORE_ENV")
    shell:
        """
        "${{CONDA_PREFIX:+$CONDA_PREFIX/bin/}}python" scripts/collect_benchmarks.py \
            --manifest {params.manifest} \
            --benchmark-dir benchmarks \
            --out {output.tsv}
        """


# Contig confidence: the last step of phase2, so that nothing downstream reads a
# contig that the annotation says is host genome.
#
# Taxonomy alone retains eukaryotic material -- a 14.7 kb contig here encodes
# titin, others collagen, dynein and laminin. MetaCerberus settles it, but not by
# presence: it reports a best hit from every database for every ORF, so every
# annotated contig has PHROG and VOG hits. The script compares the best cellular
# (KOFam) score against the best viral (PHROG/VOG) score per ORF, and an
# extracted RdRp makes a contig high confidence outright.
rule contig_confidence:
    input:
        annotation = "results/annotation/metacerberus",
        taxonomy = "results/taxonomy/master_metadata.csv",
        rdrp = get_all_rdrp_domains,
    output:
        table = "results/reporting/contig_confidence.tsv",
        keep = "results/reporting/contigs_keep.txt",
        summary = "results/reporting/contig_confidence_summary.tsv",
        families = "results/reporting/analysis_families.yaml",
    log:
        "logs/contig_confidence.log",
    params:
        # Serratus keeps the shared scripts in workflow/common; conifer keeps
        # them beside its own. Resolve rather than duplicate the file.
        script = next((c for c in ("scripts/filter_contigs_by_annotation.py",
                                   "../common/scripts/filter_contigs_by_annotation.py")
                       if os.path.exists(c)),
                      "scripts/filter_contigs_by_annotation.py"),
        keep_levels = config.get("contig_confidence", {}).get("keep", "high,medium"),
        score_ratio = config.get("contig_confidence", {}).get("score_ratio", 1.5),
        min_viral_frac = config.get("contig_confidence", {}).get("min_viral_frac", 0.5),
        min_tools = config.get("contig_confidence", {}).get("min_tools", 2),
    conda: get_env("TAX_CORE_ENV")
    threads: 1
    resources: mem_gb=16, runtime=120
    shell:
        """
        ANNOT=$(ls -d {input.annotation}/final/annotations-* 2>/dev/null | head -1)
        if [ -z "$ANNOT" ]; then
            echo "no MetaCerberus annotation directory under {input.annotation}" >&2
            exit 1
        fi
        "${{CONDA_PREFIX:+$CONDA_PREFIX/bin/}}python" -I {params.script} \
            --annotation-dir "$ANNOT" \
            --taxonomy {input.taxonomy} \
            --rdrp "data/06.rdrp/*/*/*_rdrp_domains.tsv" \
            --out-tsv {output.table} \
            --out-keep {output.keep} \
            --out-summary {output.summary} \
            --out-families-yaml {output.families} \
            --score-ratio {params.score_ratio} \
            --min-viral-frac {params.min_viral_frac} \
            --min-tools {params.min_tools} \
            --keep "{params.keep_levels}" > {log} 2>&1
        """
