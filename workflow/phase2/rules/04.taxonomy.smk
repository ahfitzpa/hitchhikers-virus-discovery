# mmseqs (nt and GPU/nr), ICCTax and the mmseqs viral/orphan partition were
# replaced by the three-tool taxonomy. Kraken2, geNomad and RdRpCATCH all run
# independently, in parallel, on the same coverage-filtered assembly -- none of
# them gates another, and VirSorter2 no longer pre-filters their input (it
# dropped RNA viruses the other tools catch). Every contig ANY of the three
# calls viral is kept; rules/05.merge_tax.smk unions them and CheckV/VirSorter2
# annotate. Output-tracking note: rules here declare only the result files,
# never directory outputs, so taxonomy_cleanup can purge tool scratch inside
# these dirs without bumping a tracked output's mtime (which, under --rerun-
# triggers mtime, would re-run everything downstream).

def get_taxonomy_contigs(wildcards):
    return f"data/03.rnaviralSPAdes/{wildcards.family}/{wildcards.sample}/filtered_contigs.fasta"

rule run_kraken2:
    conda: get_env("TAX_CORE_ENV")
    threads: 8
    resources: mem_gb=16, runtime=60
    input:
        contigs=get_taxonomy_contigs,
        db=os.path.join(config["databases"]["kraken2_db"], "hash.k2d")
    output:
        tsv="data/04.taxonomy/{family}/{sample}/kraken2/kraken2_output.tsv",
        report="data/04.taxonomy/{family}/{sample}/kraken2/kraken2_report.txt"
    benchmark:
        "benchmarks/run_kraken2/{family}_{sample}.tsv"
    shell:
        """
        kraken2 --db {config[databases][kraken2_db]} --threads {threads} --use-names \
            --report {output.report} --output {output.tsv} {input.contigs}
        """

rule run_genomad:
    # Pinned to the geNomad version the method comparison and CAMI2
    # benchmark ran (envs/genomad.yaml).
    conda: get_env("GENOMAD_ENV")
    threads: 16
    resources: mem_gb=64, runtime=1440
    input:
        contigs=get_taxonomy_contigs,
        db=os.path.join(config["databases"]["genomad_db"], "version.txt")
    output:
        # geNomad names everything after the input fasta's basename.
        summary="data/04.taxonomy/{family}/{sample}/genomad/filtered_contigs_summary/filtered_contigs_virus_summary.tsv"
    benchmark:
        "benchmarks/run_genomad/{family}_{sample}.tsv"
    shell:
        """
        OUT_DIR=data/04.taxonomy/{wildcards.family}/{wildcards.sample}/genomad
        # geNomad resumes from whatever it finds in OUT_DIR; start clean
        # so a retry can't silently mix in a half-finished earlier attempt.
        rm -rf "$OUT_DIR"
        mkdir -p "$(dirname {output.summary})"

        # geNomad hard-errors on an empty fasta; zero contigs surviving
        # filter_contigs is legitimate, not a failure.
        if [ ! -s {input.contigs} ]; then
            printf "seq_name\\tlength\\ttopology\\tcoordinates\\tn_genes\\tgenetic_code\\tvirus_score\\tfdr\\tn_hallmarks\\tmarker_enrichment\\ttaxonomy\\n" > {output.summary}
            exit 0
        fi

        genomad end-to-end --cleanup --splits 8 --threads {threads} \
            {input.contigs} "$OUT_DIR" {config[databases][genomad_db]}

        # --cleanup only removes geNomad's own tmp files. The per-module
        # dirs (annotate/, marker_classification/, ...: ~68 GB over the run)
        # are intermediates; collect_viral_positives reads only *_summary/.
        find "$OUT_DIR" -mindepth 1 -maxdepth 1 -type d \
            ! -name 'filtered_contigs_summary' -exec rm -rf {{}} +
        """

rule run_rdrpcatch:
    conda: get_env("TAX_CORE_ENV")
    threads: 8
    # 120 -> 240: now scans every coverage-filtered contig, not just the
    # mmseqs orphans it used to get.
    # 240 -> 1440 (2026-11-20): 193 samples finished inside 3.9 h, but
    # six with 1.0-4.7 M contigs TIMEOUT'd at 4 h twice (jobs 1162814-1162822);
    # four were on Zayed_HMM, the last of the 7 databases, and the two biggest
    # were still on NeoRdRp. runtime is not a rerun trigger, so this changes
    # scheduling only.
    resources: mem_gb=16, runtime=1440
    input:
        contigs=get_taxonomy_contigs,
        db=os.path.join(config["databases"]["rdrp_catch_hmm"], ".done")
    output:
        tsv="data/04.taxonomy/{family}/{sample}/rdrpcatch/rdrp_annotated.tsv"
    benchmark:
        "benchmarks/run_rdrpcatch/{family}_{sample}.tsv"
    # MAX_LEN is a shell variable, not a params: entry, on purpose
    shell:
        """
        OUT_DIR=data/04.taxonomy/{wildcards.family}/{wildcards.sample}/rdrpcatch/scan
        # scan.log lives outside OUT_DIR -- --overwrite makes rdrpcatch
        # recreate that directory mid-run, which would orphan a log file
        # tee'd to a path nested inside it.
        SCAN_LOG="$OUT_DIR.scan.log"
        mkdir -p "$(dirname "$OUT_DIR")"

        # Contigs >300 kb are dropped: RdRpCATCH hard-errors on them
        MAX_LEN=299999
        FILTERED="$(dirname "$OUT_DIR")/contigs_le_${{MAX_LEN}}.fasta"
        # On every exit: the filtered copy of the assembly and the scan/ dir
        # (~1-2 GB/sample) are scratch. Previously only removed on success /
        # by taxonomy_cleanup, so failed or killed runs left them behind.
        # SCAN_LOG is kept on failure for debugging.
        cleanup_rdrpcatch() {{
            rc=$?
            rm -rf "$FILTERED" "$OUT_DIR"
            if [ "$rc" -eq 0 ]; then rm -f "$SCAN_LOG"; fi
        }}
        trap cleanup_rdrpcatch EXIT
        # SLURM time limit / scancel send TERM: exit so the EXIT trap fires
        trap 'exit 143' TERM INT
        seqkit seq -M "$MAX_LEN" {input.contigs} -o "$FILTERED" 2>/dev/null
        N_IN=$(grep -c '^>' {input.contigs} || true)
        N_KEPT=$(grep -c '^>' "$FILTERED" || true)
        echo "[run_rdrpcatch] {wildcards.sample}: $N_KEPT of $N_IN contigs <= $MAX_LEN bp"
        if [ "$N_IN" -ne "$N_KEPT" ]; then
            echo "[run_rdrpcatch] dropped $((N_IN - N_KEPT)) contig(s) over $MAX_LEN bp"
        fi

        # Zero scannable contigs is legitimate (an empty or all-over-length
        # assembly), not a failure -- emit the same header-only placeholder
        # the no-hits branch below uses.
        if [ ! -s "$FILTERED" ]; then
            echo "Contig_name" > {output.tsv}
            rm -f "$FILTERED"
            exit 0
        fi

        set +e
        rdrpcatch scan --input "$FILTERED" --db-dir {config[databases][rdrp_catch_hmm]} \
            --seq-type nuc --cpus {threads} --overwrite --output "$OUT_DIR" 2>&1 | tee "$SCAN_LOG"
        set -e

        FILE=$(find "$OUT_DIR" -name "*_rdrpcatch_output_annotated.tsv" | head -n 1)
        # The HMM-search hits are written before rdrpcatch's final UpSet
        # plot, which can crash on an upsetplot/matplotlib mismatch
        # unrelated to the search itself (seen 2026-10-13 in the method
        # comparison) -- fall back to the pre-plot combined hits table.
        COMBINED=$(find "$OUT_DIR" -path "*tsv_files*" -name "*_combined.tsv" | head -n 1)
        # Strip ANSI escapes and collapse whitespace so rich's line-wrapped
        # "no hits" message can't dodge a plain grep.
        SCAN_LOG_FLAT=$(sed 's/\\x1b\\[[0-9;]*[a-zA-Z]//g' "$SCAN_LOG" | tr -s '[:space:]' ' ')
        if [ -f "$FILE" ]; then
            cp "$FILE" {output.tsv}
        elif [ -f "$COMBINED" ]; then
            cp "$COMBINED" {output.tsv}
        elif echo "$SCAN_LOG_FLAT" | grep -q "No hits found by RdRpCATCH"; then
            # Legitimate "nothing matched any HMM" exit, not a failure.
            echo "Contig_name" > {output.tsv}
        else
            echo "Error: RdRpCATCH output not found, and scan.log doesn't show the known no-hits exit -- check $SCAN_LOG" >&2
            exit 1
        fi
        rm -f "$SCAN_LOG" "$FILTERED"
        """

rule run_checkv:
    # Quality/completeness annotation on the union of all three tools'
    # positives (rules/05.merge_tax.smk: collect_viral_positives) -- not a
    # filter, every positive is kept regardless of CheckV's verdict.
    conda: get_env("TAX_QC_ENV")
    threads: 8
    # Not yet sacct-validated.
    resources: mem_gb=16, runtime=120
    input:
        fasta="data/05.final_selection/{family}/{sample}/viral_positives.fasta",
        db=os.path.join(config["databases"]["checkv_db"], ".done")
    output:
        tsv="data/04.taxonomy/{family}/{sample}/checkv/quality_summary.tsv"
    benchmark:
        "benchmarks/run_checkv/{family}_{sample}.tsv"
    shell:
        # checkv download_database creates a versioned subdir rather than
        # downloading straight into the configured path; resolve it here.
        """
        OUT_DIR=$(dirname {output.tsv})
        if [ ! -s {input.fasta} ]; then
            mkdir -p "$OUT_DIR"
            printf "contig_id\\tcontig_length\\tprovirus\\tproviral_length\\tgene_count\\tviral_genes\\thost_genes\\tcheckv_quality\\tmiuvig_quality\\tcompleteness\\tcompleteness_method\\tcontamination\\tkmer_freq\\twarnings\\n" > {output.tsv}
            exit 0
        fi
        DB_DIR=$(find {config[databases][checkv_db]} -maxdepth 1 -type d -name 'checkv-db-*' | head -n1)
        checkv end_to_end {input.fasta} "$OUT_DIR" -d "$DB_DIR" -t {threads} --restart
        """
