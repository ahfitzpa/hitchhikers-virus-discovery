# =========================================================================
# STEP 1 — MERGE SRA DISCOVERY METADATA MATRIX
# =========================================================================
checkpoint merge_sra_metadata:
    conda:
        "../" + config["tools"]["DISCOVERY_ENV"]
    input:
        filtered_csv = "data/02.SRA_metadata/{mode}/filtered_unique.csv",
        metadata_tsv = "data/02.SRA_metadata/{mode}/unified_metadata.tsv.gz"
    output:
        merged_csv    = "data/03.merged_metadata/{mode}/merged_metadata.csv",
        shortlist_csv = "data/03.merged_metadata/{mode}/shortlist_summary.csv",
        manifest      = "data/03.merged_metadata/{mode}/download_manifest.txt"
    log:
        "logs/metadata_merge/{mode}_merge_execution.log"
    script:
        "../scripts/04.metadata_merge.py"


# =========================================================================
# STEP 2 — NON-INTERACTIVE AUTOMATED GATEKEEPER
# =========================================================================
rule automated_review_gate:
    input:
        shortlist = "data/03.merged_metadata/{mode}/shortlist_summary.csv",
        manifest  = "data/03.merged_metadata/{mode}/download_manifest.txt"
    output:
        flag = "data/03.merged_metadata/{mode}/approved.flag"
    shell:
        """
        # Ensure the directory exists before logging
        mkdir -p data/03.merged_metadata/{wildcards.mode}
        
        echo "[REVIEW GATE] Successfully processed candidate datasets for mode: {wildcards.mode}"
        
        # Calculate total runs
        TOTAL_RUNS=$(tail -n +2 {input.manifest} | wc -l || echo 0)
        echo "[REVIEW GATE] Candidates compiled into manifest target pool: $TOTAL_RUNS runs."
        
        if [ "$TOTAL_RUNS" -eq 0 ]; then
            echo "[WARN] Manifest file is empty. Phase 2 assembly will be marked idle."
        fi
        
        # Touch flag to signal completion
        touch {output.flag}
        """