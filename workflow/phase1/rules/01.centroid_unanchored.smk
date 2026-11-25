# Unanchored Pathway: Skips MMseqs, directly transforms discovery candidates

# =========================================================
# STEP 1 — Convert palm db output to hits.tsv
# =========================================================
rule transform_db_to_hits:
    input: runs="data/01.sOTUs_centroid/{mode}/candidate_runs.txt"
    output: hits="data/01.sOTUs_centroid/{mode}/hits_novel.tsv"
    run:
        import pandas as pd  # <--- ADD THIS
        df = pd.read_csv(input.runs)
        df.assign(target='novel_discovery', pident=100.0, alnlen=100, evalue=0.0)[
            ['run_id', 'target', 'pident', 'alnlen', 'evalue']
        ].to_csv(output.hits, sep='\t', header=False, index=False)

# =========================================================
# STEP 2 — Extract SRA IDs for metadata fetching
# =========================================================
rule extract_sra_ids:
    input: 
        hits="data/01.sOTUs_centroid/{mode}/hits_novel.tsv"  # Changed from hits.tsv
    output: 
        "data/01.sOTUs_centroid/{mode}/sra_ids.txt"
    script:
        "../scripts/extract_sra_ids.py"

rule build_matrix_novel:
    # unique_novel.csv, not unique.csv: the novel matrix is built separately.
    input: hits="data/01.sOTUs_centroid/{mode}/hits_novel.tsv"
    output: "data/01.sOTUs_centroid/{mode}/unique_novel.csv"
    run:
        import pandas as pd  # <--- ADD THIS
        pd.DataFrame(columns=["sotu", "palm_id", "run_id"]).to_csv(output[0], index=False)