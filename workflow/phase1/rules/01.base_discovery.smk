# =========================================================
# STEP 1 — DB QUERY
# =========================================================

# Path definitions
DIR_CENTROID = "data/01.sOTUs_centroid/{mode}"
PATH_RUNS = os.path.join(DIR_CENTROID, "candidate_runs.txt")
PATH_UNANCHORED = os.path.join(DIR_CENTROID, "unanchored_palm_ids.txt")
PATH_REF_FASTA = os.path.join(DIR_CENTROID, "ncbi_nucleotide_reference.fasta")
PATH_COMBINED = os.path.join(DIR_CENTROID, "combined.fasta")
PATH_MATRIX = os.path.join(DIR_CENTROID, "unique.csv")  # unused

checkpoint centroid_step1_query:
    input: scenarios="configs/phase1_sra_discovery.yaml"
    output:
        runs="data/01.sOTUs_centroid/{mode}/candidate_runs.txt",
        unanchored_meta="data/01.sOTUs_centroid/{mode}/unanchored_palm_ids.txt",
        ref_fasta="data/01.sOTUs_centroid/{mode}/ncbi_nucleotide_reference.fasta"
    conda: "../" + config["tools"]["DISCOVERY_ENV"]
    script: "../scripts/01.centroid_query.py"
    