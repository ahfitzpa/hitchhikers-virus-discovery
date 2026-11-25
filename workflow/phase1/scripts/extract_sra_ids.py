"""Extract unique SRA run accessions from headerless MMseqs2 hit tables.

Invoked via Snakemake `script:` by `rules/01.centroid_anchored.smk` and
`rules/01.centroid_unanchored.smk`. Parses the `<run_id>_...` prefix from the
first column of `snakemake.input.hits` and writes one accession per line to
`snakemake.output[0]`.
"""

import pandas as pd
import sys

# Snakemake injects this object automatically
hits_path = snakemake.input.hits
output_path = snakemake.output[0]

try:
    # 1. Read the hits file (Assuming it remains a tab-delimited file without headers)
    # If this hits.tsv file structure hasn't changed, this part is fine.
    df = pd.read_csv(hits_path, sep="\t", 
                     names=["query", "target", "pident", "alnlen", "evalue"],
                     header=None) 
    
    # 2. Extract IDs (The 'query' column in hits.tsv still likely looks like 'ERR12345_...')
    # This remains the same to extract the SRA run_id from the MMseqs query
    ids = df["query"].str.split("_").str[0].dropna().unique()
    
    # 3. Write output in a clean format
    with open(output_path, "w") as f:
        # We keep this as a simple list because it is consumed by downstream logic
        for run_id in ids:
            f.write(str(run_id) + "\n")
            
except pd.errors.EmptyDataError:
    with open(output_path, "w") as f:
        pass
except Exception as e:
    print(f"Error parsing {hits_path}: {e}")
    sys.exit(1)