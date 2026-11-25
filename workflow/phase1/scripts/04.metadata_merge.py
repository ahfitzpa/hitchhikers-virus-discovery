"""Merge Serratus sOTU discovery tables with NCBI SRA run metadata.

Invoked via Snakemake `script:` by `rules/04.metadata_merge.smk`. Combines
per-scenario Serratus hit matrices with fetched SRA metadata XML/TSV records
and emits the unified annotated candidate table for the Phase 1 review gate.
"""

import os
import sys
import yaml
import pandas as pd
import gzip

# 1. SETUP SNAKEMAKE CONTEXT
if hasattr(snakemake, "log") and snakemake.log:
    log_file = open(str(snakemake.log), "w")
    sys.stdout = log_file
    sys.stderr = log_file

try:
    config = snakemake.config
    wildcards = snakemake.wildcards
    output = snakemake.output

    # DYNAMIC MODE IDENTIFICATION
    current_mode = wildcards.mode 
    scen = config.get('scenarios', {}).get(current_mode, {})

    if not scen:
        raise ValueError(f"CRITICAL: Mode '{current_mode}' not found in configuration scenarios.")

    family = scen.get('family', 'Hepeviridae')
    cfg = scen.get('filtering', {}) 

    print(f"[INFO] Initializing metadata merge for mode: [{current_mode}] (Family: {family})", file=sys.stderr)

    # 2. LOAD UNIQUE RUNS (Now expects headers: run_id, sotu, palm_id)
    filtered_input_path = snakemake.input.filtered_csv
    df_unique = pd.read_csv(filtered_input_path)
    
    # Inject family from config
    df_unique['tax_family'] = family 

    # 3. STREAM LOAD DATASETS (NCBI METADATA)
    with gzip.open(snakemake.input.metadata_tsv, 'rt') as f_delim_test:
        first_raw_line = f_delim_test.readline()
    
    chosen_sep = ',' if ',' in first_raw_line and '\t' not in first_raw_line else '\t'
    df_ncbi = pd.read_csv(snakemake.input.metadata_tsv, sep=chosen_sep, compression='gzip')
    df_ncbi.columns = df_ncbi.columns.str.lower()

    # Join Keys
    df_unique['run_key'] = df_unique['run_id'].astype(str).str.strip()
    if 'run' in df_ncbi.columns:
        df_ncbi['run_key'] = df_ncbi['run'].astype(str).str.strip()
    else:
        df_ncbi['run_key'] = df_ncbi.iloc[:, 0].astype(str).str.strip()

    # 4. RELATIONAL JOIN
    df_filtered = pd.merge(df_unique, df_ncbi, on='run_key', how='inner')
    os.makedirs(os.path.dirname(snakemake.output.merged_csv), exist_ok=True)
    df_filtered.to_csv(snakemake.output.merged_csv, index=False)
    
    # 5. DYNAMIC SHORTLIST FILTERING
    min_bioproj = int(cfg.get('min_bioprojects_per_sotu', 1))
    min_hosts = int(cfg.get('min_unique_hosts', 1))

    # We trust the 'sotu' column exists now because it comes from your 01.py output
    df_candidates = df_filtered.copy()
    
    # INDENTATION MUST BE 4 SPACES HERE
    host_col = 'host'
    if 'host' not in df_filtered.columns:
        if 'scientificname' in df_filtered.columns:
            host_col = 'scientificname'
        elif 'organism' in df_filtered.columns:
            host_col = 'organism'
        else:
            df_filtered['host'] = 'unknown'
            host_col = 'host'

    stats = df_filtered.groupby('sotu').agg(
        bioproject_count=('bioproject', 'nunique'),
        host_count=(host_col, 'nunique') 
    ).reset_index()
    
    valid = stats[(stats['bioproject_count'] >= min_bioproj) & (stats['host_count'] >= min_hosts)]['sotu']
    df_final = df_candidates[df_candidates['sotu'].isin(valid)].copy()

    df_final.to_csv(snakemake.output.shortlist_csv, index=False)

    # 6. PORTABLE DOWNLOAD MANIFEST GENERATION (Raw host names)
    use_bioproj = 'bioproject' if 'bioproject' in df_ncbi.columns else 'run'
    use_host = 'host' if 'host' in df_ncbi.columns else 'scientificname'

    tokens = []
    if not df_final.empty and use_bioproj in df_final.columns:
        target_bioprojects = df_final[use_bioproj].unique()
        df_download_pool = df_ncbi[df_ncbi[use_bioproj].isin(target_bioprojects)].copy()
        
        if not df_download_pool.empty:
            # Raw host values, deliberately not normalised here.
            df_download_pool['host_clean'] = df_download_pool[use_host].fillna('unknown').astype(str)
            df_download_pool['task_token'] = (df_download_pool['run'].astype(str) + ":" + 
                                            df_download_pool[use_bioproj].astype(str) + ":" + 
                                            df_download_pool['host_clean'].astype(str))
            
            tokens = df_download_pool['task_token'].drop_duplicates().tolist()
        
    with open(str(snakemake.output.manifest), 'w') as f_out:
        for token in tokens:
            f_out.write(f"{token}\n")

except Exception as e:
    print(f"[ERROR] Script crashed: {str(e)}", file=sys.stderr)
    raise e