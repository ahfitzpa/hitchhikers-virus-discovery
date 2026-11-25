#!/bin/bash
set -euo pipefail

mode="$1"
run_file="$2"
filtered_csv="$3"
metadata_out="$4"

mkdir -p "$(dirname "$filtered_csv")"
mkdir -p "$(dirname "$metadata_out")"

# 1. VALIDATE
if [ ! -s "$run_file" ]; then
    echo "[WARN] No run IDs found for mode: $mode"
    echo "run_id,sotu,palm_id" > "$filtered_csv"
    echo "Run,BioProject,ScientificName" | gzip > "$metadata_out"
    exit 0
fi

# 2. FETCH STRUCTURED DATA FROM PALM DB
sql_ids=$(tr '\n' ',' < "$run_file" | sed 's/,$//' | sed "s/,/','/g" | sed "s/^/'/" | sed "s/$/'/")

echo "[INFO] Querying palm_db for SRA metadata and sOTU associations..."

# Fetch from DB and save to filtered_csv
if ! psql -h "$DB_HOST" -d "$DB_NAME" -U "$DB_USER" -t -A -F ',' -c "
    SELECT ps.run_id, pt.sotu, pt.palm_id 
    FROM palm_tax pt
    JOIN palm_sra2 ps ON pt.sotu = ps.sotu
    WHERE ps.run_id IN ($sql_ids)
    AND pt.pp_cov > 50 
    AND pt.centroid = 't' 
    AND ps.coverage > 100
" > "$filtered_csv"; then
    echo "[ERROR] Database query failed."
    exit 1
fi

# Add the header
sed -i '1i run_id,sotu,palm_id' "$filtered_csv"

# 3. FETCH NCBI METADATA
# Define a standard header that matches the rest of your pipeline
HEADER="Run,ReleaseDate,LoadDate,spots,bases,spots_with_mates,avgLength,size_MB,AssemblyName,download_path,Experiment,LibraryName,LibraryStrategy,LibrarySelection,LibrarySource,LibraryLayout,InsertSize,InsertDev,Platform,Model,SRAStudy,BioProject,Study_Pubmed_id,ProjectID,Sample,BioSample,SampleType,TaxID,ScientificName,SampleName,g1k_pop_code,source,g1k_analysis_group,Subject_ID,Sex,Disease,Tumor,Affection_Status,Analyte_Type,Histological_Type,Body_Site,CenterName,Submission,dbgap_study_accession,Consent,RunHash,ReadHash"

tmp=$(mktemp)
echo "$HEADER" > "$tmp"

xargs -n 200 < "$run_file" | while read -r chunk; do
    ids=$(echo "$chunk" | tr ' ' ',')
    echo "[DEBUG] Fetching IDs: $ids" # ADD THIS
    if data=$(efetch -db sra -id "$ids" -format runinfo 2>/dev/null); then
        echo "[DEBUG] Got data for $ids" # ADD THIS
        echo "$data" | tail -n +2 >> "$tmp"
    else
        echo "[DEBUG] efetch failed for $ids" # ADD THIS
    fi
done

# Compress the final consistent file
gzip -c "$tmp" > "$metadata_out"
rm -f "$tmp"
echo "[INFO] Metadata fetch complete for mode: $mode"