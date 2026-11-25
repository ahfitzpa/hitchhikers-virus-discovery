"""Build the anchored sOTU-by-SRA alignment matrix from filtered MMseqs2 hits.

Invoked via Snakemake `script:` by `rules/01.centroid_anchored.smk`. Joins
`snakemake.input.hits` and `snakemake.input.unanchored` with Serratus metadata
to produce a downstream-compatible matrix (emitting mandatory header columns
even when zero hits pass filters).
"""

import os
import sys
import pandas as pd
import psycopg  # Migrated to Psycopg v3

# Inputs and Outputs from Snakemake mapping
hits_file = snakemake.input.hits
unanchored_file = snakemake.input.unanchored
out = snakemake.output[0]
params = snakemake.params

# Establish structure expectations for downstream scripts if things are empty
MANDATORY_COLUMNS = ["sotu", "palm_id", "run_id"]

# =========================================================
# 1. LOAD HITS
# =========================================================
column_names = ["query", "target", "pident", "alnlen", "evalue"]

if os.path.exists(hits_file) and os.path.getsize(hits_file) > 0:
    # Explicitly use the column mapping since hits.tsv has no native headers
    hits = pd.read_csv(hits_file, sep="\t", header=None, names=column_names)
else:
    hits = pd.DataFrame(columns=column_names)

# =========================================================
# 2. EXTRACT RUN IDS
# =========================================================
sra_ids = []
if not hits.empty and "query" in hits.columns:
    # FIXED LOGIC: Safely cast to string, split by underscore, and isolate SRA IDs
    hits["query"] = hits["query"].astype(str)
    sra_ids = hits["query"].str.split("_").str[0].dropna().unique().tolist()

# =========================================================
# 3. OPTIONAL: UNANCHORED PATH (fallback)
# =========================================================
palm_ids = []
if os.path.exists(unanchored_file) and os.path.getsize(unanchored_file) > 0:
    try:
        palm_ids = pd.read_csv(unanchored_file, header=None)[0].dropna().tolist()
    except Exception:
        palm_ids = []

# =========================================================
# 4. DB QUERY PREPARATION & CLEAN FALLBACKS
# =========================================================
db_config = dict(params.db)

# FIXED FOR PSYCOPG v3: Map "name" parameter keys cleanly to "dbname"
if "name" in db_config:
    db_config["dbname"] = db_config.pop("name")
elif "database" in db_config:
    db_config["dbname"] = db_config.pop("database")

db_config.setdefault("password", "serratus")

if sra_ids:
    ids = ",".join([f"'{x}'" for x in sra_ids])
    sql = f"""
        SELECT pt.sotu, pt.palm_id, ps.run_id
        FROM palm_tax pt
        JOIN palm_sra2 ps ON pt.sotu = ps.sotu
        WHERE ps.run_id IN ({ids})
    """
elif palm_ids:
    ids = ",".join([f"'{x}'" for x in palm_ids])
    sql = f"""
        SELECT pt.sotu, pt.palm_id, ps.run_id
        FROM palm_tax pt
        JOIN palm_sra2 ps ON pt.sotu = ps.sotu
        WHERE pt.palm_id IN ({ids})
    """
else:
    print("[INFO] No sequence tracking IDs found to match. Creating structural empty fallback layout.")
    pd.DataFrame(columns=MANDATORY_COLUMNS).to_csv(out, index=False)
    sys.exit(0)

# =========================================================
# 5. EXECUTE MATRIX EXTRACTION
# =========================================================
try:
    with psycopg.connect(**db_config) as conn:
        with conn.cursor() as cur:
            cur.execute(sql)
            rows = cur.fetchall()
            colnames = [desc[0] for desc in cur.description] if cur.description else []

    # Build dataframe from extracted rows
    if rows:
        df = pd.DataFrame(rows, columns=colnames)
    else:
        df = pd.DataFrame(columns=MANDATORY_COLUMNS)

    # Corrected logic to ensure missing run_id column fallback scales cleanly
    if "run_id" not in df.columns:
        if "Run" in df.columns:
            df["run_id"] = df["Run"]
        else:
            # Force empty column generation to guarantee merge script stability
            df["run_id"] = None

    df.to_csv(out, index=False)
    print(f"[INFO] Successfully exported matrix tracking frame to: {out}")

except Exception as e:
    print(f"[CRITICAL ERROR] Failed during matrix construction query: {str(e)}", file=sys.stderr)
    # Generate backup structure so cluster steps do not freeze completely
    pd.DataFrame(columns=MANDATORY_COLUMNS).to_csv(out, index=False)
    raise e