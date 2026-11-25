"""Build the anchored sOTU-by-SRA alignment matrix from filtered MMseqs2 hits.

Invoked via Snakemake `script:` by `rules/01.centroid_anchored.smk`. Joins
`snakemake.input.hits` and `snakemake.input.unanchored` with Serratus metadata
to produce a downstream-compatible matrix (emitting mandatory header columns
even when zero hits pass filters).
"""

import os
import sys
import tempfile
from pathlib import Path

import pandas as pd
import psycopg  # Psycopg v3


def _atomic_write_csv(df: pd.DataFrame, dest: str | Path, **kwargs) -> None:
    """Write a DataFrame to a sibling temp file and atomically rename into place."""
    dest_path = Path(dest)
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{dest_path.name}.", suffix=".tmp", dir=dest_path.parent)
    os.close(fd)
    try:
        df.to_csv(tmp_name, **kwargs)
        os.replace(tmp_name, dest_path)
    except BaseException:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)
        raise


hits_file = Path(snakemake.input.hits)
unanchored_file = Path(snakemake.input.unanchored)
out = Path(snakemake.output[0])
params = snakemake.params

MANDATORY_COLUMNS = ["sotu", "palm_id", "run_id"]
column_names = ["query", "target", "pident", "alnlen", "evalue"]

if hits_file.exists() and hits_file.stat().st_size > 0:
    hits = pd.read_csv(hits_file, sep="\t", header=None, names=column_names)
else:
    hits = pd.DataFrame(columns=column_names)

sra_ids: list[str] = []
if not hits.empty and "query" in hits.columns:
    hits["query"] = hits["query"].astype(str)
    sra_ids = hits["query"].str.split("_").str[0].dropna().unique().tolist()

palm_ids: list[str] = []
if unanchored_file.exists() and unanchored_file.stat().st_size > 0:
    try:
        palm_ids = pd.read_csv(unanchored_file, header=None)[0].dropna().astype(str).tolist()
    except pd.errors.EmptyDataError:
        palm_ids = []

db_config = dict(params.db)
if "name" in db_config:
    db_config["dbname"] = db_config.pop("name")
elif "database" in db_config:
    db_config["dbname"] = db_config.pop("database")
db_config.setdefault("password", os.environ.get("SERRATUS_DB_PASSWORD", "serratus"))

if sra_ids:
    sql = """
        SELECT pt.sotu, pt.palm_id, ps.run_id
        FROM palm_tax pt
        JOIN palm_sra2 ps ON pt.sotu = ps.sotu
        WHERE ps.run_id = ANY(%s)
    """
    query_param = (sra_ids,)
elif palm_ids:
    sql = """
        SELECT pt.sotu, pt.palm_id, ps.run_id
        FROM palm_tax pt
        JOIN palm_sra2 ps ON pt.sotu = ps.sotu
        WHERE pt.palm_id = ANY(%s)
    """
    query_param = (palm_ids,)
else:
    print("[INFO] No sequence tracking IDs found to match. Creating structural empty fallback layout.")
    _atomic_write_csv(pd.DataFrame(columns=MANDATORY_COLUMNS), out, index=False)
    sys.exit(0)

with psycopg.connect(**db_config) as conn:
    with conn.cursor() as cur:
        cur.execute(sql, query_param)
        rows = cur.fetchall()
        colnames = [desc[0] for desc in cur.description] if cur.description else []

if rows:
    df = pd.DataFrame(rows, columns=colnames)
else:
    df = pd.DataFrame(columns=MANDATORY_COLUMNS)

if "run_id" not in df.columns:
    df["run_id"] = df["Run"] if "Run" in df.columns else None

_atomic_write_csv(df, out, index=False)
print(f"[INFO] Successfully exported matrix tracking frame to: {out}")
