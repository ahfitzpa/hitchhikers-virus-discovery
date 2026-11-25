"""Query Serratus/PalmDB for viral family sOTU centroids matching a discovery scenario.

Invoked via Snakemake `script:` by `rules/01.base_discovery.smk`. Reads the active
scenario configuration (`snakemake.config["scenarios"][snakemake.wildcards.mode]`),
queries the Serratus PostgreSQL index and NCBI Entrez taxonomy, and writes filtered
sOTU tables and anchor sequences to `snakemake.output`.
"""

import os
import re
import sys
import tempfile
import time
import xml.etree.ElementTree as ET
from pathlib import Path

import pandas as pd
import psycopg
import requests

NCBI_TIMEOUT_SEC = 30
NCBI_MAX_RETRIES = 3


def _atomic_write_text(dest: str | Path, text: str) -> None:
    """Write text to a sibling temporary file and atomically rename into place."""
    dest_path = Path(dest)
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{dest_path.name}.", suffix=".tmp", dir=dest_path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as fh:
            fh.write(text)
        os.replace(tmp_name, dest_path)
    except BaseException:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)
        raise


def _atomic_write_csv(df: pd.DataFrame, dest: str | Path, **kwargs) -> None:
    """Write a DataFrame to a sibling temporary file and atomically rename into place."""
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


def _http_get_with_retry(url: str, params: dict[str, object]) -> requests.Response:
    """GET with explicit timeout, status check, and bounded exponential backoff."""
    last_exc: Exception | None = None
    for attempt in range(1, NCBI_MAX_RETRIES + 1):
        try:
            resp = requests.get(url, params=params, timeout=NCBI_TIMEOUT_SEC)
            resp.raise_for_status()
            return resp
        except requests.RequestException as exc:
            last_exc = exc
            if attempt < NCBI_MAX_RETRIES:
                wait_sec = 2 ** (attempt - 1)
                print(
                    f"[WARNING] NCBI request failed (attempt {attempt}/{NCBI_MAX_RETRIES}): {exc}; "
                    f"retrying in {wait_sec}s",
                    file=sys.stderr,
                )
                time.sleep(wait_sec)
    raise RuntimeError(f"NCBI request failed after {NCBI_MAX_RETRIES} attempts: {url}") from last_exc


config = snakemake.config
wildcards = snakemake.wildcards
output = snakemake.output

# STRICT INFERENCE: Pull scenario parameters directly from YAML config.
current_mode = wildcards.mode
scen = config["scenarios"][current_mode]

family = str(scen["family"])
filt = scen.get("filtering", {})
min_read_coverage = int(filt["min_sra_read_coverage"])

# Database Setup (public read-only Serratus RDS; allow env override)
db_config = dict(config["database"])
if "name" in db_config:
    db_config["dbname"] = db_config.pop("name")
db_config.setdefault("password", os.environ.get("SERRATUS_DB_PASSWORD", "serratus"))

# =========================================================
# RUN PARAMETERIZED QUERY
# =========================================================
mmseqs_params = scen.get("mmseqs_params", {})
min_seq_id = float(mmseqs_params.get("min_seq_id", 0.45)) * 100.0
min_cov = float(mmseqs_params.get("coverage_threshold", 0.5)) * 100.0

# Optional run-level filters. Absent from a scenario's config => not applied,
# so the novel/surrogate/host_expansion scenarios keep their previous behaviour.
min_run_pid = filt.get("min_run_percent_identity")
require_geo = bool(filt.get("require_geo", False))
exclude_cell_culture = bool(filt.get("exclude_cell_culture", False))
geo_include = [str(x) for x in (filt.get("geo_include") or [])]
sample_type_include = [str(x) for x in (filt.get("sample_type_include") or [])]

# Lab-passage / cell-line markers. These samples carry the LAB's address as
# geo_loc_name, which is actively misleading for phylogeography.
CELL_CULTURE_RE = (
    r"(cell[ _-]?culture|cell[ _-]?line|passage|serial[ _-]?passag|"
    r"vero|huh7|hek293|bhk|c6/36|s2 cells|in vitro|synthetic)"
)
TISSUE_CELL_CULTURE_RE = r"(cell[ _-]?culture|cell[ _-]?line|derived cell)"

# 1. Collapse palm_tax to one row per sOTU first so joining palm_sra2 does not
#    multiply rows by the number of palm_ids per sOTU.
tax_cte = """WITH tax AS (
    SELECT sotu,
           max(percent_identity) AS best_pid,
           max(pp_cov)           AS best_pp_cov
    FROM palm_tax
    WHERE tax_family = %s
    GROUP BY sotu
)"""
sql_params: list[object] = [family]

# 2. ps.percent_identity is the run's OWN contig-vs-palmprint identity, a
#    different number from the sOTU assignment confidence above. min_seq_id
#    gates the latter; min_run_percent_identity (if set) gates the former.
clauses: list[str] = [
    "ps.coverage >= %s",
    "tax.best_pid >= %s",
    "tax.best_pp_cov >= %s",
]
sql_params.extend([min_read_coverage, min_seq_id, min_cov])

# 3. Geography and sample-type filters need srarun + biosample, which live in
#    the same Serratus DB -- no separate BioSample efetch required.
needs_srarun = bool(
    require_geo or exclude_cell_culture or geo_include or sample_type_include
)
joins = [
    "FROM palm_sra2 ps",
    "JOIN tax ON tax.sotu = ps.sotu",
]
if needs_srarun:
    joins.append("JOIN srarun sr ON sr.run = ps.run_id")
    joins.append("LEFT JOIN biosample bs ON bs.run_id = sr.sample")

if min_run_pid is not None:
    clauses.append("ps.percent_identity >= %s")
    sql_params.append(float(min_run_pid))

if require_geo:
    clauses.append(
        "bs.geo_text_all->>'geo_loc_name' IS NOT NULL "
        "AND btrim(lower(bs.geo_text_all->>'geo_loc_name')) NOT IN "
        "('', 'not determined', 'missing', 'na', 'n/a', 'uncalculated', "
        "'not applicable', 'not collected', 'unknown', 'none')"
    )

if geo_include:
    geo_re = "(" + "|".join(re.escape(g) for g in geo_include) + ")"
    clauses.append("bs.geo_text_all->>'geo_loc_name' ~* %s")
    sql_params.append(geo_re)

if sample_type_include:
    st_re = "(" + "|".join(sample_type_include) + ")"
    clauses.append(
        "(sr.scientific_name ~* %s OR sr.library_name ~* %s OR sr.sample_name ~* %s "
        "OR EXISTS (SELECT 1 FROM biosample_tissue bt "
        "WHERE bt.biosample_id = bs.biosample_id AND bt.tissue ~* %s))"
    )
    sql_params.extend([st_re, st_re, st_re, st_re])

if exclude_cell_culture:
    for col in ("sr.scientific_name", "sr.library_name", "sr.sample_name"):
        clauses.append(f"({col} IS NULL OR {col} !~* %s)")
        sql_params.append(CELL_CULTURE_RE)

    clauses.append(
        "NOT EXISTS (SELECT 1 FROM biosample_tissue bt "
        "WHERE bt.biosample_id = bs.biosample_id "
        "AND bt.tissue ~* %s)"
    )
    sql_params.append(TISSUE_CELL_CULTURE_RE)

# Explicit BioProject denylist: auditable backstop for lab-derived projects.
exclude_bioprojects = [str(p) for p in (filt.get("exclude_bioprojects") or [])]
if exclude_bioprojects and not needs_srarun:
    raise ValueError(
        "exclude_bioprojects requires require_geo or "
        "exclude_cell_culture so that srarun is joined"
    )
for proj in exclude_bioprojects:
    clauses.append("sr.bio_project <> %s")
    sql_params.append(proj)

max_runs = scen.get("max_runs") or scen.get("test_limit")

query = "\n".join([
    tax_cte,
    "SELECT ps.run_id, ps.sotu, ps.palm_id, max(ps.coverage) AS coverage",
    *joins,
    "WHERE " + "\n  AND ".join(clauses),
    "GROUP BY ps.run_id, ps.sotu, ps.palm_id",
    "ORDER BY max(ps.coverage) DESC, ps.run_id ASC",
])

print(f"[INFO] Querying mode: {current_mode} | Family: {family}")
print(
    f"[INFO] filters: cov>={min_read_coverage} tax_pid>={min_seq_id} "
    f"pp_cov>={min_cov} run_pid>={min_run_pid} require_geo={require_geo} "
    f"exclude_cell_culture={exclude_cell_culture} max_runs={max_runs} "
    f"geo_include={len(geo_include)} sample_type_include={len(sample_type_include)}"
)

with psycopg.connect(**db_config) as conn:
    with conn.cursor() as cur:
        cur.execute(query, tuple(sql_params))
        rows = cur.fetchall()
        colnames = [desc[0] for desc in cur.description] if cur.description else []

print(f"[INFO] {len(rows)} hit rows returned ({len({r[0] for r in rows})} distinct runs)")

# 4. Cap on DISTINCT RUNS, not hit rows, so all (sotu, palm_id) hits of a kept run survive together.
if max_runs and rows:
    keep, seen = [], set()
    for r in rows:
        if r[0] not in seen:
            if len(seen) >= int(max_runs):
                continue
            seen.add(r[0])
        keep.append(r)
    print(f"[INFO] max_runs={max_runs}: kept {len(keep)} hit rows across {len(seen)} runs")
    rows = keep

# 5. Drop ordering-only coverage column so downstream schema stays run_id,sotu,palm_id.
if colnames and colnames[-1] == "coverage":
    rows = [r[:-1] for r in rows]
    colnames = colnames[:-1]

if rows:
    df_out = pd.DataFrame(rows, columns=colnames)
    _atomic_write_csv(df_out, output.runs, index=False)
else:
    _atomic_write_text(output.runs, "run_id,sotu,palm_id\n")

# =========================================================
# REF FASTA & METADATA
# =========================================================
taxid = scen.get("anchor", {}).get("taxid") or scen.get("reference_anchor_taxid")

if taxid:
    esearch_resp = _http_get_with_retry(
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi",
        params={"db": "nuccore", "term": f"txid{taxid}[Organism]", "retmax": 200},
    )
    ids = [x.text for x in ET.fromstring(esearch_resp.content).findall(".//Id") if x.text]
    if ids:
        efetch_resp = _http_get_with_retry(
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
            params={"db": "nuccore", "id": ",".join(ids), "rettype": "fasta", "retmode": "text"},
        )
        _atomic_write_text(output.ref_fasta, efetch_resp.text)
    else:
        print(f"[WARNING] NCBI esearch returned 0 sequences for txid{taxid}", file=sys.stderr)
        _atomic_write_text(output.ref_fasta, "")
else:
    _atomic_write_text(output.ref_fasta, "")

_atomic_write_csv(pd.DataFrame(columns=["palm_id"]), output.unanchored_meta, index=False, header=False)