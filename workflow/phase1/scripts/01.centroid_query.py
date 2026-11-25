import os
import re
import sys
import yaml
import psycopg
import pandas as pd
import requests
import xml.etree.ElementTree as ET

config = snakemake.config
wildcards = snakemake.wildcards
output = snakemake.output

# STRICT INFERENCE: No hardcoded fallbacks. 
# It pulls exactly what is defined in your YAML config.
current_mode = wildcards.mode
scen = config["scenarios"][current_mode] 

# The rest of your script follows...
family = scen['family']
min_read_coverage = int(scen['filtering']['min_sra_read_coverage'])

# Retrieve family from config instead of wildcards
family = scen.get('family', 'Hepeviridae')

# Extract Filtering Parameters
min_read_coverage = int(scen.get('filtering', {}).get('min_sra_read_coverage', 100))

# Database Setup
db_config = dict(config["database"])
if "name" in db_config: db_config["dbname"] = db_config.pop("name")
db_config.setdefault("password", "serratus")

# =========================================================
# RUN QUERY
# =========================================================
# Now uses the 'family' variable derived from your config/mode mapping
# Extract MMseqs parameters from config
mmseqs_params = scen.get("mmseqs_params", {})
min_seq_id = float(mmseqs_params.get("min_seq_id", 0.45)) * 100 # Convert 0.45 to 45
min_cov = float(mmseqs_params.get("coverage_threshold", 0.5)) * 100

# Optional run-level filters. Absent from a scenario's config => not applied,
# so the novel/surrogate/host_expansion scenarios keep their previous behaviour.
filt = scen.get("filtering", {})
min_run_pid = filt.get("min_run_percent_identity")
require_geo = bool(filt.get("require_geo", False))
exclude_cell_culture = bool(filt.get("exclude_cell_culture", False))
# 2a: geo_include is an allowlist of country/region substrings; sample_type_include
#     an allowlist of sample-type patterns (wastewater and the like). Both need
#     srarun + biosample, so both force the join below.
geo_include = [str(x) for x in (filt.get("geo_include") or [])]
sample_type_include = [str(x) for x in (filt.get("sample_type_include") or [])]

# Lab-passage / cell-line markers. These samples carry the LAB's address as
# geo_loc_name, which is actively misleading for phylogeography.
CELL_CULTURE_RE = (
    r"(cell[ _-]?culture|cell[ _-]?line|passage|serial[ _-]?passag|"
    r"vero|huh7|hek293|bhk|c6/36|s2 cells|in vitro|synthetic)"
)

# 1. palm_tax holds MANY rows per sotu (144 for the WNV sOTU u4485), because a
#    sotu is a cluster of palm_ids each with its own taxon assignment. Joining
#    it directly on sotu multiplies palm_sra2 by that factor -- 2,914 rows
#    become 419,616 for u4485 alone -- and the old SELECT DISTINCT then threw
#    the duplicates away again. Collapsing to one row per sotu first keeps the
#    identical result set and makes the filter's real meaning explicit: it is
#    "this sOTU has SOME confident family assignment", a property of the
#    reference cluster, not of the run.
tax_cte = f"""WITH tax AS (
    SELECT sotu,
           max(percent_identity) AS best_pid,
           max(pp_cov)           AS best_pp_cov
    FROM palm_tax
    WHERE tax_family = '{family}'
    GROUP BY sotu
)"""

# 2. ps.percent_identity is the run's OWN contig-vs-palmprint identity, a
#    different number from the sOTU assignment confidence above. min_seq_id
#    gates the latter; min_run_percent_identity (if set) gates the former.
clauses = [
    f"ps.coverage >= {min_read_coverage}",
    f"tax.best_pid >= {min_seq_id}",
    f"tax.best_pp_cov >= {min_cov}",
]

# 3. geography and sample-type filters need srarun + biosample, which live in
#    the same Serratus DB -- no separate BioSample efetch required.
needs_srarun = bool(require_geo or exclude_cell_culture or geo_include
                    or sample_type_include)
joins = [
    "FROM palm_sra2 ps",
    "JOIN tax ON tax.sotu = ps.sotu",
]
if needs_srarun:
    joins.append("JOIN srarun sr ON sr.run = ps.run_id")
    joins.append("LEFT JOIN biosample bs ON bs.run_id = sr.sample")

if min_run_pid is not None:
    clauses.append(f"ps.percent_identity >= {float(min_run_pid)}")

if require_geo:
    clauses.append(
        "bs.geo_text_all->>'geo_loc_name' IS NOT NULL "
        "AND btrim(lower(bs.geo_text_all->>'geo_loc_name')) NOT IN "
        "('', 'not determined', 'missing', 'na', 'n/a', 'uncalculated', "
        "'not applicable', 'not collected', 'unknown', 'none')"
    )

# 3a: geography allowlist, against the same geo_loc_name free text require_geo
#     gates. Substrings, so "Czech" covers both spellings of the country.
if geo_include:
    _geo_re = "|".join(re.escape(g) for g in geo_include)
    clauses.append(f"bs.geo_text_all->>'geo_loc_name' ~* '({_geo_re})'")

# 3b: sample-type allowlist, matched the way analysis_geography's
#     wastewater_flavi_scope.py matches it -- srarun free text plus Serratus's
#     normalised tissue call, since library_source alone says METAGENOMIC for
#     most sewage runs. Entries are regexes, not literals.
if sample_type_include:
    _st_re = "|".join(sample_type_include)
    _cols = " OR ".join(f"{c} ~* '({_st_re})'" for c in
                        ("sr.scientific_name", "sr.library_name", "sr.sample_name"))
    clauses.append(
        f"({_cols} OR EXISTS (SELECT 1 FROM biosample_tissue bt "
        f"WHERE bt.biosample_id = bs.biosample_id AND bt.tissue ~* '({_st_re})'))"
    )

if exclude_cell_culture:
    # The free-text pass catches the easy cases only. It CANNOT catch a project
    # like PRJNA686429 (Zika serial passage in Huh7.5.1 cells): there,
    # scientific_name is "Zika virus" and library_name is "DD9_D", so nothing
    # in srarun reveals the provenance -- it lives in the BioProject title,
    # which this schema does not carry. geo_text_all was no help either: for
    # those biosamples it holds the single key 'geo_loc_name'.
    for col in ("sr.scientific_name", "sr.library_name", "sr.sample_name"):
        clauses.append(f"({col} IS NULL OR {col} !~* '{CELL_CULTURE_RE}')")

    # biosample_tissue carries Serratus's NORMALISED tissue call, derived from
    # the real BioSample isolation_source: for PRJNA686429 every sample is
    # source='isolation_source', text='huh7.5.1 cell culture',
    # tissue='cell culture'. Matching the normalised column beats guessing at
    # free text. NOT EXISTS rather than a JOIN, because biosample_tissue has
    # several rows per biosample and a join would multiply palm_sra2 again --
    # the same mistake that made the palm_tax join 144x more expensive.
    clauses.append(
        "NOT EXISTS (SELECT 1 FROM biosample_tissue bt "
        "WHERE bt.biosample_id = bs.biosample_id "
        "AND bt.tissue ~* '(cell[ _-]?culture|cell[ _-]?line|derived cell)')"
    )

# Explicit BioProject denylist: the auditable backstop for lab-derived projects
# that carry no machine-readable marker at all. Each entry needs a reason in
# the config comment, so an exclusion is never silent.
for proj in (filt.get("exclude_bioprojects") or []):
    clauses.append(f"sr.bio_project <> '{proj}'")
if filt.get("exclude_bioprojects") and not needs_srarun:
    raise ValueError("exclude_bioprojects requires require_geo or "
                     "exclude_cell_culture so that srarun is joined")

# 3. deterministic ordering, so max_runs takes the best-supported runs and a
#    re-run reproduces the same shortlist. The old test_limit was rows[:N] on
#    an UNORDERED result, which made the shortlist an arbitrary slice.
#    GROUP BY rather than DISTINCT: a run can have several assembly nodes
#    hitting the same palm_id at different coverage, and DISTINCT over a
#    coverage column would emit one row per node instead of collapsing them.
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
print(f"[INFO] filters: cov>={min_read_coverage} tax_pid>={min_seq_id} "
      f"pp_cov>={min_cov} run_pid>={min_run_pid} require_geo={require_geo} "
      f"exclude_cell_culture={exclude_cell_culture} max_runs={max_runs} "
      f"geo_include={len(geo_include)} sample_type_include={len(sample_type_include)}")
print(f"[INFO] SQL:\n{query}")

with psycopg.connect(**db_config) as conn:
    with conn.cursor() as cur:
        cur.execute(query)
        rows = cur.fetchall()
        colnames = [desc[0] for desc in cur.description] if cur.description else []

print(f"[INFO] {len(rows)} hit rows returned "
      f"({len({r[0] for r in rows})} distinct runs)")

# 4. cap on DISTINCT RUNS, not hit rows: a run may contribute several
#    (sotu, palm_id) hits and all of them must survive together, or the
#    shortlist's per-sOTU statistics downstream are computed on a truncated
#    view of the same run. Rows are already ordered by coverage desc, so this
#    keeps the best-supported runs and is reproducible.
if max_runs and rows:
    keep, seen = [], set()
    for r in rows:
        if r[0] not in seen:
            if len(seen) >= int(max_runs):
                continue
            seen.add(r[0])
        keep.append(r)
    print(f"[INFO] max_runs={max_runs}: kept {len(keep)} hit rows "
          f"across {len(seen)} runs")
    rows = keep

# 5. coverage was selected only to make the ordering deterministic; drop it so
#    the downstream schema stays run_id,sotu,palm_id.
if colnames and colnames[-1] == "coverage":
    rows = [r[:-1] for r in rows]
    colnames = colnames[:-1]

# Write outputs
os.makedirs(os.path.dirname(output.runs), exist_ok=True)
# 2. Update the write logic
# We remove 'header=False' to allow pandas to write the column names
if rows:
    # Use the column names derived from the cursor
    df_out = pd.DataFrame(rows, columns=colnames)
    df_out.to_csv(output.runs, index=False) 
else:
    # Write a header even if empty
    with open(output.runs, "w") as f:
        f.write("run_id,sotu,palm_id\n")

# =========================================================
# REF FASTA & METADATA (Logic remains identical)
# =========================================================
taxid = scen.get("anchor", {}).get("taxid") or scen.get("reference_anchor_taxid")

os.makedirs(os.path.dirname(output.ref_fasta), exist_ok=True)
if taxid:
    r = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi", 
                     params={"db": "nuccore", "term": f"txid{taxid}[Organism]", "retmax": 200})
    ids = [x.text for x in ET.fromstring(r.content).findall(".//Id")]
    if ids:
        fetch = requests.get("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi",
                             params={"db": "nuccore", "id": ",".join(ids), "rettype": "fasta", "retmode": "text"})
        with open(output.ref_fasta, "w") as f: f.write(fetch.text)
    else:
        with open(output.ref_fasta, "w") as f: pass
else:
    with open(output.ref_fasta, "w") as f: pass

os.makedirs(os.path.dirname(output.unanchored_meta), exist_ok=True)
pd.DataFrame(columns=["palm_id"]).to_csv(output.unanchored_meta, index=False, header=False)