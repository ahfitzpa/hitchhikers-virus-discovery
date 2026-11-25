#!/bin/bash
# One-time provisioning of phase2's shared reference databases for a new
# checkout/fork. See SETUP.md for prerequisites, sizing, and the SLURM/
# internet-access caveats, read that before running this unattended.
#
# Provisions six databases: Kraken2 viral, geNomad, RdRpCATCH, CheckV,
# VirSorter2 and the MetaCerberus HMMs. No MMseqs2 database is needed by any
# part of this workflow -- phase2's taxonomy is Kraken2 + geNomad + RdRpCATCH,
# and phase1 runs mmseqs only against DBs it builds from its own sOTU
# centroids, so nothing downloads an MMseqs2 nt/nr database.
#
# Usage: ./setup_databases.sh /path/to/database/root [threads]
set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 /path/to/database/root [threads]" >&2
    echo "See SETUP.md for what this provisions and how long it takes." >&2
    exit 1
fi

DB_ROOT="$1"
THREADS="${2:-8}"
cd "$(dirname "$0")"

if ! command -v snakemake >/dev/null 2>&1; then
    # Best-effort for Environment-Modules clusters; harmless no-op
    # elsewhere (e.g. a plain conda install where snakemake is already
    # on PATH, or a system without `module` at all).
    module load snakemake/8.30.0 2>/dev/null || true
fi
if ! command -v snakemake >/dev/null 2>&1; then
    echo "ERROR: snakemake not found on PATH (and no 'module load snakemake' available)." >&2
    echo "Install/activate it first, see SETUP.md prerequisites." >&2
    exit 1
fi

mkdir -p "$DB_ROOT"

# Point paths.yaml's databases: block at $DB_ROOT. Edits only the known
# `key: "value"` lines in place (not a full YAML rewrite) so every comment
# and everything else in the file is left untouched.
echo "[setup_databases] Pointing configs/paths.yaml databases: block at $DB_ROOT"
sed -i \
    -e "s#^\(\s*kraken2_db:\s*\).*#\1\"$DB_ROOT/kraken2_viral_db\"#" \
    -e "s#^\(\s*genomad_db:\s*\).*#\1\"$DB_ROOT/genomad_db\"#" \
    -e "s#^\(\s*rdrp_catch_hmm:\s*\).*#\1\"$DB_ROOT/rdrpcatch_db\"#" \
    -e "s#^\(\s*virsorter2_db:\s*\).*#\1\"$DB_ROOT/virsorter2_db\"#" \
    -e "s#^\(\s*checkv_db:\s*\).*#\1\"$DB_ROOT/checkv_db\"#" \
    -e "s#^\(\s*metacerberus:\s*\).*#\1\"$DB_ROOT/metacerberus\"#" \
    configs/paths.yaml
grep -A8 "^databases:" configs/paths.yaml

SNAKE_OPTS="--use-conda --conda-prefix .snakemake/conda --conda-frontend conda --rerun-incomplete"

run_target() {
    local rule="$1" label="$2"
    echo
    echo "[setup_databases] $label"
    snakemake -s Snakefile_phase2 --until "$rule" --cores "$THREADS" $SNAKE_OPTS
}

# Cheapest first, so if something's going to fail on a missing tool or
# bad path, you find out in minutes rather than after the big downloads.
run_target setup_kraken2       "Fetching Kraken2 viral DB (~2GB)..."
run_target setup_rdrpcatch     "Fetching RdRpCATCH DB (~4GB)..."
run_target setup_checkv        "Fetching CheckV DB (~6GB)..."
run_target setup_virsorter2    "Fetching VirSorter2 DB (~12GB)..."
run_target setup_genomad       "Fetching geNomad DB (~1.5GB)..."
run_target setup_metacerberus  "Fetching MetaCerberus HMMs (KOFam/VOG/PHROG)..."

echo
echo "[setup_databases] Fetching ICTV Virus Metadata Resource (VMR MSL41.v1)..."
mkdir -p references
curl -fsSL -A "Mozilla/5.0" \
    "https://ictv.global/sites/default/files/VMR/VMR_MSL41.v1.20260729.xlsx" \
    -o "references/VMR_MSL41.v1.20260729.xlsx"

echo
echo "[setup_databases] Done. All 6 databases and ICTV VMR are in place under $DB_ROOT and references/."
