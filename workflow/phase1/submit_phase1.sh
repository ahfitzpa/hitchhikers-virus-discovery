#!/bin/bash
#SBATCH --job-name=snakemake_master
#SBATCH --output=logs/snakemake_%j.log
#SBATCH --error=logs/snakemake_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=14-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE
# Pass your own notification address at submit time (SBATCH directives
# can't read a config file): sbatch --mail-user=you@example.org submit_phase1.sh

# This job is just the Snakemake scheduler -- it submits every rule as its own
# SLURM job via --executor slurm, so the actual compute happens in those per-
# rule jobs, not here.

set -uo pipefail

module purge
module load snakemake/8.30.0
export SNAKEMAKE_CONDA_FRONTEND=conda

source "$(dirname "$(which conda)")/../etc/profile.d/conda.sh"
conda deactivate 2>/dev/null || true

cd "$SLURM_SUBMIT_DIR"
mkdir -p logs

# Scratch-backed tmp, not node-local /tmp (small, shared, fills up under
# concurrent jobs). Exported as TMPDIR and passed to Snakemake's default
# tmpdir resource so every rule picks it up.
mkdir -p tmp
export TMPDIR="$(pwd)/tmp"

SNAKE_OPTS="--executor slurm --slurm-no-account --jobs 40 --cores 64 --latency-wait 120 --use-conda --conda-prefix .snakemake/conda --rerun-incomplete --conda-frontend conda --keep-going --rerun-triggers mtime code input params --default-resources tmpdir='$TMPDIR' slurm_partition=shared runtime=180 mem_mb=0"

echo "[SLURM INFO] Starting Phase 1: SRA Discovery (fanning out across the shared partition via --executor slurm)..."

snakemake -s Snakefile_phase1 $SNAKE_OPTS

# Reporting (kept local -- trivial, no need to submit these)
snakemake -s Snakefile_phase1 --dag --jobs 1 | dot -Tsvg > logs/phase1_dag.svg
snakemake -s Snakefile_phase1 --report report_phase1.html --jobs 1

echo "[SLURM INFO] Phase 1 completed."
