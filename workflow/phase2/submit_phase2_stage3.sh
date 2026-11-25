#!/usr/bin/env bash
# Stage 3 only: master viral database + MetaCerberus over the samples that
# completed, plus results/taxonomy/excluded_samples.tsv naming the ones left out.
#
# Supplies stage 3 to a job whose submit script predates it (SLURM runs the copy
# of the script taken at submit time). A no-op if nothing failed. See
# ENGINEERING_NOTES.md § build_master_database.
#
#     sbatch --dependency=afterany:<main job id> submit_phase2_stage3.sh
#SBATCH --job-name=p2_stage3
#SBATCH --partition=shared
#SBATCH --cpus-per-task=4
#SBATCH --time=1-00:00:00
#SBATCH --output=logs/stage3_%j.out
#SBATCH --mail-type=END,FAIL
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$PWD/workflow/phase2}"
mkdir -p logs

# Caches go in this workflow dir on scratch, never ~ (2026-11-16: a full home
# quota killed run 1099433 at startup -- snakemake's source cache lives in
# $XDG_CACHE_HOME, default ~/.cache). Exported here so it holds even when the
# job's environment didn't come from a shell that sourced ~/.bashrc.
export XDG_CACHE_HOME="$PWD/.cache"
export CONDA_PKGS_DIRS="$XDG_CACHE_HOME/conda_pkgs"
mkdir -p "$XDG_CACHE_HOME" "$CONDA_PKGS_DIRS"
module purge
module load snakemake/8.30.0
export SNAKEMAKE_CONDA_FRONTEND=mamba
source "$(dirname "$(which conda)")/../etc/profile.d/conda.sh"
conda deactivate 2>/dev/null || true

SNAKE_OPTS="--executor slurm --slurm-no-account --jobs 40 --cores 64 --latency-wait 300 \
    --use-conda --conda-prefix .snakemake/conda --conda-frontend conda --rerun-incomplete \
    --keep-going --default-resources slurm_partition=shared runtime=180 mem_mb=0"

# --rerun-triggers mtime only, and an env var rather than --config: see
# ENGINEERING_NOTES.md § build_master_database.
echo "[stage3] $(date): aggregating over completed samples"
PHASE2_ALLOW_PARTIAL_MASTER=true snakemake -s Snakefile_phase2 $SNAKE_OPTS \
    --rerun-triggers mtime \
    -- results/master_viral_database.fasta results/taxonomy/excluded_samples.tsv \
    results/annotation/metacerberus
rc=$?
echo "[stage3] $(date): snakemake exit $rc"
if [ -s results/taxonomy/excluded_samples.tsv ]; then
    echo "[stage3] excluded samples:"
    cat results/taxonomy/excluded_samples.tsv
fi
exit $rc
