#!/bin/bash
#SBATCH --job-name=serratus_confidence
#SBATCH --output=logs/contig_confidence_%j.log
#SBATCH --error=logs/contig_confidence_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=04:00:00
#SBATCH --mail-type=END,FAIL
#
# 1: score contig confidence from MetaCerberus annotation + taxonomy + RdRp.
# --rerun-triggers mtime is mandatory, not stylistic: Snakemake 8's default
# triggers include `code`, so an edited rule file marks every rule in it stale
# and the cascade reaches back through reclaimed temp() files to re-downloading
# reads (2368 jobs). With mtime it is the one job it should be.
set -uo pipefail
cd "${SLURM_SUBMIT_DIR:-$(dirname "$0")}"
mkdir -p logs
module purge
module load snakemake/8.30.0
export SNAKEMAKE_CONDA_FRONTEND=conda
export XDG_CACHE_HOME=$PWD/.cache
export CONDA_PKGS_DIRS=$XDG_CACHE_HOME/conda_pkgs
source "$(dirname "$(which conda)")/../etc/profile.d/conda.sh"
conda deactivate 2>/dev/null || true
export PHASE2_MANIFEST_OVERRIDE=${PHASE2_MANIFEST_OVERRIDE:-data/test_run/manifest_novel.csv}
snakemake -s Snakefile_phase2 --cores 2 --use-conda \
    --conda-prefix .snakemake/conda --conda-frontend conda \
    --rerun-triggers mtime --rerun-incomplete \
    -- results/reporting/contig_confidence.tsv
rc=$?
echo "exit $rc"
[ -f results/reporting/contig_confidence_summary.tsv ] && cat results/reporting/contig_confidence_summary.tsv
exit $rc
