#!/bin/bash
#SBATCH --job-name=snakemake_p2
#SBATCH --output=logs/snakemake_p2_%j.log
#SBATCH --error=logs/snakemake_p2_%j.err
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=2
#SBATCH --time=14-00:00:00
#SBATCH --mail-type=BEGIN,END,FAIL,REQUEUE

# This job is just the Snakemake scheduler -- it submits every rule as its own
# SLURM job via --executor slurm, so the actual compute happens in those per-
# rule jobs, not here.

cd "$SLURM_SUBMIT_DIR"
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

# sra_fetch=4: caps concurrent SRA downloads regardless of the --jobs ceiling
# above, which is shared across every rule type, not fetch-specific.
# This starting value hasn't been through the sacct
# wall-clock-vs-N benchmark the review calls for yet.
SNAKE_OPTS="--executor slurm --slurm-no-account --jobs 40 --cores 64 --latency-wait 300 --use-conda --conda-prefix .snakemake/conda --conda-frontend conda --rerun-incomplete --rerun-triggers mtime params input --keep-going --default-resources slurm_partition=shared runtime=180 mem_mb=0 --resources sra_fetch=4"

echo "[SLURM INFO] Starting Phase 2: Assembly & Analysis (fanning out across the shared partition via --executor slurm)"

# Stage 1: finish VirSorter2 before anything consumes its score -- the score is
# a params path in annotate, so nothing else enforces that ordering. Allowed to
# fail.
# --notemp: contig_coverage_stats is outside this partial run, so without it
# Snakemake reclaims each sorted.bam once filter_contigs is done and stage 2
# has to reassemble to get it back (what happened to conifer, 2026-11). Stage 2
# runs without it and reclaims the temp files once every consumer is done.
# rules/99.retention.smk refuses an --until run without --notemp.
echo "[SLURM INFO] Stage 1/3: VirSorter2 (best-effort; failures are logged, not fatal)"
snakemake -s Snakefile_phase2 $SNAKE_OPTS --notemp --until run_virsorter2 || true

echo "[SLURM INFO] Stage 2/3: full workflow"
snakemake -s Snakefile_phase2 $SNAKE_OPTS || true

# Stage 3: aggregate whatever stage 2 produced, so one unassemblable sample
# cannot block MetaCerberus. Separate invocation (a DAG is built once) and an
# env var, not --config (this Snakefile rebuilds config from the YAML files).
echo "[SLURM INFO] Stage 3/3: master database + MetaCerberus over the samples that completed"
PHASE2_ALLOW_PARTIAL_MASTER=true snakemake -s Snakefile_phase2 $SNAKE_OPTS \
    --rerun-triggers mtime \
    -- results/master_viral_database.fasta results/taxonomy/excluded_samples.tsv \
    results/annotation/metacerberus

# Reporting (kept local -- trivial, no need to submit these)
snakemake -s Snakefile_phase2 --dag --jobs 1 | dot -Tsvg > logs/phase2_dag.svg
snakemake -s Snakefile_phase2 --report report_phase2.html --jobs 1
