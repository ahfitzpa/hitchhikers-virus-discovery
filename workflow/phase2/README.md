# Phase 2 — Viral Assembly, Three-Tool Taxonomy, and Functional Annotation

## Problem Solved
Phase 2 converts a CSV manifest of SRA run accessions and/or local FASTQ files into quality-controlled viral assemblies, consensus taxonomic assignments (`Kraken2`, `geNomad`, `RdRpCATCH`), quality/completeness summaries (`CheckV`), functional domain annotations (`MetaCerberus`), and ICTV VMR-mapped candidate tables. It runs either downstream of Phase 1 or completely standalone.

## Directory Map

| Path | Role |
| :--- | :--- |
| `Snakefile_phase2` | Workflow entry point (supports Illumina, ONT, and PacBio manifests) |
| `submit_phase2.sh` | Three-stage SLURM driver (per-sample assembly/taxonomy → reporting → master DB) |
| `submit_phase2_stage3.sh` | Stage 3 standalone driver (master viral DB + MetaCerberus) |
| `setup_databases.sh` / `SETUP.md` | Idempotent external database provisioning (`Kraken2`, `geNomad`, `RdRpCATCH`, `CheckV`, `VirSorter2`, `MetaCerberus`) |
| `configs/paths.yaml` | Manifest path and external database locations |
| `rules/00.setup.smk` – `09.novelty.smk` | Modular Snakemake rules in execution order |
| `rules/99.retention.smk` | Retention guard preventing accidental deletion of downstream-required files |
| `scripts/` | Manifest intake, read staging, coverage stats, three-tool taxonomy merge, ICTV prep, reporting |
| `envs/` | Per-rule Conda environment definitions |

## Primary Outputs
- `results/<sample>/taxonomy_merged.tsv` — Unified per-contig taxonomy and annotation table
- `results/master_viral_database.fasta` — Deduplicated viral positive contigs across the manifest
- `results/reporting/` — Stage attrition counts, escalation fractions, RdRp domains, and candidate tables
