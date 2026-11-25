# Phase 1 — Serratus/Logan SRA Discovery Scenarios

## Problem Solved
Public sequencing archives contain millions of transcriptomes, making exhaustive *de novo* assembly infeasible. Phase 1 screens SRA runs via Serratus/PalmDB and Logan unitigs under four configurable discovery scenarios (`novel`, `surrogate`, `host_expansion`, and `geographic`) to produce a compact, metadata-annotated shortlist of candidate SRA runs for Phase 2.

## Directory Map

| Path | Role |
| :--- | :--- |
| `Snakefile_phase1` | Workflow entry point loading configs and rules |
| `submit_phase1.sh` | SLURM driver (`--executor slurm`) |
| `configs/phase1_sra_discovery.yaml` | Active `current_mode` and per-scenario family/anchor/threshold settings |
| `configs/paths.yaml` | Input/output directory layout |
| `rules/01.base_discovery.smk` | Initial Serratus/PalmDB sOTU centroid query |
| `rules/01.centroid_anchored.smk` | MMseqs2 alignment against anchor references and best-hit filtering |
| `rules/01.centroid_unanchored.smk` | Unanchored sOTU expansion to SRA accessions |
| `rules/03.metadata_fetch.smk` | NCBI SRA XML/TSV metadata retrieval |
| `rules/04.metadata_merge.smk` | Metadata integration and automated review gate (`approved.flag`) |
| `scripts/` | Python and Bash helpers invoked by the rules above |
| `envs/env_phase1_discovery.yaml` | Conda environment definition |

## Primary Outputs
- `data/01.sOTUs_centroid/<mode>/` — sOTU centroid hits and MMseqs2 alignments
- `data/02.SRA_metadata/<mode>/` — Raw NCBI SRA metadata records
- `data/03.merged_metadata/<mode>/` — Annotated candidate matrices and `approved.flag` / `sra_ids.txt`
