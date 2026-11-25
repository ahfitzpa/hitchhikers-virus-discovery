# Changelog

All notable changes to **A Hitchhiker's Guide to Virus Discovery** are documented in this file.
The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).

## [1.0.0] - 2026-11-25

### Added
- Public release of **Phase 1** (Serratus/Logan SRA discovery scenarios) and **Phase 2** (standalone viral assembly, three-tool taxonomy, and functional annotation).
- Machine-readable run provenance writer (`workflow/phase2/scripts/write_run_manifest.py`) recording workflow version, git commit SHA, command, and SHA-256 checksums of inputs.
- Module-level READMEs (`workflow/phase1/README.md`, `workflow/phase2/README.md`), `CITATION.cff`, `CONTRIBUTING.md`, `pyproject.toml`, and `.pre-commit-config.yaml`.
- Automated unit, syntax, hygiene, and synthetic-fixture integration test suite (`tests/`) with GitHub Actions CI (`.github/workflows/ci.yml`).

### Changed
- Parameterized all PostgreSQL queries in `01.centroid_query.py` and `align_matrix.py` using `psycopg` v3 placeholders.
- Added atomic sibling temp-file + `os.replace()` writes across Phase 1 and core Phase 2 output writers.
- Added explicit 30 s timeout and 3-attempt exponential backoff retries for NCBI Entrez HTTP calls.
- Added type hints, `pathlib.Path`, structured `logging`, and `--log-level` / `--version` flags to core CLI scripts (`filter_best_hit.py`, `merge_taxonomy_calls.py`, `convert_manifest.py`, `phase1_intake.py`).

### Fixed
- Removed misleading error-path empty CSV write in `align_matrix.py` so failed database queries fail cleanly in Snakemake.
- Fixed colliding Conda env name (`tax_qc`) and `ruamel.yaml` package name in `workflow/phase2/envs/tax_qc.yaml`.
- Fixed unexpanded `$PWD` inside `#SBATCH --output` in `workflow/phase2/submit_phase2_stage3.sh`.
- Restored `workflow/phase2/scripts/filter_contigs_by_annotation.py` referenced by `rules/07.reporting.smk`.
