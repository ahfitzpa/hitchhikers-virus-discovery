#!/usr/bin/env python3
"""
Migrate the legacy download_manifest_parsed.csv schema into the schema
rules/01.input.smk expects. Every legacy row was SRA-fetched, so sample_id
== run_id and source == "sra" for all of them; new local/ONT rows don't
need to follow that.

Usage:
    python convert_manifest.py INPUT.csv OUTPUT.csv \
        [--search-dir DIR ...] [--dry-run]
"""
from __future__ import annotations

import argparse
import logging
import os
import sys
import tempfile
from collections.abc import Sequence
from pathlib import Path

import pandas as pd

__version__ = "1.0.0"
LOGGER = logging.getLogger(__name__)

LEGACY_COLS = [
    "run_id", "bioproject", "scientificname", "ref_genome_path",
    "scenario", "family", "host_status",
]

DEFAULT_SEARCH_DIRS = [
    "data/00.shared_sra",
    "data/01.SRA_files",
]


def _atomic_write_csv(df: pd.DataFrame, dest: str | Path) -> None:
    """Write a DataFrame as CSV to a sibling temp file and atomically rename into place."""
    dest_path = Path(dest)
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{dest_path.name}.", suffix=".tmp", dir=dest_path.parent)
    os.close(fd)
    try:
        df.to_csv(tmp_name, index=False)
        os.replace(tmp_name, dest_path)
    except BaseException:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)
        raise


def load_legacy(path: str | Path) -> pd.DataFrame:
    """Load a legacy download_manifest_parsed.csv and collapse duplicate run_id rows across scenarios."""
    df = pd.read_csv(path, dtype=str).fillna("")
    missing = [c for c in LEGACY_COLS if c not in df.columns]
    if missing:
        sys.exit(f"[convert_manifest] legacy file missing columns: {missing}")

    other_cols = [c for c in LEGACY_COLS if c not in ("run_id", "scenario")]
    inconsistent = [
        c for c in other_cols
        if df.groupby("run_id")[c].nunique().gt(1).any()
    ]
    if inconsistent:
        sys.exit(
            f"[convert_manifest] run_id duplicates disagree on {inconsistent}, "
            f"can't collapse safely, inspect the legacy file manually"
        )

    return (
        df.groupby("run_id", as_index=False)
        .agg({**{c: "first" for c in other_cols}, "scenario": lambda s: ";".join(sorted(set(s)))})
    )


def infer_layout(run_id: str, family: str, search_dirs: Sequence[str | Path]) -> str | None:
    """Infer paired vs single layout by checking existing downloaded FASTQ files on disk."""
    candidates: list[Path] = []
    for d in search_dirs:
        base = Path(d)
        candidates.extend([
            base / f"{run_id}_2.fastq.gz",
            base / family / run_id / f"{run_id}_2.fastq.gz",
            base / f"{run_id}_1.fastq.gz",
            base / family / run_id / f"{run_id}_1.fastq.gz",
        ])

    for c in candidates:
        if c.name.endswith("_2.fastq.gz") and c.exists():
            return "paired"
    for c in candidates:
        if c.name.endswith("_1.fastq.gz") and c.exists():
            return "single"
    return None


def convert(legacy_df: pd.DataFrame, search_dirs: Sequence[str | Path]) -> tuple[pd.DataFrame, list[str]]:
    """Convert legacy DataFrame rows into the Phase 2 manifest schema and list unresolved run IDs."""
    unresolved: list[str] = []
    rows: list[dict[str, str]] = []
    for _, r in legacy_df.iterrows():
        run_id = str(r["run_id"])
        family = str(r["family"])
        layout = infer_layout(run_id, family, search_dirs)
        if layout is None:
            unresolved.append(run_id)
            layout = "paired"

        rows.append({
            "sample_id": run_id,
            "family": family,
            "source": "sra",
            "run_id": run_id,
            "fq1": "",
            "fq2": "",
            "platform": "illumina",
            "layout": layout,
            "barcode": "",
            "library_id": "",
            "batch": "",
            "condition": "",
            "replicate": "",
            "host_species": "",
            "collection_date": "",
            "notes": "",
            "md5_fq1": "",
            "md5_fq2": "",
            "adapter": "",
            "genome_size": "",
            "expected_coverage": "",
            "bioproject": str(r["bioproject"]),
            "scientificname": str(r["scientificname"]),
            "ref_genome_path": str(r["ref_genome_path"]),
            "scenario": str(r["scenario"]),
            "host_status": str(r["host_status"]),
        })

    return pd.DataFrame(rows), unresolved


def main(argv: list[str] | None = None) -> int:
    """CLI entry point for legacy manifest migration."""
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_csv")
    ap.add_argument("output_csv")
    ap.add_argument(
        "--search-dir",
        action="append",
        default=None,
        help=f"directory to check for existing downloads when inferring layout (repeatable; default: {DEFAULT_SEARCH_DIRS})",
    )
    ap.add_argument("--dry-run", action="store_true", help="print the result, don't write output_csv")
    ap.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    ap.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = ap.parse_args(argv)

    logging.basicConfig(level=getattr(logging, args.log_level), format="[%(levelname)s] %(message)s")
    search_dirs = args.search_dir or DEFAULT_SEARCH_DIRS

    legacy = load_legacy(args.input_csv)
    new_df, unresolved = convert(legacy, search_dirs)

    dupes = new_df.loc[new_df["sample_id"].duplicated(), "sample_id"].unique().tolist()
    if dupes:
        sys.exit(f"[convert_manifest] duplicate sample_id after conversion: {dupes}")

    LOGGER.info("%d legacy rows -> %d manifest rows", len(legacy), len(new_df))
    LOGGER.info("families: %s", sorted(new_df["family"].unique()))
    LOGGER.info("layout counts: %s", new_df["layout"].value_counts().to_dict())
    if unresolved:
        LOGGER.warning(
            "layout could not be resolved from disk for %d run(s), defaulted to 'paired': %s",
            len(unresolved),
            unresolved,
        )

    if args.dry_run:
        print(new_df.to_string())
        return 0

    _atomic_write_csv(new_df, args.output_csv)
    LOGGER.info("wrote %s", args.output_csv)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
