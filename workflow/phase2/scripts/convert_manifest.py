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
import argparse
import os
import sys

import pandas as pd

LEGACY_COLS = [
    "run_id", "bioproject", "scientificname", "ref_genome_path",
    "scenario", "family", "host_status",
]

# Default places to look for already-downloaded reads when inferring
# paired vs. single layout (see infer_layout()).
DEFAULT_SEARCH_DIRS = [
    "data/00.shared_sra",
    "data/01.SRA_files",
]


def load_legacy(path):
    df = pd.read_csv(path, dtype=str).fillna("")
    missing = [c for c in LEGACY_COLS if c not in df.columns]
    if missing:
        sys.exit(f"[convert_manifest] legacy file missing columns: {missing}")

    # The same run_id can appear under different discovery scenarios
    # (e.g. both a "surrogate" match and a "host_expansion" candidate).
    # One physical run is one sample, so collapse to one row per run_id
    # and keep all distinct scenario values.
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

    df = (
        df.groupby("run_id", as_index=False)
        .agg({**{c: "first" for c in other_cols}, "scenario": lambda s: ";".join(sorted(set(s)))})
    )
    return df


def infer_layout(run_id, family, search_dirs):
    """Paired if a _2 mate exists on disk, else single. Checked against
    real files (most legacy rows are already downloaded) rather than
    re-derived from SRA metadata, to avoid a network round trip per row."""
    candidates = []
    for d in search_dirs:
        candidates.append(os.path.join(d, f"{run_id}_2.fastq.gz"))
        candidates.append(os.path.join(d, family, run_id, f"{run_id}_2.fastq.gz"))
        candidates.append(os.path.join(d, f"{run_id}_1.fastq.gz"))
        candidates.append(os.path.join(d, family, run_id, f"{run_id}_1.fastq.gz"))

    for c in candidates:
        if c.endswith("_2.fastq.gz") and os.path.exists(c):
            return "paired"
    for c in candidates:
        if c.endswith("_1.fastq.gz") and os.path.exists(c):
            return "single"
    return None  # unresolved, caller decides how to handle


def convert(legacy_df, search_dirs):
    unresolved = []
    rows = []
    for _, r in legacy_df.iterrows():
        layout = infer_layout(r["run_id"], r["family"], search_dirs)
        if layout is None:
            unresolved.append(r["run_id"])
            layout = "paired"  # default for this study; flagged below regardless

        rows.append({
            "sample_id": r["run_id"],
            "family": r["family"],
            "source": "sra",
            "run_id": r["run_id"],
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
            "bioproject": r["bioproject"],
            "scientificname": r["scientificname"],
            "ref_genome_path": r["ref_genome_path"],
            "scenario": r["scenario"],
            "host_status": r["host_status"],
        })

    return pd.DataFrame(rows), unresolved


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input_csv")
    ap.add_argument("output_csv")
    ap.add_argument("--search-dir", action="append", default=None,
                     help=f"directory to check for existing downloads when inferring layout "
                          f"(repeatable; default: {DEFAULT_SEARCH_DIRS})")
    ap.add_argument("--dry-run", action="store_true", help="print the result, don't write output_csv")
    args = ap.parse_args()

    search_dirs = args.search_dir or DEFAULT_SEARCH_DIRS

    legacy = load_legacy(args.input_csv)
    new_df, unresolved = convert(legacy, search_dirs)

    dupes = new_df.loc[new_df["sample_id"].duplicated(), "sample_id"].unique().tolist()
    if dupes:
        sys.exit(f"[convert_manifest] duplicate sample_id after conversion: {dupes}")

    print(f"[convert_manifest] {len(legacy)} legacy rows -> {len(new_df)} manifest rows")
    print(f"[convert_manifest] families: {sorted(new_df['family'].unique())}")
    print(f"[convert_manifest] layout counts: {new_df['layout'].value_counts().to_dict()}")
    if unresolved:
        print(f"[convert_manifest] WARNING: layout could not be resolved from disk for "
              f"{len(unresolved)} run(s), defaulted to 'paired', verify manually: {unresolved}")

    if args.dry_run:
        print(new_df.to_string())
        return

    new_df.to_csv(args.output_csv, index=False)
    print(f"[convert_manifest] wrote {args.output_csv}")


if __name__ == "__main__":
    main()
