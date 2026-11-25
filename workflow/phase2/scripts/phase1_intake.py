#!/usr/bin/env python3
"""
Bring a phase1 scenario's candidate runs into phase2, with a manual review
step in between.

  review   phase1 shortlist_summary.csv -> one-row-per-run review CSV with
           an empty `include` column, a `suggested` value and `flags`
           explaining anything phase2 can't handle or that deserves a look.
  build    reviewed CSV (rows with include = Y) -> phase2 manifest rows +
           the matching shortlist rows (build_master_database joins host/
           bioproject metadata from the shortlist). Writes new files; it
           never edits the live manifest, so it's safe while phase2 runs.

Usage:
    python scripts/phase1_intake.py review \\
        --shortlist ../phase1/data/03.merged_metadata/geographic/shortlist_summary.csv \\
        --scenario geographic --out data/intake/flaviviridae_geographic_review.csv
    # edit the `include` column (Y / N), then:
    python scripts/phase1_intake.py build \\
        --review data/intake/flaviviridae_geographic_review.csv \\
        --shortlist ../phase1/data/03.merged_metadata/geographic/shortlist_summary.csv \\
        --existing-manifest data/test_run/manifest.csv \\
        --out-manifest data/intake/flaviviridae_geographic_manifest.csv \\
        --out-shortlist data/intake/flaviviridae_geographic_shortlist.csv
"""
from __future__ import annotations

import argparse
import logging
import os
import sys
import tempfile
from collections.abc import Mapping
from pathlib import Path
from typing import Any

import pandas as pd

__version__ = "1.0.0"
LOGGER = logging.getLogger(__name__)

REVIEW_COLS = [
    "run_id", "include", "suggested", "flags", "scientificname", "bioproject",
    "librarysource", "librarystrategy", "librarylayout", "model", "spots",
    "size_mb", "avglength", "releasedate", "biosample", "sample", "centername",
    "body_site", "disease", "n_sotus", "sotus", "tax_family",
]

MANIFEST_COLS = [
    "sample_id", "family", "source", "run_id", "fq1", "fq2", "platform", "layout",
    "barcode", "library_id", "batch", "condition", "replicate", "host_species",
    "collection_date", "notes", "md5_fq1", "md5_fq2", "adapter", "genome_size",
    "expected_coverage", "bioproject", "scientificname", "ref_genome_path",
    "scenario", "host_status",
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


def flag_run(row: Mapping[str, Any]) -> tuple[str, str]:
    """Return (suggested, flags) for a candidate run; only hard Phase 2 blockers suggest 'N'."""
    flags: list[str] = []
    suggested = "Y"
    if str(row.get("librarylayout", "")).upper() == "SINGLE":
        flags.append("single-end: phase2 fastp/assembly are paired-only")
        suggested = "N"
    if str(row.get("platform", "")).upper() not in ("ILLUMINA", ""):
        flags.append(f"platform {row.get('platform')}: phase2 is short-read only")
        suggested = "N"
    name = str(row.get("scientificname", ""))
    if "virus" in name.lower():
        flags.append("organism is the virus itself: likely isolate/culture sequencing")
    if str(row.get("librarystrategy", "")).upper() in ("AMPLICON", "TARGETED-CAPTURE"):
        flags.append("amplicon/targeted library: partial genome coverage expected")
    if "metagenome" in name.lower():
        flags.append("host unresolved (metagenome)")
    return suggested, "; ".join(flags)


def review(args: argparse.Namespace) -> None:
    """Collapse a Phase 1 shortlist to one row per run with automated blocker flags."""
    sl = pd.read_csv(args.shortlist, dtype=str).fillna("")
    sotus = sl.groupby("run_id")["sotu"].agg(lambda s: ";".join(sorted(set(s))))
    runs = sl.drop_duplicates("run_id").set_index("run_id")
    runs["sotus"] = sotus
    runs["n_sotus"] = runs["sotus"].str.count(";") + 1
    runs = runs.reset_index()
    runs[["suggested", "flags"]] = runs.apply(lambda r: pd.Series(flag_run(r)), axis=1)
    runs["include"] = ""
    runs = runs.sort_values(
        ["suggested", "scientificname", "bioproject", "run_id"],
        ascending=[False, True, True, True],
    )
    _atomic_write_csv(runs[[c for c in REVIEW_COLS if c in runs.columns]], args.out)
    LOGGER.info(
        "review: %d runs -> %s (%d suggested Y, %d N)",
        len(runs),
        args.out,
        int((runs.suggested == "Y").sum()),
        int((runs.suggested == "N").sum()),
    )


def build(args: argparse.Namespace) -> None:
    """Convert reviewed include=Y rows into Phase 2 manifest and filtered shortlist CSVs."""
    rv = pd.read_csv(args.review, dtype=str).fillna("")
    undecided = rv[~rv["include"].str.upper().isin(["Y", "N"])]
    if len(undecided) and not args.allow_undecided:
        sys.exit(
            f"[build] {len(undecided)} rows have no Y/N in `include` "
            f"(e.g. {', '.join(undecided.run_id.head(5))}); finish the review "
            "or pass --allow-undecided to treat blanks as N."
        )
    keep = rv[rv["include"].str.upper() == "Y"]
    blocked = keep[keep["suggested"] == "N"]
    if len(blocked):
        sys.exit(
            f"[build] {len(blocked)} included runs have a phase2 blocker "
            f"(e.g. {blocked.run_id.iloc[0]}: {blocked.flags.iloc[0]})."
        )

    if args.existing_manifest:
        existing = set(pd.read_csv(args.existing_manifest, dtype=str)["sample_id"])
        dup = keep[keep.run_id.isin(existing)]
        if len(dup):
            LOGGER.info("build: skipping %d runs already in %s", len(dup), args.existing_manifest)
            keep = keep[~keep.run_id.isin(existing)]

    sl = pd.read_csv(args.shortlist, dtype=str).fillna("")
    meta = sl.drop_duplicates("run_id").set_index("run_id")
    scenario = args.scenario or ""
    rows: list[dict[str, str]] = []
    for rid in keep.run_id:
        m = meta.loc[rid]
        rows.append({
            "sample_id": str(rid),
            "family": str(m["tax_family"]),
            "source": "sra",
            "run_id": str(rid),
            "platform": "illumina",
            "layout": "paired",
            "bioproject": str(m["bioproject"]),
            "scientificname": str(m["scientificname"]),
            "host_species": str(m["scientificname"]),
            "scenario": scenario,
        })
    man = pd.DataFrame(rows, columns=MANIFEST_COLS)
    _atomic_write_csv(man, args.out_manifest)
    _atomic_write_csv(sl[sl.run_id.isin(keep.run_id)], args.out_shortlist)
    LOGGER.info("build: %d runs -> %s, shortlist -> %s", len(man), args.out_manifest, args.out_shortlist)


def main(argv: list[str] | None = None) -> int:
    """CLI entry point for Phase 1 intake review and manifest build."""
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    sub = p.add_subparsers(dest="cmd", required=True)

    r = sub.add_parser("review")
    r.add_argument("--shortlist", required=True)
    r.add_argument("--scenario", default="")
    r.add_argument("--out", required=True)
    r.set_defaults(func=review)

    b = sub.add_parser("build")
    b.add_argument("--review", required=True)
    b.add_argument("--shortlist", required=True)
    b.add_argument("--scenario", default="")
    b.add_argument("--existing-manifest")
    b.add_argument("--allow-undecided", action="store_true")
    b.add_argument("--out-manifest", required=True)
    b.add_argument("--out-shortlist", required=True)
    b.set_defaults(func=build)

    args = p.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level), format="[%(levelname)s] %(message)s")
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
