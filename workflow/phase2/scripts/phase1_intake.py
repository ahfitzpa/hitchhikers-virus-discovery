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
import argparse
import sys

import pandas as pd

REVIEW_COLS = [
    "run_id", "include", "suggested", "flags", "scientificname", "bioproject",
    "librarysource", "librarystrategy", "librarylayout", "model", "spots",
    "size_mb", "avglength", "releasedate", "biosample", "sample", "centername",
    "body_site", "disease", "n_sotus", "sotus", "tax_family",
]

# Phase2's manifest schema (rules/01.input.smk). Unset optional columns
# stay blank.
MANIFEST_COLS = [
    "sample_id", "family", "source", "run_id", "fq1", "fq2", "platform", "layout",
    "barcode", "library_id", "batch", "condition", "replicate", "host_species",
    "collection_date", "notes", "md5_fq1", "md5_fq2", "adapter", "genome_size",
    "expected_coverage", "bioproject", "scientificname", "ref_genome_path",
    "scenario", "host_status",
]


def flag_run(row):
    """Return (suggested, flags). Only hard blockers suggest N."""
    flags, suggested = [], "Y"
    if str(row.get("librarylayout", "")).upper() == "SINGLE":
        # fastp_cleaning and rnaviral_assembly are paired-end only.
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


def review(args):
    sl = pd.read_csv(args.shortlist, dtype=str).fillna("")
    sotus = sl.groupby("run_id")["sotu"].agg(lambda s: ";".join(sorted(set(s))))
    runs = sl.drop_duplicates("run_id").set_index("run_id")
    runs["sotus"] = sotus
    runs["n_sotus"] = runs["sotus"].str.count(";") + 1
    runs = runs.reset_index()
    runs[["suggested", "flags"]] = runs.apply(lambda r: pd.Series(flag_run(r)), axis=1)
    runs["include"] = ""
    runs = runs.sort_values(["suggested", "scientificname", "bioproject", "run_id"],
                            ascending=[False, True, True, True])
    runs[[c for c in REVIEW_COLS if c in runs.columns]].to_csv(args.out, index=False)
    print(f"[review] {len(runs)} runs -> {args.out} "
          f"({(runs.suggested == 'Y').sum()} suggested Y, {(runs.suggested == 'N').sum()} N)")


def build(args):
    rv = pd.read_csv(args.review, dtype=str).fillna("")
    undecided = rv[~rv["include"].str.upper().isin(["Y", "N"])]
    if len(undecided) and not args.allow_undecided:
        sys.exit(f"[build] {len(undecided)} rows have no Y/N in `include` "
                 f"(e.g. {', '.join(undecided.run_id.head(5))}); finish the review "
                 "or pass --allow-undecided to treat blanks as N.")
    keep = rv[rv["include"].str.upper() == "Y"]
    blocked = keep[keep["suggested"] == "N"]
    if len(blocked):
        sys.exit(f"[build] {len(blocked)} included runs have a phase2 blocker "
                 f"(e.g. {blocked.run_id.iloc[0]}: {blocked.flags.iloc[0]}).")

    if args.existing_manifest:
        existing = set(pd.read_csv(args.existing_manifest, dtype=str)["sample_id"])
        dup = keep[keep.run_id.isin(existing)]
        if len(dup):
            print(f"[build] skipping {len(dup)} runs already in {args.existing_manifest}")
            keep = keep[~keep.run_id.isin(existing)]

    sl = pd.read_csv(args.shortlist, dtype=str).fillna("")
    meta = sl.drop_duplicates("run_id").set_index("run_id")
    scenario = args.scenario or ""
    rows = []
    for rid in keep.run_id:
        m = meta.loc[rid]
        rows.append({
            "sample_id": rid, "family": m["tax_family"], "source": "sra", "run_id": rid,
            "platform": "illumina", "layout": "paired",
            "bioproject": m["bioproject"], "scientificname": m["scientificname"],
            "host_species": m["scientificname"], "scenario": scenario,
        })
    man = pd.DataFrame(rows, columns=MANIFEST_COLS)
    man.to_csv(args.out_manifest, index=False)
    sl[sl.run_id.isin(keep.run_id)].to_csv(args.out_shortlist, index=False)
    print(f"[build] {len(man)} runs -> {args.out_manifest}, shortlist -> {args.out_shortlist}")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
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
    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
