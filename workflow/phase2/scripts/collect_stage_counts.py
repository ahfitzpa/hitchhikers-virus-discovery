#!/usr/bin/env python3
"""Item 1: per-library, per-stage counts as one tidy TSV.

Emits long/tidy form (one row per library x stage) as the primary output, plus
a wide pivot for reading. Every stage name is defined with its unit in
results/reporting/metrics_definitions.tsv.

Read counts come from fastp's own JSON report (before/after filtering), so no
extra pass over the FASTQs is needed. Contig counts come from the BAM idxstats
and coverage_stats written by contig_coverage_stats, and from Kraken2's
per-contig output -- the only per-contig source that covers EVERY contig rather
than only the tool-called positives.
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import sys

# (stage_name, unit) in funnel order.
STAGES = [
    ("raw_reads", "reads"),
    ("post_fastp_reads", "reads"),
    ("contigs_assembled", "contigs"),
    ("contigs_after_coverage_filter", "contigs"),
    ("contigs_ge_min_len", "contigs"),
    ("contigs_viral_union", "contigs"),
    ("contigs_kraken2_viral", "contigs"),
    ("contigs_kraken2_orphan", "contigs"),
]


def count_fasta(path: str) -> int:
    n = 0
    try:
        with open(path) as fh:
            for line in fh:
                if line.startswith(">"):
                    n += 1
    except FileNotFoundError:
        return -1
    return n


def fastp_counts(path: str) -> tuple[int, int]:
    try:
        with open(path) as fh:
            d = json.load(fh)
    except (FileNotFoundError, json.JSONDecodeError):
        return -1, -1
    summ = d.get("summary", {})
    return (int(summ.get("before_filtering", {}).get("total_reads", -1)),
            int(summ.get("after_filtering", {}).get("total_reads", -1)))


def idxstats_contigs(path: str) -> int:
    n = 0
    try:
        with open(path) as fh:
            for line in fh:
                parts = line.split("\t")
                if parts and parts[0] != "*" and parts[0].strip():
                    n += 1
    except FileNotFoundError:
        return -1
    return n


def coverage_counts(path: str, min_len: int) -> tuple[int, int]:
    """(retained contigs, retained contigs >= min_len)."""
    retained = 0
    long_enough = 0
    try:
        with open(path) as fh:
            for rec in csv.DictReader(fh, delimiter="\t"):
                if rec.get("retained") == "1":
                    retained += 1
                    try:
                        if int(rec["contig_length"]) >= min_len:
                            long_enough += 1
                    except (KeyError, ValueError) as exc:
                        print(f"[collect_stage_counts] WARNING: skipping malformed row in {path}: {exc}", file=sys.stderr)
    except FileNotFoundError:
        return -1, -1
    return retained, long_enough


def kraken2_counts(path: str) -> tuple[int, int]:
    """(classified, unclassified). The Kraken2 DB here is viral-only, so any
    classified contig sits under TaxID 10239; unclassified carries TaxID 0."""
    viral = orphan = 0
    try:
        with open(path) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                taxid = parts[2].strip()
                if taxid in ("0", ""):
                    orphan += 1
                else:
                    viral += 1
    except FileNotFoundError:
        return -1, -1
    return viral, orphan


def count_tsv_rows(path: str) -> int:
    try:
        with open(path) as fh:
            first = fh.readline()
            if not first:
                return 0
            return sum(1 for line in fh if line.strip())
    except FileNotFoundError:
        return -1
    return -1


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--root", default=".")
    ap.add_argument("--min-len", type=int, default=1000)
    ap.add_argument("--out-long", required=True)
    ap.add_argument("--out-wide", required=True)
    args = ap.parse_args()

    with open(args.manifest) as fh:
        rows = list(csv.DictReader(fh))

    R = args.root
    long_rows = []
    wide_rows = []

    for m in rows:
        fam, s = m["family"], m["sample_id"]
        raw, post = fastp_counts(
            os.path.join(R, f"data/02.fastp/{fam}/{s}/{s}_report.json"))
        assembled = idxstats_contigs(
            os.path.join(R, f"data/03.rnaviralSPAdes/{fam}/{s}/idxstats.tsv"))
        retained, long_enough = coverage_counts(
            os.path.join(R, f"data/03.rnaviralSPAdes/{fam}/{s}/coverage_stats.tsv"),
            args.min_len)
        union = count_tsv_rows(
            os.path.join(R, f"data/05.final_selection/{fam}/{s}/taxonomy_merged.tsv"))
        k_viral, k_orphan = kraken2_counts(
            os.path.join(R, f"data/04.taxonomy/{fam}/{s}/kraken2/kraken2_output.tsv"))

        values = {
            "raw_reads": raw,
            "post_fastp_reads": post,
            "contigs_assembled": assembled,
            "contigs_after_coverage_filter": retained,
            "contigs_ge_min_len": long_enough,
            "contigs_viral_union": union,
            "contigs_kraken2_viral": k_viral,
            "contigs_kraken2_orphan": k_orphan,
        }

        meta = {
            "run": s,
            "cohort_family": fam,
            "bioproject": m.get("bioproject", ""),
            "host_organism": m.get("host_species", ""),
            "cohort": m.get("condition", ""),
        }

        for order, (stage, unit) in enumerate(STAGES, start=1):
            v = values[stage]
            long_rows.append({
                **meta,
                "stage_order": order,
                "stage": stage,
                "count": "" if v < 0 else v,
                "unit": unit,
                "status": "missing" if v < 0 else "ok",
            })
        wide_rows.append({**meta,
                          **{k: ("" if v < 0 else v) for k, v in values.items()}})

    os.makedirs(os.path.dirname(args.out_long) or ".", exist_ok=True)
    long_cols = ["run", "cohort_family", "bioproject", "host_organism", "cohort",
                 "stage_order", "stage", "count", "unit", "status"]
    with open(args.out_long, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=long_cols, delimiter="\t")
        w.writeheader()
        w.writerows(long_rows)

    wide_cols = ["run", "cohort_family", "bioproject", "host_organism", "cohort"] \
        + [s for s, _ in STAGES]
    with open(args.out_wide, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=wide_cols, delimiter="\t")
        w.writeheader()
        w.writerows(wide_rows)

    n_missing = sum(1 for r in long_rows if r["status"] == "missing")
    print(f"{args.out_long}: {len(long_rows)} rows ({len(rows)} libraries x "
          f"{len(STAGES)} stages), {n_missing} missing", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
