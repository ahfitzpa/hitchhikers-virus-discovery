#!/usr/bin/env python3
"""Item 2: every Snakemake benchmark file combined into one TSV.

Snakemake writes one benchmark per (rule, wildcard set) at
benchmarks/{rule}/{stem}.tsv. The stem is NOT uniform across rules:

    benchmarks/rnaviral_assembly/Hylobius_abietis_SRR6765939.tsv   {family}_{sample}
    benchmarks/fetch_and_dump_sra_paired/SRR6765939.tsv            {sra_id} only
    benchmarks/setup_metacerberus/setup.tsv                        no wildcards

Cohort labels contain underscores (Hylobius_abietis), so the stem cannot be
split naively -- it is resolved against the manifest's known families and
sample IDs instead.

Wall time and peak RSS are re-exposed as explicitly unit-named columns
(wall_time_s, peak_rss_mib) alongside Snakemake's own headers, which are
undocumented in the file itself.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

# Snakemake's benchmark header -> (our name, unit). Snakemake reports memory in
# MiB and io in MB; see snakemake/benchmark.py.
PASSTHROUGH_UNITS = {
    "s": "seconds", "h:m:s": "hh:mm:ss",
    "max_rss": "MiB", "max_vms": "MiB", "max_uss": "MiB", "max_pss": "MiB",
    "io_in": "MB", "io_out": "MB",
    "mean_load": "ratio_cpu_cores", "cpu_time": "seconds",
}


def resolve_stem(stem: str, families: set[str], samples: set[str]) -> tuple[str, str]:
    """(family, sample) from a benchmark filename stem."""
    if stem in samples:
        return "", stem
    for fam in families:
        prefix = fam + "_"
        if stem.startswith(prefix):
            rest = stem[len(prefix):]
            if rest in samples:
                return fam, rest
    # Global rule (setup.tsv, run.tsv, select.tsv) or an unknown stem.
    return "", ""


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--benchmark-dir", default="benchmarks")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with open(args.manifest) as fh:
        man = list(csv.DictReader(fh))
    families = {m["family"] for m in man}
    samples = {m["sample_id"] for m in man}
    sample_meta = {m["sample_id"]: m for m in man}

    out_cols = ["rule", "cohort_family", "run", "bioproject", "host_organism",
                "wall_time_s", "peak_rss_mib", "cpu_time_s", "mean_load",
                "io_in_mb", "io_out_mb", "max_vms_mib", "max_uss_mib",
                "max_pss_mib", "wall_clock_hms", "benchmark_file"]

    rows = []
    if not os.path.isdir(args.benchmark_dir):
        sys.exit(f"no benchmark dir: {args.benchmark_dir}")

    for rule_name in sorted(os.listdir(args.benchmark_dir)):
        rule_dir = os.path.join(args.benchmark_dir, rule_name)
        if not os.path.isdir(rule_dir):
            continue
        for fname in sorted(os.listdir(rule_dir)):
            if not fname.endswith(".tsv"):
                continue
            path = os.path.join(rule_dir, fname)
            stem = fname[:-4]
            fam, samp = resolve_stem(stem, families, samples)
            meta = sample_meta.get(samp, {})
            # A rule with retries appends one row per attempt; keep them all.
            with open(path) as fh:
                for rec in csv.DictReader(fh, delimiter="\t"):
                    if not rec or not (rec.get("s") or "").strip():
                        continue
                    rows.append({
                        "rule": rule_name,
                        "cohort_family": fam or meta.get("family", ""),
                        "run": samp,
                        "bioproject": meta.get("bioproject", ""),
                        "host_organism": meta.get("host_species", ""),
                        "wall_time_s": rec.get("s", ""),
                        "peak_rss_mib": rec.get("max_rss", ""),
                        "cpu_time_s": rec.get("cpu_time", ""),
                        "mean_load": rec.get("mean_load", ""),
                        "io_in_mb": rec.get("io_in", ""),
                        "io_out_mb": rec.get("io_out", ""),
                        "max_vms_mib": rec.get("max_vms", ""),
                        "max_uss_mib": rec.get("max_uss", ""),
                        "max_pss_mib": rec.get("max_pss", ""),
                        "wall_clock_hms": rec.get("h:m:s", ""),
                        "benchmark_file": path,
                    })

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    unresolved = sum(1 for r in rows if not r["run"])
    print(f"{args.out}: {len(rows)} benchmark rows across "
          f"{len({r['rule'] for r in rows})} rules "
          f"({unresolved} not tied to a library)", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
