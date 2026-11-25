#!/usr/bin/env python3
"""Item 3: per-library fraction of contigs a cheap k-mer prescreen would have
to escalate to the full three-tool taxonomy.

The denominator is the contigs a prescreen would actually see: every contig
surviving the coverage filter. Against that:

  unassigned  no tool called the contig viral -- absent from taxonomy_merged.tsv
  ambiguous   called, but the call is not corroborated. Taken from the
              merge's own family_agreement column:
                single_tool  exactly one of the three tools called it
                unresolved   called by >1 tool, but no family could be resolved
                conflict     tools resolved to DIFFERENT families
              `agree` is the only value treated as confident.
  escalate    unassigned + ambiguous

A prescreen can only safely skip the confident fraction, so escalate_fraction
is the share of contigs that still needs the expensive path.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

AMBIGUOUS_AGREEMENT = {"single_tool", "unresolved", "conflict"}
# agree_rank is one taxon at two ranks, so it is as confident as agree.
# rdrp_resolved is a disagreement the RdRp settled: counted confident, and
# auditable through the family_disagreement column.
CONFIDENT_AGREEMENT = {"agree", "agree_rank", "rdrp_resolved"}


def read_retained(path: str) -> int:
    n = 0
    try:
        with open(path) as fh:
            for rec in csv.DictReader(fh, delimiter="\t"):
                if rec.get("retained") == "1":
                    n += 1
    except FileNotFoundError:
        return -1
    return n


def read_positives(path: str) -> tuple[int, dict[str, int]]:
    buckets: dict[str, int] = {}
    total = 0
    try:
        with open(path) as fh:
            for rec in csv.DictReader(fh, delimiter="\t"):
                if not (rec.get("contig_id") or "").strip():
                    continue
                total += 1
                agr = (rec.get("family_agreement") or "").strip() or "missing"
                buckets[agr] = buckets.get(agr, 0) + 1
    except FileNotFoundError:
        return -1, {}
    return total, buckets


def kraken2_orphans(path: str) -> tuple[int, int]:
    viral = orphan = 0
    try:
        with open(path) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 3:
                    continue
                if parts[2].strip() in ("0", ""):
                    orphan += 1
                else:
                    viral += 1
    except FileNotFoundError:
        return -1, -1
    return viral, orphan


def frac(num: int, den: int) -> str:
    if den <= 0 or num < 0:
        return ""
    return f"{num / den:.6f}"


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--root", default=".")
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    with open(args.manifest) as fh:
        man = list(csv.DictReader(fh))

    cols = ["run", "cohort_family", "bioproject", "host_organism", "cohort",
            "contigs_retained", "contigs_called_viral",
            "n_unassigned", "n_ambiguous", "n_confident",
            "n_single_tool", "n_unresolved", "n_conflict", "n_agree",
            "frac_unassigned", "frac_ambiguous", "frac_confident",
            "n_escalate", "frac_escalate",
            "kraken2_orphan_contigs", "frac_kraken2_orphan", "status"]

    rows = []
    R = args.root
    for m in man:
        fam, s = m["family"], m["sample_id"]
        retained = read_retained(
            os.path.join(R, f"data/03.rnaviralSPAdes/{fam}/{s}/coverage_stats.tsv"))
        positives, buckets = read_positives(
            os.path.join(R, f"data/05.final_selection/{fam}/{s}/taxonomy_merged.tsv"))
        _, k_orphan = kraken2_orphans(
            os.path.join(R, f"data/04.taxonomy/{fam}/{s}/kraken2/kraken2_output.tsv"))

        status = "ok"
        if retained < 0 or positives < 0:
            status = "missing"

        n_amb = sum(v for k, v in buckets.items() if k in AMBIGUOUS_AGREEMENT)
        n_conf = sum(v for k, v in buckets.items() if k in CONFIDENT_AGREEMENT)
        n_other = positives - n_amb - n_conf if positives >= 0 else 0
        if n_other > 0:
            # family_agreement value we don't recognise -- count as ambiguous
            # rather than silently confident.
            n_amb += n_other
        unassigned = (retained - positives) if (retained >= 0 and positives >= 0) else -1
        if unassigned is not None and unassigned < 0 and retained >= 0 and positives >= 0:
            # positives should never exceed retained; flag rather than hide.
            status = "positives_exceed_retained"
            unassigned = 0
        escalate = (unassigned + n_amb) if unassigned >= 0 else -1

        rows.append({
            "run": s,
            "cohort_family": fam,
            "bioproject": m.get("bioproject", ""),
            "host_organism": m.get("host_species", ""),
            "cohort": m.get("condition", ""),
            "contigs_retained": "" if retained < 0 else retained,
            "contigs_called_viral": "" if positives < 0 else positives,
            "n_unassigned": "" if unassigned < 0 else unassigned,
            "n_ambiguous": n_amb if positives >= 0 else "",
            "n_confident": n_conf if positives >= 0 else "",
            "n_single_tool": buckets.get("single_tool", 0) if positives >= 0 else "",
            "n_unresolved": buckets.get("unresolved", 0) if positives >= 0 else "",
            "n_conflict": buckets.get("conflict", 0) if positives >= 0 else "",
            "n_agree": buckets.get("agree", 0) if positives >= 0 else "",
            "frac_unassigned": frac(unassigned, retained),
            "frac_ambiguous": frac(n_amb, retained) if positives >= 0 else "",
            "frac_confident": frac(n_conf, retained) if positives >= 0 else "",
            "n_escalate": "" if escalate < 0 else escalate,
            "frac_escalate": frac(escalate, retained),
            "kraken2_orphan_contigs": "" if k_orphan < 0 else k_orphan,
            "frac_kraken2_orphan": frac(k_orphan, retained),
            "status": status,
        })

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"{args.out}: {len(rows)} libraries", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
