#!/usr/bin/env python3
"""One-off: re-derive the family consensus in place, without re-running phase2.

Why this exists
---------------
scripts/family_taxonomy.py changed how the three taxonomy callers are
reconciled (mixed ranks and superseded names are no longer "conflict", and a
sister-family disagreement is settled by RdRpCATCH). The pipeline applies that
from rules/05.merge_tax.smk, but only on a re-run, and re-running the merge
drags CheckV and MetaCerberus back in through mtime staleness -- ~9.7 h on
Serratus -- even though the file that makes them stale, viral_positives.fasta,
is byte-identical (its content is `keep = set(merged)`, the union of contigs
any tool called viral, which does not depend on the consensus).

Nothing else has to move, because the confidence decision never reads the
family: in filter_contigs_by_annotation.py `consensus_family` is used at one
place only, to count how many kept contigs each family has. So contigs_keep.txt
is unaffected and MetaCerberus is not needed at all. This script rewrites the
three derived products directly:

  results/taxonomy/master_metadata.csv        consensus columns recomputed
  results/reporting/contig_confidence.tsv     same, joined on run_id + contig
  results/reporting/analysis_families.yaml    families with >= --min-kept-contigs
                                              kept contigs, recounted

Every file is copied to <name>.pre_retax.<timestamp> before being replaced.

This is a migration, not part of the workflow. The workflow carries the same
logic for every future run; this only spares one expensive re-run.
"""
from __future__ import annotations

import argparse
import collections
import csv
import os
import shutil
import sys
import time

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import family_taxonomy  # noqa: E402

TOOL_PRIORITY = ["genomad", "rdrpcatch", "kraken2"]
CONSENSUS_COLS = ["consensus_family", "family_agreement", "family_disagreement"]


def backup(path: str) -> str:
    dst = f"{path}.pre_retax.{time.strftime('%Y%m%dT%H%M%S')}"
    shutil.copy2(path, dst)
    return dst


def reconcile_row(row):
    fams = {t: (row.get(f"{t}_family") or "").strip() for t in TOOL_PRIORITY}
    fams = {t: f for t, f in fams.items() if f}
    return family_taxonomy.reconcile(fams, TOOL_PRIORITY)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--metadata", default="results/taxonomy/master_metadata.csv")
    ap.add_argument("--confidence", default="results/reporting/contig_confidence.tsv")
    ap.add_argument("--families-yaml", default="results/reporting/analysis_families.yaml")
    ap.add_argument("--keep", default="high,medium",
                    help="confidence levels counted as kept")
    ap.add_argument("--min-kept-contigs", type=int, default=1)
    ap.add_argument("--dry-run", action="store_true")
    a = ap.parse_args()

    table = family_taxonomy.load_lineage()
    if not table:
        print("[retax] ERROR: family_lineage.tsv not found; refusing to run "
              "without the ICTV lineage", file=sys.stderr)
        return 1
    print(f"[retax] lineage: {len(table)} taxa")

    # 1: master_metadata.csv -- the only file that carries the per-tool calls
    if not os.path.exists(a.metadata):
        print(f"[retax] ERROR: {a.metadata} not found", file=sys.stderr)
        return 1
    with open(a.metadata, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)
        fields = list(reader.fieldnames or [])
    for c in CONSENSUS_COLS:
        if c not in fields:
            fields.append(c)

    before = collections.Counter(r.get("family_agreement", "") for r in rows)
    moved = collections.Counter()
    # key the confidence table by the same identity the pipeline joins on
    by_key = {}
    for r in rows:
        old_fam = (r.get("consensus_family") or "").strip()
        fam, status, note = reconcile_row(r)
        if status != r.get("family_agreement") or fam != old_fam:
            moved[(r.get("family_agreement", ""), status, old_fam, fam)] += 1
        r["consensus_family"], r["family_agreement"], r["family_disagreement"] = fam, status, note
        # contig_confidence.tsv names a contig "<run_id>|<contig>", master
        # metadata keeps the bare contig in contig_id, so index both spellings.
        run, cid = r.get("run_id", ""), r.get("contig_id", "")
        by_key[(run, cid)] = (fam, status, note)
        by_key[(run, f"{run}|{cid}")] = (fam, status, note)
    after = collections.Counter(r["family_agreement"] for r in rows)

    print(f"[retax] {len(rows)} contigs")
    print(f"[retax] before: {dict(before)}")
    print(f"[retax] after:  {dict(after)}")
    if moved:
        print("[retax] reassignments:")
        for (os_, ns, of, nf), n in moved.most_common(20):
            arrow = f"{of} -> {nf}" if of != nf else of
            print(f"   {n:4}  {os_:10} -> {ns:14} {arrow}")

    if a.dry_run:
        print("[retax] dry run: nothing written")
        return 0

    print(f"[retax] backup: {backup(a.metadata)}")
    tmp = a.metadata + ".tmp"
    with open(tmp, "w", newline="", encoding="utf-8") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)
    os.replace(tmp, a.metadata)
    print(f"[retax] wrote {a.metadata}")

    # 2: contig_confidence.tsv -- carries the consensus columns but not the
    # per-tool calls, so the new values are joined in by run_id + contig
    if os.path.exists(a.confidence):
        with open(a.confidence, newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            crows = list(reader)
            cfields = list(reader.fieldnames or [])
        for c in CONSENSUS_COLS:
            if c not in cfields:
                cfields.append(c)
        hit = miss = 0
        for r in crows:
            got = by_key.get((r.get("run_id", ""), r.get("contig", "")))
            if got is None:
                miss += 1
                r.setdefault("family_disagreement", "")
                continue
            hit += 1
            r["consensus_family"], r["family_agreement"], r["family_disagreement"] = got
        print(f"[retax] confidence rows matched {hit}, unmatched {miss}")
        print(f"[retax] backup: {backup(a.confidence)}")
        tmp = a.confidence + ".tmp"
        with open(tmp, "w", newline="", encoding="utf-8") as fh:
            w = csv.DictWriter(fh, fieldnames=cfields, delimiter="\t",
                               extrasaction="ignore", lineterminator="\n")
            w.writeheader()
            w.writerows(crows)
        os.replace(tmp, a.confidence)
        print(f"[retax] wrote {a.confidence}")

        # 3: the analysed family list, recounted over the kept rows only
        levels = {s.strip() for s in a.keep.split(",") if s.strip()}
        kept = [r for r in crows
                if (r.get("confidence") or "").strip() in levels
                and (r.get("consensus_family") or "").strip()]
        fams = collections.Counter(r["consensus_family"].strip() for r in kept)
        chosen = sorted(f for f, n in fams.items() if n >= a.min_kept_contigs)
        if a.families_yaml:
            if os.path.exists(a.families_yaml):
                print(f"[retax] backup: {backup(a.families_yaml)}")
            with open(a.families_yaml, "w", encoding="utf-8") as fh:
                fh.write("# Generated by scripts/retax_metadata_oneoff.py -- "
                         "do not hand-edit.\n")
                fh.write(f"# Families with >= {a.min_kept_contigs} contig(s) of "
                         f"{a.keep} confidence.\n")
                # Both blocks, as filter_contigs_by_annotation.py writes them:
                # rules 09-12 read phylogenetics.families and the pathway reads
                # novel.families, and a file with only one of them would change
                # which list a pathway falls back to.
                for block in ("phylogenetics", "novel"):
                    fh.write(f"{block}:\n  families:\n")
                    for f in chosen:
                        fh.write(f"    - {f}\n")
            print(f"[retax] wrote {a.families_yaml}: {len(chosen)} families "
                  f"from {len(kept)} kept contigs")
    else:
        print(f"[retax] {a.confidence} absent; skipped confidence + families")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
