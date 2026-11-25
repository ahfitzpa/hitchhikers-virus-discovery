#!/usr/bin/env python3
"""Item 7 SCAFFOLD: prepare per-virus-family bundles for an ICTV phylogenetic
comparison.

This script deliberately implements NO phylogenetic method. It establishes the
interface the method will plug into, and nothing else:

  INPUT   results/reporting/candidate_contigs.tsv   retained viral contigs
          results/master_viral_database.fasta       their nucleotide sequences
          data/06.rdrp/**/**_rdrp_domains.faa       their RdRp protein domains

  OUTPUT  results/ictv/families.tsv                 families found + counts
          results/ictv/<family>/contigs.tsv         member contigs + metadata
          results/ictv/<family>/contigs.fasta       member nucleotide sequences
          results/ictv/<family>/rdrp_domains.faa    member RdRp domains
          results/ictv/<family>/REFERENCES.todo     what the method must supply
          results/ictv/<family>/.prepared           completion sentinel

The family grouping key is `assigned_family` from the candidate table (the
majority vote across Kraken2 / geNomad / RdRpCATCH). Contigs whose family could
not be resolved are grouped under `_unresolved` so they are visible rather than
dropped -- a later method may choose to place them by phylogeny instead.

TO BE SPECIFIED (later prompt): which ICTV release and reference set, the
alignment and tree method, the substitution model, how an ICTV
genus/family demarcation threshold is applied, and what distance the comparison
should report (in named units -- see the metric note in
scripts/build_candidate_table.py).
"""
from __future__ import annotations

import argparse
import csv
import os
import re
import sys
from collections import defaultdict

UNRESOLVED = "_unresolved"
SAFE = re.compile(r"[^A-Za-z0-9._-]+")


def safe_name(name: str) -> str:
    s = SAFE.sub("_", (name or "").strip())
    return s or UNRESOLVED


def iter_fasta(path: str):
    name, chunks = None, []
    try:
        fh = open(path)
    except FileNotFoundError:
        return
    with fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                hdr = line[1:].strip()
                name = hdr.split()[0] if hdr else ""
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def write_fasta(fh, name: str, seq: str) -> None:
    fh.write(f">{name}\n")
    for i in range(0, len(seq), 60):
        fh.write(seq[i:i + 60] + "\n")


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--candidates", required=True)
    ap.add_argument("--master-fasta", required=True)
    ap.add_argument("--rdrp-faa", nargs="*", default=[])
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--families-tsv", required=True)
    ap.add_argument("--min-contigs", type=int, default=1,
                    help="families with fewer members still get a bundle, but "
                         "are flagged in families.tsv")
    args = ap.parse_args()

    with open(args.candidates) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))

    by_family: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        fam = (r.get("assigned_family") or "").strip()
        if not fam or fam.lower() in ("na", "nan", "none", "unknown", "unclassified"):
            fam = UNRESOLVED
        by_family[fam].append(r)

    # The master FASTA holds every retained contig across all libraries. Its
    # headers are rewritten by collect_viral_positives as
    # sample|family_context|scenario|tax|tools|host|contig_id, so index on the
    # trailing contig id as well as the whole header.
    wanted_ids: set[str] = {(r.get("contig_id") or "").strip() for r in rows}
    seqs: dict[str, str] = {}
    for name, seq in iter_fasta(args.master_fasta):
        if name in wanted_ids:
            seqs.setdefault(name, seq)
        tail = name.rsplit("|", 1)[-1]
        if tail in wanted_ids:
            seqs.setdefault(tail, seq)

    # RdRp domain headers start with "<sample>|<contig_id>|frame=..."
    rdrp_by_contig: dict[str, list[tuple[str, str]]] = defaultdict(list)
    for path in args.rdrp_faa:
        for name, seq in iter_fasta(path):
            parts = name.split("|")
            if len(parts) >= 2:
                rdrp_by_contig[parts[1]].append((name, seq))

    os.makedirs(args.out_dir, exist_ok=True)
    fam_rows = []

    for fam, members in sorted(by_family.items()):
        fdir = os.path.join(args.out_dir, safe_name(fam))
        os.makedirs(fdir, exist_ok=True)

        cols = ["contig_id", "run", "bioproject", "host_organism", "cohort",
                "contig_length", "coverage_mean", "coverage_evenness",
                "checkv_completeness", "assigned_family", "family_agreement",
                "nearest_reference_protein", "similarity_to_nearest_reference"]
        with open(os.path.join(fdir, "contigs.tsv"), "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t",
                               extrasaction="ignore")
            w.writeheader()
            w.writerows(members)

        n_seq = 0
        with open(os.path.join(fdir, "contigs.fasta"), "w") as fh:
            for r in members:
                cid = (r.get("contig_id") or "").strip()
                seq = seqs.get(cid)
                if seq:
                    write_fasta(fh, f"{r.get('run','')}|{cid}", seq)
                    n_seq += 1

        n_dom = 0
        with open(os.path.join(fdir, "rdrp_domains.faa"), "w") as fh:
            for r in members:
                cid = (r.get("contig_id") or "").strip()
                for name, seq in rdrp_by_contig.get(cid, []):
                    write_fasta(fh, name, seq)
                    n_dom += 1

        with open(os.path.join(fdir, "REFERENCES.todo"), "w") as fh:
            fh.write(
                f"# ICTV comparison inputs still to be specified for: {fam}\n"
                f"#\n"
                f"# members                 {len(members)} contigs\n"
                f"# nucleotide sequences    {n_seq} written to contigs.fasta\n"
                f"# RdRp protein domains    {n_dom} written to rdrp_domains.faa\n"
                f"#\n"
                f"# Needed from the later prompt:\n"
                f"#   1. ICTV release (e.g. VMR MSL40) and the reference set to\n"
                f"#      compare against -- family exemplars, or all members\n"
                f"#   2. marker: whole genome, RdRp domain, or a named protein\n"
                f"#   3. alignment method and trimming\n"
                f"#   4. tree method and substitution model\n"
                f"#   5. the demarcation criterion to apply, and\n"
                f"#   6. the distance to report, WITH UNITS -- see\n"
                f"#      scripts/build_candidate_table.py on why the previous\n"
                f"#      run's unnamed 3.0-3.5 value was uninterpretable\n"
            )

        open(os.path.join(fdir, ".prepared"), "w").close()

        fam_rows.append({
            "assigned_family": fam,
            "bundle_dir": fdir,
            "n_contigs": len(members),
            "n_sequences_written": n_seq,
            "n_rdrp_domains": n_dom,
            "n_runs": len({r.get("run", "") for r in members}),
            "n_bioprojects": len({r.get("bioproject", "") for r in members}),
            "below_min_contigs": "1" if len(members) < args.min_contigs else "0",
            "method_status": "awaiting_specification",
        })

    os.makedirs(os.path.dirname(args.families_tsv) or ".", exist_ok=True)
    with open(args.families_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, delimiter="\t", fieldnames=[
            "assigned_family", "bundle_dir", "n_contigs", "n_sequences_written",
            "n_rdrp_domains", "n_runs", "n_bioprojects", "below_min_contigs",
            "method_status"])
        w.writeheader()
        w.writerows(sorted(fam_rows, key=lambda r: -r["n_contigs"]))

    print(f"{args.families_tsv}: {len(fam_rows)} families prepared "
          f"({sum(r['n_contigs'] for r in fam_rows)} contigs). "
          f"Method not implemented -- scaffold only.", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
