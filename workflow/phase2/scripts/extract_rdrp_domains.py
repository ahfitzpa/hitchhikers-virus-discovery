#!/usr/bin/env python3
"""Item 6 (RdRp half): extract RdRpCATCH's RdRp domains as protein FASTA.

RdRpCATCH reports each hit against a six-frame translation of the contigs, with
the domain bounds in amino-acid coordinates on that translated sequence:

    Translated_contig_name (frame)   NODE_22495_length_929_cov_6.826291_frame=-3
    RdRp_from(AA) / RdRp_to(AA)      186 / 309      (1-based, inclusive)

Those translated sequences live in the scan working directory
(scan/tmp/seqkit_translate_output/*.fasta), which taxonomy_cleanup otherwise
purges -- so this rule is wired to run BEFORE the purge.

The translation holds six frames for every contig (millions of records), so the
needed IDs are collected first and the FASTA is streamed once.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys
from collections import defaultdict


def find_col(header: list[str], *candidates: str) -> str | None:
    for cand in candidates:
        for col in header:
            if col.strip() == cand:
                return col
    # Fall back to a substring match (RdRpCATCH column names carry
    # parenthesised units that have moved between versions).
    for cand in candidates:
        for col in header:
            if cand.split("(")[0].strip().lower() in col.strip().lower():
                return col
    return None


def iter_fasta(path: str):
    name, chunks = None, []
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    yield name, "".join(chunks)
                # Headers carry a trailing space: ">NODE_1_..._frame=1 "
                name = line[1:].strip().split()[0] if line[1:].strip() else ""
                chunks = []
            else:
                chunks.append(line.strip())
    if name is not None:
        yield name, "".join(chunks)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rdrp-tsv", required=True)
    ap.add_argument("--translated", required=True,
                    help="seqkit_translate_output FASTA (six-frame translation)")
    ap.add_argument("--family", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out-faa", required=True)
    ap.add_argument("--out-tsv", required=True,
                    help="per-domain coordinate/score table")
    args = ap.parse_args()

    os.makedirs(os.path.dirname(args.out_faa) or ".", exist_ok=True)

    # --- read the hit table -------------------------------------------------
    wanted: dict[str, list[dict]] = defaultdict(list)
    n_rows = 0
    if os.path.exists(args.rdrp_tsv) and os.path.getsize(args.rdrp_tsv) > 0:
        with open(args.rdrp_tsv) as fh:
            rdr = csv.DictReader(fh, delimiter="\t")
            hdr = list(rdr.fieldnames or [])
            c_contig = find_col(hdr, "Contig_name")
            c_trans = find_col(hdr, "Translated_contig_name (frame)",
                               "Translated_contig_name")
            c_from = find_col(hdr, "RdRp_from(AA)", "RdRp_from")
            c_to = find_col(hdr, "RdRp_to(AA)", "RdRp_to")
            c_prof = find_col(hdr, "Profile_name")
            c_db = find_col(hdr, "db_name")
            c_ev = find_col(hdr, "E-value")
            c_nbp = find_col(hdr, "norm_bitscore_profile")
            c_score = find_col(hdr, "score")
            if not all([c_trans, c_from, c_to]):
                sys.exit(f"{args.rdrp_tsv}: missing translated-name/coordinate "
                         f"columns; found {hdr}")
            for rec in rdr:
                tname = (rec.get(c_trans) or "").strip()
                if not tname:
                    continue
                try:
                    a = int(float(rec.get(c_from) or 0))
                    b = int(float(rec.get(c_to) or 0))
                except ValueError:
                    continue
                if a <= 0 or b < a:
                    continue
                n_rows += 1
                wanted[tname].append({
                    "contig_id": (rec.get(c_contig) or "").strip() if c_contig else "",
                    "translated_name": tname,
                    "aa_from": a, "aa_to": b,
                    "profile": (rec.get(c_prof) or "").strip() if c_prof else "",
                    "db": (rec.get(c_db) or "").strip() if c_db else "",
                    "evalue": (rec.get(c_ev) or "").strip() if c_ev else "",
                    "norm_bitscore_profile": (rec.get(c_nbp) or "").strip() if c_nbp else "",
                    "bitscore": (rec.get(c_score) or "").strip() if c_score else "",
                })

    # --- stream the translation, slicing the wanted ranges -------------------
    out_cols = ["family", "sample", "contig_id", "translated_name", "frame",
                "aa_from", "aa_to", "domain_aa_length", "profile", "db",
                "evalue", "norm_bitscore_profile", "bitscore", "extracted"]
    written = 0
    rows_out = []
    seen: set[str] = set()

    faa = open(args.out_faa, "w")
    try:
        if wanted and os.path.exists(args.translated):
            for name, seq in iter_fasta(args.translated):
                hits = wanted.get(name)
                if not hits:
                    continue
                seen.add(name)
                frame = name.rsplit("_frame=", 1)[-1] if "_frame=" in name else ""
                for h in hits:
                    a, b = h["aa_from"], h["aa_to"]
                    sub = seq[a - 1:b]
                    ok = bool(sub)
                    if ok:
                        header = (f"{args.sample}|{h['contig_id']}|frame={frame}"
                                  f"|aa={a}-{b}|profile={h['profile']}"
                                  f"|db={h['db']}|evalue={h['evalue']}"
                                  f"|norm_bitscore_profile={h['norm_bitscore_profile']}")
                        faa.write(f">{header}\n")
                        for i in range(0, len(sub), 60):
                            faa.write(sub[i:i + 60] + "\n")
                        written += 1
                    rows_out.append({
                        "family": args.family, "sample": args.sample,
                        "contig_id": h["contig_id"],
                        "translated_name": name, "frame": frame,
                        "aa_from": a, "aa_to": b,
                        "domain_aa_length": len(sub),
                        "profile": h["profile"], "db": h["db"],
                        "evalue": h["evalue"],
                        "norm_bitscore_profile": h["norm_bitscore_profile"],
                        "bitscore": h["bitscore"],
                        "extracted": "1" if ok else "0",
                    })
    finally:
        faa.close()

    # Hits whose translated sequence never appeared get a row with extracted=0,
    # so a silent mismatch between the two files is visible rather than lost.
    for tname, hits in wanted.items():
        if tname in seen:
            continue
        frame = tname.rsplit("_frame=", 1)[-1] if "_frame=" in tname else ""
        for h in hits:
            rows_out.append({
                "family": args.family, "sample": args.sample,
                "contig_id": h["contig_id"], "translated_name": tname,
                "frame": frame, "aa_from": h["aa_from"], "aa_to": h["aa_to"],
                "domain_aa_length": 0, "profile": h["profile"], "db": h["db"],
                "evalue": h["evalue"],
                "norm_bitscore_profile": h["norm_bitscore_profile"],
                "bitscore": h["bitscore"], "extracted": "0",
            })

    with open(args.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=out_cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows_out)

    missed = len(rows_out) - written
    print(f"{args.sample}: {n_rows} RdRp hits, {written} domains written, "
          f"{missed} not extracted -> {args.out_faa}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
