#!/usr/bin/env python3
"""Best-hit max_seq_id/min_seq_id filter for Phase1 anchored modes.

See ENGINEERING_NOTES.md / docs/SPEEDUP_PLAN.md step 15 for why this exists:
mmseqs convertalis can report several partial alignments per query/target
pair (one toy case: a 39-nt fragment at 100% identity, scored far below the
real best alignment by evalue), so "best hit" has to mean the top-scoring
alignment, not the maximum fident seen anywhere for that query. max_seq_id
is enforced only against that one best-scoring alignment per query.

Input: convertalis output with columns query,target,pident,alnlen,evalue,fident
(fident is a fraction, pident is a percentage -- only fident is used for the
min_seq_id/max_seq_id comparison, since those config values are fractions).

Output: same 5 columns as before this step existed (query,target,pident,
alnlen,evalue) -- fident is dropped after use so build_matrix/extract_sra_ids
and expand_cluster_hits need no schema changes.
"""
import argparse
import csv
import sys


def best_hit_per_query(rows):
    """rows: iterable of dicts with query, evalue (float), alnlen (int).
    Returns dict query -> best row (lowest evalue, ties broken by longest alnlen)."""
    best = {}
    for row in rows:
        q = row["query"]
        current = best.get(q)
        if current is None:
            best[q] = row
            continue
        if row["evalue"] < current["evalue"] or (
            row["evalue"] == current["evalue"] and row["alnlen"] > current["alnlen"]
        ):
            best[q] = row
    return best


def filter_by_seq_id(best, min_seq_id, max_seq_id):
    return {
        q: row
        for q, row in best.items()
        if min_seq_id <= row["fident"] <= max_seq_id
    }


def read_hits_raw(path):
    rows = []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            query, target, pident, alnlen, evalue, fident = line.split("\t")
            rows.append(
                {
                    "query": query,
                    "target": target,
                    "pident": pident,
                    "alnlen": int(alnlen),
                    "evalue": float(evalue),
                    "fident": float(fident),
                }
            )
    return rows


def write_hits(path, best):
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        for row in best.values():
            writer.writerow(
                [row["query"], row["target"], row["pident"], row["alnlen"], row["evalue"]]
            )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hits-raw", required=True)
    parser.add_argument("--min-seq-id", required=True, type=float)
    parser.add_argument("--max-seq-id", required=True, type=float)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    rows = read_hits_raw(args.hits_raw)
    best = best_hit_per_query(rows)
    kept = filter_by_seq_id(best, args.min_seq_id, args.max_seq_id)
    write_hits(args.out, kept)

    dropped = len(best) - len(kept)
    print(
        f"filter_best_hit: {len(rows)} raw alignments, {len(best)} queries, "
        f"{len(kept)} kept in [{args.min_seq_id}, {args.max_seq_id}], "
        f"{dropped} dropped outside window",
        file=sys.stderr,
    )


if __name__ == "__main__":
    main()
