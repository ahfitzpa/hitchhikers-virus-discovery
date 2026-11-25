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
from __future__ import annotations

import argparse
import csv
import logging
import os
import tempfile
from collections.abc import Iterable
from pathlib import Path
from typing import Any

__version__ = "1.0.0"
LOGGER = logging.getLogger(__name__)


def best_hit_per_query(rows: Iterable[dict[str, Any]]) -> dict[str, dict[str, Any]]:
    """Select the single best alignment per query (lowest evalue, ties broken by longest alnlen)."""
    best: dict[str, dict[str, Any]] = {}
    for row in rows:
        q = str(row["query"])
        current = best.get(q)
        if current is None:
            best[q] = row
            continue
        if row["evalue"] < current["evalue"] or (
            row["evalue"] == current["evalue"] and row["alnlen"] > current["alnlen"]
        ):
            best[q] = row
    return best


def filter_by_seq_id(
    best: dict[str, dict[str, Any]],
    min_seq_id: float,
    max_seq_id: float,
) -> dict[str, dict[str, Any]]:
    """Retain queries whose best alignment fident falls within [min_seq_id, max_seq_id]."""
    return {
        q: row
        for q, row in best.items()
        if min_seq_id <= row["fident"] <= max_seq_id
    }


def read_hits_raw(path: str | Path) -> list[dict[str, Any]]:
    """Read 6-column MMseqs2 convertalis TSV (query, target, pident, alnlen, evalue, fident)."""
    src = Path(path)
    rows: list[dict[str, Any]] = []
    with src.open(encoding="utf-8") as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t")
            if len(parts) < 6:
                raise ValueError(f"{src}:{lineno}: expected 6 tab-separated columns, found {len(parts)}")
            query, target, pident, alnlen, evalue, fident = parts[:6]
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


def write_hits(path: str | Path, best: dict[str, dict[str, Any]]) -> None:
    """Atomically write 5-column filtered hits TSV (dropping fident)."""
    dest = Path(path)
    dest.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{dest.name}.", suffix=".tmp", dir=dest.parent)
    try:
        with os.fdopen(fd, "w", newline="", encoding="utf-8") as fh:
            writer = csv.writer(fh, delimiter="\t")
            for row in best.values():
                writer.writerow(
                    [row["query"], row["target"], row["pident"], row["alnlen"], row["evalue"]]
                )
        os.replace(tmp_name, dest)
    except BaseException:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)
        raise


def main(argv: list[str] | None = None) -> int:
    """CLI entry point for best-hit sequence identity filtering."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--hits-raw", required=True)
    parser.add_argument("--min-seq-id", required=True, type=float)
    parser.add_argument("--max-seq-id", required=True, type=float)
    parser.add_argument("--out", required=True)
    parser.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
    args = parser.parse_args(argv)

    logging.basicConfig(level=getattr(logging, args.log_level), format="[%(levelname)s] %(message)s")

    if args.min_seq_id > args.max_seq_id:
        parser.error(f"--min-seq-id ({args.min_seq_id}) cannot exceed --max-seq-id ({args.max_seq_id})")

    rows = read_hits_raw(args.hits_raw)
    best = best_hit_per_query(rows)
    kept = filter_by_seq_id(best, args.min_seq_id, args.max_seq_id)
    write_hits(args.out, kept)

    dropped = len(best) - len(kept)
    LOGGER.info(
        "filter_best_hit: %d raw alignments, %d queries, %d kept in [%s, %s], %d dropped outside window",
        len(rows),
        len(best),
        len(kept),
        args.min_seq_id,
        args.max_seq_id,
        dropped,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
