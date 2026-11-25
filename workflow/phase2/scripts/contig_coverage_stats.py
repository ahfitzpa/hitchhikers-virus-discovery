#!/usr/bin/env python3
"""Per-contig coverage statistics from the assembly BAM.

The BAM is bwa-mem3-derived (bowtie2 was dropped from the assembly env over a
zlib conflict with hostile's pin -- see ENGINEERING_NOTES.md § rnaviral_assembly).

Every metric emitted here is defined with units in
`results/reporting/metrics_definitions.tsv`; nothing is reported unnamed.

Reads `samtools depth -a` on stdin (one pass, streaming -- depth histograms per
contig rather than full depth vectors, so memory stays bounded) plus
`samtools idxstats` for contig lengths and mapped-read counts.
"""
from __future__ import annotations

import argparse
import sys
from collections import Counter, defaultdict


def evenness(hist: Counter[int], length: int, mean_depth: float) -> float:
    """Coverage evenness E, after Mortazavi et al. 2010 (Nat Methods 7:e1000).

        D2 = { i : d_i <= C }           C = mean depth
        E  = 1 - ( |D2| - (sum_{i in D2} d_i) / C ) / L

    Dimensionless, 0-1. E = 1 is perfectly uniform depth; E -> 0 as coverage
    concentrates into a shrinking fraction of the contig. Returns 0.0 for a
    contig with no coverage at all.
    """
    if length <= 0 or mean_depth <= 0:
        return 0.0
    n_le = 0
    sum_le = 0
    for depth, count in hist.items():
        if depth <= mean_depth:
            n_le += count
            sum_le += depth * count
    return 1.0 - (n_le - (sum_le / mean_depth)) / length


def median_from_hist(hist: Counter[int], length: int) -> float:
    """Median depth from a depth histogram, including zero-depth positions."""
    if length <= 0:
        return 0.0
    # Positions absent from `samtools depth -a` output are depth 0.
    counted = sum(hist.values())
    if counted < length:
        hist = Counter(hist)
        hist[0] += length - counted
    target_lo = (length - 1) // 2
    target_hi = length // 2
    seen = 0
    lo = hi = None
    for depth in sorted(hist):
        seen += hist[depth]
        if lo is None and seen > target_lo:
            lo = depth
        if seen > target_hi:
            hi = depth
            break
    if lo is None:
        lo = 0
    if hi is None:
        hi = lo
    return (lo + hi) / 2.0


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--idxstats", required=True,
                    help="samtools idxstats output (contig lengths, mapped reads)")
    ap.add_argument("--depth", default="-",
                    help="samtools depth -a output, or '-' for stdin")
    ap.add_argument("--retained", default=None,
                    help="optional: filtered_contigs.fasta -- adds a 'retained' flag")
    ap.add_argument("--family", required=True)
    ap.add_argument("--sample", required=True)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    lengths: dict[str, int] = {}
    mapped: dict[str, int] = {}
    with open(args.idxstats) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4 or parts[0] == "*":
                continue
            contig, ln, m = parts[0], int(parts[1]), int(parts[2])
            lengths[contig] = ln
            mapped[contig] = m

    retained: set[str] | None = None
    if args.retained:
        retained = set()
        with open(args.retained) as fh:
            for line in fh:
                if line.startswith(">"):
                    retained.add(line[1:].split()[0])

    # Streaming accumulation: histogram + running sums per contig.
    hists: dict[str, Counter[int]] = defaultdict(Counter)
    sum_d: dict[str, int] = defaultdict(int)
    sum_d2: dict[str, int] = defaultdict(int)
    covered: dict[str, int] = defaultdict(int)

    depth_fh = sys.stdin if args.depth == "-" else open(args.depth)
    try:
        for line in depth_fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            contig = parts[0]
            d = int(parts[2])
            hists[contig][d] += 1
            sum_d[contig] += d
            sum_d2[contig] += d * d
            if d > 0:
                covered[contig] += 1
    finally:
        if depth_fh is not sys.stdin:
            depth_fh.close()

    cols = ["family", "sample", "contig_id", "contig_length", "mapped_reads",
            "coverage_mean", "coverage_median", "coverage_sd", "coverage_cv",
            "coverage_breadth", "coverage_evenness", "retained"]

    with open(args.out, "w") as out:
        out.write("\t".join(cols) + "\n")
        for contig in sorted(lengths):
            L = lengths[contig]
            total = sum_d.get(contig, 0)
            mean = total / L if L else 0.0
            # Variance over all L positions (uncovered positions contribute 0).
            if L > 0:
                var = max(sum_d2.get(contig, 0) / L - mean * mean, 0.0)
            else:
                var = 0.0
            sd = var ** 0.5
            cv = (sd / mean) if mean > 0 else 0.0
            breadth = (covered.get(contig, 0) / L) if L else 0.0
            med = median_from_hist(hists.get(contig, Counter()), L)
            even = evenness(hists.get(contig, Counter()), L, mean)
            flag = "" if retained is None else ("1" if contig in retained else "0")
            out.write("\t".join([
                args.family, args.sample, contig, str(L), str(mapped.get(contig, 0)),
                f"{mean:.4f}", f"{med:.4f}", f"{sd:.4f}", f"{cv:.4f}",
                f"{breadth:.4f}", f"{even:.4f}", flag,
            ]) + "\n")
    return 0


if __name__ == "__main__":
    sys.exit(main())
