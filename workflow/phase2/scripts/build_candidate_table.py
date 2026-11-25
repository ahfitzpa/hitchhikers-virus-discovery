#!/usr/bin/env python3
"""Build the candidate contig table: one row per retained viral contig.

Joins, per library:
  data/05.final_selection/{family}/{sample}/taxonomy_merged.tsv  tool calls + CheckV
  data/03.rnaviralSPAdes/{family}/{sample}/coverage_stats.tsv    coverage metrics
  data/04.taxonomy/{family}/{sample}/rdrpcatch/rdrp_annotated.tsv RdRp profile scores
  data/intake/manifest.csv                                        run -> BioProject/host

On the reference-similarity columns
-----------------------------------
The previous run reported an unnamed "distance" that saturated near 3.0-3.5.
That value was RdRpCATCH's `norm_bitscore_profile` -- the HMMER bitscore divided
by profile length, i.e. BITS PER PROFILE POSITION. Two things were wrong with
reporting it as a distance:

  1. direction -- it is a SIMILARITY; higher means a closer match, not further
  2. units     -- it is bits/position, with a ceiling set by the HMM's maximum
                  per-position information content (observed max 3.430 over
                  126,992 hits in the Hepeviridae run), which is why it
                  appeared to "saturate"

It is therefore carried here as `similarity_to_nearest_reference`, with its
name, units and direction in dedicated columns. No distance is derived from it:
there is no principled monotone conversion from an HMM bitscore to an
evolutionary or identity distance, so `distance_to_nearest_reference` is left
empty and `distance_metric_name` records why. Populating a real distance needs
an alignment against sequence references (e.g. DIAMOND blastp against ICTV
exemplar or RefSeq viral proteins, giving 100 - %AAI); that step is not in this
pipeline.

Likewise `nearest_reference_accession` is empty by construction: RdRpCATCH
matches profile HMMs (built from many sequences), not individual accessions, so
`nearest_reference_protein` holds the profile name and
`nearest_reference_type` records that it is an HMM profile.
"""
from __future__ import annotations

import argparse
import csv
import os
import sys

SIM_NAME = "rdrpcatch_norm_bitscore_profile"
SIM_UNITS = "bits per profile position (HMMER bitscore / profile length)"
SIM_DIRECTION = "higher_is_more_similar"
SIM_NOTE = ("similarity, not a distance; ceiling set by the profile's maximum "
            "per-position information content (observed max 3.430)")
DIST_NAME = "not_computed"
DIST_NOTE = ("no sequence-level reference alignment in this pipeline; an HMM "
             "bitscore has no principled conversion to a distance. Add a "
             "DIAMOND blastp step against ICTV/RefSeq viral proteins to report "
             "100 - %AAI in named units.")

OUT_COLS = [
    "contig_id", "run", "bioproject", "host_organism", "cohort", "cohort_family",
    "contig_length",
    "coverage_mean", "coverage_median", "coverage_sd", "coverage_cv",
    "coverage_breadth", "coverage_evenness", "mapped_reads",
    "checkv_quality", "checkv_completeness", "checkv_contamination",
    "checkv_viral_genes", "checkv_host_genes",
    "assigned_family", "family_agreement", "family_disagreement",
    "detected_by", "n_tools",
    "genomad_family", "kraken2_family", "kraken2_taxid", "rdrpcatch_family",
    "nearest_reference_type", "nearest_reference_accession",
    "nearest_reference_protein", "nearest_reference_db",
    "similarity_to_nearest_reference", "similarity_metric_name",
    "similarity_metric_units", "similarity_metric_direction",
    "similarity_metric_note",
    "rdrpcatch_evalue", "rdrpcatch_bitscore", "rdrpcatch_id_score",
    "rdrpcatch_profile_hits",
    "distance_to_nearest_reference", "distance_metric_name",
    "distance_metric_units", "distance_metric_note",
    "vs2_max_score", "virsorter2_viral",
]

METRIC_DEFS = [
    # name, definition, units, range, direction, source
    ("contig_length", "Contig length from the assembly BAM header",
     "nucleotides", ">0", "n/a", "samtools idxstats"),
    ("coverage_mean", "Sum of per-base read depth divided by contig length; "
     "uncovered positions counted as zero", "reads per base", ">=0",
     "higher_is_deeper", "samtools depth -a"),
    ("coverage_median", "Median per-base depth over all contig positions, "
     "including zero-depth positions", "reads per base", ">=0",
     "higher_is_deeper", "samtools depth -a"),
    ("coverage_sd", "Population standard deviation of per-base depth",
     "reads per base", ">=0", "higher_is_more_variable", "samtools depth -a"),
    ("coverage_cv", "Coefficient of variation of depth: coverage_sd / "
     "coverage_mean", "dimensionless", ">=0", "higher_is_more_variable",
     "derived"),
    ("coverage_breadth", "Fraction of contig positions with depth > 0",
     "fraction", "0-1", "higher_is_better", "samtools depth -a"),
    ("coverage_evenness", "Coverage evenness after Mortazavi et al. 2010: "
     "E = 1 - (|D2| - sum_{i in D2} d_i / C) / L, where C is mean depth and "
     "D2 = {i : d_i <= C}", "dimensionless", "0-1",
     "higher_is_more_uniform", "derived"),
    ("mapped_reads", "Reads mapped to the contig (bwa-mem3, unfiltered)",
     "reads", ">=0", "n/a", "samtools idxstats"),
    ("checkv_completeness", "CheckV estimated genome completeness",
     "percent", "0-100", "higher_is_more_complete", "CheckV quality_summary"),
    ("checkv_contamination", "CheckV estimated host contamination",
     "percent", "0-100", "lower_is_better", "CheckV quality_summary"),
    ("assigned_family", "Majority-vote virus family across the three taxonomy "
     "tools (consensus_family)", "taxon name", "n/a", "n/a",
     "merge_taxonomy_calls.py"),
    ("family_agreement", "agree | single_tool | unresolved | conflict",
     "category", "n/a", "n/a", "merge_taxonomy_calls.py"),
    (SIM_NAME, "HMMER bitscore divided by profile length for the best-scoring "
     "RdRp profile hit on the contig", SIM_UNITS, "0 to ~3.43",
     SIM_DIRECTION, "RdRpCATCH rdrp_annotated.tsv"),
    ("rdrpcatch_evalue", "HMMER E-value of the best-scoring RdRp profile hit",
     "expected false hits per search", ">0", "lower_is_more_significant",
     "RdRpCATCH rdrp_annotated.tsv"),
    ("distance_to_nearest_reference", DIST_NOTE, "n/a", "n/a", "n/a",
     "not computed"),
]


def read_tsv(path: str) -> list[dict]:
    try:
        with open(path) as fh:
            return [r for r in csv.DictReader(fh, delimiter="\t")]
    except FileNotFoundError:
        return []


def to_float(v: str | None):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def best_rdrp_hits(path: str) -> dict[str, dict]:
    """Best RdRp profile hit per contig, ranked by norm_bitscore_profile then
    bitscore, with E-value as the final tie-break."""
    best: dict[str, dict] = {}
    counts: dict[str, int] = {}
    rows = read_tsv(path)
    if not rows:
        return {}
    hdr = list(rows[0].keys())

    def col(*cands):
        for c in cands:
            if c in hdr:
                return c
        for c in cands:
            key = c.split("(")[0].strip().lower()
            for h in hdr:
                if key and key in h.strip().lower():
                    return h
        return None

    c_contig = col("Contig_name")
    c_prof = col("Profile_name")
    c_db = col("db_name")
    c_ev = col("E-value")
    c_nbp = col("norm_bitscore_profile")
    c_score = col("score")
    c_id = col("ID_score")
    if not c_contig:
        return {}

    for rec in rows:
        contig = (rec.get(c_contig) or "").strip()
        if not contig:
            continue
        counts[contig] = counts.get(contig, 0) + 1
        nbp = to_float(rec.get(c_nbp)) if c_nbp else None
        sc = to_float(rec.get(c_score)) if c_score else None
        ev = to_float(rec.get(c_ev)) if c_ev else None
        key = (nbp if nbp is not None else -1.0,
               sc if sc is not None else -1.0,
               -(ev if ev is not None else 1e308))
        cur = best.get(contig)
        if cur is None or key > cur["_key"]:
            best[contig] = {
                "_key": key,
                "profile": (rec.get(c_prof) or "").strip() if c_prof else "",
                "db": (rec.get(c_db) or "").strip() if c_db else "",
                "evalue": (rec.get(c_ev) or "").strip() if c_ev else "",
                "nbp": (rec.get(c_nbp) or "").strip() if c_nbp else "",
                "score": (rec.get(c_score) or "").strip() if c_score else "",
                "id_score": (rec.get(c_id) or "").strip() if c_id else "",
            }
    for contig, d in best.items():
        d["n_hits"] = counts.get(contig, 0)
    return best


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--root", default=".")
    ap.add_argument("--out", required=True)
    ap.add_argument("--out-metrics", required=True)
    args = ap.parse_args()

    with open(args.manifest) as fh:
        man = list(csv.DictReader(fh))

    R = args.root
    out_rows = []
    for m in man:
        fam, s = m["family"], m["sample_id"]
        merged = read_tsv(os.path.join(
            R, f"data/05.final_selection/{fam}/{s}/taxonomy_merged.tsv"))
        if not merged:
            continue
        cov = {r["contig_id"]: r for r in read_tsv(os.path.join(
            R, f"data/03.rnaviralSPAdes/{fam}/{s}/coverage_stats.tsv"))}
        rdrp = best_rdrp_hits(os.path.join(
            R, f"data/04.taxonomy/{fam}/{s}/rdrpcatch/rdrp_annotated.tsv"))

        for rec in merged:
            contig = (rec.get("contig_id") or "").strip()
            if not contig:
                continue
            c = cov.get(contig, {})
            h = rdrp.get(contig, {})
            out_rows.append({
                "contig_id": contig,
                "run": s,
                "bioproject": m.get("bioproject", ""),
                "host_organism": m.get("host_species", ""),
                "cohort": m.get("condition", ""),
                "cohort_family": fam,
                "contig_length": c.get("contig_length", ""),
                "coverage_mean": c.get("coverage_mean", ""),
                "coverage_median": c.get("coverage_median", ""),
                "coverage_sd": c.get("coverage_sd", ""),
                "coverage_cv": c.get("coverage_cv", ""),
                "coverage_breadth": c.get("coverage_breadth", ""),
                "coverage_evenness": c.get("coverage_evenness", ""),
                "mapped_reads": c.get("mapped_reads", ""),
                "checkv_quality": rec.get("checkv_quality", ""),
                "checkv_completeness": rec.get("checkv_completeness", ""),
                "checkv_contamination": rec.get("checkv_contamination", ""),
                "checkv_viral_genes": rec.get("checkv_viral_genes", ""),
                "checkv_host_genes": rec.get("checkv_host_genes", ""),
                "assigned_family": rec.get("consensus_family", ""),
                "family_agreement": rec.get("family_agreement", ""),
                "family_disagreement": rec.get("family_disagreement", ""),
                "detected_by": rec.get("detected_by", ""),
                "n_tools": rec.get("n_tools", ""),
                "genomad_family": rec.get("genomad_family", ""),
                "kraken2_family": rec.get("kraken2_family", ""),
                "kraken2_taxid": rec.get("kraken2_taxid", ""),
                "rdrpcatch_family": rec.get("rdrpcatch_family", ""),
                # RdRpCATCH matches profile HMMs, not single accessions.
                "nearest_reference_type": "hmm_profile" if h else "",
                "nearest_reference_accession": "",
                "nearest_reference_protein": h.get("profile", "")
                                             or rec.get("rdrpcatch_profile", ""),
                "nearest_reference_db": h.get("db", "")
                                        or rec.get("rdrpcatch_db", ""),
                "similarity_to_nearest_reference": h.get("nbp", ""),
                "similarity_metric_name": SIM_NAME if h else "",
                "similarity_metric_units": SIM_UNITS if h else "",
                "similarity_metric_direction": SIM_DIRECTION if h else "",
                "similarity_metric_note": SIM_NOTE if h else "",
                "rdrpcatch_evalue": h.get("evalue", "")
                                    or rec.get("rdrpcatch_evalue", ""),
                "rdrpcatch_bitscore": h.get("score", ""),
                "rdrpcatch_id_score": h.get("id_score", ""),
                "rdrpcatch_profile_hits": h.get("n_hits", ""),
                "distance_to_nearest_reference": "",
                "distance_metric_name": DIST_NAME,
                "distance_metric_units": "n/a",
                "distance_metric_note": DIST_NOTE,
                "vs2_max_score": rec.get("vs2_max_score", ""),
                "virsorter2_viral": rec.get("virsorter2_viral", ""),
            })

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=OUT_COLS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(out_rows)

    with open(args.out_metrics, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["metric", "definition", "units", "range", "direction", "source"])
        w.writerows(METRIC_DEFS)

    n_missing_cov = sum(1 for r in out_rows if r["coverage_mean"] == "")
    print(f"{args.out}: {len(out_rows)} retained viral contigs "
          f"({n_missing_cov} without coverage stats)", file=sys.stderr)
    print(f"{args.out_metrics}: {len(METRIC_DEFS)} metric definitions",
          file=sys.stderr)
    return 0


if __name__ == "__main__":
    sys.exit(main())
