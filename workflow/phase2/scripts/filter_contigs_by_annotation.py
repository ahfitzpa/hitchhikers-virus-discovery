#!/usr/bin/env python3
"""Score each candidate contig's confidence from annotation plus taxonomy.

Taxonomy alone keeps a lot of eukaryotic genomic material: a 14.7 kb contig in
this dataset is called viral and encodes titin, and others encode collagen,
dynein, laminin, ubiquitin and Down syndrome cell adhesion molecule. Those are
insect genome, not virus.

MetaCerberus annotation settles it, but NOT by presence: it reports a best hit
from every database for every ORF, so all 5991 annotated contigs here have
PHROG and VOG hits. The signal is in WHICH database scores better.

  cellular databases  KOFam (KEGG, FOAM) -- cellular and universal functions
  viral databases     PHROG, VOG -- virus orthologous groups only

Per ORF the best cellular score is compared with the best viral score; a side
wins if it beats the other by `--score-ratio`. Measured on conifer: 59% of ORFs
are clearly cellular, 30% clearly viral, 11% ambiguous.

Confidence per contig:

  high    an RdRp domain was extracted -- a replicase is definitive, whatever
          else is on the contig -- or the ORFs are predominantly viral with no
          eukaryotic marker and at least `--min-tools` taxonomy tools agreeing
  medium  viral evidence is real but mixed: a viral majority resting on one
          tool, or viral ORFs alongside cellular ones. Chimeric assemblies and
          integrated elements land here, so these are worth inspecting.
  low     cellular ORFs dominate, or a named eukaryotic marker is present and
          no RdRp. This is the host genome.

Use `--keep high,medium` to carry the defensible set forward; `low` contigs stay
in the table with their reason, so nothing disappears silently.
"""
from __future__ import annotations

import argparse
import csv
import glob
import os
import re
import sys
from collections import Counter, defaultdict

CELLULAR_FILES = ("annotation_summary_KOFam_all_KEGG.tsv",
                  "annotation_summary_KOFam_all_FOAM.tsv")
VIRAL_FILES = ("annotation_summary_PHROG.tsv", "annotation_summary_VOG.tsv")

# Unambiguously eukaryotic or organellar products. A viral genome does not
# encode these, so one good hit is enough to distrust the contig.
EUK_MARKERS = (
    r"titin", r"\bactin\b", r"tubulin", r"myosin", r"collagen", r"histone",
    r"ribosomal protein", r"spliceosome", r"kinesin", r"dynein", r"laminin",
    r"cadherin", r"ubiquitin", r"cytochrome c oxidase", r"NADH dehydrogenase",
    r"ATP synthase", r"cell adhesion molecule", r"TATA-binding",
    r"SWI/SNF", r"zinc finger", r"ankyrin repeat", r"WD40", r"dehydrogenase",
    r"proteasome", r"chaperonin", r"elongation factor", r"aminoacyl",
    r"glycosyltransferase", r"cytoskelet",
)

ORF_SUFFIX = re.compile(r"_\d+_\d+_[+-]$")


def log(*a):
    print(*a, file=sys.stderr, flush=True)


def clean(v) -> str:
    v = str(v or "").strip()
    return "" if v.lower() in {"", "nan", "na", "none", "-"} else v


def as_float(v, default=0.0) -> float:
    try:
        return float(clean(v))
    except ValueError:
        return default


def contig_of(orf: str) -> str:
    # 1: MetaCerberus names an ORF <run>|<contig>_<start>_<end>_<strand>
    return ORF_SUFFIX.sub("", orf)


def read_annotations(annotation_dir: str, files: tuple[str, ...]):
    """orf -> (best score, product) over the given databases."""
    best: dict[str, tuple[float, str]] = {}
    for fn in files:
        path = os.path.join(annotation_dir, fn)
        if not os.path.exists(path):
            log(f"[filter] missing annotation file: {fn}")
            continue
        n = 0
        with open(path, newline="") as fh:
            for r in csv.DictReader(fh, delimiter="\t"):
                orf = clean(r.get("target"))
                if not orf:
                    continue
                score = as_float(r.get("score"))
                product = clean(r.get("product"))
                if score > best.get(orf, (0.0, ""))[0]:
                    best[orf] = (score, product)
                n += 1
        log(f"[filter] {fn}: {n} annotation(s)")
    return best


def load_rdrp(paths: list[str]) -> set[str]:
    """Contigs with an extracted RdRp domain."""
    contigs: set[str] = set()
    for pattern in paths:
        for path in glob.glob(pattern):
            if not os.path.isfile(path):
                continue
            if path.endswith((".faa", ".fasta")):
                # headers are pipe-delimited: <run>|<contig>|frame=N|aa=...
                with open(path) as fh:
                    for line in fh:
                        if not line.startswith(">"):
                            continue
                        parts = line[1:].split()[0].split("|")
                        if len(parts) >= 2:
                            contigs.add(parts[1])
                            contigs.add(f"{parts[0]}|{parts[1]}")
                        elif parts:
                            contigs.add(contig_of(parts[0]))
            else:
                with open(path, newline="") as fh:
                    for r in csv.DictReader(fh, delimiter="\t"):
                        # the extractor writes contig_id plus sample
                        contig = ""
                        for key in ("contig_id", "contig_name", "Contig_name",
                                    "contig", "target"):
                            if clean(r.get(key)):
                                contig = clean(r[key]).split("|")[-1]
                                break
                        if not contig:
                            continue
                        contig = re.sub(r"\|?frame=.*$", "", contig)
                        contigs.add(contig)
                        run = clean(r.get("sample") or r.get("run_id"))
                        if run:
                            contigs.add(f"{run}|{contig}")
    return contigs


def load_taxonomy(path: str) -> dict[str, dict]:
    """contig key -> taxonomy evidence, keyed both bare and run-qualified."""
    out: dict[str, dict] = {}
    if not path or not os.path.exists(path):
        return out
    with open(path, newline="") as fh:
        sample = fh.readline()
        fh.seek(0)
        delim = "\t" if "\t" in sample else ","
        for r in csv.DictReader(fh, delimiter=delim):
            contig = clean(r.get("contig_id") or r.get("contig_name"))
            run = clean(r.get("run_id") or r.get("sample"))
            if not contig:
                continue
            rec = {
                "run_id": run,
                "consensus_family": clean(r.get("consensus_family")),
                "family_agreement": clean(r.get("family_agreement")),
                "n_tools": as_float(r.get("n_tools")),
                "checkv_completeness": as_float(r.get("checkv_completeness"), -1),
                "checkv_quality": clean(r.get("checkv_quality")),
            }
            out[contig] = rec
            if run:
                out[f"{run}|{contig}"] = rec
    return out


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--annotation-dir", required=True,
                   help="MetaCerberus final/annotations-* directory")
    p.add_argument("--taxonomy", help="master_metadata.csv or candidate_contigs.tsv")
    p.add_argument("--rdrp", action="append", default=[],
                   help="glob for rdrp_domains files (repeatable)")
    p.add_argument("--out-tsv", required=True)
    p.add_argument("--out-keep", help="contig ids of the kept set")
    p.add_argument("--out-summary")
    p.add_argument("--out-families-yaml",
                   help="YAML fragment naming the families that survive in the "
                        "kept set, for the analysis pathway's config")
    p.add_argument("--min-kept-contigs", type=int, default=2,
                   help="a family needs this many kept contigs to be analysed")
    p.add_argument("--score-ratio", type=float, default=1.5,
                   help="a side must beat the other by this factor to win")
    p.add_argument("--min-viral-frac", type=float, default=0.5,
                   help="fraction of decided ORFs that must be viral for high")
    p.add_argument("--min-tools", type=float, default=2,
                   help="taxonomy tools that must agree for high without RdRp")
    p.add_argument("--keep", default="high,medium")
    args = p.parse_args()

    cellular = read_annotations(args.annotation_dir, CELLULAR_FILES)
    viral = read_annotations(args.annotation_dir, VIRAL_FILES)
    rdrp = load_rdrp(args.rdrp)
    taxonomy = load_taxonomy(args.taxonomy)
    log(f"[filter] RdRp contigs: {len(rdrp)}; taxonomy rows: {len(taxonomy)}")

    euk_re = re.compile("|".join(EUK_MARKERS), re.I)
    keep_levels = {k.strip() for k in args.keep.split(",") if k.strip()}

    # 2: decide each ORF, then aggregate to its contig
    per_contig = defaultdict(lambda: {
        "n_orfs": 0, "viral": 0, "cellular": 0, "ambiguous": 0,
        "euk_hits": 0, "euk_products": [], "best_viral": 0.0, "best_cell": 0.0})
    for orf in set(cellular) | set(viral):
        c_score, c_prod = cellular.get(orf, (0.0, ""))
        v_score, v_prod = viral.get(orf, (0.0, ""))
        rec = per_contig[contig_of(orf)]
        rec["n_orfs"] += 1
        rec["best_viral"] = max(rec["best_viral"], v_score)
        rec["best_cell"] = max(rec["best_cell"], c_score)
        if v_score > c_score * args.score_ratio:
            rec["viral"] += 1
        elif c_score > v_score * args.score_ratio:
            rec["cellular"] += 1
        else:
            rec["ambiguous"] += 1
        if c_prod and euk_re.search(c_prod):
            rec["euk_hits"] += 1
            if c_prod not in rec["euk_products"]:
                rec["euk_products"].append(c_prod)

    rows, counts = [], Counter()
    for contig, rec in per_contig.items():
        decided = rec["viral"] + rec["cellular"]
        viral_frac = (rec["viral"] / decided) if decided else 0.0
        tax = taxonomy.get(contig, {})
        n_tools = tax.get("n_tools", 0)
        has_rdrp = contig in rdrp or contig.split("|")[-1] in rdrp

        # 3: a replicase is definitive; otherwise require a viral majority,
        # no eukaryotic marker, and corroborating taxonomy
        if has_rdrp:
            level, why = "high", "RdRp domain extracted"
        elif (viral_frac >= args.min_viral_frac and rec["euk_hits"] == 0
              and n_tools >= args.min_tools):
            level, why = "high", (f"{viral_frac:.0%} of decided ORFs viral, no "
                                  f"eukaryotic marker, {n_tools:.0f} tools")
        elif viral_frac >= args.min_viral_frac and rec["euk_hits"] == 0:
            level, why = "medium", (f"{viral_frac:.0%} viral but only "
                                    f"{n_tools:.0f} taxonomy tool(s)")
        elif rec["viral"] > 0 and rec["euk_hits"] > 0:
            level, why = "medium", (f"{rec['viral']} viral ORF(s) alongside "
                                    f"{rec['euk_hits']} eukaryotic marker(s): "
                                    f"possible chimera or integrated element")
        elif rec["viral"] > 0:
            level, why = "medium", f"{viral_frac:.0%} of decided ORFs viral"
        elif rec["euk_hits"] > 0:
            level, why = "low", (f"eukaryotic marker(s) and no RdRp: "
                                 f"{'; '.join(rec['euk_products'][:3])}")
        else:
            level, why = "low", (f"{rec['cellular']} cellular vs "
                                 f"{rec['viral']} viral ORF(s), no RdRp")
        counts[level] += 1

        rows.append({
            "contig": contig,
            "run_id": tax.get("run_id", contig.split("|")[0] if "|" in contig else ""),
            "confidence": level,
            "keep": "1" if level in keep_levels else "0",
            "has_rdrp": "1" if has_rdrp else "0",
            "n_orfs": rec["n_orfs"],
            "n_orfs_viral": rec["viral"],
            "n_orfs_cellular": rec["cellular"],
            "n_orfs_ambiguous": rec["ambiguous"],
            "viral_fraction": f"{viral_frac:.3f}",
            "best_viral_score": f"{rec['best_viral']:.1f}",
            "best_cellular_score": f"{rec['best_cell']:.1f}",
            "n_eukaryotic_markers": rec["euk_hits"],
            "eukaryotic_products": "; ".join(rec["euk_products"][:5]),
            "consensus_family": tax.get("consensus_family", ""),
            "family_agreement": tax.get("family_agreement", ""),
            "n_tools": f"{n_tools:.0f}" if n_tools else "",
            "checkv_completeness": ("" if tax.get("checkv_completeness", -1) < 0
                                    else f"{tax['checkv_completeness']:.1f}"),
            "reason": why,
        })

    order = {"high": 0, "medium": 1, "low": 2}
    rows.sort(key=lambda r: (order[r["confidence"]], -int(r["n_orfs_viral"]),
                             r["contig"]))
    cols = list(rows[0].keys()) if rows else ["contig"]
    os.makedirs(os.path.dirname(os.path.abspath(args.out_tsv)) or ".", exist_ok=True)
    with open(args.out_tsv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", lineterminator="\n")
        w.writeheader()
        w.writerows(rows)

    kept = [r for r in rows if r["keep"] == "1"]
    if args.out_keep:
        with open(args.out_keep, "w") as fh:
            fh.write("".join(f"{r['contig']}\n" for r in kept))

    fams = Counter(r["consensus_family"] for r in kept if r["consensus_family"])
    if args.out_summary:
        with open(args.out_summary, "w", newline="") as fh:
            w = csv.writer(fh, delimiter="\t", lineterminator="\n")
            w.writerow(["metric", "value"])
            w.writerow(["n_contigs", len(rows)])
            for lvl in ("high", "medium", "low"):
                w.writerow([f"n_{lvl}", counts[lvl]])
            w.writerow(["n_kept", len(kept)])
            w.writerow(["n_with_rdrp", sum(1 for r in rows if r["has_rdrp"] == "1")])
            w.writerow(["n_with_eukaryotic_marker",
                        sum(1 for r in rows if int(r["n_eukaryotic_markers"]) > 0)])
            w.writerow(["n_families_in_kept_set", len(fams)])
            w.writerow(["score_ratio", args.score_ratio])
            w.writerow(["min_viral_fraction", args.min_viral_frac])
            w.writerow(["min_tools_for_high", args.min_tools])
            w.writerow(["kept_levels", args.keep])

    log(f"[filter] {len(rows)} contig(s): high {counts['high']}, "
        f"medium {counts['medium']}, low {counts['low']}; "
        f"{len(kept)} kept ({args.keep})")
    log(f"[filter] with an RdRp: {sum(1 for r in rows if r['has_rdrp']=='1')}; "
        f"with a eukaryotic marker: "
        f"{sum(1 for r in rows if int(r['n_eukaryotic_markers'])>0)}")
    # 4: the analysis family list follows from the contigs that survived, so
    # the two filters cannot disagree about what is real
    if args.out_families_yaml:
        chosen = sorted(f for f, n in fams.items() if n >= args.min_kept_contigs)
        with open(args.out_families_yaml, "w") as fh:
            fh.write("# Generated by filter_contigs_by_annotation.py -- do not "
                     "hand-edit.\n")
            fh.write(f"# Families with >= {args.min_kept_contigs} contig(s) of "
                     f"{args.keep} confidence.\n")
            fh.write("phylogenetics:\n  families:\n")
            for f in chosen:
                fh.write(f"    - {f}\n")
            fh.write("novel:\n  families:\n")
            for f in chosen:
                fh.write(f"    - {f}\n")
        log(f"[filter] {len(chosen)} family(ies) with >="
            f"{args.min_kept_contigs} kept contig(s) -> {args.out_families_yaml}")
    log(f"[filter] families surviving in the kept set: {len(fams)}")
    log(f"[filter] -> {args.out_tsv}")
    return 0 if kept else 1


if __name__ == "__main__":
    raise SystemExit(main())
