#!/usr/bin/env python3
"""Merge per-sample viral calls from Kraken2, geNomad and RdRpCATCH.

Two subcommands, one per stage of rules/05.merge_tax.smk:

  collect   union of every contig called viral by ANY of the three tools
            -> tool_calls.tsv (one row per positive contig, one column
            block per tool) + viral_positives.fasta
  annotate  joins CheckV quality (run on viral_positives.fasta) and the
            upstream VirSorter2 score onto tool_calls.tsv as annotation
            columns -- nothing is dropped here, every positive is kept.

consensus_family is a convenience column: majority vote across the tools
that resolved a family, ties broken geNomad > RdRpCATCH > Kraken2 (geNomad
is protein/marker-based and scored best on the CAMI2 benchmark).
"""
import argparse
import csv
import os
from collections import Counter, defaultdict

import pandas as pd
from Bio import SeqIO

import family_taxonomy

TOOL_PRIORITY = ["genomad", "rdrpcatch", "kraken2"]


# --------------------------------------------------------------------------
# Kraken2
# --------------------------------------------------------------------------
def load_ktaxonomy(path):
    """ktaxonomy.tsv: taxid | parent | rank | depth | name"""
    nodes = {}
    with open(path) as f:
        for line in f:
            parts = [p.strip() for p in line.rstrip("\n").split("|")]
            if len(parts) < 5:
                continue
            nodes[parts[0]] = (parts[1], parts[2], parts[4])
    return nodes


def kraken2_lineage(taxid, nodes):
    names, family = [], ""
    seen = set()
    while taxid in nodes and taxid not in seen:
        seen.add(taxid)
        parent, rank, name = nodes[taxid]
        if rank != "R":  # skip "root"
            names.append(name)
        if rank == "F":
            family = name
        if parent == taxid:
            break
        taxid = parent
    return ";".join(reversed(names)), family


def parse_kraken2(output_tsv, ktaxonomy):
    nodes = load_ktaxonomy(ktaxonomy)
    calls = {}
    if not os.path.exists(output_tsv):
        return calls
    with open(output_tsv) as f:
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3 or parts[0] != "C":
                continue
            # --use-names writes "Name (taxid N)"; plain mode writes "N".
            tax_field = parts[2]
            taxid = tax_field.rsplit("taxid", 1)[-1].strip(" )") if "taxid" in tax_field else tax_field.strip()
            lineage, family = kraken2_lineage(taxid, nodes)
            calls[parts[1]] = {
                "kraken2_taxid": taxid,
                "kraken2_lineage": lineage,
                "kraken2_family": family,
            }
    return calls


# --------------------------------------------------------------------------
# geNomad
# --------------------------------------------------------------------------
def parse_genomad(virus_summary):
    calls = {}
    if not os.path.exists(virus_summary) or os.path.getsize(virus_summary) == 0:
        return calls
    df = pd.read_csv(virus_summary, sep="\t", dtype=str)
    for _, row in df.iterrows():
        # Proviruses are reported as "<contig>|provirus_<start>_<end>".
        seq_name = row["seq_name"]
        contig = seq_name.split("|", 1)[0]
        lineage = row.get("taxonomy", "")
        lineage = "" if pd.isna(lineage) else lineage
        family = next((r for r in lineage.split(";") if r.endswith("viridae")), "")
        prev = calls.get(contig)
        score = float(row.get("virus_score", 0) or 0)
        if prev and prev["genomad_score"] >= score:
            continue
        calls[contig] = {
            "genomad_seq_name": seq_name,
            "genomad_score": score,
            "genomad_n_hallmarks": row.get("n_hallmarks", ""),
            "genomad_lineage": lineage,
            "genomad_family": family,
        }
    return calls


# --------------------------------------------------------------------------
# RdRpCATCH
# --------------------------------------------------------------------------
def rdrp_family_name(db, profile):
    # Same family-level profile selection as
    # virus_discovery_comparison/scripts/extract_sra_family_data.py.
    if db == "RVMT":
        name = profile.replace("_RVMT", "")
        if not (name.startswith("f.") or name.startswith("c.")):
            return name
    elif db == "TSA_Olendraite_fam":
        return profile.replace(".aln", "")
    return ""


def parse_rdrpcatch(tsv):
    calls = {}
    if not os.path.exists(tsv) or os.path.getsize(tsv) == 0:
        return calls
    df = pd.read_csv(tsv, sep="\t", dtype=str)
    if df.empty:
        return calls
    contig_col = "Contig_name" if "Contig_name" in df.columns else df.columns[0]
    has_hits = {"Profile_name", "db_name"}.issubset(df.columns)
    for contig, grp in df.groupby(contig_col):
        call = {"rdrpcatch_profile": "", "rdrpcatch_db": "",
                "rdrpcatch_evalue": "", "rdrpcatch_family": ""}
        if has_hits:
            grp = grp.copy()
            if "score" in grp.columns:
                grp["_s"] = pd.to_numeric(grp["score"], errors="coerce")
                best = grp.sort_values("_s", ascending=False).iloc[0]
            else:
                best = grp.iloc[0]
            fams = [rdrp_family_name(d, p) for d, p in zip(grp["db_name"], grp["Profile_name"])]
            fams = [f for f in fams if f]
            call.update({
                "rdrpcatch_profile": best["Profile_name"],
                "rdrpcatch_db": best["db_name"],
                "rdrpcatch_evalue": best.get("E-value", ""),
                "rdrpcatch_family": Counter(fams).most_common(1)[0][0] if fams else "",
            })
        calls[contig] = call
    return calls


# --------------------------------------------------------------------------
# collect / annotate
# --------------------------------------------------------------------------
def consensus(row):
    """consensus_family, family_agreement, family_disagreement.

    Delegated to family_taxonomy.reconcile: comparing the names as plain
    strings called Flaviviridae vs Amarillovirales a conflict and dropped the
    contig, when the two are the same taxon at different ranks. See that
    module for the statuses.
    """
    fams = {t: row.get(f"{t}_family", "") for t in TOOL_PRIORITY}
    fams = {t: f for t, f in fams.items() if isinstance(f, str) and f}
    return family_taxonomy.reconcile(fams, TOOL_PRIORITY)


def collect(args):
    per_tool = {
        "kraken2": parse_kraken2(args.kraken2, args.ktaxonomy),
        "genomad": parse_genomad(args.genomad),
        "rdrpcatch": parse_rdrpcatch(args.rdrpcatch),
    }
    merged = defaultdict(dict)
    for tool, calls in per_tool.items():
        for contig, fields in calls.items():
            merged[contig].update(fields)
            merged[contig].setdefault("_tools", []).append(tool)

    columns = [
        "contig_id", "detected_by", "n_tools", "consensus_family", "family_agreement",
        "family_disagreement",
        "genomad_seq_name", "genomad_score", "genomad_n_hallmarks", "genomad_lineage", "genomad_family",
        "rdrpcatch_profile", "rdrpcatch_db", "rdrpcatch_evalue", "rdrpcatch_family",
        "kraken2_taxid", "kraken2_lineage", "kraken2_family",
    ]
    rows = []
    for contig, fields in merged.items():
        tools = [t for t in TOOL_PRIORITY if t in fields["_tools"]]
        row = {k: v for k, v in fields.items() if k != "_tools"}
        row.update({"contig_id": contig, "detected_by": ",".join(tools), "n_tools": len(tools)})
        (row["consensus_family"], row["family_agreement"],
         row["family_disagreement"]) = consensus(row)
        rows.append(row)
    df = pd.DataFrame(rows, columns=columns).sort_values(["n_tools", "contig_id"], ascending=[False, True])
    df.to_csv(args.out_tsv, sep="\t", index=False)

    keep = set(merged)
    with open(args.out_fasta, "w") as out:
        SeqIO.write((r for r in SeqIO.parse(args.contigs, "fasta") if r.id in keep), out, "fasta")


def annotate(args):
    calls = pd.read_csv(args.calls, sep="\t", dtype=str)

    checkv_cols = ["contig_id", "checkv_quality", "completeness", "contamination", "viral_genes", "host_genes", "warnings"]
    if os.path.getsize(args.checkv) > 0:
        checkv = pd.read_csv(args.checkv, sep="\t", dtype=str)
        checkv = checkv[[c for c in checkv_cols if c in checkv.columns]]
        checkv = checkv.rename(columns={c: f"checkv_{c}" for c in checkv.columns
                                        if c not in ("contig_id", "checkv_quality")})
        calls = calls.merge(checkv, on="contig_id", how="left")

    # VirSorter2 seqnames carry a "||full"/"||lt2gene"/"||<n>_partial"
    # suffix; keep the best-scoring call per original contig.
    vs2_cols = ["vs2_max_score", "vs2_max_score_group"]
    if os.path.exists(args.virsorter2) and os.path.getsize(args.virsorter2) > 0:
        vs2 = pd.read_csv(args.virsorter2, sep="\t", dtype=str)
        if not vs2.empty:
            vs2["contig_id"] = vs2["seqname"].str.split("||", regex=False).str[0]
            vs2["vs2_max_score"] = pd.to_numeric(vs2["max_score"], errors="coerce")
            vs2["vs2_max_score_group"] = vs2["max_score_group"]
            vs2 = vs2.sort_values("vs2_max_score", ascending=False).drop_duplicates("contig_id")
            calls = calls.merge(vs2[["contig_id"] + vs2_cols], on="contig_id", how="left")
    for c in vs2_cols:
        if c not in calls.columns:
            calls[c] = ""
    calls["virsorter2_viral"] = calls["vs2_max_score"].notna() & (calls["vs2_max_score"].astype(str) != "")

    calls.to_csv(args.out, sep="\t", index=False)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = p.add_subparsers(dest="cmd", required=True)

    c = sub.add_parser("collect")
    c.add_argument("--contigs", required=True)
    c.add_argument("--kraken2", required=True, help="kraken2 --output tsv")
    c.add_argument("--ktaxonomy", required=True, help="kraken2 DB ktaxonomy.tsv")
    c.add_argument("--genomad", required=True, help="geNomad *_virus_summary.tsv")
    c.add_argument("--rdrpcatch", required=True, help="RdRpCATCH annotated/combined tsv")
    c.add_argument("--out-tsv", required=True)
    c.add_argument("--out-fasta", required=True)
    c.set_defaults(func=collect)

    a = sub.add_parser("annotate")
    a.add_argument("--calls", required=True)
    a.add_argument("--checkv", required=True)
    a.add_argument("--virsorter2", required=True, help="VirSorter2 final-viral-score.tsv")
    a.add_argument("--out", required=True)
    a.set_defaults(func=annotate)

    args = p.parse_args()
    args.func(args)


if __name__ == "__main__":
    main()
