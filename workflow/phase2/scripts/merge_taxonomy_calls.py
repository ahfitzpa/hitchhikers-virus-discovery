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
from __future__ import annotations

import argparse
import logging
import os
import tempfile
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

import pandas as pd
from Bio import SeqIO

import family_taxonomy

__version__ = "1.0.0"
LOGGER = logging.getLogger(__name__)
TOOL_PRIORITY = ["genomad", "rdrpcatch", "kraken2"]


def _atomic_write_tsv(df: pd.DataFrame, dest: str | Path) -> None:
    """Write a DataFrame as TSV to a sibling temp file and atomically rename into place."""
    dest_path = Path(dest)
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{dest_path.name}.", suffix=".tmp", dir=dest_path.parent)
    os.close(fd)
    try:
        df.to_csv(tmp_name, sep="\t", index=False)
        os.replace(tmp_name, dest_path)
    except BaseException:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)
        raise


def _atomic_write_fasta(records: Any, dest: str | Path) -> None:
    """Write SeqIO records to a sibling temp file and atomically rename into place."""
    dest_path = Path(dest)
    dest_path.parent.mkdir(parents=True, exist_ok=True)
    fd, tmp_name = tempfile.mkstemp(prefix=f".{dest_path.name}.", suffix=".tmp", dir=dest_path.parent)
    try:
        with os.fdopen(fd, "w", encoding="utf-8") as out:
            SeqIO.write(records, out, "fasta")
        os.replace(tmp_name, dest_path)
    except BaseException:
        if os.path.exists(tmp_name):
            os.unlink(tmp_name)
        raise


# --------------------------------------------------------------------------
# Kraken2
# --------------------------------------------------------------------------
def load_ktaxonomy(path: str | Path) -> dict[str, tuple[str, str, str]]:
    """Parse Kraken2 ktaxonomy.tsv (taxid | parent | rank | depth | name)."""
    nodes: dict[str, tuple[str, str, str]] = {}
    with Path(path).open(encoding="utf-8") as fh:
        for line in fh:
            parts = [p.strip() for p in line.rstrip("\n").split("|")]
            if len(parts) < 5:
                continue
            nodes[parts[0]] = (parts[1], parts[2], parts[4])
    return nodes


def kraken2_lineage(taxid: str, nodes: dict[str, tuple[str, str, str]]) -> tuple[str, str]:
    """Walk ktaxonomy parent links to return (semicolon lineage, family name)."""
    names: list[str] = []
    family = ""
    seen: set[str] = set()
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


def parse_kraken2(output_tsv: str | Path, ktaxonomy: str | Path) -> dict[str, dict[str, Any]]:
    """Parse classified Kraken2 output rows and resolve taxonomy lineage and family."""
    nodes = load_ktaxonomy(ktaxonomy)
    calls: dict[str, dict[str, Any]] = {}
    out_path = Path(output_tsv)
    if not out_path.exists():
        return calls
    with out_path.open(encoding="utf-8") as fh:
        for line in fh:
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
def parse_genomad(virus_summary: str | Path) -> dict[str, dict[str, Any]]:
    """Parse geNomad virus_summary.tsv keeping the highest-scoring call per contig."""
    calls: dict[str, dict[str, Any]] = {}
    summary_path = Path(virus_summary)
    if not summary_path.exists() or summary_path.stat().st_size == 0:
        return calls
    df = pd.read_csv(summary_path, sep="\t", dtype=str)
    for _, row in df.iterrows():
        # Proviruses are reported as "<contig>|provirus_<start>_<end>".
        seq_name = str(row["seq_name"])
        contig = seq_name.split("|", 1)[0]
        raw_lineage = row.get("taxonomy", "")
        lineage = "" if pd.isna(raw_lineage) else str(raw_lineage)
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
def rdrp_family_name(db: str, profile: str) -> str:
    """Extract a family-level name from RVMT or TSA_Olendraite_fam profile IDs."""
    if db == "RVMT":
        name = profile.replace("_RVMT", "")
        if not (name.startswith("f.") or name.startswith("c.")):
            return name
    elif db == "TSA_Olendraite_fam":
        return profile.replace(".aln", "")
    return ""


def parse_rdrpcatch(tsv: str | Path) -> dict[str, dict[str, Any]]:
    """Parse RdRpCATCH hits TSV and retain the top-scoring hit and consensus family per contig."""
    calls: dict[str, dict[str, Any]] = {}
    tsv_path = Path(tsv)
    if not tsv_path.exists() or tsv_path.stat().st_size == 0:
        return calls
    df = pd.read_csv(tsv_path, sep="\t", dtype=str)
    if df.empty:
        return calls
    contig_col = "Contig_name" if "Contig_name" in df.columns else df.columns[0]
    has_hits = {"Profile_name", "db_name"}.issubset(df.columns)
    for contig, grp in df.groupby(contig_col):
        call: dict[str, Any] = {
            "rdrpcatch_profile": "",
            "rdrpcatch_db": "",
            "rdrpcatch_evalue": "",
            "rdrpcatch_family": "",
        }
        if has_hits:
            grp = grp.copy()
            if "score" in grp.columns:
                grp["_s"] = pd.to_numeric(grp["score"], errors="coerce")
                best = grp.sort_values("_s", ascending=False).iloc[0]
            else:
                best = grp.iloc[0]
            fams = [rdrp_family_name(str(d), str(p)) for d, p in zip(grp["db_name"], grp["Profile_name"])]
            fams = [f for f in fams if f]
            call.update({
                "rdrpcatch_profile": best["Profile_name"],
                "rdrpcatch_db": best["db_name"],
                "rdrpcatch_evalue": best.get("E-value", ""),
                "rdrpcatch_family": Counter(fams).most_common(1)[0][0] if fams else "",
            })
        calls[str(contig)] = call
    return calls


# --------------------------------------------------------------------------
# collect / annotate
# --------------------------------------------------------------------------
def consensus(row: dict[str, Any]) -> tuple[str, str, str]:
    """Reconcile per-tool family calls via family_taxonomy.reconcile."""
    fams = {t: row.get(f"{t}_family", "") for t in TOOL_PRIORITY}
    fams = {t: f for t, f in fams.items() if isinstance(f, str) and f}
    return family_taxonomy.reconcile(fams, TOOL_PRIORITY)


def collect(args: argparse.Namespace) -> None:
    """Union viral calls from Kraken2, geNomad, and RdRpCATCH and emit TSV + FASTA."""
    per_tool = {
        "kraken2": parse_kraken2(args.kraken2, args.ktaxonomy),
        "genomad": parse_genomad(args.genomad),
        "rdrpcatch": parse_rdrpcatch(args.rdrpcatch),
    }
    merged: dict[str, dict[str, Any]] = defaultdict(dict)
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
    rows: list[dict[str, Any]] = []
    for contig, fields in merged.items():
        tools = [t for t in TOOL_PRIORITY if t in fields["_tools"]]
        row = {k: v for k, v in fields.items() if k != "_tools"}
        row.update({"contig_id": contig, "detected_by": ",".join(tools), "n_tools": len(tools)})
        (
            row["consensus_family"],
            row["family_agreement"],
            row["family_disagreement"],
        ) = consensus(row)
        rows.append(row)
    df = pd.DataFrame(rows, columns=columns).sort_values(["n_tools", "contig_id"], ascending=[False, True])
    _atomic_write_tsv(df, args.out_tsv)

    keep = set(merged)
    _atomic_write_fasta((r for r in SeqIO.parse(args.contigs, "fasta") if r.id in keep), args.out_fasta)
    LOGGER.info("collect: wrote %d viral positive contigs to %s", len(df), args.out_tsv)


def annotate(args: argparse.Namespace) -> None:
    """Left-join CheckV quality and VirSorter2 scores onto tool_calls.tsv."""
    calls = pd.read_csv(args.calls, sep="\t", dtype=str)

    checkv_cols = [
        "contig_id", "checkv_quality", "completeness", "contamination",
        "viral_genes", "host_genes", "warnings",
    ]
    checkv_path = Path(args.checkv)
    if checkv_path.exists() and checkv_path.stat().st_size > 0:
        checkv = pd.read_csv(checkv_path, sep="\t", dtype=str)
        checkv = checkv[[c for c in checkv_cols if c in checkv.columns]]
        checkv = checkv.rename(
            columns={c: f"checkv_{c}" for c in checkv.columns if c not in ("contig_id", "checkv_quality")}
        )
        calls = calls.merge(checkv, on="contig_id", how="left")

    # VirSorter2 seqnames carry a "||full"/"||lt2gene"/"||<n>_partial" suffix; keep best score per contig.
    vs2_cols = ["vs2_max_score", "vs2_max_score_group"]
    vs2_path = Path(args.virsorter2)
    if vs2_path.exists() and vs2_path.stat().st_size > 0:
        vs2 = pd.read_csv(vs2_path, sep="\t", dtype=str)
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

    _atomic_write_tsv(calls, args.out)
    LOGGER.info("annotate: wrote %d annotated rows to %s", len(calls), args.out)


def main(argv: list[str] | None = None) -> int:
    """CLI entry point for collect and annotate subcommands."""
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])
    p.add_argument("--version", action="version", version=f"%(prog)s {__version__}")
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

    args = p.parse_args(argv)
    logging.basicConfig(level=getattr(logging, args.log_level), format="[%(levelname)s] %(message)s")
    args.func(args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
