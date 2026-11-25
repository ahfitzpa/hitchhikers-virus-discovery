#!/usr/bin/env python3
"""Build the master viral database and its metadata from the per-sample tables.

Run by rules/05.merge_tax.smk as a `script:`, NOT a `run:` block. Snakemake
ignores the rule's `conda:` directive for `run:` and executes the body in its
own interpreter, so `from Bio import SeqIO` raised ModuleNotFoundError and the
rule could never succeed -- which is what stopped job 1145155. A `script:` runs
inside the rule's conda environment, where biopython is installed.

Inputs come from the snakemake object: `tables` (one taxonomy_merged.tsv per
completed sample), `fastas` (the matching viral_positives.fasta), `meta_dict`
(the shortlist metadata) and `manifest` (so the excluded-sample report can name
every sample, which the Snakefile's MANIFEST global used to supply).
"""
import os
import sys

import pandas as pd
from Bio import SeqIO

snakemake = globals()["snakemake"]                      # injected by Snakemake

tables = list(snakemake.input.tables)
fastas = list(snakemake.input.fastas)
sep = "\t" if str(snakemake.input.manifest).endswith((".tsv", ".tab")) else ","
manifest = pd.read_csv(snakemake.input.manifest, sep=sep, dtype=str, comment="#").fillna("")
manifest.columns = [c.strip().lower() for c in manifest.columns]

# Record the exclusions in the results, not just in a log.
included = {t.split("/")[-2] for t in tables}
os.makedirs("results/taxonomy", exist_ok=True)
rows_ex = []
for _, row in manifest.iterrows():
    s_id, fam = row["sample_id"], row["family"]
    if s_id in included:
        continue
    # Last stage that produced output: where this sample stopped.
    trail = [
        ("reads", f"data/01.reads/{fam}/{s_id}/{s_id}_1.fastq.gz"),
        ("fastp", f"data/02.fastp/{fam}/{s_id}/{s_id}_cleaned_1.fastq.gz"),
        ("assembly", f"data/03.rnaviralSPAdes/{fam}/{s_id}/contigs.fasta"),
        ("filtered_contigs", f"data/03.rnaviralSPAdes/{fam}/{s_id}/filtered_contigs.fasta"),
        ("kraken2", f"data/04.taxonomy/{fam}/{s_id}/kraken2/kraken2_output.tsv"),
        ("rdrpcatch", f"data/04.taxonomy/{fam}/{s_id}/rdrpcatch/rdrp_annotated.tsv"),
        ("tool_calls", f"data/05.final_selection/{fam}/{s_id}/tool_calls.tsv"),
    ]
    reached = [name for name, pth in trail if os.path.exists(pth)]
    rows_ex.append(dict(sample_id=s_id, family=fam,
                        last_stage_completed=(reached[-1] if reached else "none"),
                        missing="taxonomy_merged.tsv",
                        excluded_from="master_viral_database.fasta, MetaCerberus"))
pd.DataFrame(rows_ex, columns=["sample_id", "family", "last_stage_completed",
                               "missing", "excluded_from"]).to_csv(
    snakemake.output.excluded, sep="\t", index=False)
if rows_ex:
    print(f"[build_master_database] WARNING: {len(rows_ex)} of {len(manifest)} samples "
          f"excluded (no taxonomy_merged.tsv): "
          f"{', '.join(r['sample_id'] for r in rows_ex)}. See {snakemake.output.excluded}.",
          file=sys.stderr)

meta_df = pd.read_csv(snakemake.input.meta_dict)
keep_cols = ["run_id", "sotu", "palm_id", "tax_family", "scientificname", "bioproject", "sample"]
meta_df = meta_df[[c for c in keep_cols if c in meta_df.columns]].drop_duplicates(subset=["run_id"])
# tax_family here is the shortlist's expected (Serratus palmprint) family, not a
# call from this pipeline -- keep it distinguishable from the per-tool
# *_family / consensus_family columns.
meta_df = meta_df.rename(columns={"tax_family": "shortlist_tax_family"})
meta_df.set_index("run_id", inplace=True)

per_sample = []
for table in tables:
    # sample_id is the 2nd-to-last path component: .../{sample}/taxonomy_merged.tsv
    sample = table.split("/")[-2]
    df = pd.read_csv(table, sep="\t", dtype=str)
    # run_id here is the shortlist's own join key; equals sample_id for
    # SRA-sourced rows, blank for local/ONT (degrades gracefully in the
    # left-merge below).
    df.insert(0, "run_id", sample)
    row = manifest[manifest["sample_id"] == sample].iloc[0]
    df["family_context"] = row["family"]
    df["scenario"] = row["bioproject"]
    per_sample.append(df)

full_df = pd.concat(per_sample, ignore_index=True).merge(meta_df, on="run_id", how="left")
full_df.to_csv(snakemake.output.master_csv, index=False)

# Keyed on (sample, contig_id): SPAdes NODE_* names repeat across samples, so
# contig_id alone would cross-wire metadata.
meta_map = full_df.set_index(["run_id", "contig_id"]).to_dict("index")
with open(snakemake.output.master_fa, "w") as fa_out:
    for fa_file in fastas:
        sample = fa_file.split("/")[-2]
        for rec in SeqIO.parse(fa_file, "fasta"):
            m = meta_map.get((sample, rec.id), {})
            tax = m.get("consensus_family") or "Unknown"
            tools = str(m.get("detected_by") or "Unknown").replace(",", "+")
            host = str(m.get("scientificname", "Unknown")).replace(" ", "_")
            fam = m.get("family_context", "Unknown")
            scen = m.get("scenario", "Unknown")
            rec.id = f"{sample}|{fam}|{scen}|{tax}|{tools}|{host}|{rec.id}"
            rec.description = ""
            SeqIO.write(rec, fa_out, "fasta")
print(f"[build_master_database] wrote {snakemake.output.master_fa} from {len(tables)} samples",
      file=sys.stderr)
