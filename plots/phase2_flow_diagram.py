"""
Phase 2 (assembly & annotation) pipeline flow diagram.

Draws the actual current rule structure of workflow/phase2/Snakefile_phase2
and its included rules/*.smk files, one box per Snakemake rule, in the
order data actually flows through them. Regenerate this whenever phase2's
rule structure changes materially -- update this script's STAGES list to
match, don't just re-run it against stale assumptions. Last redrawn for
the three-tool taxonomy + MetaCerberus structure (2026-11-02).

TOOL_VERSIONS below are a snapshot queried directly from each rule's
actual conda env in workflow/phase2/.snakemake/conda/ on 2026-11-02 (e.g.
`mmseqs version`, `virsorter --version`), not hand-typed from memory or
from the env yaml pins (which can drift from what's actually installed).
Re-query and update if envs get rebuilt.

Usage:
    python phase2_flow_diagram.py
Output:
    phase2_flow_diagram.png (written next to this script)
"""

import os

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "phase2_flow_diagram.png")

# Queried 2026-11-02 from each rule's actual .snakemake/conda/<hash>_ env
# (`<tool> --version`/`version`), not from env yaml pins or memory -- see
# module docstring. None = no single version string applies (internal
# python/pandas scripts).
TOOL_VERSIONS = {
    "input": "sra-tools 3.1.1",
    "fastp": "fastp 1.4.0",
    "assembly": "SPAdes 4.3.0 / bwa-mem3 0.14.0 / Samtools 1.24",
    "filter": "Samtools 1.24 / Seqkit",
    "virsorter2": "VirSorter2 2.2.4",
    "kraken2": "Kraken2 2.17.2 (k2_viral 2026-06-26)",
    "genomad": "geNomad 1.12.0 (DB v1.9)",
    "rdrpcatch": "RdRpCATCH 0.0.9",
    "checkv": "CheckV 1.1.1",
    "metacerberus": "MetaCerberus 1.4.0",
}

# Each stage: (label, rule file, column, row). Column groups stages that
# run in parallel at the same point in the pipeline; row orders them
# top-to-bottom.
STAGES = {
    "input":      ("Stage reads\n(SRA fetch or local)",                 "01.input.smk",      0, 0),
    "fastp":      ("fastp_cleaning",                                    "02.fastp.smk",      0, 1),
    "assembly":   ("rnaviral_assembly\n(rnaviralSPAdes)",               "03.assembly.smk",   0, 2),
    "filter":     ("filter_contigs\n(drop zero-coverage)",              "03.assembly.smk",   0, 3),
    "kraken2":    ("run_kraken2\n(k-mer, known viruses)",               "04.taxonomy.smk",   -1.5, 4.3),
    "genomad":    ("run_genomad\n(markers + NN, ICTV taxonomy)",        "04.taxonomy.smk",   -0.5, 4.3),
    "rdrpcatch":  ("run_rdrpcatch\n(RdRp HMMs)",                        "04.taxonomy.smk",   0.5, 4.3),
    "virsorter2": ("run_virsorter2\n(viral score, annotation only)",    "03.assembly.smk",   1.5, 4.3),
    "collect":    ("collect_viral_positives\n(union: any tool positive)", "05.merge_tax.smk", -0.5, 5.6),
    "checkv":     ("run_checkv\n(quality, annotation only)",            "04.taxonomy.smk",   -0.5, 6.6),
    "annotate":   ("annotate_viral_positives\n(+ CheckV, VirSorter2 columns)", "05.merge_tax.smk", 0, 7.6),
    "merge":      ("build_master_database\n(all samples)",              "05.merge_tax.smk",  0, 8.6),
    "select":     ("select_contigs_for_annotation\n(positives >= 1 kb)", "06.annotation.smk", 0, 9.6),
    "metacerberus": ("run_metacerberus\n(KOFam / VOG / PHROG)",         "06.annotation.smk", 0, 10.6),
    "final":      ("results/master_viral_database.fasta\nresults/annotation/metacerberus/", None, 0, 11.6),
}

EDGES = [
    ("input", "fastp"),
    ("fastp", "assembly"),
    ("assembly", "filter"),
    ("filter", "kraken2"),
    ("filter", "genomad"),
    ("filter", "rdrpcatch"),
    ("filter", "virsorter2"),
    ("kraken2", "collect"),
    ("genomad", "collect"),
    ("rdrpcatch", "collect"),
    ("collect", "checkv"),
    ("checkv", "annotate"),
    ("virsorter2", "annotate"),
    ("annotate", "merge"),
    ("merge", "select"),
    ("select", "metacerberus"),
    ("metacerberus", "final"),
]

# Colour by rule file so it's visually obvious which .smk each stage lives
# in -- not a measured/data palette, just a fixed lookup for this diagram.
FILE_COLORS = {
    "01.input.smk": "#8ecae6",
    "02.fastp.smk": "#fb8500",
    "03.assembly.smk": "#8338ec",
    "04.taxonomy.smk": "#06a77d",
    "05.merge_tax.smk": "#d62828",
    "06.annotation.smk": "#ffb703",
    None: "#e9ecef",
}

BOX_W, BOX_H = 2.6, 1.05
COL_SPACING = 2.9
ROW_SPACING = 1.55


def stage_xy(col, row):
    return col * COL_SPACING, -row * ROW_SPACING


def draw():
    fig, ax = plt.subplots(figsize=(11, 17))
    ax.set_xlim(-6, 6)
    ax.set_ylim(-19.5, 2)
    ax.axis("off")
    ax.set_title(
        "Phase 2 data flow (three-tool taxonomy + MetaCerberus, 2026-11-02)\n"
        "workflow/phase2/Snakefile_phase2 + rules/*.smk",
        fontsize=13, fontweight="bold", pad=20,
    )

    centers = {}
    for key, (label, rule_file, col, row) in STAGES.items():
        x, y = stage_xy(col, row)
        centers[key] = (x, y)
        color = FILE_COLORS.get(rule_file, "#e9ecef")
        box = FancyBboxPatch(
            (x - BOX_W / 2, y - BOX_H / 2), BOX_W, BOX_H,
            boxstyle="round,pad=0.08,rounding_size=0.12",
            linewidth=1.2, edgecolor="#343a40", facecolor=color, alpha=0.85,
            zorder=2,
        )
        ax.add_patch(box)
        version = TOOL_VERSIONS.get(key)
        if version:
            ax.text(x, y + 0.12, label, ha="center", va="center", fontsize=8.5,
                     fontweight="medium", color="#212529", zorder=3, wrap=True)
            ax.text(x, y - BOX_H / 2 + 0.16, version, ha="center", va="center",
                     fontsize=5.8, style="italic", color="#212529", zorder=3)
        else:
            ax.text(x, y, label, ha="center", va="center", fontsize=8.5,
                     fontweight="medium", color="#212529", zorder=3, wrap=True)

    for src, dst in EDGES:
        x1, y1 = centers[src]
        x2, y2 = centers[dst]
        arrow = FancyArrowPatch(
            (x1, y1 - BOX_H / 2 - 0.02), (x2, y2 + BOX_H / 2 + 0.02),
            arrowstyle="-|>", mutation_scale=14, linewidth=1.1,
            color="#495057", zorder=1,
            connectionstyle="arc3,rad=0.0" if x1 == x2 else "arc3,rad=0.15",
        )
        ax.add_patch(arrow)

    # Legend: which colour = which rule file
    legend_x, legend_y = 3.3, 1.0
    ax.text(legend_x, legend_y + 0.3, "rules/", fontsize=9, fontweight="bold")
    seen = []
    for rule_file, color in FILE_COLORS.items():
        if rule_file is None or rule_file in seen:
            continue
        seen.append(rule_file)
        legend_y -= 0.35
        ax.add_patch(plt.Rectangle((legend_x, legend_y - 0.1), 0.25, 0.2, facecolor=color, edgecolor="#343a40"))
        ax.text(legend_x + 0.35, legend_y, rule_file, fontsize=7.5, va="center")

    fig.tight_layout()
    fig.savefig(OUT_PATH, dpi=150, bbox_inches="tight")
    print(f"Wrote {OUT_PATH}")


if __name__ == "__main__":
    draw()
