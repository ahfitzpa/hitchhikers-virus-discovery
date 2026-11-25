"""
Phase 1 (SRA discovery) pipeline flow diagram.

Draws the actual current rule structure of workflow/phase1/Snakefile_phase1
and its included rules/*.smk files, one box per Snakemake rule, in the
order data actually flows through them. Mirrors phase2_flow_diagram.py's
style. Regenerate whenever phase1's rule structure changes materially.

Usage:
    python phase1_flow_diagram.py
Output:
    phase1_flow_diagram.png (written next to this script)
"""

import os

import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

OUT_PATH = os.path.join(os.path.dirname(os.path.abspath(__file__)), "phase1_flow_diagram.png")

# Each stage: (label, rule file, column, row).
STAGES = {
    "query":         ("centroid_step1_query\n(checkpoint)",            "01.base_discovery.smk",      0, 0),
    "novel_label":   ("discovery_type = novel",                        None,                          -2, 1),
    "known_label":   ("discovery_type = known",                        None,                           2, 1),

    # Unanchored (novel) branch
    "transform":     ("transform_db_to_hits",                          "01.centroid_unanchored.smk",  -2, 2),

    # Anchored (known) branch
    "combine":       ("combine_fasta",                                 "01.centroid_anchored.smk",     2, 2),
    "dedup_split":   ("centroid_dedup_and_split",                      "01.centroid_anchored.smk",     2, 3),
    "cluster_chunk": ("centroid_cluster_chunk\n(x16, parallel)",       "01.centroid_anchored.smk",     2, 4),
    "merge_querydb": ("centroid_merge_querydb",                        "01.centroid_anchored.smk",     2, 5),
    "build_refdb":   ("mmseqs_build_refdb",                            "01.centroid_anchored.smk",     3.3, 5),
    "search":        ("mmseqs_search",                                 "01.centroid_anchored.smk",     2, 6),
    "expand":        ("expand_cluster_hits",                           "01.centroid_anchored.smk",     2, 7),

    # Convergence
    "extract_sra":   ("extract_sra_ids",                               "01.centroid_anchored.smk /\n01.centroid_unanchored.smk", 0, 8),
    "build_matrix":  ("build_matrix /\nbuild_matrix_novel",            "01.centroid_anchored.smk /\n01.centroid_unanchored.smk", -1.3, 8),

    "fetch_meta":    ("fetch_sra_metadata",                             "03.metadata_fetch.smk",       0, 9),
    "merge_meta":    ("merge_sra_metadata\n(checkpoint)",               "04.metadata_merge.smk",       0, 10),
    "review_gate":   ("automated_review_gate\n(known only)",            "04.metadata_merge.smk",       1.3, 11),
    "final":         ("sra_ids.txt (novel) /\nunique.csv + approved.flag (known)", None, 0, 12),
}

EDGES = [
    ("query", "novel_label"), ("query", "known_label"),
    ("novel_label", "transform"),
    ("known_label", "combine"),
    ("combine", "dedup_split"),
    ("dedup_split", "cluster_chunk"),
    ("cluster_chunk", "merge_querydb"),
    ("build_refdb", "search"),
    ("merge_querydb", "search"),
    ("search", "expand"),
    ("expand", "extract_sra"),
    ("expand", "build_matrix"),
    ("transform", "extract_sra"),
    ("transform", "build_matrix"),
    ("extract_sra", "fetch_meta"),
    ("build_matrix", "fetch_meta"),
    ("fetch_meta", "merge_meta"),
    ("merge_meta", "review_gate"),
    ("merge_meta", "final"),
    ("review_gate", "final"),
]

# Colour by rule file so it's visually obvious which .smk each stage lives in.
FILE_COLORS = {
    "01.base_discovery.smk": "#8ecae6",
    "01.centroid_unanchored.smk": "#ffb703",
    "01.centroid_anchored.smk": "#8338ec",
    "01.centroid_anchored.smk /\n01.centroid_unanchored.smk": "#fb8500",
    "03.metadata_fetch.smk": "#219ebc",
    "04.metadata_merge.smk": "#06a77d",
    None: "#e9ecef",
}

BOX_W, BOX_H = 2.7, 1.0
COL_SPACING = 2.5
ROW_SPACING = 1.5


def stage_xy(col, row):
    return col * COL_SPACING, -row * ROW_SPACING


def draw():
    fig, ax = plt.subplots(figsize=(13, 16))
    ax.set_xlim(-8, 10)
    ax.set_ylim(-19.5, 2)
    ax.axis("off")
    ax.set_title(
        "Phase 1 data flow\n"
        "workflow/phase1/Snakefile_phase1 + rules/*.smk",
        fontsize=13, fontweight="bold", pad=20,
    )

    centers = {}
    for key, (label, rule_file, col, row) in STAGES.items():
        x, y = stage_xy(col, row)
        centers[key] = (x, y)
        color = FILE_COLORS.get(rule_file, "#e9ecef")
        is_label_only = key in ("novel_label", "known_label")
        box = FancyBboxPatch(
            (x - BOX_W / 2, y - BOX_H / 2), BOX_W, BOX_H,
            boxstyle="round,pad=0.08,rounding_size=0.12",
            linewidth=0 if is_label_only else 1.2,
            edgecolor="#343a40",
            facecolor="none" if is_label_only else color,
            alpha=0.85,
            zorder=2,
        )
        ax.add_patch(box)
        ax.text(x, y, label, ha="center", va="center",
                 fontsize=9.5 if is_label_only else 8.5,
                 fontweight="bold" if is_label_only else "medium",
                 style="italic" if is_label_only else "normal",
                 color="#495057" if is_label_only else "#212529",
                 zorder=3, wrap=True)

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
    legend_x, legend_y = 5.2, 1.0
    ax.text(legend_x, legend_y + 0.3, "rules/", fontsize=9, fontweight="bold")
    seen = []
    for rule_file, color in FILE_COLORS.items():
        if rule_file is None or rule_file in seen:
            continue
        seen.append(rule_file)
        legend_y -= 0.45
        ax.add_patch(plt.Rectangle((legend_x, legend_y - 0.1), 0.25, 0.2, facecolor=color, edgecolor="#343a40"))
        ax.text(legend_x + 0.35, legend_y, rule_file, fontsize=7, va="center")

    fig.tight_layout()
    fig.savefig(OUT_PATH, dpi=150, bbox_inches="tight")
    print(f"Wrote {OUT_PATH}")


if __name__ == "__main__":
    draw()
