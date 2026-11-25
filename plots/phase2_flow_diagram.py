#!/usr/bin/env python3
"""Generate a clean, publication-grade SVG workflow diagram for Phase 2.

Run:
    python3 plots/phase2_flow_diagram.py
Writes:
    plots/phase2_flow_diagram.svg
"""
from __future__ import annotations

from html import escape
from pathlib import Path

OUT_SVG = Path(__file__).resolve().with_suffix(".svg")

COLORS = {
    "input": ("#e0f2fe", "#0284c7", "#0c4a6e"),      # 01.input.smk
    "qc": ("#ffedd5", "#ea580c", "#7c2d12"),         # 02.fastp.smk
    "asm": ("#ede9fe", "#6366f1", "#312e81"),        # 03.assembly.smk
    "tax": ("#d1fae5", "#059669", "#064e3b"),        # 04.taxonomy.smk
    "vs2": ("#f3e8ff", "#7c3aed", "#4c1d95"),        # VirSorter2 (annotation only)
    "merge": ("#ffe4e6", "#e11d48", "#881337"),      # 05.merge_tax.smk
    "annot": ("#fef3c7", "#d97706", "#78350f"),      # 06.annotation.smk
    "output": ("#f1f5f9", "#334155", "#0f172a"),     # final outputs
}


def card(
    x: int,
    y: int,
    w: int,
    h: int,
    title: str,
    subtitle: str,
    tool: str,
    style: str,
    dashed: bool = False,
) -> str:
    """Render a crisp rounded card with left accent bar and optional tool badge."""
    bg, accent, fg = COLORS[style]
    dash_attr = ' stroke-dasharray="5,3"' if dashed else ""
    cx = x + w / 2 + 3
    lines = [
        f'<g filter="url(#shadow)">',
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="8" fill="{bg}" '
        f'stroke="{accent}" stroke-width="1.5"{dash_attr}/>',
        f'<rect x="{x}" y="{y}" width="6" height="{h}" rx="3" fill="{accent}"/>',
    ]
    if tool:
        lines.append(
            f'<text x="{cx}" y="{y + 21}" text-anchor="middle" font-size="13" font-weight="600" fill="{fg}">{escape(title)}</text>'
            f'<text x="{cx}" y="{y + 37}" text-anchor="middle" font-size="11" fill="#475569">{escape(subtitle)}</text>'
            f'<text x="{cx}" y="{y + 53}" text-anchor="middle" font-size="10" font-weight="600" fill="{accent}">{escape(tool)}</text>'
        )
    elif subtitle:
        lines.append(
            f'<text x="{cx}" y="{y + h / 2 - 4}" text-anchor="middle" font-size="13" font-weight="600" fill="{fg}">{escape(title)}</text>'
            f'<text x="{cx}" y="{y + h / 2 + 13}" text-anchor="middle" font-size="11" fill="#475569">{escape(subtitle)}</text>'
        )
    else:
        lines.append(
            f'<text x="{cx}" y="{y + h / 2 + 4}" text-anchor="middle" font-size="13" font-weight="600" fill="{fg}">{escape(title)}</text>'
        )
    lines.append("</g>")
    return "\n".join(lines)


def arrow(points: list[tuple[float, float]], dashed: bool = False) -> str:
    """Render a crisp orthogonal connector with arrowhead."""
    d = f"M {points[0][0]},{points[0][1]}"
    for px, py in points[1:]:
        d += f" L {px},{py}"
    dash_attr = ' stroke-dasharray="4,3"' if dashed else ""
    return (
        f'<path d="{d}" fill="none" stroke="#64748b" stroke-width="1.8" '
        f'stroke-linejoin="round" stroke-linecap="round"{dash_attr} marker-end="url(#ah)"/>'
    )


def build_svg() -> str:
    """Assemble the complete Phase 2 SVG diagram."""
    w, h = 1040, 670
    parts: list[str] = [
        f'<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 {w} {h}" width="{w}" height="{h}" '
        f'font-family="Inter, -apple-system, BlinkMacSystemFont, Segoe UI, Helvetica, Arial, sans-serif">',
        "<defs>",
        '  <marker id="ah" markerWidth="8" markerHeight="6" refX="7" refY="3" orient="auto">',
        '    <polygon points="0 0, 8 3, 0 6" fill="#64748b"/>',
        "  </marker>",
        '  <filter id="shadow" x="-5%" y="-5%" width="110%" height="115%">',
        '    <feDropShadow dx="0" dy="1.5" stdDeviation="2" flood-color="#0f172a" flood-opacity="0.08"/>',
        "  </filter>",
        "</defs>",
        f'<rect x="0" y="0" width="{w}" height="{h}" rx="12" fill="#ffffff" stroke="#e2e8f0" stroke-width="1.5"/>',
        # Header
        '<text x="36" y="40" font-size="18" font-weight="700" fill="#0f172a">'
        "Phase 2 — Viral Assembly, Three-Tool Taxonomy &amp; Functional Annotation</text>",
        '<text x="36" y="60" font-size="12" fill="#64748b">'
        "workflow/phase2/Snakefile_phase2 · standalone from CSV manifest (SRA or local FASTQ)</text>",
        '<line x1="36" y1="74" x2="1004" y2="74" stroke="#e2e8f0" stroke-width="1"/>',
    ]

    # Row 1: Stage 1 Intake & Assembly (4 cards horizontally)
    parts.append('<text x="36" y="96" font-size="11" font-weight="700" fill="#64748b">1 · READ INTAKE, QC &amp; ASSEMBLY (PER SAMPLE)</text>')
    parts.append(card(36, 106, 215, 66, "Stage reads", "SRA fetch or local FASTQ", "sra-tools 3.1.1", "input"))
    parts.append(card(287, 106, 215, 66, "fastp_cleaning", "adapter & quality trim", "fastp 1.4.0", "qc"))
    parts.append(card(538, 106, 215, 66, "rnaviral_assembly", "de novo viral assembly", "rnaviralSPAdes · bwa-mem3", "asm"))
    parts.append(card(789, 106, 215, 66, "filter_contigs", "drop zero-coverage contigs", "samtools · seqkit", "asm"))

    parts.append(arrow([(251, 139), (285, 139)]))
    parts.append(arrow([(502, 139), (536, 139)]))
    parts.append(arrow([(753, 139), (787, 139)]))

    # Orthogonal fan-out bus from filter_contigs to Row 2
    parts.append('<text x="36" y="196" font-size="11" font-weight="700" fill="#64748b">2 · PARALLEL VIRAL SCREENING &amp; SCORING</text>')
    parts.append(arrow([(896, 172), (896, 206), (143, 206), (143, 224)]))
    parts.append(arrow([(896, 206), (394, 206), (394, 224)]))
    parts.append(arrow([(896, 206), (645, 206), (645, 224)]))
    parts.append(arrow([(896, 206), (896, 224)], dashed=True))

    # Row 2: Parallel Viral Screening (3 screening tools + 1 annotation-only)
    parts.append(card(36, 226, 215, 66, "run_kraken2", "k-mer known viruses", "Kraken2 2.17.2", "tax"))
    parts.append(card(287, 226, 215, 66, "run_genomad", "markers + NN, ICTV tax", "geNomad 1.12.0 (DB v1.9)", "tax"))
    parts.append(card(538, 226, 215, 66, "run_rdrpcatch", "RdRp HMM profiles", "RdRpCATCH 0.0.9", "tax"))
    parts.append(card(789, 226, 215, 66, "run_virsorter2", "viral score (annotation only)", "VirSorter2 2.2.4", "vs2", dashed=True))

    # Orthogonal merge from 3 screening tools into collect_viral_positives
    parts.append(arrow([(143, 292), (143, 314), (394, 314), (394, 342)]))
    parts.append(arrow([(394, 292), (394, 342)]))
    parts.append(arrow([(645, 292), (645, 314), (394, 314), (394, 342)]))

    # Row 3: Consensus Merge & Quality Annotation
    parts.append('<text x="36" y="334" font-size="11" font-weight="700" fill="#64748b">3 · CONSENSUS UNION &amp; QUALITY ANNOTATION</text>')
    parts.append(card(236, 344, 316, 66, "collect_viral_positives", "union: any screening tool positive", "tool_calls.tsv + viral_positives.fasta", "merge"))
    parts.append(card(36, 344, 172, 66, "run_checkv", "quality (annot. only)", "CheckV 1.1.1", "tax", dashed=True))
    parts.append(card(586, 344, 260, 66, "annotate_viral_positives", "+ CheckV & VirSorter2 cols", "taxonomy_merged.tsv", "merge"))

    # Clean connectors in Row 3 (separated incoming vs outgoing x-coords)
    parts.append(arrow([(236, 377), (210, 377)]))
    parts.append(arrow([(122, 410), (122, 430), (676, 430), (676, 412)], dashed=True))
    parts.append(arrow([(552, 377), (584, 377)]))
    parts.append(arrow([(896, 292), (896, 377), (848, 377)], dashed=True))

    # Row 4: Stage 3 Cohort Aggregation & Functional Annotation
    parts.append('<text x="36" y="454" font-size="11" font-weight="700" fill="#64748b">4 · COHORT AGGREGATION &amp; FUNCTIONAL ANNOTATION (ALL SAMPLES)</text>')
    parts.append(arrow([(756, 410), (756, 464), (143, 464), (143, 482)]))
    parts.append(card(36, 484, 215, 66, "build_master_database", "deduplicated cohort DB", "master_metadata.csv", "merge"))
    parts.append(card(287, 484, 215, 66, "select_contigs", "positives ≥ 1 kb", "filtered FASTA", "annot"))
    parts.append(card(538, 484, 215, 66, "run_metacerberus", "KOFam / VOG / PHROG", "MetaCerberus 1.4.0", "annot"))
    parts.append(card(789, 484, 215, 66, "Phase 2 Outputs", "master_viral_database.fasta", "results/annotation/metacerberus/", "output"))

    parts.append(arrow([(251, 517), (285, 517)]))
    parts.append(arrow([(502, 517), (536, 517)]))
    parts.append(arrow([(753, 517), (787, 517)]))

    # Bottom Legend
    legend_items = [
        ("input", "01.input.smk"),
        ("qc", "02.fastp.smk"),
        ("asm", "03.assembly.smk"),
        ("tax", "04.taxonomy.smk"),
        ("merge", "05.merge_tax.smk"),
        ("annot", "06.annotation.smk"),
        ("vs2", "dashed = annotation only"),
    ]
    lx = 36
    ly = 620
    for key, label in legend_items:
        bg, accent, _ = COLORS[key]
        dash_attr = ' stroke-dasharray="3,2"' if key == "vs2" else ""
        parts.append(
            f'<rect x="{lx}" y="{ly - 10}" width="12" height="12" rx="3" fill="{bg}" stroke="{accent}" stroke-width="1.5"{dash_attr}/>'
            f'<text x="{lx + 18}" y="{ly}" font-size="11" fill="#475569">{escape(label)}</text>'
        )
        lx += 140 if key != "vs2" else 160

    parts.append("</svg>")
    return "\n".join(parts) + "\n"


def main() -> int:
    OUT_SVG.write_text(build_svg(), encoding="utf-8")
    print(f"Wrote {OUT_SVG}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
