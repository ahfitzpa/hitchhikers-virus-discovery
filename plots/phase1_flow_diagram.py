#!/usr/bin/env python3
"""Generate a clean, publication-grade SVG workflow diagram for Phase 1.

Run:
    python3 plots/phase1_flow_diagram.py
Writes:
    plots/phase1_flow_diagram.svg
"""
from __future__ import annotations

from html import escape
from pathlib import Path

OUT_SVG = Path(__file__).resolve().with_suffix(".svg")

# Modern high-contrast palette (light surface cards with crisp colored accent bars)
COLORS = {
    "base": ("#e0f2fe", "#0284c7", "#0c4a6e"),       # 01.base_discovery.smk
    "unanch": ("#fef3c7", "#d97706", "#78350f"),     # 01.centroid_unanchored.smk
    "anch": ("#ede9fe", "#6366f1", "#312e81"),       # 01.centroid_anchored.smk
    "shared": ("#ffedd5", "#ea580c", "#7c2d12"),     # shared matrix / ID extraction
    "meta": ("#ccfbf1", "#0d9488", "#134e4a"),       # 03.metadata_fetch.smk
    "merge": ("#d1fae5", "#059669", "#064e3b"),      # 04.metadata_merge.smk
    "output": ("#f1f5f9", "#334155", "#0f172a"),     # final outputs
}


def card(
    x: int,
    y: int,
    w: int,
    h: int,
    title: str,
    subtitle: str,
    style: str,
    dashed: bool = False,
) -> str:
    """Render a rounded stage card with a crisp left accent bar."""
    bg, accent, fg = COLORS[style]
    dash_attr = ' stroke-dasharray="5,3"' if dashed else ""
    cx = x + w / 2 + 3
    if subtitle:
        text_svg = (
            f'<text x="{cx}" y="{y + h / 2 - 4}" text-anchor="middle" '
            f'font-size="13" font-weight="600" fill="{fg}">{escape(title)}</text>'
            f'<text x="{cx}" y="{y + h / 2 + 13}" text-anchor="middle" '
            f'font-size="11" fill="#475569">{escape(subtitle)}</text>'
        )
    else:
        text_svg = (
            f'<text x="{cx}" y="{y + h / 2 + 4}" text-anchor="middle" '
            f'font-size="13" font-weight="600" fill="{fg}">{escape(title)}</text>'
        )
    return (
        f'<g filter="url(#shadow)">'
        f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="8" fill="{bg}" '
        f'stroke="{accent}" stroke-width="1.5"{dash_attr}/>'
        f'<rect x="{x}" y="{y}" width="6" height="{h}" rx="3" fill="{accent}"/>'
        f"{text_svg}"
        f"</g>"
    )


def arrow(points: list[tuple[float, float]], dashed: bool = False) -> str:
    """Render a crisp orthogonal connector with rounded corners and arrowhead."""
    d = f"M {points[0][0]},{points[0][1]}"
    for px, py in points[1:]:
        d += f" L {px},{py}"
    dash_attr = ' stroke-dasharray="4,3"' if dashed else ""
    return (
        f'<path d="{d}" fill="none" stroke="#64748b" stroke-width="1.8" '
        f'stroke-linejoin="round" stroke-linecap="round"{dash_attr} marker-end="url(#ah)"/>'
    )


def build_svg() -> str:
    """Assemble the complete Phase 1 SVG diagram."""
    w, h = 1040, 690
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
        "Phase 1 — Serratus / Logan SRA Discovery Scenarios</text>",
        '<text x="36" y="60" font-size="12" fill="#64748b">'
        "workflow/phase1/Snakefile_phase1 · discovery_type = novel vs. known (surrogate, host_expansion, geographic)</text>",
        '<line x1="36" y1="74" x2="1004" y2="74" stroke="#e2e8f0" stroke-width="1"/>',
    ]

    # Top checkpoint
    parts.append(card(360, 92, 320, 52, "centroid_step1_query", "checkpoint · Serratus RDS + NCBI anchor", "base"))

    # Branch labels & containers
    parts.append(
        '<rect x="36" y="176" width="276" height="190" rx="10" fill="#fffbeb" '
        'stroke="#fde68a" stroke-width="1.2" stroke-dasharray="5,3"/>'
    )
    parts.append(
        '<text x="52" y="198" font-size="11.5" font-weight="700" fill="#b45309">discovery_type = novel</text>'
    )

    parts.append(
        '<rect x="348" y="176" width="656" height="190" rx="10" fill="#f5f3ff" '
        'stroke="#ddd6fe" stroke-width="1.2" stroke-dasharray="5,3"/>'
    )
    parts.append(
        '<text x="986" y="198" text-anchor="end" font-size="11.5" font-weight="700" fill="#4338ca">'
        "discovery_type = known (anchored MMseqs2)</text>"
    )

    # Arrows from top checkpoint to both branches (enter cleanly to the right of branch labels)
    parts.append(arrow([(450, 144), (450, 160), (220, 160), (220, 224)]))
    parts.append(arrow([(590, 144), (590, 160), (508, 160), (508, 212)]))

    # Left branch (novel)
    parts.append(card(64, 226, 220, 56, "transform_db_to_hits", "unanchored sOTU hits", "unanch"))

    # Right branch (known) - clean 2-row flow inside container
    parts.append(card(366, 214, 184, 46, "combine_fasta", "Logan unitigs", "anch"))
    parts.append(card(578, 214, 196, 46, "centroid_dedup_and_split", "16 balanced chunks", "anch"))
    parts.append(card(802, 214, 184, 46, "centroid_cluster_chunk", "×16 parallel", "anch"))

    parts.append(card(802, 294, 184, 48, "centroid_merge_querydb", "merged query DB", "anch"))
    parts.append(card(578, 294, 196, 48, "mmseqs_search", "vs. mmseqs_build_refdb", "anch"))
    parts.append(card(366, 294, 184, 48, "expand_cluster_hits", "best-hit fident filter", "anch"))

    # Connectors inside known branch
    parts.append(arrow([(550, 237), (576, 237)]))
    parts.append(arrow([(774, 237), (800, 237)]))
    parts.append(arrow([(894, 260), (894, 292)]))
    parts.append(arrow([(802, 318), (776, 318)]))
    parts.append(arrow([(578, 318), (552, 318)]))

    # Shared convergence row
    parts.append(card(140, 404, 260, 52, "build_matrix / build_matrix_novel", "sOTU × SRA matrix", "shared"))
    parts.append(card(448, 404, 210, 52, "extract_sra_ids", "unique run accessions", "shared"))
    parts.append(card(706, 404, 230, 52, "fetch_sra_metadata", "NCBI SRA XML / TSV", "meta"))

    # Connectors from novel & known to shared row
    parts.append(arrow([(174, 282), (174, 384), (240, 384), (240, 402)]))
    parts.append(arrow([(458, 342), (458, 384), (310, 384), (310, 402)]))
    parts.append(arrow([(458, 342), (458, 384), (553, 384), (553, 402)]))
    parts.append(arrow([(658, 430), (704, 430)]))

    # Merge + optional review gate row (zero box crossings!)
    parts.append(card(350, 492, 260, 50, "merge_sra_metadata", "checkpoint · shortlist CSV", "merge"))
    parts.append(card(706, 492, 230, 50, "automated_review_gate", "known scenarios only", "merge", dashed=True))

    parts.append(arrow([(270, 456), (270, 517), (348, 517)]))
    parts.append(arrow([(821, 456), (821, 474), (480, 474), (480, 490)]))
    parts.append(arrow([(610, 517), (704, 517)], dashed=True))
    parts.append(arrow([(480, 542), (480, 568)]))
    parts.append(arrow([(821, 542), (821, 595), (687, 595)], dashed=True))

    # Output box
    parts.append(
        card(
            275,
            570,
            410,
            52,
            "Phase 1 Outputs → Phase 2 Intake",
            "sra_ids.txt (novel)  ·  unique.csv + approved.flag (known)",
            "output",
        )
    )

    # Clean bottom legend (generous even spacing, zero overlap)
    legend_items = [
        ("base", "01.base_discovery.smk", 158),
        ("unanch", "01.centroid_unanchored.smk", 184),
        ("anch", "01.centroid_anchored.smk", 174),
        ("shared", "matrix & ID extraction", 156),
        ("meta", "03.metadata_fetch.smk", 156),
        ("merge", "04.metadata_merge.smk", 150),
    ]
    lx = 36
    ly = 656
    for key, label, step in legend_items:
        bg, accent, _ = COLORS[key]
        parts.append(
            f'<rect x="{lx}" y="{ly - 10}" width="12" height="12" rx="3" fill="{bg}" stroke="{accent}" stroke-width="1.5"/>'
            f'<text x="{lx + 18}" y="{ly}" font-size="11" fill="#475569">{escape(label)}</text>'
        )
        lx += step

    parts.append("</svg>")
    return "\n".join(parts) + "\n"


def main() -> int:
    OUT_SVG.write_text(build_svg(), encoding="utf-8")
    print(f"Wrote {OUT_SVG}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
