#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fig. 1  Cannabinoid acid biosynthesis and the materials analysed in this study.

v11 (2026-10-22): the chemotype panels now show only the materials analysed
(cultivar and haplotypes). Classification results were removed from this figure
-- they are given in Fig. 2, Table 4 and Additional file 1. The CBCAS
copy-number note (stale counts) was also removed and the layout tightened.

Output: 08_results/figures/Fig1_cannabinoid_pathway_v12.{pdf,png}

    python3 Fig1_pathway_v11.py
"""
import os
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "DejaVu Sans"

MM = 1 / 25.4
FIG_W, FIG_H = 180 * MM, 92 * MM
OUTDIR = os.path.dirname(os.path.abspath(__file__))
OUTBASE = os.path.join(OUTDIR, "Fig1_cannabinoid_pathway_v12")

C = {
    "cbga":   "#2E7D32",
    "enzyme": "#F5A623", "enzyme_e": "#D98806",
    "thca":   "#C62828", "cbca": "#7B1FA2", "cbda": "#1565C0",
    "p1":     "#FADCE2", "p1e": "#C62828",
    "p2":     "#FDF3C8", "p2e": "#D4A017",
    "p3":     "#D6E8F7", "p3e": "#1565C0",
    "p4":     "#FBE4BF", "p4e": "#E07B00",
    "p5":     "#F7F7F7", "p5e": "#AAAAAA",
    "ink":    "#1A2B4C",
}


def box(ax, x, y, w, h, fc, ec, lw=1.2, r=0.6, ls="solid", z=2):
    p = FancyBboxPatch((x, y), w, h,
                       boxstyle=f"round,pad=0,rounding_size={r}",
                       facecolor=fc, edgecolor=ec, linewidth=lw,
                       linestyle=ls, zorder=z)
    ax.add_patch(p)
    return p


def arrow(ax, x0, y0, x1, y1, color, ls="solid", lw=1.3, z=3,
          rad=0.0, head=7):
    ax.add_patch(FancyArrowPatch(
        (x0, y0), (x1, y1), arrowstyle=f"-|>,head_width=2.6,head_length={head*0.55}",
        connectionstyle=f"arc3,rad={rad}", color=color, linewidth=lw,
        linestyle=ls, mutation_scale=1.0, zorder=z,
        shrinkA=0, shrinkB=0))


def txt(ax, x, y, s, size=6.4, color="#1A2B4C", weight="normal",
        ha="center", va="center", style="normal", ls=1.35, z=5):
    ax.text(x, y, s, fontsize=size, color=color, fontweight=weight,
            ha=ha, va=va, style=style, linespacing=ls, zorder=z)


def main():
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))

    # --- common precursor ---------------------------------------------------
    box(ax, 3, 87.5, 76, 6.5, C["cbga"], "#1B5E20")
    txt(ax, 41, 90.7, "CBGA  (cannabigerolic acid)   —   common precursor of "
        "all cannabinoid acids", size=7.4, weight="bold", color="white")

    # --- enzymes ------------------------------------------------------------
    enz = [(14, "THCAS", C["thca"]),
           (41, "CBCAS", C["cbca"]),
           (68, "CBDAS", C["cbda"])]
    for cx, name, acc in enz:
        arrow(ax, cx, 87.3, cx, 80.6, acc)
        box(ax, cx - 9, 73.5, 18, 7, C["enzyme"], C["enzyme_e"])
        txt(ax, cx, 77.0, name, size=7.4, weight="bold", color="#7A4A00")

    # --- products -----------------------------------------------------------
    prod = [(14, "THCA", "tetrahydrocannabinolic acid", C["thca"]),
            (41, "CBCA", "cannabichromenic acid", C["cbca"]),
            (68, "CBDA", "cannabidiolic acid", C["cbda"])]
    for cx, name, full, acc in prod:
        arrow(ax, cx, 73.3, cx, 67.0, acc)
        box(ax, cx - 11, 59.5, 22, 7.2, acc, acc)
        txt(ax, cx, 64.6, name, size=7.4, weight="bold", color="white")
        txt(ax, cx, 61.6, full, size=5.5, color="white")

    # --- chemotype panels: materials analysed (no classification results) ---
    PY, PH = 28.0, 23.0          # panel bottom / height
    panels = [
        (1.5, 21.0, C["p1"], C["p1e"], "Chemotype I", "(high-THCA)",
         "Sour Diesel", "SODLa Hap1 / SODLb Hap2"),
        (23.5, 21.0, C["p2"], C["p2e"], "Chemotype II", "(mixed)",
         "ERBxHO40_23\n(HO40 × ERB F1)", "EH23a Hap1 / EH23b Hap2"),
        (45.5, 21.0, C["p3"], C["p3e"], "Chemotype III", "(high-CBDA)",
         "Kompolti", "KOMPa Hap1 / KOMPb Hap2"),
        (67.5, 21.0, C["p4"], C["p4e"], "Chemotype IV", "(high-CBGA)",
         "Santhica 27", "SAN2a Hap1 / SAN2b Hap2"),
    ]
    for x, w, fc, ec, t1, t2, body, foot in panels:
        box(ax, x, PY, w, PH, fc, ec, lw=1.1, r=0.5)
        cx = x + w / 2
        txt(ax, cx, PY + 19.6, t1, size=7.0, weight="bold", color=ec)
        txt(ax, cx, PY + 16.8, t2, size=5.8, color=ec)
        ax.plot([cx - 6, cx + 6], [PY + 14.6, PY + 14.6], color=ec, lw=0.6, zorder=4)
        txt(ax, cx, PY + 9.8, body, size=5.6, color="#333333", ls=1.65)
        ax.plot([cx - 6, cx + 6], [PY + 5.4, PY + 5.4], color=ec, lw=0.6, zorder=4)
        txt(ax, cx, PY + 2.8, foot, size=5.4, weight="bold", color=ec, ls=1.5)

    PT = PY + PH                 # panel top (= 51.0)

    # --- links from products to chemotypes ---------------------------------
    arrow(ax, 12.0, 59.3, 12.0, PT + 0.3, C["p1e"], ls=(0, (3, 2)), lw=1.0)
    arrow(ax, 65.0, 59.3, 56.0, PT + 0.3, C["p3e"], ls=(0, (3, 2)), lw=1.0)
    # chemotype II receives both THCA and CBDA routes
    arrow(ax, 20.0, 59.3, 30.0, PT + 0.3, C["p1e"], ls=(0, (3, 2)), lw=0.9)
    arrow(ax, 57.5, 59.3, 38.0, PT + 0.3, C["p3e"], ls=(0, (3, 2)), lw=0.9)
    # chemotype IV: CBGA accumulates, no conversion
    arrow(ax, 79.2, 90.7, 84.0, 90.7, C["p4e"], ls=(0, (3, 2)), lw=1.0,
          head=0)
    ax.plot([84, 84], [90.7, PT + 1.5], color=C["p4e"], lw=1.0,
            ls=(0, (3, 2)), zorder=3)
    arrow(ax, 84, PT + 2.0, 84, PT + 0.3, C["p4e"], ls="solid", lw=1.0)
    txt(ax, 84, PT + 4.6, "CBGA accumulates", size=5.2, weight="bold",
        color=C["p4e"])

    # --- chemotype V (not analysed) ----------------------------------------
    box(ax, 87.5, 59.0, 12.0, 32.0, C["p5"], C["p5e"], lw=0.9, r=0.4,
        ls=(0, (3, 2)))
    txt(ax, 93.5, 87.5, "Chemotype V", size=6.2, weight="bold",
        color="#666666")
    txt(ax, 93.5, 85.0, "(cannabinoid-free)", size=5.2, color="#666666")
    ax.plot([89.5, 97.5], [82.6, 82.6], color=C["p5e"], lw=0.6, zorder=4)
    txt(ax, 93.5, 77.5, "candidates\n8 accessions (CLR)\nLynch et al. 2025",
        size=5.1, color="#666666", ls=1.6)
    ax.plot([89.5, 97.5], [72.6, 72.6], color=C["p5e"], lw=0.6, zorder=4)
    txt(ax, 93.5, 69.8, "not analysed", size=5.6, weight="bold",
        color="#777777")
    txt(ax, 93.5, 64.0, "contig-level assembly;\nno haplotype phasing\n"
        "(contig N50 4.0–12.5 Mb)", size=4.8, color="#888888", style="italic", ls=1.6)

    ax.set_xlim(0, 100)
    ax.set_ylim(25.0, 95.5)
    ax.set_aspect("auto")
    ax.axis("off")
    fig.tight_layout()

    os.makedirs(OUTDIR, exist_ok=True)
    fig.savefig(OUTBASE + ".pdf", bbox_inches="tight")
    fig.savefig(OUTBASE + ".png", dpi=400, bbox_inches="tight")
    print("written:", OUTBASE + ".pdf", "/", OUTBASE + ".png")


if __name__ == "__main__":
    main()
