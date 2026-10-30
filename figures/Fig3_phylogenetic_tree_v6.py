#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fig. 3  Maximum-likelihood tree of the functional candidates and the three
        copies discussed as Cases 1-3.

v6 (2026-10-23): built from the tree produced by 07_cross_haplotype_tree.sh,
which uses the same settings as the per-haplotype trees (MAFFT L-INS-i,
IQ-TREE -m MFP -B 1000 -seed 12345, rooted on BBE). The previous version was
drawn from a separate 2026-10-08 run that used -m TEST and an unfixed seed and
so could not be reproduced; it also lacked SAN2a g539950 and SAN2b g313950.

The outgroup branch is far longer than any ingroup branch, so it is drawn
shortened, with a break mark and its true length given.

Input : Fig3_cross_haplotype.contree  (same directory, or argv[1])
Output: Fig3_phylogenetic_tree_v6.{pdf,png}

    python3 Fig3_phylogenetic_tree_v6.py [tree.contree]
"""
import os
import re
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import matplotlib.lines as mlines

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "DejaVu Sans"

MM = 1 / 25.4
FIG_W, FIG_H = 180 * MM, 112 * MM

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.join(
    HERE, "Fig3_cross_haplotype.contree")
OUTBASE = os.path.join(HERE, "Fig3_phylogenetic_tree_v6")

OUTGROUP = "sp|P30986|BBE_ESCCA"
BBE_DRAW = 0.030          # substitutions/site drawn for the outgroup branch
MIN_SUPPORT = 0           # show every support value

# series of each tip, for colour
SERIES = {
    "SODLa.chr7.v1.g255560.t1": "THCAS",
    "SODLb.chr7.v1.g250920.t1": "THCAS",
    "EH23a.chr7.v1.g200100.t1": "THCAS",
    "SAN2a.chr7.v1.g547390.t1": "THCAS",
    "SAN2b.chr7.v1.g313950.t1": "THCAS",
    "sp|Q8GTB6.1|THCAS_CANSA": "THCAS",
    "EH23b.chr7.v1.g195030.t1": "CBDAS",
    "KOMPa.chr7.v1.g259030.t1": "CBDAS",
    "KOMPb.chr7.v1.g226860.t1": "CBDAS",
    "sp|A6P6V9|CBDAS_CANSA": "CBDAS",
    "SAN2a.chr7.v1.g539950.t1": "CBCAS",
    "sp|Q33DQ2|CBCAS_CANSA": "CBCAS",
    OUTGROUP: "outgroup",
}
COLOUR = {"THCAS": "#B71C1C", "CBDAS": "#1565C0",
          "CBCAS": "#6A1B9A", "outgroup": "#555555"}

# classification, for the tip marker
TIER = {
    "SODLa.chr7.v1.g255560.t1": "functional",
    "SODLb.chr7.v1.g250920.t1": "functional",
    "EH23a.chr7.v1.g200100.t1": "functional",
    "EH23b.chr7.v1.g195030.t1": "functional",
    "KOMPa.chr7.v1.g259030.t1": "functional",
    "KOMPb.chr7.v1.g226860.t1": "functional",
    "SAN2a.chr7.v1.g539950.t1": "divergent",
    "SAN2a.chr7.v1.g547390.t1": "divergent",
    "SAN2b.chr7.v1.g313950.t1": "pseudo",
}
MARKER = {"functional": ("o", True), "divergent": ("s", False),
          "pseudo": ("X", True)}

# display name, and the haplotype prefix shown in grey
DISPLAY = {
    "sp|Q8GTB6.1|THCAS_CANSA": ("THCAS", "reference (Q8GTB6)"),
    "sp|A6P6V9|CBDAS_CANSA": ("CBDAS", "reference (A6P6V9)"),
    "sp|Q33DQ2|CBCAS_CANSA": ("CBCAS", "reference (Q33DQ2)"),
    OUTGROUP: ("BBE", "outgroup, E. californica (P30986)"),
}
CASE = {"SAN2a.chr7.v1.g539950.t1": "Case 1",
        "SAN2a.chr7.v1.g547390.t1": "Case 2",
        "SAN2b.chr7.v1.g313950.t1": "Case 3"}


class Node:
    __slots__ = ("name", "length", "support", "children", "x", "y")

    def __init__(self):
        self.name = ""
        self.length = 0.0
        self.support = None
        self.children = []
        self.x = self.y = 0.0

    @property
    def is_tip(self):
        return not self.children


def parse_newick(text):
    """Minimal newick parser; tip names may contain '|' and '.'."""
    s = text.strip().rstrip(";")
    pos = 0

    def read_node():
        nonlocal pos
        node = Node()
        if s[pos] == "(":
            pos += 1
            while True:
                node.children.append(read_node())
                if s[pos] == ",":
                    pos += 1
                    continue
                if s[pos] == ")":
                    pos += 1
                    break
                raise ValueError(f"unexpected {s[pos]!r} at {pos}")
        # label (support for internal nodes, name for tips)
        start = pos
        while pos < len(s) and s[pos] not in ",():":
            pos += 1
        label = s[start:pos]
        if node.children:
            node.support = int(label) if label.isdigit() else None
        else:
            node.name = label
        if pos < len(s) and s[pos] == ":":
            pos += 1
            start = pos
            while pos < len(s) and s[pos] not in ",()":
                pos += 1
            node.length = float(s[start:pos])
        return node

    root = read_node()
    return root


def tips_in_order(node, out=None):
    if out is None:
        out = []
    if node.is_tip:
        out.append(node)
    for c in node.children:
        tips_in_order(c, out)
    return out


def reorder_outgroup_last(root):
    """Draw the outgroup at the bottom of the figure.

    Tip y values are assigned in traversal order from 0 upwards and y = 0 is
    the bottom of the axes, so the outgroup has to come first in the child
    order to end up at the bottom.
    """
    def has_og(n):
        return n.name == OUTGROUP or any(has_og(c) for c in n.children)
    root.children.sort(key=lambda n: not has_og(n))
    return root


def layout(root):
    """Assign x (cumulative length) and y (tip order)."""
    counter = [0]

    def walk(node, x0):
        # the outgroup branch is drawn shortened
        blen = BBE_DRAW if node.name == OUTGROUP else node.length
        node.x = x0 + blen
        if node.is_tip:
            node.y = counter[0]
            counter[0] += 1
        else:
            for c in node.children:
                walk(c, node.x)
            node.y = (node.children[0].y + node.children[-1].y) / 2

    root.x = 0.0
    for c in root.children:
        walk(c, 0.0)
    root.y = (root.children[0].y + root.children[-1].y) / 2
    return root


def draw(ax, node):
    for c in node.children:
        colour = "#555555"
        ax.plot([node.x, node.x], [c.y, node.y], color=colour, lw=1.0,
                solid_capstyle="round", zorder=2)
        ax.plot([node.x, c.x], [c.y, c.y], color=colour, lw=1.0,
                solid_capstyle="round", zorder=2)
        if c.name == OUTGROUP:
            # break mark on the shortened outgroup branch
            xm = (node.x + c.x) / 2
            for dx in (-0.0016, 0.0016):
                ax.plot([xm + dx - 0.0013, xm + dx + 0.0013],
                        [c.y - 0.17, c.y + 0.17], color="white", lw=2.2,
                        zorder=3)
                ax.plot([xm + dx - 0.0013, xm + dx + 0.0013],
                        [c.y - 0.17, c.y + 0.17], color="#555555", lw=0.9,
                        zorder=4)
        if c.support is not None and c.support >= MIN_SUPPORT:
            ax.text(c.x - 0.0012, c.y + 0.17, str(c.support), ha="right",
                    va="bottom", fontsize=7.0, color="#333333", zorder=5)
        draw(ax, c)


def main():
    if not os.path.exists(DATA):
        sys.exit(f"tree not found: {DATA}")
    with open(DATA) as fh:
        root = parse_newick(fh.read())
    reorder_outgroup_last(root)
    layout(root)
    tips = tips_in_order(root)

    unknown = [t.name for t in tips if t.name not in SERIES]
    if unknown:
        print("[warn] tips not in SERIES:", unknown, file=sys.stderr)

    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    fig.patch.set_facecolor("white")
    draw(ax, root)

    xmax = max(t.x for t in tips)
    gap = xmax * 0.022

    # A generous x range to start with: label positions are measured in pixels
    # under this range and stored in data coordinates, and the range is only
    # narrowed afterwards, which can widen the gaps but never close them.
    ax.set_xlim(-0.004, xmax * 1.95)
    ax.set_ylim(-2.10, len(tips) - 0.4)
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()

    def right_edge(artist, pad_px=7.0):
        """Data x just right of a drawn text, so the next item cannot overlap."""
        bb = artist.get_window_extent(renderer=renderer)
        return ax.transData.inverted().transform((bb.x1 + pad_px, 0))[0]

    right_most = xmax
    case_pending = []          # badges are drawn afterwards, in one column
    for t in tips:
        ser = SERIES.get(t.name, "outgroup")
        col = COLOUR[ser]
        tier = TIER.get(t.name)
        if tier:
            m, filled = MARKER[tier]
            ax.plot(t.x, t.y, marker=m, markersize=4.6,
                    markerfacecolor=col if filled else "white",
                    markeredgecolor=col, markeredgewidth=1.0, zorder=6)

        if t.name in DISPLAY:
            main_txt, sub = DISPLAY[t.name]
            a = ax.text(t.x + gap, t.y, main_txt, ha="left", va="center",
                        fontsize=7.0, color=col, fontweight="bold", zorder=6)
            b = ax.text(right_edge(a), t.y, sub, ha="left", va="center",
                        fontsize=7.0, color="#777777", style="italic",
                        zorder=6)
            right_most = max(right_most, right_edge(b, 0.0))
        else:
            hap, gene = t.name.split(".")[0], t.name.split(".")[-2]
            a = ax.text(t.x + gap, t.y, f"{hap}  {gene}", ha="left",
                        va="center", fontsize=7.0, color=col, zorder=6)
            right_most = max(right_most, right_edge(a, 0.0))
            if t.name in CASE:
                case_pending.append((right_edge(a, 9.0), t.y, CASE[t.name]))

    # all Case badges share one column, so they read as a set
    if case_pending:
        col = max(x for x, _, _ in case_pending)
        for _, y, text in case_pending:
            b = ax.text(col, y, text, ha="left", va="center", fontsize=7.0,
                        color="#333333", fontweight="bold", zorder=6,
                        bbox=dict(boxstyle="round,pad=0.20",
                                  facecolor="#F0F0F0", edgecolor="#AAAAAA",
                                  linewidth=0.6))
            right_most = max(right_most, right_edge(b, 2.0))

    # scale bar
    bar = 0.01
    y0 = -1.35
    ax.plot([0.004, 0.004 + bar], [y0, y0], color="#333333", lw=1.2)
    for x in (0.004, 0.004 + bar):
        ax.plot([x, x], [y0 - 0.12, y0 + 0.12], color="#333333", lw=1.2)
    ax.text(0.004 + bar / 2, y0 - 0.30, f"{bar} substitutions per site",
            ha="center", va="top", fontsize=7.0, color="#333333")
    ax.text(0.004 + bar * 2.4, y0 - 0.02,
            "outgroup branch shortened (true length 1.187)",
            ha="left", va="center", fontsize=7.0, color="#777777",
            style="italic")

    handles = [
        mlines.Line2D([], [], color="#555555", marker="o", linestyle="none",
                      markersize=4.6, markerfacecolor="#555555",
                      label="functional candidate"),
        mlines.Line2D([], [], color="#555555", marker="s", linestyle="none",
                      markersize=4.6, markerfacecolor="white",
                      label="divergent copy"),
        mlines.Line2D([], [], color="#555555", marker="X", linestyle="none",
                      markersize=4.6, markerfacecolor="#555555",
                      label="pseudogenized copy"),
        mpatches.Patch(facecolor=COLOUR["THCAS"], label="THCAS series"),
        mpatches.Patch(facecolor=COLOUR["CBDAS"], label="CBDAS series"),
        mpatches.Patch(facecolor=COLOUR["CBCAS"], label="CBCAS series"),
    ]
    ax.legend(handles=handles, loc="lower right", fontsize=7.0, ncol=2,
              framealpha=0.95, handlelength=1.4, borderpad=0.6,
              bbox_to_anchor=(1.0, -0.02))

    ax.set_xlim(-0.004, right_most + xmax * 0.03)
    ax.set_ylim(y0 - 0.75, len(tips) - 0.4)
    ax.axis("off")
    fig.tight_layout()

    fig.savefig(OUTBASE + ".pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(OUTBASE + ".png", dpi=400, bbox_inches="tight",
                facecolor="white")
    print("written:", OUTBASE + ".pdf", "/", OUTBASE + ".png")


if __name__ == "__main__":
    main()
