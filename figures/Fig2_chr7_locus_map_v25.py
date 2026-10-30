#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fig. 2  Chromosome 7 synthase locus map across chemotype I-IV haplotypes.

v25 (2026-10-27): hit regions with no overlapping gene model are excluded from
the copy counts, from the per-series functional counts and from the coordinate
range, and are reported on a line of their own. The range therefore covers the
classified copies only, matching the ranges quoted in the main text. Functional
candidates are counted per series rather than as one total, because a bare
"functional xN" in a mixed cluster reads as a claim about the wrong series.

v24 (2026-10-23): keeps the v21 layout and the v21 habit of naming the series
on the chromosome. Each cluster is labelled with the series names and copy
counts, written in the series colours, instead of a bare "n =" total. Content
comes from Additional file 1, the workflow output, so the figure, the main text
and Table 4 all come from the same source.

Copies are not drawn individually. On chromosome 7 the synthase copies fall into
one or two tight clusters per haplotype, each spanning less than 2 Mb, so
drawing all of them on a 65-70 Mb axis produces an unreadable band. Each cluster
is drawn as a single block, segmented by series in proportion to the number of
copies it contains.

Chromosome 7 is drawn to each haplotype's own length, measured from the
assemblies in 00_genomes (see CHR7_MB). Cluster positions are each haplotype's
own assembly coordinates and are not directly comparable between haplotypes.

Encoding
    block         one synthase cluster (copies separated by <= CLUSTER_GAP Mb)
    segments      series of the nearest reference, segment width proportional
                  to the number of copies
    dark outline  the cluster contains at least one functional candidate
    label         coordinate range, then one line per series in its own colour

Only chromosome 7 copies with >= MIN_IDENTITY identity to their nearest
reference are included (Methods section 7); the more distant BBE-like copies are
listed in Additional file 1.

Input : Additional_file_1_synthase_copy_detail.tsv
Output: Fig2_chr7_locus_map_v25.{pdf,png}

    python3 Fig2_chr7_locus_map_v25.py [Additional_file_1_synthase_copy_detail.tsv]
"""
import os
import sys
import csv
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from matplotlib.offsetbox import TextArea, VPacker, AnnotationBbox

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "DejaVu Sans"

MM = 1 / 25.4
FIG_W, FIG_H = 180 * MM, 262 * MM

HERE = os.path.dirname(os.path.abspath(__file__))
DEFAULT_DATA = os.path.join(HERE, "..", "..",
                            "Additional_file_1_synthase_copy_detail.tsv")
DATA = sys.argv[1] if len(sys.argv) > 1 else os.path.normpath(DEFAULT_DATA)
OUTBASE = os.path.join(HERE, "Fig2_chr7_locus_map_v25")

MIN_IDENTITY = 84.0      # Methods section 7: synthase paralogues
CLUSTER_GAP = 1.0        # Mb; copies further apart start a new cluster
BLOCK_W = 2.6            # Mb drawn per cluster, so that a 1-2 Mb cluster reads
BLOCK_H = 0.62

# chromosome 7 length per haplotype, measured from the assembly FASTA in
# 00_genomes/<set>/<hap>/<hap>.unmasked.fasta.gz  (2026-10-23)
CHR7_MB = {"SODLa": 66.7, "SODLb": 70.7, "EH23a": 65.5, "EH23b": 64.2,
           "KOMPa": 65.1, "KOMPb": 65.6, "SAN2a": 66.6, "SAN2b": 66.5}
XMIN, XMAX = 0, 71       # axis covers the longest chromosome 7

HAPS = ["SODLa", "SODLb", "EH23a", "EH23b",
        "KOMPa", "KOMPb", "SAN2a", "SAN2b"]

LABELS = {
    "SODLa": "SODLa\n(Type I Hap1)",
    "SODLb": "SODLb\n(Type I Hap2)",
    "EH23a": "EH23a\n(Type II Hap1, HO40-bg)",
    "EH23b": "EH23b\n(Type II Hap2, ERB-bg)",
    "KOMPa": "KOMPa\n(Type III Hap1)",
    "KOMPb": "KOMPb\n(Type III Hap2)",
    "SAN2a": "SAN2a\n(Type IV Hap1)",
    "SAN2b": "SAN2b\n(Type IV Hap2)",
}
LABEL_BG = {
    "SODLa": "#FFCDD2", "SODLb": "#FFCDD2",
    "EH23a": "#FFF9C4", "EH23b": "#FFF9C4",
    "KOMPa": "#BBDEFB", "KOMPb": "#BBDEFB",
    "SAN2a": "#FFE0B2", "SAN2b": "#FFE0B2",
}

SERIES = ["THCAS", "CBDAS", "CBCAS", "unassigned"]
SERIES_COLOUR = {"THCAS": "#B71C1C", "CBDAS": "#1565C0",
                 "CBCAS": "#6A1B9A", "unassigned": "#757575"}

LAB_W = 11.6             # Mb, label box width
LAB_LINE = 0.345         # y units per label line
LAB_BOTTOM = 0.52        # y of the label box bottom
GREY = "#555555"


def load(path):
    rows = defaultdict(list)
    skipped = 0
    with open(path, newline="") as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["Chromosome"] != "chr7":
                continue
            try:
                pid = float(r["Nearest_reference_identity_pct"])
                start, end = float(r["Start_Mb"]), float(r["End_Mb"])
            except ValueError:
                skipped += 1
                continue
            if pid < MIN_IDENTITY:
                skipped += 1
                continue
            ser = r["Series"] if r["Series"] in SERIES_COLOUR else "unassigned"
            rows[r["Haplotype"]].append(dict(
                pos=(start + end) / 2, series=ser,
                nomodel=(r["Gene_ID"] == "(no gene model)"),
                functional=r["Classification"].startswith("Functional")))
    for h in rows:
        rows[h].sort(key=lambda c: c["pos"])
    print(f"[Fig2] chromosome 7 copies at >= {MIN_IDENTITY}% identity: "
          f"{sum(len(v) for v in rows.values())}   excluded: {skipped}",
          file=sys.stderr)
    return rows


def cluster(copies):
    out = []
    for c in copies:
        if out and c["pos"] - out[-1][-1]["pos"] <= CLUSTER_GAP:
            out[-1].append(c)
        else:
            out.append([c])
    return out


def describe(group):
    """Label lines as (text, colour), plus the series counts for the block."""
    classified = [c for c in group if not c["nomodel"]]
    nomodel = len(group) - len(classified)
    counts = {s: sum(1 for c in classified if c["series"] == s)
              for s in SERIES}
    # functional candidates per series: a bare "functional xN" cannot say which
    # series they belong to, and in a mixed cluster that reads as a claim about
    # the wrong series.
    func = {s: sum(1 for c in classified
                   if c["series"] == s and c["functional"]) for s in SERIES}
    nfunc = sum(func.values())

    # The coordinate range covers the classified copies only, so that it matches
    # the ranges quoted in the text. Hit regions with no gene model are counted
    # on a line of their own and are not classified.
    span = classified or group
    lo = min(c["pos"] for c in span)
    hi = max(c["pos"] for c in span)

    lines = [(f"{lo:.1f}–{hi:.1f} Mb", "#222222")]
    for s in SERIES:
        if not counts[s]:
            continue
        text = f"{s} ×{counts[s]}"
        if func[s]:
            text += f"  ({func[s]} functional)"
        lines.append((text, SERIES_COLOUR[s]))
    if nomodel:
        lines.append((f"+{nomodel} hit{'s' if nomodel > 1 else ''}, "
                      f"no gene model", GREY))
    return lines, counts, nfunc, (lo + hi) / 2


def main():
    if not os.path.exists(DATA):
        sys.exit(f"data file not found: {DATA}")
    data = load(DATA)

    fig, axes = plt.subplots(len(HAPS), 1, figsize=(FIG_W, FIG_H))
    fig.patch.set_facecolor("white")

    for i, (hap, ax) in enumerate(zip(HAPS, axes)):
        ax.set_xlim(XMIN - 1, XMAX + 1)
        ax.set_ylim(-2.60, 2.60)
        ax.axis("off")

        ax.text(-2.2, 0, LABELS[hap], ha="right", va="center",
                fontsize=7.0, fontweight="bold", color="#222222",
                linespacing=1.4,
                bbox=dict(boxstyle="round,pad=0.30",
                          facecolor=LABEL_BG[hap],
                          edgecolor="#888888", linewidth=0.9))

        # chromosome backbone, drawn to this haplotype's own length
        length = CHR7_MB[hap]
        ax.add_patch(mpatches.FancyBboxPatch(
            (XMIN, -0.16), length, 0.32,
            boxstyle="round,pad=0.02,rounding_size=0.3",
            facecolor="#EEEEEE", edgecolor="#AAAAAA",
            linewidth=0.8, zorder=1))
        ax.text(length + 1.1, 0, f"{length:.1f} Mb", ha="left", va="center",
                fontsize=7.0, color="#777777")

        for t in range(XMIN, int(length) + 1, 5):
            ax.plot([t, t], [-0.16, -0.30], color="#888888",
                    lw=0.8 if t % 10 == 0 else 0.4, zorder=2)
            if i == len(HAPS) - 1 and t % 10 == 0:
                ax.text(t, -0.46, str(t), ha="center", va="top",
                        fontsize=7.0, color="#555555")

        placed = []        # labels already drawn in this row
        for group in cluster(data.get(hap, [])):
            lines, counts, nfunc, centre = describe(group)
            total = sum(counts.values())
            if not total:
                continue

            # cluster block, segmented by series
            x0 = min(max(centre - BLOCK_W / 2, XMIN), length - BLOCK_W)
            x = x0
            for s in SERIES:
                if not counts[s]:
                    continue
                seg = BLOCK_W * counts[s] / total
                ax.add_patch(mpatches.Rectangle(
                    (x, -BLOCK_H / 2), seg, BLOCK_H,
                    facecolor=SERIES_COLOUR[s], edgecolor="none", zorder=4))
                x += seg
            ax.add_patch(mpatches.Rectangle(
                (x0, -BLOCK_H / 2), BLOCK_W, BLOCK_H, facecolor="none",
                edgecolor="#1A1A1A" if nfunc else "#FFFFFF",
                linewidth=1.0 if nfunc else 0.8, zorder=5))

            # label box: one line per series, each in its own colour.
            # Packed by matplotlib so that the box is sized in points and the
            # lines cannot overlap, whatever the subplot height.
            areas = [TextArea(text, textprops=dict(
                color=colour, fontsize=7.0,
                fontweight="bold" if k == 0 else "normal"))
                for k, (text, colour) in enumerate(lines)]
            pack = VPacker(children=areas, pad=1.0, sep=1.6, align="center")
            lx = min(max(centre, LAB_W / 2), length - LAB_W / 2)

            # A haplotype has at most two clusters, and the labels are too wide
            # to sit side by side. The first goes above the chromosome bar and
            # the second below it, which needs no measurement and cannot
            # collide however long the labels become.
            above = (len(placed) == 0)
            ab = AnnotationBbox(
                pack,
                (centre, (BLOCK_H / 2 + 0.03) * (1 if above else -1)),
                xybox=(lx, LAB_BOTTOM if above else -LAB_BOTTOM),
                xycoords="data", boxcoords="data",
                box_alignment=(0.5, 0.0 if above else 1.0),
                pad=0.0, zorder=8,
                bboxprops=dict(boxstyle="round,pad=0.30", facecolor="white",
                               edgecolor="#BBBBBB", linewidth=0.7),
                arrowprops=dict(arrowstyle="-", color="#999999",
                                linewidth=0.7, linestyle=":",
                                shrinkA=0, shrinkB=0))
            ab.set_clip_on(False)
            ax.add_artist(ab)
            placed.append(ab)

        if i == 0:
            ax.text(XMAX + 1.6, 1.30, "Chr7", ha="left", va="center",
                    fontsize=7.0, color="#555555", style="italic")

    axes[-1].text(XMAX / 2, -2.20, "Chromosome 7 position (Mb)",
                  ha="center", va="top", fontsize=8.0, color="#333333")

    handles = [mpatches.Patch(facecolor=SERIES_COLOUR[s], label=f"{s} series")
               for s in ("THCAS", "CBDAS", "CBCAS")]
    handles += [
        mpatches.Patch(facecolor=SERIES_COLOUR["unassigned"],
                       label="series unassigned"),
        mpatches.Patch(facecolor="#CCCCCC", edgecolor="#1A1A1A",
                       linewidth=1.0,
                       label="cluster with a functional candidate"),
        mpatches.Patch(facecolor="#EEEEEE", edgecolor="#AAAAAA",
                       linewidth=0.8, label="chromosome 7 backbone"),
    ]
    fig.legend(handles=handles, loc="lower center", ncol=3, fontsize=7.0,
               framealpha=0.95, bbox_to_anchor=(0.5, 0.010),
               handlelength=1.5, columnspacing=1.6)

    plt.subplots_adjust(left=0.215, right=0.925, top=0.990,
                        bottom=0.085, hspace=0.95)

    fig.savefig(OUTBASE + ".pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(OUTBASE + ".png", dpi=400, bbox_inches="tight",
                facecolor="white")
    print("written:", OUTBASE + ".pdf", "/", OUTBASE + ".png")


if __name__ == "__main__":
    main()
