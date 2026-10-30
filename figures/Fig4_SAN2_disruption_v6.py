#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Fig. 4  Domain architecture and marker residues of the three cannabinoid
        synthase copies of the two Santhica 27 haplotypes, compared with a
        functional synthase.

v5 (2026-10-23). Changes from v4, all forced by the data or by a source:

  * Content now comes from the workflow output (Additional file 1 and
    04_functional_residues.tsv) instead of the hand-maintained Fig4_SAN2_data.tsv.
  * Domain boundaries follow the crystal structures of the three synthases:
    an FAD-binding domain of residues 29-251 and 477-544 and a substrate-binding
    domain of residues 252-476. v4 drew a single FAD-binding block, labelled
    477-544 as "C-terminal", and drew "active site" as a separate domain, which
    it is not: the active-site residues lie in both domains.
  * Residue numbering is reconciled. FAD-related and active-site residues are in
    THCAS numbering, specificity-related residues in CBDAS numbering (as in
    Table 3), and CBDAS numbering runs one behind THCAS numbering from residue
    253 onwards because of the deletion of THCAS Ser253 in CBDAS. v4 plotted
    both sets on one axis without converting, so A414, W443 and I445 were each
    drawn one residue off, and it also mislabelled H292 as H291 and Y484 as Y483.
  * SAN2a g539950: residues 229-301 (THCAS numbering) of the predicted protein
    have no counterpart in the reference, so His292 is drawn as having no
    counterpart (hatched segment). The locus carries a 7-nt deletion that the
    GFF3 annotation models as an intron (manuscript Section 2.7).
  * Conserved residues and series-defining residues are drawn with different
    symbols. Specificity-related residues differ between functional synthases
    (THCAS/CBCAS carry A/V/W/T, CBDAS carries P/A/W/I), so a difference there
    is not evidence of lost function; v4 drew them alongside the loss markers.
  * The claims "CDS fragmentation", "CDS split", "structurally compromised" and
    "catalytically inactive" are gone. g539950 is classified as a divergent
    copy with no internal stop, and in silico analysis cannot establish that an
    enzyme is inactive.

v6 (2026-10-30): rows are ordered as Cases 1-3 of the manuscript, and SAN2b
g313950 lists L5 among its unmet levels, as in Table 4.

Input : Additional_file_1_synthase_copy_detail.tsv, and
        step4_rooted_20261012/<hap>/04_functional_residues.tsv
Output: Fig4_SAN2_disruption_v6.{pdf,png}

    python3 Fig4_SAN2_disruption_v6.py [project_root]
"""
import os
import sys
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle
from matplotlib.lines import Line2D

matplotlib.rcParams["pdf.fonttype"] = 42
matplotlib.rcParams["ps.fonttype"] = 42
matplotlib.rcParams["font.family"] = "DejaVu Sans"

MM = 1 / 25.4
FIG_W, FIG_H = 180 * MM, 150 * MM

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = sys.argv[1] if len(sys.argv) > 1 else os.path.normpath(
    os.path.join(HERE, "..", ".."))
AF1 = os.path.join(ROOT, "Additional_file_1_synthase_copy_detail.tsv")
RESID = os.path.join(ROOT, "step4_rooted_20261012", "{hap}",
                     "04_functional_residues.tsv")
OUTBASE = os.path.join(HERE, "Fig4_SAN2_disruption_v6")

REF_LEN = 545        # THCAS reference precursor length
BAR_H = 0.46
TRACK_DY = 1.65      # vertical spacing between tracks

# Domain boundaries, THCAS/CBCAS numbering, from the crystal structures of
# CBCAS, CBDAS and THCAS: "an FAD binding domain (residues N29-V251 and
# 477-544) and a substrate binding domain (residues P252-N476)". Residues
# 1-28 precede the modelled N-terminus and are the cleaved signal peptide.
DOMAINS = [
    (1, 28, "Signal peptide", "#C8E6C9"),
    (29, 251, "FAD-binding domain", "#66BB6A"),
    (252, 476, "Substrate-binding domain", "#43A047"),
    (477, 544, "FAD-binding domain\n(C-terminal part)", "#66BB6A"),
]

# Marker residues. "ref" is the numbering the manuscript uses in Table 3;
# "axis" is the position on the THCAS axis used for drawing.
#   FAD-related and active-site: THCAS numbering, expected H/C/H and H/Y/Y
#   specificity-related:         CBDAS numbering, expected P/A/W/I
# CBDAS numbering is one behind THCAS numbering from residue 253 onwards.
def to_thcas(cbdas_pos):
    return cbdas_pos if cbdas_pos <= 252 else cbdas_pos + 1


FAD = [(114, "H"), (176, "C"), (184, "H")]
ACTIVE = [(292, "H"), (417, "Y"), (484, "Y")]
SPECIFICITY = [(192, "P"), (414, "A"), (443, "W"), (445, "I")]

GENES = [
    ("SAN2a", "SAN2a.chr7.v1.g539950"),
    ("SAN2a", "SAN2a.chr7.v1.g547390"),
    ("SAN2b", "SAN2b.chr7.v1.g313950"),
]

# The criterion each copy fails, restating Table 4. Kept here because L5 is a
# phylogenetic result that this figure cannot show.
# Regions of the predicted protein that have no counterpart in the reference,
# so that the residue sitting in a given alignment column there is not a
# homologous residue and must not be read as a substitution. Given in THCAS
# reference numbering, from the four-sequence alignment of g539950 with the
# three references (09_san2a_mapping/g539950_validation/aligned.fasta):
# identity with CBCAS is exact outside 229-301 and absent within it, and the
# whole 31-residue length deficit falls inside this interval. The underlying
# cause is a 7 nt deletion in the locus near codon 231/232 (Results section 5,
# Case 1), which the GFF3 models as an 86 nt intron.
NON_HOMOLOGOUS = {
    "SAN2a.chr7.v1.g539950": [(229, 301)],
}

REASON = {
    "SAN2a.chr7.v1.g539950":
        "L4 unmet: no counterpart to active-site His292\n"
        "(7 nt deletion in the locus; frame lost)",
    "SAN2b.chr7.v1.g313950":
        "L2, L3, L5 and L6 unmet: premature stop codon;\n"
        "N-terminal region not encoded",
    "SAN2a.chr7.v1.g547390":
        "L5 unmet: not assigned to a functional clade\n"
        "(Fig. 3); all residue criteria met",
}
SERIES_COLOUR = {"THCAS": "#B71C1C", "CBDAS": "#1565C0", "CBCAS": "#6A1B9A"}

C_OK = "#333333"
C_SUB = "#D32F2F"
C_ABS = "#999999"
C_SPEC = "#E65100"


def read_tsv(path):
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def load():
    af1 = {r["Gene_ID"]: r for r in read_tsv(AF1)}
    out = []
    for hap, gid in GENES:
        if gid not in af1:
            sys.exit(f"{gid} not in {AF1}")
        res = {r["Sequence"].rsplit(".t", 1)[0]: r
               for r in read_tsv(RESID.format(hap=hap))}
        if gid not in res:
            sys.exit(f"{gid} not in {RESID.format(hap=hap)}")
        a, r = af1[gid], res[gid]
        out.append(dict(
            gene=gid, hap=hap,
            short=gid.split(".")[-1],
            series=a["Series"],
            identity=float(a["Nearest_reference_identity_pct"]),
            length=int(a["Protein_length_aa"]),
            coverage=float(a["Reference_coverage"]),
            classification=a["Classification"],
            fad=r["FAD_related_residues"].split("/"),
            active=r["Active_site_residues"].split("/"),
            spec=r["Specificity_related_residues"].split("/"),
            pattern=r["Specificity_pattern"],
        ))
    return out


def main():
    for p in (AF1,):
        if not os.path.exists(p):
            sys.exit(f"not found: {p}")
    copies = load()
    ntrack = len(copies) + 1

    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    fig.patch.set_facecolor("white")

    def ytop(t):
        return (ntrack - 1 - t) * TRACK_DY

    # ---- reference track ---------------------------------------------------
    y = ytop(0)
    for a, b, name, colour in DOMAINS:
        ax.add_patch(Rectangle((a - 1, y - BAR_H / 2), b - a + 1, BAR_H,
                               facecolor=colour, edgecolor="#2E7D32",
                               linewidth=0.7, zorder=2))
    lift = [0.50, 0.14, 0.14, 0.50]
    for (a, b, name, _), lf in zip(DOMAINS, lift):
        xc = (a + b) / 2
        ax.plot([xc, xc], [y + BAR_H / 2, y + BAR_H / 2 + lf],
                color="#9E9E9E", lw=0.5, zorder=1)
        ax.text(xc, y + BAR_H / 2 + lf + 0.02, name, ha="center", va="bottom",
                fontsize=7.0, color="#1B5E20", fontweight="bold",
                linespacing=1.3, zorder=5)

    # marker residues on the reference, drawn under the bar
    ruler = ([(to_thcas(p), f"{aa}{p}", C_SPEC) for p, aa in SPECIFICITY]
             + [(p, f"{aa}{p}", C_OK) for p, aa in FAD + ACTIVE])
    for i, (pos, name, colour) in enumerate(sorted(ruler)):
        drop = (0.12, 0.31, 0.50, 0.69)[i % 4]
        ax.plot([pos, pos], [y - BAR_H / 2, y - BAR_H / 2 - drop],
                color="#AAAAAA", lw=0.6, zorder=3)
        ax.text(pos, y - BAR_H / 2 - drop - 0.02, name, ha="center", va="top",
                fontsize=7.0, color=colour, zorder=5)

    # ---- one track per copy ------------------------------------------------
    for t, c in enumerate(copies, start=1):
        y = ytop(t)
        col = SERIES_COLOUR.get(c["series"], "#555555")

        # a leading region that is not encoded, when the leading markers are
        # absent; its extent follows from the protein length
        lead_absent = all(x == "-" for x in c["fad"])
        start = (REF_LEN - c["length"] + 1) if lead_absent else 1
        if lead_absent:
            ax.add_patch(Rectangle((0, y - BAR_H / 2), start - 1, BAR_H,
                                   facecolor="white", edgecolor="#BBBBBB",
                                   linewidth=0.7, linestyle=(0, (2, 2)),
                                   zorder=2))
            ax.text((start - 1) / 2, y - BAR_H / 2 - 0.08,
                    f"residues 1–{start - 1} not encoded",
                    ha="center", va="top", fontsize=7.0, color="#777777",
                    style="italic", zorder=5)
        ax.add_patch(Rectangle((start - 1, y - BAR_H / 2),
                               REF_LEN - start + 1, BAR_H,
                               facecolor="#E8E8E8", edgecolor=col,
                               linewidth=0.9, zorder=2))

        # regions with no counterpart in the reference
        nonhom = NON_HOMOLOGOUS.get(c["gene"], [])
        for a, b in nonhom:
            ax.add_patch(Rectangle((a - 1, y - BAR_H / 2), b - a + 1, BAR_H,
                                   facecolor="#FFFFFF", edgecolor=col,
                                   hatch="///", linewidth=0.9, zorder=3))
            ax.text((a + b) / 2, y - BAR_H / 2 - 0.08,
                    f"residues {a}–{b}: no counterpart\nin the reference",
                    ha="center", va="top", fontsize=7.0, color="#777777",
                    style="italic", linespacing=1.4, zorder=5)

        def in_nonhom(p):
            return any(a <= p <= b for a, b in nonhom)

        # conserved-or-not markers: FAD-related and active-site
        for (pos, exp), obs in zip(FAD + ACTIVE, c["fad"] + c["active"]):
            if obs == "-" or in_nonhom(pos):
                # inside a non-homologous region the column occupant is not a
                # homologous residue, so it is shown as having no counterpart
                ax.plot([pos], [y], marker="x", ms=4.6, color=C_ABS, mew=1.3,
                        zorder=6)
            elif obs == exp:
                ax.plot([pos], [y], marker="o", ms=3.4, mfc="white",
                        mec=C_OK, mew=0.9, zorder=6)
            else:
                ax.plot([pos], [y], marker="v", ms=5.2, color=C_SUB, zorder=6)
                ax.text(pos, y + BAR_H / 2 + 0.06, f"{exp}{pos}{obs}",
                        ha="center", va="bottom", fontsize=7.0, color=C_SUB,
                        fontweight="bold", zorder=6)

        # series-defining markers: specificity-related
        for (pos, _), obs in zip(SPECIFICITY, c["spec"]):
            xp = to_thcas(pos)
            if obs == "-" or in_nonhom(xp):
                ax.plot([xp], [y], marker="x", ms=4.6, color=C_ABS, mew=1.3,
                        zorder=6)
            else:
                ax.plot([xp], [y], marker="D", ms=3.2, mfc=C_SPEC,
                        mec=C_SPEC, zorder=6)

        # row labels
        head = f"{c['hap']}  {c['short']}  —  {c['length']} aa"
        ax.text(-18, y + 0.46, head, ha="right", va="center", fontsize=7.0,
                fontweight="bold", color=col, zorder=5)
        sub = (f"{c['classification']}\n"
               f"{c['identity']:.1f}% identity, "
               f"{c['coverage'] * 100:.1f}% of reference length\n"
               f"specificity pattern: {c['pattern']}\n"
               f"{REASON[c['gene']]}")
        ax.text(-18, y - 0.26, sub, ha="right", va="center", fontsize=7.0,
                color="#555555", linespacing=1.5, zorder=5)

    ax.text(-18, ytop(0) + 0.40, f"Functional synthase (reference, {REF_LEN} aa)",
            ha="right", va="center", fontsize=7.0, fontweight="bold",
            color="#1B5E20", zorder=5)
    ax.text(-18, ytop(0) - 0.14,
            "THCAS Q8GTB6 / CBDAS A6P6V9 / CBCAS Q33DQ2",
            ha="right", va="center", fontsize=7.0, color="#555555",
            style="italic", linespacing=1.45, zorder=5)

    # ---- amino acid axis ---------------------------------------------------
    yax = ytop(ntrack - 1) - 0.92
    ax.plot([0, REF_LEN], [yax, yax], color="#555555", lw=0.8)
    for xx in range(0, REF_LEN + 1, 100):
        ax.plot([xx, xx], [yax, yax - 0.07], color="#555555", lw=0.8)
        ax.text(xx, yax - 0.10, str(xx), ha="center", va="top", fontsize=7.0,
                color="#555555")
    ax.text(REF_LEN / 2, yax - 0.34,
            "Amino acid position (THCAS reference numbering)",
            ha="center", va="top", fontsize=7.0, color="#444444")

    handles = [
        Rectangle((0, 0), 1, 1, fc="#66BB6A", ec="#2E7D32",
                  label="FAD-binding domain"),
        Rectangle((0, 0), 1, 1, fc="#43A047", ec="#2E7D32",
                  label="Substrate-binding domain"),
        Rectangle((0, 0), 1, 1, fc="#C8E6C9", ec="#2E7D32",
                  label="Signal peptide"),
        Rectangle((0, 0), 1, 1, fc="white", ec="#BBBBBB", ls="--",
                  label="Not encoded"),
        Rectangle((0, 0), 1, 1, fc="white", ec="#555555", hatch="///",
                  label="No counterpart in the reference"),
        Line2D([], [], ls="none", marker="o", mfc="white", mec=C_OK, ms=5,
               label="FAD/active-site residue conserved"),
        Line2D([], [], ls="none", marker="v", color=C_SUB, ms=6,
               label="FAD/active-site residue substituted"),
        Line2D([], [], ls="none", marker="x", color=C_ABS, ms=6, mew=1.3,
               label="Residue not present in the protein"),
        Line2D([], [], ls="none", marker="D", color=C_SPEC, ms=5,
               label="Specificity-related residue (series-defining)"),
    ]
    ax.legend(handles=handles, loc="upper left", bbox_to_anchor=(0.0, -0.015),
              ncol=2, frameon=True, fontsize=7.0, handlelength=1.5,
              columnspacing=1.4, borderpad=0.7, labelspacing=0.45)

    ax.set_xlim(-265, REF_LEN + 20)
    ax.set_ylim(yax - 0.60, ytop(0) + 1.20)
    ax.axis("off")
    fig.tight_layout()

    fig.savefig(OUTBASE + ".pdf", bbox_inches="tight", facecolor="white")
    fig.savefig(OUTBASE + ".png", dpi=400, bbox_inches="tight",
                facecolor="white")
    print("written:", OUTBASE + ".pdf", "/", OUTBASE + ".png")


if __name__ == "__main__":
    main()
