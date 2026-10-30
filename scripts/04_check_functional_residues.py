#!/usr/bin/env python3
"""Extract function-related residues from an aligned protein FASTA.

Study defaults:
  FAD-related and active-site residues: THCAS reference numbering
  Specificity-related residues: CBDAS reference numbering
"""
from __future__ import annotations
import argparse, csv
from Bio import AlignIO

GROUPS = {
    "FAD-related": [114,176,184],
    "Active-site": [292,417,484],
    "Specificity-related": [192,414,443,445],
}
EXPECTED = {
    "FAD-related": list("HCH"),
    "Active-site": list("HYY"),
    "Specificity-related": list("PAWI"),
}


def posmap(seq):
    d={}; n=0
    for i,c in enumerate(str(seq)):
        if c not in "-.":
            n+=1; d[n]=i
    return d


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--alignment", required=True)
    ap.add_argument("--thcas-id", required=True, help="THCAS anchor sequence ID in the alignment")
    ap.add_argument("--cbdas-id", required=True, help="CBDAS anchor sequence ID in the alignment")
    ap.add_argument("--out", required=True)
    args=ap.parse_args()

    aln=AlignIO.read(args.alignment,"fasta")
    seqs={r.id:r for r in aln}
    for x in (args.thcas_id,args.cbdas_id):
        if x not in seqs: raise SystemExit(f"Anchor ID not found: {x}")
    maps={args.thcas_id:posmap(seqs[args.thcas_id].seq), args.cbdas_id:posmap(seqs[args.cbdas_id].seq)}

    fields=["Sequence","FAD_related_residues","FAD_related_conserved","Active_site_residues",
            "Active_site_conserved","Specificity_related_residues","Specificity_pattern",
            "L3_FAD_evidence","L4_active_site_evidence"]
    with open(args.out,"w",newline="") as oh:
        w=csv.DictWriter(oh,fieldnames=fields,delimiter="\t",lineterminator="\n"); w.writeheader()
        for r in aln:
            s=str(r.seq)
            vals={}
            for group, positions in GROUPS.items():
                anchor=args.cbdas_id if group=="Specificity-related" else args.thcas_id
                aa=[]
                for p in positions:
                    col=maps[anchor].get(p)
                    aa.append(s[col] if col is not None and col < len(s) and s[col] not in "." else "-")
                vals[group]=aa
            fad="/".join(vals["FAD-related"]); act="/".join(vals["Active-site"]); spec="/".join(vals["Specificity-related"])
            fad_ok=vals["FAD-related"]==EXPECTED["FAD-related"]
            act_ok=vals["Active-site"]==EXPECTED["Active-site"]
            if vals["Specificity-related"] == list("PAWI"):
                pattern="CBDAS-type"
            elif vals["Specificity-related"] == list("AVWT"):
                pattern="THCAS/CBCAS-type"
            elif "-" in vals["Specificity-related"]:
                pattern="Not applicable"
            else:
                pattern="Mixed"
            w.writerow({"Sequence":r.id,"FAD_related_residues":fad,"FAD_related_conserved":"yes" if fad_ok else "no",
                        "Active_site_residues":act,"Active_site_conserved":"yes" if act_ok else "no",
                        "Specificity_related_residues":spec,"Specificity_pattern":pattern,
                        "L3_FAD_evidence":"yes" if fad_ok else "no","L4_active_site_evidence":"yes" if act_ok else "no"})

if __name__=="__main__": main()
