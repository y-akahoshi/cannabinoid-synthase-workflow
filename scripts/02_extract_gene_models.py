#!/usr/bin/env python3
"""Map BLAST/tBLASTn genomic hits to overlapping GFF3 gene models.

The script deliberately reports *all* overlapping gene IDs. A single genomic
hit may overlap two predicted genes when an annotation has split one biological
locus; such cases are reported, not merged.
"""
from __future__ import annotations
import argparse, csv, gzip, re
from pathlib import Path


def op(path, mode="rt"):
    return gzip.open(path, mode) if str(path).endswith(".gz") else open(path, mode)


def parse_attrs(s: str) -> dict[str, str]:
    d = {}
    for token in s.strip().split(";"):
        if not token:
            continue
        if "=" in token:
            k, v = token.split("=", 1)
            d[k] = v
        elif " " in token:  # tolerant of GTF-like attributes
            k, v = token.split(" ", 1)
            d[k] = v.strip('"')
    return d


def load_genes(gff: str, feature: str = "gene"):
    by_chrom = {}
    with op(gff) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9 or p[2] != feature:
                continue
            attrs = parse_attrs(p[8])
            gid = attrs.get("ID") or attrs.get("gene_id") or attrs.get("Name") or "NA"
            rec = (int(p[3]), int(p[4]), p[6], gid)
            by_chrom.setdefault(p[0], []).append(rec)
    for chrom in by_chrom:
        by_chrom[chrom].sort()
    return by_chrom


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--hits", required=True, help="BLAST outfmt 6 from script 01")
    ap.add_argument("--gff", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--feature", default="gene")
    args = ap.parse_args()

    genes = load_genes(args.gff, args.feature)
    cols = ["qseqid","sseqid","pident","length","qstart","qend","sstart","send","evalue","bitscore","qlen","slen","qcovs"]
    outcols = cols + ["hit_start","hit_end","hit_strand","overlapping_gene_ids","n_overlapping_genes"]

    with open(args.hits, newline="") as ih, open(args.out, "w", newline="") as oh:
        rd = csv.reader(ih, delimiter="\t")
        wr = csv.DictWriter(oh, fieldnames=outcols, delimiter="\t", lineterminator="\n")
        wr.writeheader()
        for row in rd:
            if not row or row[0].startswith("#"):
                continue
            if len(row) < len(cols):
                raise SystemExit(f"Expected {len(cols)} BLAST columns, got {len(row)}: {row}")
            r = dict(zip(cols, row[:len(cols)]))
            s1, s2 = int(r["sstart"]), int(r["send"])
            start, end = sorted((s1, s2))
            strand = "+" if s1 <= s2 else "-"
            overlap = []
            for gs, ge, gstrand, gid in genes.get(r["sseqid"], []):
                if ge < start:
                    continue
                if gs > end:
                    break
                if not (ge < start or gs > end):
                    overlap.append(gid)
            r.update(hit_start=start, hit_end=end, hit_strand=strand,
                     overlapping_gene_ids=";".join(overlap) if overlap else "NO_GENE_MODEL",
                     n_overlapping_genes=len(overlap))
            wr.writerow(r)

if __name__ == "__main__":
    main()
