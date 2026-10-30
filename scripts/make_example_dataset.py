#!/usr/bin/env python3
"""make_example_dataset.py - build the small example dataset shipped with the workflow.

The full assemblies are several hundred megabytes each, so they cannot be
distributed with the source code. This script extracts the window around the
synthase cluster of one haplotype and re-bases every coordinate to 1, so that
the genome FASTA, the GFF3 annotation and the protein FASTA form a
self-consistent miniature input. The whole workflow then runs on it in a few
minutes, exercising the same code path as the full analysis.

Defaults: SAN2b, chromosome 7, 54.0-55.2 Mb. That window holds the entire
THCAS-series cluster of this haplotype (54.5-54.9 Mb), including g313950, the
copy discussed as Case 3, so the example output exercises the premature-stop
rule as well as the ordinary path.

Sequence names and gene identifiers are left unchanged, so the example output
can be compared directly with the published results. Only the coordinates
differ: subtract OFFSET (written to example_data/README_offset.txt) to recover
the coordinates of the published assembly.

Run from the project directory ($ROOT):

    python cannabinoid-synthase-workflow/scripts/make_example_dataset.py \
        --out example_data
"""
import argparse
import gzip
import os
import re
import sys


def open_maybe_gz(path, mode="rt"):
    return gzip.open(path, mode) if path.endswith(".gz") else open(path, mode)


def find_one(*patterns):
    import glob
    for pat in patterns:
        hits = sorted(glob.glob(pat))
        if hits:
            return hits[0]
    return None


def extract_window(genome_path, seqid, start, end):
    """Return the 1-based inclusive window [start, end] of the named sequence."""
    want, chunks, taken = False, [], 0
    with open_maybe_gz(genome_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if want:
                    break
                want = line[1:].split()[0] == seqid
                continue
            if want:
                chunks.append(line.strip())
                taken += len(line.strip())
                if taken >= end:
                    break
    if not chunks:
        sys.exit(f"sequence {seqid!r} not found in {genome_path}")
    seq = "".join(chunks)
    if len(seq) < end:
        sys.exit(f"{seqid} is only {len(seq)} bp; window end {end} is past it")
    return seq[start - 1:end]


def write_fasta(path, seqid, seq, width=60):
    with gzip.open(path, "wt") as fh:
        fh.write(f">{seqid}\n")
        for i in range(0, len(seq), width):
            fh.write(seq[i:i + width] + "\n")


def subset_gff3(gff_path, out_path, seqid, start, end):
    """Keep features of seqid lying wholly inside the window; shift to 1-based."""
    kept, dropped = [], 0
    with open_maybe_gz(gff_path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9 or f[0] != seqid:
                continue
            s, e = int(f[3]), int(f[4])
            if s < start or e > end:
                dropped += 1
                continue
            f[3], f[4] = str(s - start + 1), str(e - start + 1)
            kept.append("\t".join(f))
    with gzip.open(out_path, "wt") as fh:
        fh.write("##gff-version 3\n")
        fh.write(f"##sequence-region {seqid} 1 {end - start + 1}\n")
        for line in kept:
            fh.write(line + "\n")
    return kept, dropped


def mrna_ids(gff_lines):
    """Identifiers of the transcripts retained, to select protein records."""
    ids = set()
    for line in gff_lines:
        f = line.split("\t")
        if f[2] not in ("mRNA", "transcript"):
            continue
        m = re.search(r"(?:^|;)ID=([^;]+)", f[8])
        if m:
            ids.add(m.group(1))
    return ids


def subset_proteins(prot_path, out_path, wanted):
    """Keep protein records whose first token matches a retained transcript."""
    kept = 0
    with open_maybe_gz(prot_path) as fh, gzip.open(out_path, "wt") as out:
        emit = False
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
                emit = name in wanted
                if emit:
                    kept += 1
            if emit:
                out.write(line)
    return kept


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", default=".", help="project root (default: .)")
    ap.add_argument("--out", default="example_data")
    ap.add_argument("--hap", default="SAN2b")
    ap.add_argument("--chrom", default=None,
                    help="sequence name; default <hap>.chr7")
    ap.add_argument("--start-mb", type=float, default=54.0)
    ap.add_argument("--end-mb", type=float, default=55.2)
    args = ap.parse_args()

    root = args.root
    hap = args.hap
    seqid = args.chrom or f"{hap}.chr7"
    start = int(args.start_mb * 1e6)
    end = int(args.end_mb * 1e6)

    genome = find_one(f"{root}/00_genomes/*/{hap}/{hap}.unmasked.fasta.gz",
                      f"{root}/00_genomes/*/{hap}/{hap}.softmasked.fasta.gz")
    gff = find_one(f"{root}/00_genomes/*/{hap}/{hap}.high_confidence.gff3.gz",
                   f"{root}/00_genomes/*/{hap}/{hap}.all.gff3.gz")
    prot = find_one(f"{root}/00_genomes/*/{hap}/{hap}.all.proteins.fasta.gz")
    for label, path in (("genome", genome), ("gff3", gff), ("proteins", prot)):
        if not path:
            sys.exit(f"{label} for {hap} not found under {root}/00_genomes")
        print(f"[in ] {label:9} {path}", file=sys.stderr)

    outdir = os.path.join(args.out, "00_genomes", "example", hap)
    os.makedirs(outdir, exist_ok=True)
    os.makedirs(os.path.join(args.out, "refs"), exist_ok=True)

    seq = extract_window(genome, seqid, start, end)
    write_fasta(os.path.join(outdir, f"{hap}.unmasked.fasta.gz"), seqid, seq)
    print(f"[out] genome    {len(seq):,} bp", file=sys.stderr)

    kept, dropped = subset_gff3(
        gff, os.path.join(outdir, f"{hap}.high_confidence.gff3.gz"),
        seqid, start, end)
    print(f"[out] gff3      {len(kept):,} features "
          f"({dropped:,} outside the window dropped)", file=sys.stderr)

    wanted = mrna_ids(kept)
    n = subset_proteins(
        prot, os.path.join(outdir, f"{hap}.all.proteins.fasta.gz"), wanted)
    print(f"[out] proteins  {n:,} of {len(wanted):,} transcripts", file=sys.stderr)
    if n < len(wanted):
        print("[warn] some transcripts had no protein record", file=sys.stderr)

    # reference sequences are small and are copied verbatim
    import shutil
    for name in ("synthase_query.fasta", "BBE_EcaBBE_P30986.fasta"):
        src = os.path.join(root, "refs", name)
        if os.path.exists(src):
            shutil.copy2(src, os.path.join(args.out, "refs", name))
            print(f"[out] ref       {name}", file=sys.stderr)
        else:
            print(f"[warn] {src} not found; copy it by hand", file=sys.stderr)

    with open(os.path.join(args.out, "README_offset.txt"), "w") as fh:
        fh.write(
            "Example dataset for the cannabinoid synthase workflow.\n\n"
            f"Haplotype          : {hap}\n"
            f"Sequence           : {seqid}\n"
            f"Window in assembly : {start:,}-{end:,} bp "
            f"({args.start_mb}-{args.end_mb} Mb)\n"
            f"Length             : {len(seq):,} bp\n\n"
            "Coordinates in this dataset start at 1. To recover the coordinates\n"
            f"of the published assembly, add {start - 1:,} to every position.\n"
            "Sequence names and gene identifiers are unchanged, so the output can\n"
            "be compared directly with the published results.\n")
    print(f"[out] {args.out}/README_offset.txt", file=sys.stderr)
    print(f"\n[done] example dataset written to {args.out}/", file=sys.stderr)


if __name__ == "__main__":
    main()
