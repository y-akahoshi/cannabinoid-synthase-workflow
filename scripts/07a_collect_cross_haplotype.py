#!/usr/bin/env python3
"""07a_collect_cross_haplotype.py - gather the sequences shown in Fig. 3.

Fig. 3 is a cross-haplotype tree: it shows the functional candidates of the
THCAS and CBDAS series together with the three copies discussed as Cases 1-3
in Results section 5, so that their clade assignment (L5) can be read from one
tree. The functional CBCAS candidates are not shown; they are many and their
clade assignment is not in question.

Protein sequences come from each haplotype's proteins.fasta, written by
run_haplotype.sh, and the references from refs/. Run from the project directory ($ROOT):

    python cannabinoid-synthase-workflow/scripts/07a_collect_cross_haplotype.py \
        --out 08_results/figures/Fig3_cross_haplotype.fasta
"""
import argparse
import os
import sys

# functional candidates, THCAS and CBDAS series (Results section 7)
FUNCTIONAL = [
    "SODLa.chr7.v1.g255560.t1",
    "SODLb.chr7.v1.g250920.t1",
    "EH23a.chr7.v1.g200100.t1",
    "EH23b.chr7.v1.g195030.t1",
    "KOMPa.chr7.v1.g259030.t1",
    "KOMPb.chr7.v1.g226860.t1",
]
# the three copies discussed as Cases 1-3 in Results section 5
CASES = [
    "SAN2a.chr7.v1.g539950.t1",
    "SAN2a.chr7.v1.g547390.t1",
    "SAN2b.chr7.v1.g313950.t1",
]
REF_FILES = ["refs/synthase_query.fasta", "refs/BBE_EcaBBE_P30986.fasta"]
WORKDIR = "results"


def read_fasta(path):
    """{first token of header: (full header, sequence)}"""
    out, name, head, buf = {}, None, None, []
    with open(path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name:
                    out[name] = (head, "".join(buf))
                head = line[1:]
                name = head.split()[0] if head.split() else head
                buf = []
            elif name:
                buf.append(line.strip())
    if name:
        out[name] = (head, "".join(buf))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True)
    ap.add_argument("--workdir", default=WORKDIR)
    args = ap.parse_args()

    wanted = FUNCTIONAL + CASES
    found, missing = [], []

    for gid in wanted:
        hap = gid.split(".")[0]
        prot = os.path.join(args.workdir, hap, "proteins.fasta")
        if not os.path.exists(prot):
            missing.append(f"{gid}  ({prot} not found)")
            continue
        recs = read_fasta(prot)
        if gid in recs:
            found.append((gid, recs[gid][1]))
        else:
            # the annotation may use a different transcript suffix
            stem = gid.rsplit(".t", 1)[0]
            alt = [k for k in recs if k.rsplit(".t", 1)[0] == stem]
            if alt:
                found.append((alt[0], recs[alt[0]][1]))
                print(f"[note] {gid} not found; used {alt[0]}", file=sys.stderr)
            else:
                missing.append(f"{gid}  (not in {prot})")

    for rf in REF_FILES:
        if not os.path.exists(rf):
            missing.append(f"{rf} not found")
            continue
        for name, (head, seq) in read_fasta(rf).items():
            found.append((head, seq))     # keep the full reference header

    if missing:
        print("[error] missing sequences:", file=sys.stderr)
        for m in missing:
            print("   ", m, file=sys.stderr)
        sys.exit(1)

    os.makedirs(os.path.dirname(args.out) or ".", exist_ok=True)
    with open(args.out, "w") as fh:
        for head, seq in found:
            fh.write(f">{head}\n")
            for i in range(0, len(seq), 60):
                fh.write(seq[i:i + 60] + "\n")

    print(f"[07a] wrote {len(found)} sequences to {args.out}", file=sys.stderr)
    for head, seq in found:
        print(f"    {head.split()[0]:34} {len(seq):4} aa", file=sys.stderr)


if __name__ == "__main__":
    main()
