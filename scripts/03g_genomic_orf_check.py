#!/usr/bin/env python3
"""03g_genomic_orf_check.py - genome-side evidence for L2 / L6 (decision F, 2026-10-12).

Rule (and nothing more):
  For each gene, take its best genomic tBLASTn hit (highest bitscore among the reference
  queries) from 02_gene_models.tsv. The genome is "cleanly full-length" only if ALL of:
    1. one HSP covers the reference from start to end
       (qstart <= 1 + start_tol and qend >= qlen - end_tol)
    2. the HSP is gapless and in frame
       (length == qend - qstart + 1 and genomic span == 3 * (qend - qstart + 1))
    3. translating that genomic interval gives no internal stop codon
  -> L2_genomic = yes, L6_genomic = yes, but ONLY if that HSP overlaps a single gene model
     (manuscript: "single overlapping GFF3 gene model"; decision 2026-10-13). A clean full-length
     HSP shared by several gene models gives no genome-side claim (annotation values apply).
  If 1 and 2 hold but 3 fails (full-length, gapless, in frame, internal stop):
  -> L6_genomic = no, Genomic_disruption = premature_stop  (decision 2026-10-12 23:55;
     the manuscript's L6 definition applied to direct genomic evidence). L2_genomic stays '—'.
     This applies to every gene model the HSP overlaps (fragments of one disrupted locus).
  Otherwise both are '—' (no genome-side claim; the annotation-based values are used).

Deliberately NOT done (decision F, "やらないこと"): joining several HSPs, searching for
start codons, inferring introns, merging neighbouring gene models. The genome side produces
'no' only for a premature stop inside a full-length, gapless, in-frame HSP.
Copies are counted per gene model; n_gene_models_on_hsp records how many models share the HSP.
"""
import argparse
import csv
import gzip
import re
import sys

MISSING = '—'
CODON = {}
_b = 'TCAG'
_aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
_i = 0
for a in _b:
    for b in _b:
        for c in _b:
            CODON[a + b + c] = _aa[_i]
            _i += 1


def translate(nt):
    nt = nt.upper()
    return ''.join(CODON.get(nt[i:i + 3], 'X') for i in range(0, len(nt) - 2, 3))


def revcomp(nt):
    return nt.upper().translate(str.maketrans('ACGTN', 'TGCAN'))[::-1]


def norm(x):
    s = str(x or '').strip()
    return re.sub(r'\.t\d+$', '', s.split()[0]) if s else ''


def opener(path):
    return gzip.open(path, 'rt') if path.endswith('.gz') else open(path)


def load_contigs(path, wanted):
    """Stream the genome once and keep only the contigs that carry hits."""
    seqs, cur, buf = {}, None, []
    with opener(path) as f:
        for line in f:
            if line.startswith('>'):
                if cur is not None:
                    seqs[cur] = ''.join(buf)
                name = line[1:].split()[0]
                cur = name if name in wanted else None
                buf = []
            elif cur is not None:
                buf.append(line.strip())
        if cur is not None:
            seqs[cur] = ''.join(buf)
    return seqs


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--gene-models', required=True, help='02_gene_models.tsv (tBLASTn hits mapped to genes)')
    ap.add_argument('--genome', required=True, help='genome FASTA (.fa/.fasta, optionally .gz)')
    ap.add_argument('--out', required=True)
    ap.add_argument('--start-tol', type=int, default=3, help='max missing residues at reference N-terminus')
    ap.add_argument('--end-tol', type=int, default=3, help='max missing residues at reference C-terminus')
    a = ap.parse_args()

    with open(a.gene_models, newline='') as f:
        rows = list(csv.DictReader(f, delimiter='\t'))
    need = ('qseqid', 'sseqid', 'length', 'qstart', 'qend', 'sstart', 'send', 'bitscore', 'qlen',
            'overlapping_gene_ids')
    if rows and any(k not in rows[0] for k in need):
        sys.exit(f'[03g] {a.gene_models}: missing columns; need {need}')

    best = {}
    for r in rows:
        for g in re.split(r'[,;\s]+', r.get('overlapping_gene_ids') or ''):
            g = norm(g)
            if g and (g not in best or float(r['bitscore']) > float(best[g]['bitscore'])):
                best[g] = r

    contigs = load_contigs(a.genome, {r['sseqid'] for r in best.values()})
    missing_contigs = {r['sseqid'] for r in best.values()} - set(contigs)
    if missing_contigs:
        print(f'[03g] WARNING: contigs not found in genome FASTA: {sorted(missing_contigs)[:5]}', file=sys.stderr)

    out_cols = ['Gene_ID', 'Genomic_query', 'qstart', 'qend', 'qlen', 'aln_length', 'genomic_span_nt',
                'c1_full_coverage', 'c2_gapless_in_frame', 'c3_no_internal_stop',
                'n_gene_models_on_hsp', 'L2_genomic', 'L6_genomic', 'Genomic_disruption', 'Genomic_evidence']
    n_yes = n_stop = n_shared = 0
    with open(a.out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=out_cols, delimiter='\t')
        w.writeheader()
        for g, r in sorted(best.items()):
            qs, qe, ql, ln = int(r['qstart']), int(r['qend']), int(r['qlen']), int(r['length'])
            ss, se = int(r['sstart']), int(r['send'])
            lo, hi = min(ss, se), max(ss, se)
            span = hi - lo + 1
            aa_len = qe - qs + 1
            c1 = qs <= 1 + a.start_tol and qe >= ql - a.end_tol
            c2 = ln == aa_len and span == 3 * aa_len
            c3 = None
            if c1 and c2 and r['sseqid'] in contigs:
                nt = contigs[r['sseqid']][lo - 1:hi]
                if ss > se:
                    nt = revcomp(nt)
                c3 = '*' not in translate(nt)
            n_models = len({norm(x) for x in re.split(r'[,;\s]+', r.get('overlapping_gene_ids') or '') if norm(x)})
            clean = bool(c1 and c2 and c3)
            ok = clean and n_models == 1
            stop = bool(c1 and c2 and c3 is False)
            n_shared += clean and not ok
            n_yes += ok
            n_stop += stop
            yn = lambda v: 'yes' if v else ('no' if v is not None else MISSING)
            w.writerow({
                'Gene_ID': g, 'Genomic_query': r['qseqid'], 'qstart': qs, 'qend': qe, 'qlen': ql,
                'aln_length': ln, 'genomic_span_nt': span,
                'c1_full_coverage': yn(c1), 'c2_gapless_in_frame': yn(c2), 'c3_no_internal_stop': yn(c3),
                'n_gene_models_on_hsp': n_models,
                'L2_genomic': 'yes' if ok else MISSING,
                'L6_genomic': 'yes' if ok else ('no' if stop else MISSING),
                'Genomic_disruption': 'premature_stop' if stop else '',
                'Genomic_evidence': (f'single gapless tBLASTn HSP {qs}-{qe}/{ql} aa, {span} nt, no internal stop' if ok else
                                     f'clean full-length HSP shared by {n_models} gene models; annotation values apply' if clean else
                                     f'single gapless tBLASTn HSP {qs}-{qe}/{ql} aa, {span} nt, INTERNAL STOP' if stop else
                                     'not cleanly full-length on genome; annotation values apply'),
            })
    print(f'[03g] genes: {len(best)}  cleanly full-length: {n_yes}  clean but shared HSP (not applied): {n_shared}  '
          f'full-length with internal stop: {n_stop}  -> {a.out}',
          file=sys.stderr)


if __name__ == '__main__':
    main()
