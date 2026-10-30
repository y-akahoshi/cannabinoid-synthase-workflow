#!/usr/bin/env python3
"""summarize_no_model_hits.py - tBLASTn hits with no overlapping GFF3 gene model, one row per region.

Decision 2026-10-13: such hits are not classified (no gene-model inference); they are listed in
Additional file 1 with coordinates, nearest reference (highest identity, BLAST rule) and identities.

Hits from the three reference queries that overlap on the same chromosome and strand are merged into
one region. Run from the project directory ($ROOT):
    python summarize_no_model_hits.py > no_gene_model_regions.tsv
"""
import csv
import re
import sys

HAPS = ['SODLa', 'SODLb', 'EH23a', 'EH23b', 'KOMPa', 'KOMPb', 'SAN2a', 'SAN2b']
GM = 'work/{h}/02_gene_models.tsv'


def ref_name(q):
    m = re.search(r'(THCAS|CBDAS|CBCAS)', q)
    return m.group(1) if m else q


out = csv.writer(sys.stdout, delimiter='\t', lineterminator='\n')
out.writerow(['Haplotype', 'Chromosome', 'Start', 'End', 'Strand', 'Nearest_reference', 'Identity_pct',
              'Ref_range_aa', 'Aln_length_aa', 'Identity_THCAS', 'Identity_CBDAS', 'Identity_CBCAS',
              'High_identity_long', 'Note'])
for h in HAPS:
    rows = [r for r in csv.DictReader(open(GM.format(h=h)), delimiter='\t')
            if r['overlapping_gene_ids'] == 'NO_GENE_MODEL']
    hits = []
    for r in rows:
        lo, hi = sorted((int(r['sstart']), int(r['send'])))
        hits.append(dict(chrom=r['sseqid'], strand=r['hit_strand'], lo=lo, hi=hi, ref=ref_name(r['qseqid']),
                         pid=float(r['pident']), qs=int(r['qstart']), qe=int(r['qend']), ln=int(r['length'])))
    hits.sort(key=lambda x: (x['chrom'], x['strand'], x['lo']))
    regions = []
    for x in hits:
        if regions and regions[-1]['chrom'] == x['chrom'] and regions[-1]['strand'] == x['strand'] \
                and x['lo'] <= regions[-1]['hi']:
            reg = regions[-1]
            reg['hi'] = max(reg['hi'], x['hi'])
            reg['hits'].append(x)
        else:
            regions.append(dict(chrom=x['chrom'], strand=x['strand'], lo=x['lo'], hi=x['hi'], hits=[x]))
    for reg in regions:
        best_by_ref = {}
        for x in reg['hits']:
            if x['ref'] not in best_by_ref or x['ln'] * x['pid'] > best_by_ref[x['ref']]['ln'] * best_by_ref[x['ref']]['pid']:
                best_by_ref[x['ref']] = x
        top = max(best_by_ref.values(), key=lambda x: (x['pid'], x['ln']))
        flag = 'yes' if top['pid'] >= 90 and top['ln'] >= 200 else 'no'
        ident = lambda r: f"{best_by_ref[r]['pid']:.1f}" if r in best_by_ref else ''
        out.writerow([h, reg['chrom'], reg['lo'], reg['hi'], reg['strand'], top['ref'], f"{top['pid']:.1f}",
                      f"{top['qs']}-{top['qe']}", top['ln'], ident('THCAS'), ident('CBDAS'), ident('CBCAS'),
                      flag, 'no overlapping GFF3 gene model; not classified'])
