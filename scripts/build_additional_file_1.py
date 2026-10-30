#!/usr/bin/env python3
"""build_additional_file_1.py - Additional file 1 from the workflow output (one row per copy).

Replaces the hand-made table: every row comes from 06_classification.tsv (classified gene models)
or no_gene_model_regions.tsv (tBLASTn hits with no overlapping gene model, not classified).
Coordinates and identity come from 02_gene_models.tsv (the best tBLASTn hit for that gene,
same source as Fig. 2).

Run from the project directory ($ROOT):
    python cannabinoid-synthase-workflow/scripts/build_additional_file_1.py \
        > Additional_file_1_synthase_copy_detail.tsv
"""
import csv
import re
import sys

HAPS = [('SODLa', 'Type I'), ('SODLb', 'Type I'), ('EH23a', 'Type II'), ('EH23b', 'Type II'),
        ('KOMPa', 'Type III'), ('KOMPb', 'Type III'), ('SAN2a', 'Type IV'), ('SAN2b', 'Type IV')]
CLS = 'results/{h}/06_classification.tsv'
ORF = 'results/{h}/03_orf.tsv'
GEN = 'results/{h}/03g_genomic_orf.tsv'
GM = 'work/{h}/02_gene_models.tsv'
NOMODEL = 'no_gene_model_regions.tsv'
MARK = {'yes': '✓', 'no': '✗'}
COLS = ['Haplotype', 'Chemotype', 'Gene_ID', 'Chromosome', 'Start_Mb', 'End_Mb', 'Strand',
        'Series', 'Nearest_reference_identity_pct', 'Protein_length_aa', 'Reference_coverage',
        'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7', 'Classification',
        'L2_L6_evidence_source', 'Clade_assignment', 'Disruption', 'Gene_models_on_same_hit', 'Note']


def series(ref):
    hits = [x for x in ('THCAS', 'CBDAS', 'CBCAS')
            if re.search(r'(?<![A-Za-z])' + x + r'(?![A-Za-z])', str(ref or ''), re.I)]
    return hits[0] if len(hits) == 1 else 'unassigned'


def key(x):
    s = str(x or '').strip()
    return re.sub(r'\.t\d+$', '', s.split()[0]) if s else ''


def read(path):
    with open(path, newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))


def best_hits(path):
    """Per gene: the tBLASTn hit with the highest bitscore (same choice as 03g)."""
    best = {}
    for r in read(path):
        if r['overlapping_gene_ids'] == 'NO_GENE_MODEL':
            continue
        for g in re.split(r'[,;\s]+', r['overlapping_gene_ids']):
            g = key(g)
            if g and (g not in best or float(r['bitscore']) > float(best[g]['bitscore'])):
                best[g] = r
    return best


def main():
    w = csv.DictWriter(sys.stdout, fieldnames=COLS, delimiter='\t', lineterminator='\n',
                       extrasaction='ignore')
    w.writeheader()
    nomodel = {}
    try:
        for r in read(NOMODEL):
            nomodel.setdefault(r['Haplotype'], []).append(r)
    except FileNotFoundError:
        print(f'[warn] {NOMODEL} not found; the unclassified regions are omitted', file=sys.stderr)

    n_cls = n_nm = 0
    for h, chemo in HAPS:
        hits = best_hits(GM.format(h=h))
        orf = {key(r['Gene_ID']): r for r in read(ORF.format(h=h))}
        gen = {key(r['Gene_ID']): r for r in read(GEN.format(h=h))}
        rows = read(CLS.format(h=h))
        rows.sort(key=lambda r: (r['Gene_ID'].split('.')[1],
                                 int(hits[key(r['Gene_ID'])]['hit_start']) if key(r['Gene_ID']) in hits else 0))
        for r in rows:
            g = key(r['Gene_ID'])
            hit, o, ge = hits.get(g), orf.get(g, {}), gen.get(g, {})
            src = r.get('L2_L6_source', '')
            note = []
            if src == 'genome_tblastn':
                note.append('L2/L6 from genomic tBLASTn hit')
            elif src == 'genome_tblastn_stop':
                note.append('L6 from premature stop in genomic hit')
            nm = int(ge['n_gene_models_on_hsp']) if ge.get('n_gene_models_on_hsp', '').isdigit() else ''
            if isinstance(nm, int) and nm > 1:
                note.append(f'best tBLASTn hit shared by {nm} gene models')
            w.writerow({
                'Haplotype': h, 'Chemotype': chemo, 'Gene_ID': r['Gene_ID'],
                'Chromosome': hit['sseqid'].split('.')[-1] if hit else '',
                'Start_Mb': f"{int(hit['hit_start'])/1e6:.3f}" if hit else '',
                'End_Mb': f"{int(hit['hit_end'])/1e6:.3f}" if hit else '',
                'Strand': hit['hit_strand'] if hit else '',
                'Series': series(r['Best_reference']),
                'Nearest_reference_identity_pct': f"{float(hit['pident']):.1f}" if hit else '',
                'Protein_length_aa': o.get('Ungapped_length_aa', ''),
                'Reference_coverage': o.get('Reference_coverage', ''),
                **{f'L{i}': MARK.get(r[f'L{i}'], r[f'L{i}']) for i in range(1, 8)},
                'Classification': r['Classification'] or r['Review_status'],
                'L2_L6_evidence_source': src, 'Clade_assignment': r.get('Clade_assignment', ''),
                'Disruption': r.get('Disruption', ''), 'Gene_models_on_same_hit': nm,
                'Note': '; '.join(note),
            })
            n_cls += 1
        for r in sorted(nomodel.get(h, []), key=lambda x: (x['Chromosome'], int(x['Start']))):
            w.writerow({
                'Haplotype': h, 'Chemotype': chemo, 'Gene_ID': '(no gene model)',
                'Chromosome': r['Chromosome'].split('.')[-1],
                'Start_Mb': f"{int(r['Start'])/1e6:.3f}", 'End_Mb': f"{int(r['End'])/1e6:.3f}",
                'Strand': r['Strand'], 'Series': r['Nearest_reference'],
                'Nearest_reference_identity_pct': r['Identity_pct'],
                'Protein_length_aa': '', 'Reference_coverage': '',
                **{f'L{i}': '—' for i in range(1, 8)},
                'Classification': 'not classified', 'L2_L6_evidence_source': '',
                'Clade_assignment': '', 'Disruption': '', 'Gene_models_on_same_hit': '',
                'Note': (f"no overlapping GFF3 gene model; tBLASTn {r['Ref_range_aa']} aa "
                         f"({r['Aln_length_aa']} aa); THCAS {r['Identity_THCAS'] or '-'}%, "
                         f"CBDAS {r['Identity_CBDAS'] or '-'}%, CBCAS {r['Identity_CBCAS'] or '-'}%"),
            })
            n_nm += 1
    print(f'[AF1] classified copies: {n_cls}  unclassified regions: {n_nm}  total rows: {n_cls + n_nm}',
          file=sys.stderr)


if __name__ == '__main__':
    main()
