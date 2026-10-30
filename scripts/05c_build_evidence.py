#!/usr/bin/env python3
"""05c_build_evidence.py - assemble the evidence table consumed by 06_assign_L1_L7.py.

Scope: pure join + column renaming. NO classification decisions are made here.
  L1  : 'yes' if the sequence reached Step 3 (present in 03_orf.tsv)
  L2  : 03_orf.tsv  L2_full_length_ORF_evidence
  L3  : 04_functional_residues.tsv  L3_FAD_evidence
  L4  : 04_functional_residues.tsv  L4_active_site_evidence
  L5  : 05b output (l5.tsv)
  L6  : 03_orf.tsv  L6_ORF_integrity_evidence
  L7  : 02_gene_models.tsv has this gene with coordinates on the expected haplotype -> yes, else '—'
Missing inputs never raise: the value becomes '—' and is counted in the diagnostics.

No derived values. Partial-length rows keep L6='—'; their tier is decided in 06
(definition 2026-10-12: L2=no, L6='—', no confirmed disruption -> tier 2 '-like copy').

L2/L6 source (decision F, 2026-10-12): if --genomic-orf is given and the genome is cleanly
full-length (03g: L2_genomic = L6_genomic = yes), L2 and L6 are taken from the genome;
if 03g reports a premature stop inside a clean full-length HSP (L6_genomic = no,
Genomic_disruption = premature_stop), L6 = no and Disruption = premature_stop
(decision 2026-10-12 23:55); otherwise the annotation values from 03 are used.
Recorded in L2_L6_source. This is the
priority rule written in the memo (section H), not a new judgement.

Series naming (definition 2026-10-12): the series label is ALWAYS the BLAST nearest
reference (Best_reference), derived by 06. L5 is used only as a tier criterion.
The 05b phylogenetic assignment is recorded in Clade_assignment for reporting only;
Series_assignment is NOT emitted, so 06's override mechanism is never triggered.
"""
import argparse
import csv
import re
import sys
from collections import Counter

MISSING = '—'
OUT_COLS = ['Gene_ID', 'Sequence', 'Best_reference', 'L1', 'L2', 'L3', 'L4', 'L5', 'L6', 'L7',
            'Disruption', 'Step2_hit_count', 'L5_evidence', 'L6_scope',
            'Clade_assignment', 'L2_L6_source']

ALIASES = {
    'id':     ['Sequence', 'Gene_ID', 'Taxon', 'Tip', 'Leaf', 'Query', 'qseqid'],
    'L3':     ['L3', 'L3_FAD_evidence', 'L3_catalytic_residue_evidence'],
    'L4':     ['L4', 'L4_active_site_evidence', 'L4_functional_residue_evidence'],
    'L5':     ['L5', 'L5_call', 'L5_clade_evidence'],
    'L5_ev':  ['L5_evidence', 'Clade_evidence', 'Evidence', 'Note', 'Reason'],
    'series': ['Assigned_reference', 'Assigned_series', 'Series_assignment', 'Series',
               'Reference_series', 'Clade_reference', 'Assignment'],
}


def norm(x):
    """Gene-level key: first token, drop trailing transcript suffix (.t1)."""
    s = str(x or '').strip().split()[0] if str(x or '').strip() else ''
    return re.sub(r'\.t\d+$', '', s)


def tri(x):
    s = str(x or '').strip().lower()
    if s in ('yes', 'y', 'true', '1'):
        return 'yes'
    if s in ('no', 'n', 'false', '0'):
        return 'no'
    return MISSING


def read_tsv(path):
    with open(path, newline='') as f:
        r = csv.DictReader(f, delimiter='\t')
        return r.fieldnames or [], list(r)


def pick(fields, key, path, required=True):
    for a in ALIASES[key]:
        if a in fields:
            return a
    if required:
        sys.exit(f'[05c] {path}: no column for {key}; tried {ALIASES[key]}; have {fields}')
    return None


def index_by_id(path, label):
    fields, rows = read_tsv(path)
    idc = pick(fields, 'id', path)
    idx = {}
    for r in rows:
        k = norm(r.get(idc))
        if k:
            idx.setdefault(k, r)
    print(f'[05c] {label}: {len(idx)} ids (key column: {idc})', file=sys.stderr)
    return fields, idx


def blastp_counts(path):
    c = Counter()
    with open(path) as f:
        for line in f:
            p = line.rstrip('\n').split('\t')
            if len(p) < 2 or p[0] == 'qseqid':
                continue
            c[norm(p[1])] += 1
    return c


def blast_series(ref):
    hits = [x for x in ('THCAS', 'CBDAS', 'CBCAS') if re.search(r'(?<![a-z])' + x + r'(?![a-z])', str(ref or ''), re.I)]
    return hits[0] if len(hits) == 1 else 'unassigned'


def gene_models(path, haplotype):
    fields, rows = read_tsv(path)
    for col in ('overlapping_gene_ids', 'hit_start', 'hit_end'):
        if col not in fields:
            sys.exit(f'[05c] {path}: missing column {col}; have {fields}')
    located = set()
    for r in rows:
        try:
            int(float(r['hit_start'])); int(float(r['hit_end']))
        except (TypeError, ValueError):
            continue
        for g in re.split(r'[,;\s]+', r.get('overlapping_gene_ids') or ''):
            k = norm(g)
            if k and (not haplotype or k.split('.')[0] == haplotype):
                located.add(k)
    return located


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--orf', required=True, help='03_orf.tsv')
    ap.add_argument('--residues', required=True, help='04_functional_residues.tsv')
    ap.add_argument('--l5', required=True, help='05b output (l5.tsv)')
    ap.add_argument('--gene-models', required=True, help='02_gene_models.tsv')
    ap.add_argument('--blastp', required=True, help='01_search/blastp_hits.tsv (outfmt 6)')
    ap.add_argument('--haplotype', default='', help='e.g. SAN2a; L7 requires the gene ID prefix to match')
    ap.add_argument('--genomic-orf', default='', help='03g output (optional; decision F)')
    ap.add_argument('--out', required=True)
    a = ap.parse_args()

    orf_fields, orf_rows = read_tsv(a.orf)
    for col in ('Gene_ID', 'Sequence', 'Best_reference', 'L2_full_length_ORF_evidence',
                'L6_ORF_integrity_evidence', 'Disruption'):
        if col not in orf_fields:
            sys.exit(f'[05c] {a.orf}: missing column {col}')

    res_fields, res = index_by_id(a.residues, '04 residues')
    c3, c4 = pick(res_fields, 'L3', a.residues), pick(res_fields, 'L4', a.residues)
    l5_fields, l5 = index_by_id(a.l5, '05b L5')
    c5 = pick(l5_fields, 'L5', a.l5)
    c5ev = pick(l5_fields, 'L5_ev', a.l5, required=False)
    cser = pick(l5_fields, 'series', a.l5, required=False)
    print(f'[05c] 05b columns: L5={c5} evidence={c5ev} series={cser}', file=sys.stderr)

    genomic, gstop = {}, {}
    if a.genomic_orf:
        with open(a.genomic_orf, newline='') as f:
            for g in csv.DictReader(f, delimiter='\t'):
                if g.get('L2_genomic') == 'yes' and g.get('L6_genomic') == 'yes':
                    genomic[norm(g['Gene_ID'])] = g
                elif g.get('L6_genomic') == 'no' and g.get('Genomic_disruption') == 'premature_stop':
                    gstop[norm(g['Gene_ID'])] = g

    hits = blastp_counts(a.blastp)
    located = gene_models(a.gene_models, a.haplotype)

    out, miss, mismatch = [], Counter(), []
    for r in orf_rows:
        keys = [norm(r['Sequence']), norm(r['Gene_ID'])]
        rr = next((res[k] for k in keys if k in res), None)
        ll = next((l5[k] for k in keys if k in l5), None)
        if rr is None: miss['04'] += 1
        if ll is None: miss['05b'] += 1

        row = {
            'Gene_ID': r['Gene_ID'], 'Sequence': r['Sequence'], 'Best_reference': r['Best_reference'],
            'L1': 'yes',
            'L2': tri(r['L2_full_length_ORF_evidence']),
            'L3': tri(rr.get(c3)) if rr else MISSING,
            'L4': tri(rr.get(c4)) if rr else MISSING,
            'L5': tri(ll.get(c5)) if ll else MISSING,
            'L6': tri(r['L6_ORF_integrity_evidence']),
            'L7': 'yes' if any(k in located for k in keys) else MISSING,
            'Disruption': r.get('Disruption') or 'not_evaluated',
            'Step2_hit_count': max(hits.get(k, 0) for k in keys),
            'L5_evidence': (ll.get(c5ev) if (ll and c5ev) else '') or '',
            'L6_scope': 'predicted_protein_only',
            'Clade_assignment': '', 'L2_L6_source': 'annotation',
        }
        if row['L7'] == MISSING: miss['L7'] += 1
        if any(k in genomic for k in keys):
            row['L2'], row['L6'], row['L2_L6_source'] = 'yes', 'yes', 'genome_tblastn'
        elif any(k in gstop for k in keys):
            row['L6'], row['Disruption'], row['L2_L6_source'] = 'no', 'premature_stop', 'genome_tblastn_stop'

        if ll and cser:
            row['Clade_assignment'] = (ll.get(cser) or '').strip()
            c, b = row['Clade_assignment'], blast_series(r['Best_reference'])
            if c and c != 'unassigned' and c != b:
                mismatch.append((r['Gene_ID'], b, c))

        out.append(row)

    with open(a.out, 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=OUT_COLS, delimiter='\t')
        w.writeheader()
        w.writerows(out)

    print(f'[05c] wrote {len(out)} rows -> {a.out}', file=sys.stderr)
    print(f'[05c] unmatched: 04={miss["04"]} 05b={miss["05b"]} L7_unlocated={miss["L7"]}', file=sys.stderr)
    for L in ('L2', 'L3', 'L4', 'L5', 'L6', 'L7'):
        print(f'[05c] {L}: {dict(Counter(x[L] for x in out))}', file=sys.stderr)
    print(f'[05c] L2_L6_source: {dict(Counter(x["L2_L6_source"] for x in out))}', file=sys.stderr)
    print(f'[05c] Clade_assignment: {dict(Counter(x["Clade_assignment"] or MISSING for x in out))}', file=sys.stderr)
    print(f'[05c] BLAST series vs clade mismatch (L5=yes only): {len(mismatch)}', file=sys.stderr)
    for g, b, c in mismatch:
        print(f'[05c]   {g}\tBLAST={b}\tclade={c}', file=sys.stderr)


if __name__ == '__main__':
    main()
