#!/usr/bin/env python3
"""Operational 5-category classifier + reproducible reference-series assignment.

Series is derived from the best BLAST reference supplied by Step 3, not from
chromosomal position. Optional Series_assignment='unassigned' with
Series_evidence documents ambiguous/conflicting phylogenetic findings.
L6 from protein-only Step 3 is not a genome-validated lesion; use an explicit
confirmed Disruption to activate pseudogenization rule A.
"""
import argparse
import csv
import re
from pathlib import Path

REFS=('THCAS','CBDAS','CBCAS')
YES={'yes','y','1','true','pass','✓'}
NO={'no','n','0','false','fail','✗'}
NA={'na','n/a','.','-','—',''}
CONFIRMED={'premature_stop','internal_stop','confirmed_frameshift','frameshift','confirmed_orf_disruption','confirmed_coding_frame_disruption','confirmed_coding_capacity_loss'}
NEUTRAL={'split_cds','multi_exon','inframe_deletion','active_site_substitution','phase_warning','phase_inconsistent','no_lesion_detected','no_internal_stop_in_full_length_like_protein','not_evaluated'}
UNCERTAIN={'orf_integrity_unresolved','structural_lesion_unresolved','large_deletion_unresolved','splice_structure_unresolved','coding_capacity_unresolved'}

def norm(v):
    v=str(v or '').strip().lower()
    if v in YES:return 'yes'
    if v in NO:return 'no'
    if v in NA:return '—'
    raise ValueError('Invalid criterion state: '+repr(v))

def tokens(raw):
    return set(filter(None,re.split(r'[,;|\s]+',str(raw or '').strip().lower())))

def reference_series(ref):
    matches=[x for x in REFS if re.search(r'(?<![a-z])'+x+r'(?![a-z])',str(ref or ''),re.I)]
    return matches[0] if len(matches)==1 else 'unassigned'

def classify(r):
    v={f'L{i}':norm(r.get(f'L{i}')) for i in range(1,8)}
    flags=tokens(r.get('Disruption'))
    unknown=flags-CONFIRMED-NEUTRAL-UNCERTAIN
    series=reference_series(r.get('Best_reference'))
    requested=str(r.get('Series_assignment') or '').strip()
    reason=str(r.get('Series_evidence') or '').strip()
    if requested:
        if requested not in (*REFS,'unassigned'):
            raise ValueError('Invalid Series_assignment '+repr(requested))
        if requested != series and not reason:
            return series,'', 'review_series_override_lacks_evidence','none'
        series=requested
    if unknown:return series,'','review_unknown_disruption_token:'+','.join(sorted(unknown)),'none'
    if v['L1']!='yes':return series,'','review_not_qualified_homology_hit','none'
    has_direct=bool(flags&CONFIRMED)
    if (has_direct and v['L6']!='no') or (v['L6']=='no' and not has_direct):
        return series,'','review_L6_disruption_evidence_conflict','none'
    a=(v['L6']=='no' and has_direct)
    b=(v['L2']==v['L3']==v['L4']=='no')
    rule='A+B' if a and b else ('A' if a else 'B' if b else 'none')
    if a or b:
        return series,f'Pseudogenized {series if series!="unassigned" else "THCAS/CBDAS/CBCAS"}-like copy','classified',rule
    if flags&UNCERTAIN:return series,'','review_structural_evidence_unresolved',rule
    if series=='unassigned':
        # Keep non-series-specific category when functional cannot be resolved.
        if all(x=='yes' for x in v.values()):return series,'','review_series_unassigned',rule
    if all(x=='yes' for x in v.values()):
        return series,f'Functional {series} candidate','classified',rule
    if v['L2']=='yes' and v['L6']=='yes' and (v['L4']=='no' or v['L5']=='no'):
        label=series if series!='unassigned' else 'THCAS/CBDAS/CBCAS'
        return series,f'Divergent {label}-like copy','classified',rule
    # 2026-10-12 definition: partial-length (L2=no), L6 not evaluable,
    # no confirmed disruption -> tier 2 '-like copy' (Pseudogenized requires demonstrated disruption)
    if v['L2']=='no' and v['L6']=='—' and not has_direct:
        label=series if series!='unassigned' else 'THCAS/CBDAS/CBCAS'
        return series,f'{label}-like copy','classified',rule
    if v['L6']=='—':
        return series,'','review_L6_not_evaluable',rule
    label=series if series!='unassigned' else 'THCAS/CBDAS/CBCAS'
    return series,f'{label}-like copy','classified',rule

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--evidence',required=True)
    p.add_argument('--out',required=True)
    a=p.parse_args()
    with open(a.evidence,encoding='utf-8',newline='') as fh:
        rd=csv.DictReader(fh,delimiter='\t')
        needed={'Sequence','Best_reference'}|{f'L{i}' for i in range(1,8)}
        if not needed.issubset(rd.fieldnames or []):p.error('Missing columns: '+','.join(sorted(needed-set(rd.fieldnames or []))))
        fields=[x for x in rd.fieldnames if x not in {'Reference_series','Classification','Review_status','Pseudogenized_rule'}]+['Reference_series','Classification','Review_status','Pseudogenized_rule']
        rows=[]
        for row in rd:
            for i in range(1,8):row[f'L{i}']=norm(row[f'L{i}'])
            series,classification,status,rule=classify(row)
            row.update(Reference_series=series,Classification=classification,Review_status=status,Pseudogenized_rule=rule)
            rows.append(row)
    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    with open(a.out,'w',encoding='utf-8',newline='') as fh:
        wr=csv.DictWriter(fh,delimiter='\t',fieldnames=fields,lineterminator='\n')
        wr.writeheader();wr.writerows(rows)
    print('rows:',len(rows),'classified:',sum(x['Review_status']=='classified' for x in rows),'review:',sum(x['Review_status']!='classified' for x in rows))

if __name__=='__main__':main()
