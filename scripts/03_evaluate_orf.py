#!/usr/bin/env python3
"""Simple, annotation-only Step 3. No genome rescue, locus merging or HSP ORF calls.
L2 is full-length-like sequence coverage, not verified biochemical function.
L6 reports only the selected predicted protein: confirmed internal '*' => no;
full-length-like and no internal '*' => yes (limited annotated-protein screening);
otherwise unresolved. This is NOT genomic locus-wide ORF verification.
"""
import argparse
import csv
import gzip
import subprocess
import sys
import tempfile
from collections import Counter
from pathlib import Path
from Bio import Align, SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord


def fasta(path):
    with (gzip.open(path, 'rt') if str(path).endswith('.gz') else open(path)) as f:
        return list(SeqIO.parse(f, 'fasta'))


def candidates(step2, min_ident, min_len, min_qcov):
    best = {}
    with open(step2, newline='') as f:
        rows = csv.DictReader(f, delimiter='\t')
        needed={'qseqid','pident','length','bitscore','qcovs','overlapping_gene_ids'}
        if not needed.issubset(rows.fieldnames or []):
            raise ValueError('Step 2 lacks fields: ' + str(sorted(needed-set(rows.fieldnames or []))))
        for r in rows:
            try:
                ident=float(r['pident']); length=float(r['length']); qcov=float(r['qcovs']); bits=float(r['bitscore'])
            except (ValueError,TypeError):
                continue
            if ident<min_ident or length<min_len or qcov<min_qcov:
                continue
            for gene in r['overlapping_gene_ids'].split(';'):
                gene=gene.strip()
                if not gene or gene=='NO_GENE_MODEL': continue
                if gene not in best or (bits, r['qseqid'])>(best[gene][0], best[gene][1]):
                    best[gene]=(bits,r['qseqid'])
    return best


def score(seq, ref, window, cov_cut):
    aligner=Align.PairwiseAligner(mode='global')
    aligner.match_score=2; aligner.mismatch_score=-1
    aligner.open_gap_score=-8; aligner.extend_gap_score=-0.5
    a=aligner.align(ref,seq)[0]
    present=set()
    for (s,e),_ in zip(*a.aligned):
        present.update(range(int(s),int(e)))
    n=len(ref); w=min(window,n)
    cov=len(present)/n
    ncv=sum(x in present for x in range(w))/w
    ccv=sum(x in present for x in range(n-w,n))/w
    full=(cov>=cov_cut and ncv>=0.5 and ccv>=0.5 and bool(seq))
    return full,cov,ncv,ccv


def run(args):
    best=candidates(args.gene_models,args.min_pident,args.min_aln_length,args.min_qcov)
    refs={r.id:str(r.seq).rstrip('*') for r in fasta(args.references)}
    proteins=fasta(args.proteins)
    for _,refid in best.values():
        if refid not in refs: raise ValueError('Missing reference ID: '+refid)
    selected=[]; outrows=[]; missing=[]
    for gene in sorted(best):
        refid=best[gene][1]
        opts=[r for r in proteins if r.id==gene or r.id.startswith((gene+'.',gene+'-',gene+'_'))]
        if not opts:
            missing.append(gene)
            outrows.append(dict(Gene_ID=gene, Sequence='NOT_FOUND', Best_reference=refid,
                Isoform_count=0, Selection_rule='not_found', Ungapped_length_aa='',
                Reference_coverage='',N_terminal_coverage='',C_terminal_coverage='',
                Internal_stop='—', L2_full_length_ORF_evidence='—',
                L6_ORF_integrity_evidence='—', Disruption='not_evaluated', QC_note='protein_not_found'))
            continue
        candidates_scored=[]
        for rec in opts:
            raw=str(rec.seq).upper().strip()
            internal='*' in raw.rstrip('*')
            clean=raw.rstrip('*').replace('*','X')
            full,cov,ncv,ccv=score(clean,refs[refid],args.terminal_window,args.min_reference_coverage)
            full=full and not internal
            candidates_scored.append((rec,clean,internal,full,cov,ncv,ccv))
        # Predeclared stable selection rule; no gene-specific exceptions.
        candidates_scored.sort(key=lambda x:(-int(x[3]),-x[4],-min(x[5],x[6]),-len(x[1]),x[0].id))
        rec,clean,internal,full,cov,ncv,ccv=candidates_scored[0]
        selected.append(SeqRecord(Seq(clean),id=rec.id,description=''))
        l6='no' if internal else ('yes' if full else '—')
        outrows.append(dict(Gene_ID=gene,Sequence=rec.id,Best_reference=refid,
            Isoform_count=len(opts),Selection_rule='L2_then_coverage_then_terminals_then_length_then_ID',
            Ungapped_length_aa=len(clean),Reference_coverage=f'{cov:.4f}',
            N_terminal_coverage=f'{ncv:.4f}',C_terminal_coverage=f'{ccv:.4f}',
            Internal_stop='yes' if internal else 'no',L2_full_length_ORF_evidence='yes' if full else 'no',
            L6_ORF_integrity_evidence=l6,
            Disruption='internal_stop' if internal else ('no_internal_stop_in_full_length_like_protein' if full else 'not_evaluated'),
            QC_note='annotation_only;no_locus_level_ORF_inference'))
    if not selected: raise ValueError('No matched proteins')
    for filename in [args.out_fasta,args.alignment,args.orf_out]:
        Path(filename).parent.mkdir(parents=True,exist_ok=True)
    SeqIO.write(selected,args.out_fasta,'fasta')
    with open(args.orf_out,'w',newline='') as f:
        wr=csv.DictWriter(f,fieldnames=list(outrows[0]),delimiter='\t',lineterminator='\n')
        wr.writeheader();wr.writerows(outrows)
    if args.skip_alignment:
        # Downstream Step4 requires a real alignment: users must explicitly opt out.
        print('MAFFT alignment intentionally skipped')
    else:
        combined=[SeqRecord(Seq(s),id=rid,description='') for rid,s in refs.items()]+selected
        if args.extra_fasta: combined+=fasta(args.extra_fasta)
        ids=[r.id for r in combined]
        if len(ids)!=len(set(ids)): raise ValueError('Duplicate FASTA IDs in MAFFT input')
        with tempfile.TemporaryDirectory() as d:
            inputfile=Path(d)/'input.fasta'
            SeqIO.write(combined,str(inputfile),'fasta')
            with open(args.alignment,'w') as output:
                subprocess.run([args.mafft,'--auto','--thread',str(args.threads),str(inputfile)],check=True,stdout=output)
    print('genes:',len(outrows),'proteins:',len(selected),'missing:',len(missing))
    for field in ['L2_full_length_ORF_evidence','L6_ORF_integrity_evidence']:
        print(field,dict(Counter(r[field] for r in outrows)))
    if missing: print('Missing Gene IDs:',','.join(missing[:20]),file=sys.stderr)


def main():
    p=argparse.ArgumentParser(description='Minimal Step 3: representative predicted protein, L2 and limited L6 QC')
    for arg in ['gene-models','proteins','references','out-fasta','alignment','orf-out']:
        p.add_argument('--'+arg,required=True)
    p.add_argument('--extra-fasta')
    p.add_argument('--min-pident',type=float,default=0)
    p.add_argument('--min-aln-length',type=float,default=0)
    p.add_argument('--min-qcov',type=float,default=0)
    p.add_argument('--min-reference-coverage',type=float,default=.90)
    p.add_argument('--terminal-window',type=int,default=20)
    p.add_argument('--threads',type=int,default=4)
    p.add_argument('--mafft',default='mafft')
    p.add_argument('--skip-alignment',action='store_true',help='For lightweight QC tests only')
    run(p.parse_args())

if __name__=='__main__': main()
