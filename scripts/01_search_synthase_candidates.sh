#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'USAGE'
Usage:
  01_search_synthase_candidates.sh \
    --refs refs.fasta \
    --proteins proteins.fasta[.gz] \
    --genome genome.fasta[.gz] \
    --outdir results/search \
    [--evalue 1e-10] [--threads 4]

Purpose:
  Search predicted proteins (BLASTP) and the genome assembly (tBLASTn)
  with THCAS/CBDAS/CBCAS reference proteins. Genome-level searching is
  retained because complete loci may be absent or split in protein annotations.

Outputs:
  blastp_hits.tsv
  tblastn_hits.tsv
  metadata.txt

Output columns:
  qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen qcovs
USAGE
}

REFS=""; PROTEINS=""; GENOME=""; OUTDIR=""; EVALUE="1e-10"; THREADS="4"
while [[ $# -gt 0 ]]; do
  case "$1" in
    --refs) REFS="$2"; shift 2 ;;
    --proteins) PROTEINS="$2"; shift 2 ;;
    --genome) GENOME="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --evalue) EVALUE="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "Unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

for x in "$REFS" "$PROTEINS" "$GENOME" "$OUTDIR"; do
  [[ -n "$x" ]] || { echo "Missing required argument." >&2; usage >&2; exit 2; }
done
for cmd in makeblastdb blastp tblastn; do
  command -v "$cmd" >/dev/null 2>&1 || { echo "Required command not found: $cmd" >&2; exit 127; }
done

mkdir -p "$OUTDIR"
TMPDIR_LOCAL="$(mktemp -d)"
trap 'rm -rf "$TMPDIR_LOCAL"' EXIT

materialize_fasta() {
  local src="$1" dst="$2"
  if [[ "$src" == *.gz ]]; then
    gzip -dc "$src" > "$dst"
  else
    cp "$src" "$dst"
  fi
}

PROT_FA="$TMPDIR_LOCAL/proteins.fa"
GENOME_FA="$TMPDIR_LOCAL/genome.fa"
materialize_fasta "$PROTEINS" "$PROT_FA"
materialize_fasta "$GENOME" "$GENOME_FA"

makeblastdb -in "$PROT_FA" -dbtype prot -out "$TMPDIR_LOCAL/protdb" >/dev/null
makeblastdb -in "$GENOME_FA" -dbtype nucl -out "$TMPDIR_LOCAL/nucdb" >/dev/null

FMT='6 qseqid sseqid pident length qstart qend sstart send evalue bitscore qlen slen qcovs'
blastp \
  -query "$REFS" -db "$TMPDIR_LOCAL/protdb" \
  -evalue "$EVALUE" -max_target_seqs 500 -num_threads "$THREADS" \
  -outfmt "$FMT" > "$OUTDIR/blastp_hits.tsv"

tblastn \
  -query "$REFS" -db "$TMPDIR_LOCAL/nucdb" \
  -evalue "$EVALUE" -max_target_seqs 1000 -num_threads "$THREADS" \
  -outfmt "$FMT" > "$OUTDIR/tblastn_hits.tsv"

{
  echo -e "parameter\tvalue"
  echo -e "reference_fasta\t$REFS"
  echo -e "protein_fasta\t$PROTEINS"
  echo -e "genome_fasta\t$GENOME"
  echo -e "evalue\t$EVALUE"
  echo -e "threads\t$THREADS"
  echo -e "blastp_version\t$(blastp -version | head -n1)"
  echo -e "tblastn_version\t$(tblastn -version | head -n1)"
} > "$OUTDIR/metadata.txt"

echo "Wrote $OUTDIR/blastp_hits.tsv and $OUTDIR/tblastn_hits.tsv" >&2
