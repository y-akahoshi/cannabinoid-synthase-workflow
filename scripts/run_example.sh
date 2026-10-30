#!/usr/bin/env bash
# run_example.sh - run the whole workflow on the small example dataset.
#
# The example dataset is a 1.2 Mb window around the synthase cluster of one
# haplotype, with every coordinate re-based to 1 (see example_data/README_offset.txt).
# It exercises the same code path as the full analysis and finishes in a few
# minutes on a laptop.
#
#     conda env create -f environment.yml
#     conda activate plant_genomics
#     bash scripts/run_example.sh
#
# Expected result: 06_classification.tsv under example_out/SAN2b/, in which
# g313950 is classified as a Pseudogenized THCAS-like copy on the evidence of a
# premature stop codon in the genomic hit.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$HERE/.." && pwd)"

H="${1:-SAN2b}"
ROOT="${ROOT:-$REPO/example_data}"
OUTBASE="${OUTBASE:-$REPO/example_out}"
THREADS="${THREADS:-2}"

G="$ROOT/00_genomes/example/$H"
W="$OUTBASE/work/$H"
REFS="$ROOT/refs/synthase_query.fasta"

die() { echo "ERROR: $*" >&2; exit 1; }
log() { echo "[example] $*"; }

for t in blastp tblastn makeblastdb mafft python3; do
  command -v "$t" >/dev/null || die "$t not found (conda activate plant_genomics)"
done
command -v iqtree2 >/dev/null || command -v iqtree >/dev/null \
  || die "iqtree not found (conda activate plant_genomics)"

[[ -d "$G" ]] || die "example dataset not found at $G
Build it first with:
  python3 $HERE/make_example_dataset.py --out $ROOT"

# Resolved with a plain loop: a failing `ls` inside $( ) would be fatal here,
# because the script runs under `set -e` with `pipefail`.
GENOME=""
for f in "$G/$H.unmasked.fasta.gz" "$G/$H.softmasked.fasta.gz"; do
  if [[ -s "$f" ]]; then GENOME="$f"; break; fi
done
[[ -n "$GENOME" ]] || die "no genome FASTA (unmasked or softmasked) in $G"

# The example always starts from a clean state, so that re-running it is safe:
# IQ-TREE refuses to overwrite the results of a previous run.
if [[ -d "$OUTBASE" ]]; then
  case "$OUTBASE" in
    */example_out) log "removing previous output in $OUTBASE"; rm -rf "$OUTBASE" ;;
    *) die "refusing to delete $OUTBASE (expected a path ending in example_out)" ;;
  esac
fi

mkdir -p "$W"

# --- 01: BLASTP against the predicted proteins, tBLASTn against the genome ---
log "01  searching with the three reference synthases"
bash "$HERE/01_search_synthase_candidates.sh" \
  --refs "$REFS" \
  --proteins "$G/$H.all.proteins.fasta.gz" \
  --genome "$GENOME" \
  --outdir "$W/01_search" \
  --threads "$THREADS"

# --- 02: map the genomic hits onto the GFF3 gene models -----------------------
log "02  matching genomic hits to GFF3 gene models"
python3 "$HERE/02_extract_gene_models.py" \
  --hits "$W/01_search/tblastn_hits.tsv" \
  --gff "$G/$H.high_confidence.gff3.gz" \
  --out "$W/02_gene_models.tsv"

# --- 03-06: evaluation and classification ------------------------------------
log "03-06  evaluating L1-L7 and classifying"
ROOT="$ROOT" OUTBASE="$OUTBASE" WORKDIR="$W" THREADS="$THREADS" \
  bash "$HERE/run_haplotype.sh" "$H"

OUT="$OUTBASE/$H/06_classification.tsv"
[[ -s "$OUT" ]] || die "no classification written to $OUT"

echo
log "done: $OUT"
echo
column -t -s$'\t' "$OUT" | head -15
echo
log "classification counts:"
awk -F'\t' 'NR==1 {for (i=1;i<=NF;i++) if ($i=="Classification") c=i; next}
           c {n[$c]++} END {for (k in n) printf "  %-34s %s\n", k, n[k]}' "$OUT"
