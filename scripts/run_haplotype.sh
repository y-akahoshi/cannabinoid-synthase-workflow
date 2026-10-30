#!/usr/bin/env bash
# run_haplotype.sh - Steps 03 -> 04 -> 05 (IQ-TREE, rooted on BBE) -> 05b -> 05c -> 06 for one haplotype.
#
# Usage:  conda activate plant_genomics
#         bash run_haplotype.sh SAN2a
#
# Tree rooting follows Methods 6: E. californica BBE (P30986) is added to the alignment
# and passed to IQ-TREE with -o. 05b does not reroot, so the tree must be written rooted.
# Output: $OUTBASE/<haplotype>/ (default ~/cbdas_project/results/<haplotype>/)
set -euo pipefail

H="${1:?usage: run_haplotype.sh <haplotype, e.g. SAN2a>}"
ROOT="${ROOT:-$HOME/cbdas_project}"
S="${SCRIPTS:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"
W="${WORKDIR:-$ROOT/work/$H}"
OUT="${OUTBASE:-$ROOT/results}/$H"
REFS="${REFS:-$ROOT/refs/synthase_query.fasta}"
BBE="${BBE:-$ROOT/refs/BBE_EcaBBE_P30986.fasta}"
THCAS_ID="${THCAS_ID:-sp|Q8GTB6.1|THCAS_CANSA}"
CBDAS_ID="${CBDAS_ID:-sp|A6P6V9|CBDAS_CANSA}"
REFMAP="${REFMAP:-THCAS=Q8GTB6,CBDAS=A6P6V9,CBCAS=Q33DQ2}"
SEED="${SEED:-12345}"

die() { echo "[$H] ERROR: $*" >&2; exit 1; }
log() { echo "[$H] $*"; }

for f in "$W/02_gene_models.tsv" "$W/01_search/blastp_hits.tsv" "$REFS" "$BBE"; do
  [[ -s "$f" ]] || die "missing input: $f"
done
PROT_GZ=$(find "$ROOT/00_genomes" -name "$H.all.proteins.fasta.gz" | head -1)
[[ -n "$PROT_GZ" ]] || die "$H.all.proteins.fasta.gz not found under 00_genomes"

mkdir -p "$OUT"; cd "$OUT"
[[ -s proteins.fasta ]] || zcat "$PROT_GZ" > proteins.fasta
BBE_ID=$(grep -m1 '^>' "$BBE" | sed 's/^>//; s/[[:space:]].*//')
log "proteins: $PROT_GZ"
log "references: $REFS ($(grep -c '^>' "$REFS") seqs)  outgroup: $BBE_ID"

# --- Step 03: representative protein, L2/L6, alignment (BBE added via --extra-fasta)
python "$S/03_evaluate_orf.py" \
  --gene-models "$W/02_gene_models.tsv" --proteins proteins.fasta \
  --references "$REFS" --extra-fasta "$BBE" \
  --out-fasta 03_candidates.fasta --alignment 03_candidates.aligned.fasta \
  --orf-out 03_orf.tsv
grep -qF ">$BBE_ID" 03_candidates.aligned.fasta || die "BBE not in alignment"
for id in "$THCAS_ID" "$CBDAS_ID"; do
  grep -qF ">$id" 03_candidates.aligned.fasta || {
    echo "non-$H headers in alignment:"; grep '^>' 03_candidates.aligned.fasta | grep -v "$H" || true
    die "anchor $id not in alignment"; }
done
if grep -q "BBE" 03_orf.tsv; then log "WARNING: BBE appears in 03_orf.tsv (will be classified as a query)"; fi

# --- Step 04: FAD (L3) and active-site (L4) residues
python "$S/04_check_functional_residues.py" --alignment 03_candidates.aligned.fasta \
  --thcas-id "$THCAS_ID" --cbdas-id "$CBDAS_ID" --out 04_functional_residues.tsv

# --- Step 05: IQ-TREE rooted on BBE (Methods 6)
bash "$S/05_run_phylogeny.sh" --alignment 03_candidates.aligned.fasta --outdir 05_phylogeny \
  --prefix "${H}_synthase" --seed "$SEED" --outgroup "$BBE_ID"
TREE="05_phylogeny/${H}_synthase.contree"
[[ -s "$TREE" ]] || die "tree not produced: $TREE"
log "tree root (first leaf in Newick): $(grep -o '[(,][^():,;]*' "$TREE" | head -1 | tr -d '(,')"

# --- Step 05b: L5 (smallest clade with exactly one reference, UFBoot >= 95)
python "$S/05b_assign_clade.py" --tree "$TREE" --out 05b_l5.tsv \
  --reference-map "$REFMAP" --outgroup BBE --min-support 95.0 > 05b.log
head -1 05b.log

# --- Step 03g: genome-side L2/L6 (decision F: trust the genome only when cleanly full-length)
GENOME=""
for _g in "$ROOT"/00_genomes/*/"$H"/"$H".unmasked.fasta.gz "$ROOT"/00_genomes/*/"$H"/"$H".softmasked.fasta.gz; do
  if [[ -s "$_g" ]]; then GENOME="$_g"; break; fi
done
[[ -n "$GENOME" ]] || die "genome FASTA not found for $H"
python "$S/03g_genomic_orf_check.py" --gene-models "$W/02_gene_models.tsv" --genome "$GENOME" \
  --out 03g_genomic_orf.tsv 2> 03g.log
cat 03g.log

# --- Step 05c: evidence table (join only)
python "$S/05c_build_evidence.py" --orf 03_orf.tsv --residues 04_functional_residues.tsv \
  --l5 05b_l5.tsv --gene-models "$W/02_gene_models.tsv" --blastp "$W/01_search/blastp_hits.tsv" \
  --haplotype "$H" --genomic-orf 03g_genomic_orf.tsv --out 06_evidence.tsv 2> 05c.log
grep -E "unmatched|L5:|mismatch|L2_L6_source" 05c.log

# --- Step 06: classification
python "$S/06_assign_L1_L7.py" --evidence 06_evidence.tsv --out 06_classification.tsv

python3 - "$H" <<'PY'
import csv, sys
from collections import Counter
rows = list(csv.DictReader(open('06_classification.tsv'), delimiter='\t'))
c = Counter(r['Classification'] or r['Review_status'] for r in rows)
print(f'[{sys.argv[1]}] ' + '; '.join(f'{k}: {v}' for k, v in c.most_common()))
PY
log "done -> $OUT"
