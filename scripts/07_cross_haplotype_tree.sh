#!/usr/bin/env bash
# 07_cross_haplotype_tree.sh - build the cross-haplotype tree shown in Fig. 3.
#
# Fig. 3 was previously built by a separate, undocumented run (2026-10-08) that
# used "-m TEST" and an unfixed seed, so it could not be reproduced. This script
# rebuilds it with the same settings as the per-haplotype trees produced by
# 05_run_phylogeny.sh, so that every tree in the paper comes from one procedure.
#
# Run from the project directory ($ROOT) with the plant_genomics environment active:
#     bash cannabinoid-synthase-workflow/scripts/07_cross_haplotype_tree.sh
set -euo pipefail

ROOT="${ROOT:-$PWD}"
S="$ROOT/cannabinoid-synthase-workflow/scripts"
OUT="${OUT:-$ROOT/08_results/figures}"
PREFIX="Fig3_cross_haplotype"
BBE_ID='sp|P30986|BBE_ESCCA'
SEED="${SEED:-12345}"
THREADS="${THREADS:-AUTO}"

die() { echo "ERROR: $*" >&2; exit 1; }

command -v mafft   >/dev/null || die "mafft not found (conda activate plant_genomics)"
command -v iqtree2 >/dev/null || command -v iqtree >/dev/null \
  || die "iqtree not found (conda activate plant_genomics)"
IQTREE=$(command -v iqtree2 || command -v iqtree)

mkdir -p "$OUT"
cd "$ROOT"

# --- 1. collect the 13 sequences --------------------------------------------
python "$S/07a_collect_cross_haplotype.py" --out "$OUT/$PREFIX.fasta"

# --- 2. align (same as the per-haplotype step: MAFFT L-INS-i) ---------------
mafft --localpair --maxiterate 1000 --quiet \
  "$OUT/$PREFIX.fasta" > "$OUT/$PREFIX.aligned.fasta"
grep -qF ">$BBE_ID" "$OUT/$PREFIX.aligned.fasta" || die "BBE not in alignment"
echo "[07] aligned: $(grep -c '^>' "$OUT/$PREFIX.aligned.fasta") sequences"

# --- 3. tree (same settings as 05_run_phylogeny.sh) -------------------------
"$IQTREE" -s "$OUT/$PREFIX.aligned.fasta" \
  -m MFP -B 1000 -T "$THREADS" -seed "$SEED" \
  -o "$BBE_ID" --prefix "$OUT/$PREFIX" -redo

echo
echo "[07] model:  $(grep -m1 'Model of substitution' "$OUT/$PREFIX.iqtree" || true)"
echo "[07] tree:   $OUT/$PREFIX.contree"
echo "[07] settings: -m MFP -B 1000 -T $THREADS -seed $SEED -o $BBE_ID"
