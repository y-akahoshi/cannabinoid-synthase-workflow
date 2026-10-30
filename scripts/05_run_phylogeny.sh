#!/usr/bin/env bash
set -euo pipefail

usage() {
  cat <<'EOF'
Step 5: infer an IQ-TREE maximum-likelihood phylogeny from the shared
MAFFT alignment produced in Step 3.

Usage:
  05_run_phylogeny.sh \
    --alignment FILE \
    --outdir DIR \
    [--prefix NAME] \
    [--include-ids FILE] \
    [--model MFP] \
    [--ufboot 1000] \
    [--threads AUTO] \
    [--seed INT] \
    [--outgroup ID]

Notes:
- --alignment is already aligned; this script does NOT run MAFFT again.
- --include-ids optionally restricts the tree to a manuscript/test subset.
  The file may contain one sequence ID per line or a TSV whose first column
  contains IDs. Blank lines and lines beginning with # are ignored.
- If --include-ids is omitted, all sequences in the Step 3 alignment are used.
- The default model search is ModelFinder Plus (MFP), with 1000 UFBoot2
  replicates. Verify these settings against the manuscript/previous IQ-TREE
  log before final submission.
EOF
}

ALIGNMENT=""
OUTDIR=""
PREFIX="synthase_tree"
INCLUDE_IDS=""
MODEL="MFP"
UFBOOT="1000"
THREADS="AUTO"
SEED=""
OUTGROUP=""

while [[ $# -gt 0 ]]; do
  case "$1" in
    --alignment) ALIGNMENT="$2"; shift 2 ;;
    --outdir) OUTDIR="$2"; shift 2 ;;
    --prefix) PREFIX="$2"; shift 2 ;;
    --include-ids) INCLUDE_IDS="$2"; shift 2 ;;
    --model) MODEL="$2"; shift 2 ;;
    --ufboot) UFBOOT="$2"; shift 2 ;;
    --threads) THREADS="$2"; shift 2 ;;
    --seed) SEED="$2"; shift 2 ;;
    --outgroup) OUTGROUP="$2"; shift 2 ;;
    -h|--help) usage; exit 0 ;;
    *) echo "ERROR: unknown argument: $1" >&2; usage >&2; exit 2 ;;
  esac
done

[[ -n "$ALIGNMENT" ]] || { echo "ERROR: --alignment is required" >&2; exit 2; }
[[ -n "$OUTDIR" ]] || { echo "ERROR: --outdir is required" >&2; exit 2; }
[[ -s "$ALIGNMENT" ]] || { echo "ERROR: alignment not found or empty: $ALIGNMENT" >&2; exit 2; }

if command -v iqtree2 >/dev/null 2>&1; then
  IQTREE="iqtree2"
elif command -v iqtree >/dev/null 2>&1; then
  IQTREE="iqtree"
else
  echo "ERROR: iqtree2/iqtree not found in PATH" >&2
  exit 127
fi

mkdir -p "$OUTDIR"
TREE_ALIGNMENT="$ALIGNMENT"

if [[ -n "$INCLUDE_IDS" ]]; then
  [[ -s "$INCLUDE_IDS" ]] || { echo "ERROR: --include-ids file not found/empty: $INCLUDE_IDS" >&2; exit 2; }
  TREE_ALIGNMENT="$OUTDIR/${PREFIX}.subset_alignment.fasta"

  python - "$ALIGNMENT" "$INCLUDE_IDS" "$TREE_ALIGNMENT" <<'PY'
import sys
from Bio import AlignIO
from Bio.Align import MultipleSeqAlignment

aln_path, ids_path, out_path = sys.argv[1:4]
aln = AlignIO.read(aln_path, "fasta")
by_id = {r.id: r for r in aln}

wanted = []
with open(ids_path) as fh:
    for line in fh:
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        token = line.split("\t", 1)[0].strip()
        # Allow a conventional header.
        if token.lower() in {"sequence", "sequence_id", "id"}:
            continue
        wanted.append(token)

missing = [x for x in wanted if x not in by_id]
if missing:
    raise SystemExit("IDs not found in alignment: " + ", ".join(missing))

if len(wanted) < 3:
    raise SystemExit("Need at least three sequences for a phylogenetic tree.")

subset = MultipleSeqAlignment([by_id[x] for x in wanted])
AlignIO.write(subset, out_path, "fasta")
print(f"Wrote subset alignment with {len(subset)} sequences: {out_path}")
PY
fi

PRE="$OUTDIR/$PREFIX"
CMD=("$IQTREE" -s "$TREE_ALIGNMENT" -pre "$PRE" -m "$MODEL" -B "$UFBOOT" -T "$THREADS")

if [[ -n "$SEED" ]]; then
  CMD+=(-seed "$SEED")
fi
if [[ -n "$OUTGROUP" ]]; then
  CMD+=(-o "$OUTGROUP")
fi

printf 'Command:' > "$OUTDIR/${PREFIX}.command.txt"
printf ' %q' "${CMD[@]}" >> "$OUTDIR/${PREFIX}.command.txt"
printf '\n' >> "$OUTDIR/${PREFIX}.command.txt"

echo "Running:"
cat "$OUTDIR/${PREFIX}.command.txt"
"${CMD[@]}"

echo
echo "IQ-TREE completed."
echo "Key outputs:"
for f in \
  "$PRE.treefile" \
  "$PRE.iqtree" \
  "$PRE.log" \
  "$PRE.contree"; do
  [[ -e "$f" ]] && echo "  $f"
done
