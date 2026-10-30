#!/usr/bin/env bash
# Runs the example dataset and checks that every copy receives the expected class.
set -euo pipefail
REPO="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
bash "$REPO/scripts/run_example.sh" > "$REPO/example_out.log" 2>&1 || { echo "FAIL: run_example.sh stopped (see example_out.log)"; tail -20 "$REPO/example_out.log"; exit 1; }
awk -F'\t' 'NR==1{for(i=1;i<=NF;i++){if($i=="Sequence")s=i;if($i=="Classification")c=i}} {print $s"\t"$c}' \
  "$REPO/example_out/SAN2b/06_classification.tsv" > "$REPO/example_out/observed.tsv"
if diff -u "$REPO/tests/expected_SAN2b_classification.tsv" "$REPO/example_out/observed.tsv"; then
  echo "PASS: all $(($(wc -l < "$REPO/tests/expected_SAN2b_classification.tsv") - 1)) copies classified as expected"
else
  echo "FAIL: classification differs from expected"; exit 1
fi
