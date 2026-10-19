#!/usr/bin/env bash
# Copy the machine-readable benchmark results (no search outputs or indexes)
# from a run_all.sh output directory into bench/results/.
#   bench/scripts/collect_results.sh /path/to/OUT bench/results
set -euo pipefail
SRC=$1; DST=$2
mkdir -p "$DST"
DST=$(cd "$DST" && pwd)   # absolute: the copy below runs inside SRC
cd "$SRC"
# times.prev.jsonl (1.2 suite) holds superseded records of earlier
# invocations, so it is skipped.
find . -type f \( -name '*.jsonl' -o -name '*.metrics.json' -o -name '*.stats.json' -o -name '*.roc.tsv' \
  -o -name 'environment*.txt' -o -name 'final_binary.txt' -o -name 'sizes.txt' -o -name 'dataset.sha256' -o -name 'rerun_check.txt' -o -name 'outputs_check.txt' \
  -o -name 'summary.json' -o -name 'summary.md' -o -name '*.png' -o -name 'run.log' \) \
  ! -name '*.prev.jsonl' ! -name '*.prev_*.jsonl' -print0 | while IFS= read -r -d '' f; do
  mkdir -p "$DST/$(dirname "$f")"
  cp "$f" "$DST/$f"
done
# keep only the test-split ROC curves (the plots use them); drop per-run stderr logs
find "$DST" -name '*.tune.roc.tsv' -delete -o -name '*.all.roc.tsv' -delete
du -sh "$DST"
