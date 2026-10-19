#!/usr/bin/env bash
#
# validate_scaling.sh — close the one gap the development sandbox could not.
#
# The parallel query loop was written and its CORRECTNESS validated (byte-
# identical, run-to-run stable output at 1/2/4/8/16 threads), but its SPEEDUP
# is unmeasured: the sandbox had a single core, where rayon overhead makes the
# parallel path look marginally *worse* than sequential.
#
# Run this on real multi-core hardware to get the real number.
#
# Usage:
#   ./validate_scaling.sh <query.fa> <db.fa> [saber_binary]
#
set -euo pipefail

QUERY="${1:?usage: $0 <query.fa> <db.fa> [saber_binary]}"
DB="${2:?usage: $0 <query.fa> <db.fa> [saber_binary]}"
SABER="${3:-./target/release/saber}"

if [[ ! -x "$SABER" ]]; then
    echo "error: $SABER not found or not executable. Run 'cargo build --release' first." >&2
    exit 1
fi

CORES=$(nproc 2>/dev/null || sysctl -n hw.ncpu 2>/dev/null || echo 1)
echo "host: ${CORES} logical cores"
echo "query: ${QUERY}"
echo "db:    ${DB}"
echo

# Thread counts to test: powers of two up to the core count, plus the core count.
THREADS=()
t=1
while (( t <= CORES )); do THREADS+=("$t"); t=$(( t * 2 )); done
[[ "${THREADS[-1]}" != "$CORES" ]] && THREADS+=("$CORES")

TMPDIR_LOCAL=$(mktemp -d)
trap 'rm -rf "$TMPDIR_LOCAL"' EXIT

printf "%-10s %10s %10s %10s %12s\n" "THREADS" "wall(s)" "speedup" "efficiency" "peakRSS(MB)"

BASE_WALL=""
for t in "${THREADS[@]}"; do
    OUT="$TMPDIR_LOCAL/out_${t}.m8"
    # Best of 3, measuring peak RSS of the child via the shell's time builtin
    # is unreliable, so prefer /usr/bin/time -v when present.
    BEST=""
    RSS="n/a"
    for _ in 1 2 3; do
        START=$(date +%s.%N)
        if command -v /usr/bin/time >/dev/null 2>&1; then
            /usr/bin/time -v "$SABER" -q "$QUERY" -d "$DB" --threads "$t" --quiet \
                -o "$OUT" 2> "$TMPDIR_LOCAL/time_${t}.txt" || true
            R=$(awk '/Maximum resident set size/ {print $NF}' "$TMPDIR_LOCAL/time_${t}.txt" 2>/dev/null || echo "")
            [[ -n "$R" ]] && RSS=$(( R / 1024 ))
        else
            "$SABER" -q "$QUERY" -d "$DB" --threads "$t" --quiet -o "$OUT" 2>/dev/null
        fi
        END=$(date +%s.%N)
        W=$(echo "$END - $START" | bc)
        if [[ -z "$BEST" ]] || (( $(echo "$W < $BEST" | bc -l) )); then BEST="$W"; fi
    done

    [[ -z "$BASE_WALL" ]] && BASE_WALL="$BEST"
    SPEEDUP=$(echo "scale=2; $BASE_WALL / $BEST" | bc)
    EFF=$(echo "scale=1; 100 * $BASE_WALL / $BEST / $t" | bc)
    printf "%-10s %10.2f %9sx %9s%% %12s\n" "$t" "$BEST" "$SPEEDUP" "$EFF" "$RSS"
done

echo
echo "--- correctness: every thread count must produce identical output ---"
REF="$TMPDIR_LOCAL/out_${THREADS[0]}.m8"
FAIL=0
for t in "${THREADS[@]}"; do
    if diff -q "$REF" "$TMPDIR_LOCAL/out_${t}.m8" >/dev/null 2>&1; then
        echo "  threads=$t  IDENTICAL"
    else
        echo "  threads=$t  *** DIFFERS FROM threads=${THREADS[0]} ***"
        FAIL=1
    fi
done

echo
if (( FAIL )); then
    echo "RESULT: FAIL — parallel output diverges. Do not ship."
    exit 1
fi
echo "RESULT: PASS — output is thread-count invariant."
echo
echo "Interpreting the scaling column: the query loop is embarrassingly"
echo "parallel, so efficiency should stay high until memory bandwidth or the"
echo "sequential index build starts to dominate. If efficiency collapses early,"
echo "the index build (which is only partly parallel — the postings fill pass"
echo "is still sequential) is the likely ceiling; profile that next."
