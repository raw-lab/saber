#!/usr/bin/env bash
# GPU scoring of candidate windows (experimental).
#
# The GPU backend is compiled only with `cargo build --release --features gpu-cuda`
# and needs a CUDA device at run time. Without either, `--gpu` prints a notice and
# uses the CPU kernels, so this script runs on any machine. The CUDA kernel was
# checked in emulation only; it has not been run on GPU hardware for this release
# (docs/LIMITATIONS.md), so compare its output with the CPU run (step 3) and run
# `saber --doctor --gpu` on a GPU machine.
#
# Run from anywhere after building.
#   SABER=/path/to/saber   use another binary (default: ./target/release/saber)
#   OUT=/some/dir          keep the outputs there (default: a new temporary directory)
set -euo pipefail
cd "$(dirname "$0")/.."
SABER="${SABER:-./target/release/saber}"
OUT="${OUT:-$(mktemp -d)}"
mkdir -p "$OUT"
Q=examples/protein_query.fa
DB=examples/protein_db.fa
[ -x "$SABER" ] || { echo "saber not found at $SABER (run: cargo build --release [--features gpu-cuda])" >&2; exit 1; }
echo "outputs: $OUT"

echo "== 1. diagnostics: SIMD support, kernel self-tests and GPU backend status"
"$SABER" --doctor | grep -E 'AVX2|GPU|Summary|Result' || true

echo "== 2. search with --gpu (falls back to the CPU kernels when no backend is available)"
"$SABER" -q "$Q" -d "$DB" --gpu -o "$OUT/pp_gpu.m8"

echo "== 3. the same search on the CPU; outputs must match"
"$SABER" -q "$Q" -d "$DB" -o "$OUT/pp_cpu.m8"
cmp "$OUT/pp_gpu.m8" "$OUT/pp_cpu.m8" && echo "GPU-path and CPU outputs are identical"
cat "$OUT/pp_cpu.m8"

echo "done: $OUT"
