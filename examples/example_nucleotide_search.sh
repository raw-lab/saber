#!/usr/bin/env bash
# Nucleotide-vs-nucleotide search (--mode nn, blastn-like) on the bundled example data.
#
# Run from anywhere after `cargo build --release`.
#   SABER=/path/to/saber   use another binary (default: ./target/release/saber)
#   OUT=/some/dir          keep the outputs there (default: a new temporary directory)
set -euo pipefail
cd "$(dirname "$0")/.."
SABER="${SABER:-./target/release/saber}"
OUT="${OUT:-$(mktemp -d)}"
mkdir -p "$OUT"
Q=examples/nucleotide_query.fa
DB=examples/nucleotide_db.fa
[ -x "$SABER" ] || { echo "saber not found at $SABER (run: cargo build --release)" >&2; exit 1; }
echo "outputs: $OUT"

echo "== 1. default nn search: match 2, mismatch -3, gaps 7/2, contiguous 11-mer seeds"
"$SABER" --mode nn -q "$Q" -d "$DB" -o "$OUT/nn.m8"
cat "$OUT/nn.m8"

echo "== 2. high-identity hits only (>= 98 % identity, E <= 1e-50)"
"$SABER" --mode nn -q "$Q" -d "$DB" --percent-identity 98 -E 1e-50 -o "$OUT/nn_strict.m8"
cat "$OUT/nn_strict.m8"

echo "== 3. scoring for diverged sequences: match 1, mismatch -2, gaps 2/2 (BLAST convention)"
# Combinations without published gapped constants are rejected unless
# --allow-uncalibrated-evalues is given (E-values are then marked UNCALIBRATED).
"$SABER" --mode nn -q "$Q" -d "$DB" --match 1 --mismatch -2 -g 2 -e 2 --gap-convention blast -o "$OUT/nn_diverged.m8"
cat "$OUT/nn_diverged.m8"

echo "== 4. SAM output (see example_sam_output.sh for samtools post-processing)"
"$SABER" --mode nn -q "$Q" -d "$DB" --outfmt sam -o "$OUT/nn.sam"
head -n 4 "$OUT/nn.sam" | cut -c1-120

echo "done: $OUT"
