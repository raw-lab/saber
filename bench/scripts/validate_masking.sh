#!/usr/bin/env bash
# Compare SABER's tantan implementation (src/masking.rs) with the reference
# tantan binary (Frith 2011; tantan 49 was used for bench/results).
#
#   bench/scripts/validate_masking.sh OUT.txt PROTEINS.fa [MORE.fa ...]
#
# The published comparison (bench/results/masking_vs_tantan.txt) used the
# first 3,000 domains of scop40.fa plus the 500 MMseqs2 example queries:
#   awk '/^>/{n++} n<=3000' DATA/scop40.fa > maskset.fa
#   cat DATA/mmseqs_QUERY.fasta >> maskset.fa
#
# Needs `tantan` on PATH and this SABER source tree. Both programs print one
# repeat probability per residue (tantan -p -f 1 layout); compare_masking.py
# reports agreement of the >= 0.5 mask calls (the masking threshold used for
# seeding), Pearson r and the mean absolute difference.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "$HERE/../.." && pwd)"
OUTF=$1
shift
TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT
cat "$@" > "$TMP/in.fa"
(cd "$ROOT" && cargo build --release --quiet --example mask_probs)
BIN="${CARGO_TARGET_DIR:-$ROOT/target}/release/examples/mask_probs"
tantan -p -f 1 "$TMP/in.fa" > "$TMP/tantan.probs"
"$BIN" "$TMP/in.fa" > "$TMP/saber.probs" 2> "$TMP/saber.log"
{
  echo "input: $* ($(grep -c '>' "$TMP/in.fa") sequences)"
  echo "reference: $(tantan --version 2>&1 | head -1) -p -f 1"
  echo "saber: examples/mask_probs ($(cat "$TMP/saber.log"))"
  python3 -I "$HERE/compare_masking.py" "$TMP/tantan.probs" "$TMP/saber.probs"
} > "$OUTF"
cat "$OUTF"
