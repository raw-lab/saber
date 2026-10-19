#!/usr/bin/env bash
# Protein-vs-protein search (--mode pp, blastp-like) on the bundled example data.
#
# Run from anywhere after `cargo build --release`.
#   SABER=/path/to/saber   use another binary (default: ./target/release/saber)
#   OUT=/some/dir          keep the outputs there (default: a new temporary directory)
set -euo pipefail
cd "$(dirname "$0")/.."
SABER="${SABER:-./target/release/saber}"
OUT="${OUT:-$(mktemp -d)}"
mkdir -p "$OUT"
Q=examples/protein_query.fa
DB=examples/protein_db.fa
[ -x "$SABER" ] || { echo "saber not found at $SABER (run: cargo build --release)" >&2; exit 1; }
echo "outputs: $OUT"

echo "== 1. default search: BLOSUM62, gaps 11/1, E <= 10, BLAST tabular output (m8)"
"$SABER" -q "$Q" -d "$DB" -o "$OUT/pp.m8"
cat "$OUT/pp.m8"

echo "== 2. sensitive preset, E <= 1e-3, at most 5 subjects per query"
"$SABER" -q "$Q" -d "$DB" --preset sensitive -E 1e-3 -k 5 -o "$OUT/pp_sensitive.m8"
cat "$OUT/pp_sensitive.m8"

echo "== 3. BLOSUM45 with NCBI's default gap costs for it (15/2 in the BLAST convention)"
# Gap costs without published gapped Karlin-Altschul constants are rejected
# unless --allow-uncalibrated-evalues is given (E-values are then marked).
"$SABER" -q "$Q" -d "$DB" -m BLOSUM45 -g 15 -e 2 --gap-convention blast -o "$OUT/pp_blosum45.m8"
cat "$OUT/pp_blosum45.m8"

echo "== 4. pairwise text (m0) and JSON output"
"$SABER" -q "$Q" -d "$DB" -E 1e-5 -k 3 --outfmt m0 -o "$OUT/pp.m0"
"$SABER" -q "$Q" -d "$DB" --outfmt json -o "$OUT/pp.json"
head -n 20 "$OUT/pp.m0"

echo "== 5. persistent index: build once, then search it memory-mapped"
"$SABER" --makedb "$OUT/protein_db.sdx" -d "$DB"
"$SABER" -q "$Q" --db "$OUT/protein_db.sdx" --verify-index --stats "$OUT/run.json" -o "$OUT/pp_db.m8"
cmp "$OUT/pp.m8" "$OUT/pp_db.m8" && echo "indexed search output is identical to the in-memory search"

echo "done: $OUT"
