#!/usr/bin/env bash
# Translated searches on the bundled example data:
#   --mode nx  nucleotide queries vs a protein database (blastx-like, six frames)
#   --mode pn  protein queries vs a nucleotide database (tblastn-like, six frames)
#
# Run from anywhere after `cargo build --release`.
#   SABER=/path/to/saber   use another binary (default: ./target/release/saber)
#   OUT=/some/dir          keep the outputs there (default: a new temporary directory)
set -euo pipefail
cd "$(dirname "$0")/.."
SABER="${SABER:-./target/release/saber}"
OUT="${OUT:-$(mktemp -d)}"
mkdir -p "$OUT"
[ -x "$SABER" ] || { echo "saber not found at $SABER (run: cargo build --release)" >&2; exit 1; }
echo "outputs: $OUT"

echo "== 1. blastx-like: DNA queries vs proteins (standard genetic code)"
"$SABER" --mode nx -q examples/nucleotide_query.fa -d examples/protein_db.fa -o "$OUT/nx.m8"
cat "$OUT/nx.m8"
# The query name carries the frame (|frame+1 ... |frame-3); query coordinates
# are in nucleotides, subject coordinates in residues.

echo "== 2. the same with the bacterial/plastid code (NCBI table 11); names or numbers work"
"$SABER" --mode nx --genetic-code 11 -q examples/nucleotide_query.fa -d examples/protein_db.fa -o "$OUT/nx_table11.m8"
cat "$OUT/nx_table11.m8"

echo "== 3. vertebrate mitochondrial code (table 2): TGA = W, AGA/AGG = stop, ATA = M"
"$SABER" --mode nx --genetic-code vertebrate -q examples/nucleotide_query.fa -d examples/protein_db.fa -o "$OUT/nx_table2.m8"
cat "$OUT/nx_table2.m8"

echo "== 4. tblastn-like: proteins vs a DNA database"
"$SABER" --mode pn -q examples/protein_query.fa -d examples/nucleotide_db.fa -o "$OUT/pn.m8"
cat "$OUT/pn.m8"

echo "== 5. translated persistent index (built with the same --mode)"
"$SABER" --mode pn --makedb "$OUT/nucleotide_db.pn.sdx" -d examples/nucleotide_db.fa
"$SABER" --mode pn -q examples/protein_query.fa --db "$OUT/nucleotide_db.pn.sdx" -o "$OUT/pn_db.m8"
cmp "$OUT/pn.m8" "$OUT/pn_db.m8" && echo "indexed search output is identical to the in-memory search"

echo "done: $OUT"
