#!/usr/bin/env bash
# SAM output on the bundled example data, with optional samtools post-processing.
# SABER writes SAM only; BAM conversion, sorting and indexing are left to samtools.
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

echo "== 1. nucleotide alignments as SAM (AS = raw score, NM = edit distance, EV = E-value)"
"$SABER" --mode nn -q examples/nucleotide_query.fa -d examples/nucleotide_db.fa --outfmt sam -o "$OUT/nn.sam"
cut -f1-9,12- "$OUT/nn.sam"

echo "== 2. protein alignments as SAM (sequence field holds the protein query)"
"$SABER" -q examples/protein_query.fa -d examples/protein_db.fa -E 1e-5 --outfmt sam -o "$OUT/pp.sam"
cut -f1-9,12- "$OUT/pp.sam"

echo "== 3. high-identity alignments only"
"$SABER" --mode nn -q examples/nucleotide_query.fa -d examples/nucleotide_db.fa --percent-identity 99 --outfmt sam -o "$OUT/nn_strict.sam"
grep -vc '^@' "$OUT/nn_strict.sam" | sed 's/$/ alignment(s) at >= 99 % identity/'

if command -v samtools >/dev/null 2>&1; then
    echo "== 4. samtools: SAM -> sorted, indexed BAM"
    # The SAM header has no @SQ lines (reference lengths are not written), so
    # give samtools the reference index explicitly.
    samtools faidx examples/nucleotide_db.fa --fai-idx "$OUT/nucleotide_db.fa.fai"
    samtools view -b -t "$OUT/nucleotide_db.fa.fai" "$OUT/nn.sam" | samtools sort -o "$OUT/nn.sorted.bam"
    samtools index "$OUT/nn.sorted.bam"
    samtools flagstat "$OUT/nn.sorted.bam"
else
    echo "== 4. samtools not installed; skipping BAM conversion"
fi

echo "done: $OUT"
