#!/bin/bash
# Example: Nucleotide-to-Nucleotide Search (N-N mode)
# DNA/RNA sequence similarity searching

set -e

SABER="./target/release/saber"
QUERIES="examples/data/dna_queries.fasta"
DATABASE="examples/data/genome_ref.fasta"

echo "================================================"
echo "SABER Example: Nucleotide-to-Nucleotide Search"
echo "================================================"

# Standard DNA search
echo "Running DNA sequence search..."
$SABER \
    -i "$QUERIES" \
    -d "$DATABASE" \
    --mode nn \
    --match 2 \
    --mismatch -3 \
    -g 5 \
    -e 2 \
    -v 1e-10 \
    -o results_nn.m8 \
    --outfmt m8 \
    -t auto \
    --verbose

echo "Results saved to: results_nn.m8"
head -5 results_nn.m8

# Relaxed search for more distant matches
echo ""
echo "Running relaxed DNA search (longer gaps allowed)..."
$SABER \
    -i "$QUERIES" \
    -d "$DATABASE" \
    --mode nn \
    --match 1 \
    --mismatch -4 \
    -g 8 \
    -e 1 \
    -v 1e-5 \
    -o results_nn_relaxed.m8 \
    --outfmt m8

echo "Relaxed search results saved to: results_nn_relaxed.m8"

# Search for very high identity matches only
echo ""
echo "Running strict high-identity search..."
$SABER \
    -i "$QUERIES" \
    -d "$DATABASE" \
    --mode nn \
    --match 5 \
    --mismatch -10 \
    -g 15 \
    -e 2 \
    -v 1e-50 \
    --percent-identity 98.0 \
    -o results_nn_strict.m8 \
    --outfmt m8

echo "Strict search results saved to: results_nn_strict.m8"

# SAM output for read mapping workflows
echo ""
echo "Generating SAM format output..."
$SABER \
    -i "$QUERIES" \
    -d "$DATABASE" \
    --mode nn \
    --match 2 \
    --mismatch -3 \
    -o alignments.sam \
    --outfmt sam \
    -t auto

echo "SAM output saved to: alignments.sam"

# BAM output (indexed binary format, requires htslib)
echo ""
echo "Generating BAM format output..."
$SABER \
    -i "$QUERIES" \
    -d "$DATABASE" \
    --mode nn \
    -o alignments.bam \
    --outfmt bam \
    --bam-index

echo "BAM output saved to: alignments.bam"
echo "BAM index saved to: alignments.bam.bai"

echo ""
echo "================================================"
echo "Examples completed successfully!"
echo "================================================"
