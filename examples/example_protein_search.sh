#!/bin/bash
# Example: Protein-to-Protein Search (P-P mode)
# This example demonstrates a standard protein sequence search

set -e

# Setup
SABER="./target/release/saber"
QUERIES="examples/data/protein_queries.fasta"
DATABASE="examples/data/nr_small.fasta"
OUTPUT="results_pp.m8"

echo "================================================"
echo "SABER Example: Protein-to-Protein Search"
echo "================================================"

# Standard protein search with default BLOSUM62
echo "Running standard protein search..."
$SABER \
    -i "$QUERIES" \
    -d "$DATABASE" \
    -m BLOSUM_62 \
    -g 10 \
    -e 1 \
    -v 10.0 \
    -a 10 \
    -o "$OUTPUT" \
    --outfmt m8 \
    -t auto \
    --verbose

echo ""
echo "Results saved to: $OUTPUT"
echo "First 5 results:"
head -5 "$OUTPUT"

# Sensitive search with stricter filtering
echo ""
echo "Running sensitive search (E-value < 1e-10)..."
$SABER \
    -i "$QUERIES" \
    -d "$DATABASE" \
    -m BLOSUM_62 \
    -g 11 \
    -e 1 \
    -v 0.000000001 \
    -a 5 \
    -o results_pp_sensitive.m8 \
    --outfmt m8 \
    -t auto

echo "Sensitive results saved to: results_pp_sensitive.m8"

# Fast search with larger candidates
echo ""
echo "Running fast search (more candidates)..."
$SABER \
    -i "$QUERIES" \
    -d "$DATABASE" \
    -m BLOSUM_62 \
    -c 50000 \
    -T 13 \
    -o results_pp_fast.m8 \
    --outfmt m8 \
    -t auto

echo "Fast search results saved to: results_pp_fast.m8"

# Search with different scoring matrix
echo ""
echo "Running with BLOSUM45 (more sensitive)..."
$SABER \
    -i "$QUERIES" \
    -d "$DATABASE" \
    -m BLOSUM_45 \
    -v 1.0 \
    -a 10 \
    -o results_pp_blosum45.m8 \
    --outfmt m8

echo "BLOSUM45 results saved to: results_pp_blosum45.m8"

# Search with formatted output
echo ""
echo "Running with formatted output (M0 format)..."
$SABER \
    -i "$QUERIES" \
    -d "$DATABASE" \
    -m BLOSUM_62 \
    -v 1e-5 \
    -a 3 \
    -o results_pp_formatted.txt \
    --outfmt m0

echo "Formatted results saved to: results_pp_formatted.txt"

# JSON output for programmatic processing
echo ""
echo "Running with JSON output..."
$SABER \
    -i "$QUERIES" \
    -d "$DATABASE" \
    -m BLOSUM_62 \
    -o results_pp.json \
    --outfmt json \
    -a 10

echo "JSON results saved to: results_pp.json"

echo ""
echo "================================================"
echo "Examples completed successfully!"
echo "================================================"
