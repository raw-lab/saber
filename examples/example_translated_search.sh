#!/bin/bash
# Example: Translated Nucleotide-to-Protein Search (N-P mode)
# Search DNA/RNA sequences against protein databases
# with automatic 6-frame translation

set -e

SABER="./target/release/saber"
QUERIES="examples/data/cds.fasta"  # Coding sequences
DATABASE="examples/data/proteins.fasta"  # Protein database

echo "================================================"
echo "SABER Example: Translated Nucleotide-to-Protein"
echo "================================================"

# Standard translated search with standard genetic code
echo "Running translated search (standard genetic code)..."
$SABER \
    -i "$QUERIES" \
    -d "$DATABASE" \
    --mode nx \
    --genetic-code standard \
    -m BLOSUM_62 \
    -g 10 \
    -e 1 \
    -v 1e-10 \
    -a 10 \
    -o results_nx_standard.m8 \
    --outfmt m8 \
    -t auto \
    --verbose

echo "Results saved to: results_nx_standard.m8"
head -5 results_nx_standard.m8

# Search with bacterial genetic code
echo ""
echo "Running with bacterial genetic code..."
$SABER \
    -i "$QUERIES" \
    -d "$DATABASE" \
    --mode nx \
    --genetic-code bacterial \
    -m BLOSUM_62 \
    -v 1e-10 \
    -o results_nx_bacterial.m8 \
    --outfmt m8

echo "Results saved to: results_nx_bacterial.m8"

# Mitochondrial genetic code (vertebrate)
echo ""
echo "Running with vertebrate mitochondrial genetic code..."
$SABER \
    -i "$QUERIES" \
    -d "$DATABASE" \
    --mode nx \
    --genetic-code vertebrate \
    -m BLOSUM_62 \
    -v 1e-8 \
    -a 5 \
    -o results_nx_vertebrate.m8 \
    --outfmt m8

echo "Results saved to: results_nx_vertebrate.m8"

# Ciliate genetic code (different stop codon usage)
echo ""
echo "Running with ciliate nuclear genetic code..."
$SABER \
    -i "$QUERIES" \
    -d "$DATABASE" \
    --mode nx \
    --genetic-code ciliate \
    -m BLOSUM_62 \
    -o results_nx_ciliate.m8 \
    --outfmt m8

echo "Results saved to: results_nx_ciliate.m8"

# Sensitive search with multiple genetic codes comparison
echo ""
echo "Running sensitive comparison search..."
$SABER \
    -i "$QUERIES" \
    -d "$DATABASE" \
    --mode nx \
    --genetic-code archaeal \
    -m BLOSUM_45 \
    -v 1e-15 \
    -a 3 \
    -o results_nx_archaeal.m8 \
    --outfmt m8

echo "Results saved to: results_nx_archaeal.m8"

# JSON output for programmatic frame-wise analysis
echo ""
echo "Running with JSON output for detailed analysis..."
$SABER \
    -i "$QUERIES" \
    -d "$DATABASE" \
    --mode nx \
    --genetic-code standard \
    -o results_nx.json \
    --outfmt json \
    -a 10

echo "JSON results saved to: results_nx.json"

# Protein-to-Translated Nucleotide (reverse search)
echo ""
echo "Running reverse search (protein to translated DNA)..."
$SABER \
    -i "$DATABASE" \
    -d "$QUERIES" \
    --mode pn \
    --genetic-code standard \
    -m BLOSUM_62 \
    -v 1e-10 \
    -o results_pn.m8 \
    --outfmt m8

echo "Reverse search results saved to: results_pn.m8"

echo ""
echo "================================================"
echo "Examples completed successfully!"
echo "================================================"
echo ""
echo "Note: The 6-frame translation is performed automatically"
echo "for all sequences in the query file. Results will show"
echo "the highest-scoring frame for each query-subject pair."
