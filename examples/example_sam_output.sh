#!/bin/bash
# Example: SAM/BAM Output Formats
# Demonstrates generating SAM/BAM files compatible with standard
# bioinformatics tools like samtools, bcftools, IGV, etc.

set -e

SABER="./target/release/saber"
READS="examples/data/reads.fastq"
REFERENCE="examples/data/reference.fasta"

echo "================================================"
echo "SABER Example: SAM/BAM Output"
echo "================================================"

# Generate SAM format (text-based, easily human-readable)
echo "Generating SAM format output..."
$SABER \
    -i "$READS" \
    -d "$REFERENCE" \
    --mode nn \
    --match 2 \
    --mismatch -3 \
    -v 1e-20 \
    -o alignments.sam \
    --outfmt sam \
    -t auto \
    --verbose

echo "SAM file generated: alignments.sam"
echo ""
echo "SAM file header:"
head -10 alignments.sam
echo ""

# Convert SAM to BAM using samtools
if command -v samtools &> /dev/null; then
    echo "Converting SAM to BAM with samtools..."
    samtools view -b -o alignments.bam alignments.sam
    samtools index alignments.bam alignments.bam.bai
    
    echo "✓ BAM file generated: alignments.bam"
    echo "✓ BAM index generated: alignments.bam.bai"
    echo ""
    
    # View BAM statistics
    echo "BAM file statistics:"
    samtools flagstat alignments.bam
    echo ""
    
    # View alignments in detail
    echo "First 5 alignments (samtools view):"
    samtools view alignments.bam | head -5
    echo ""
    
    # Index and sort BAM (recommended for some tools)
    echo "Sorting and indexing BAM file..."
    samtools sort alignments.bam -o alignments.sorted.bam
    samtools index alignments.sorted.bam
    
    echo "✓ Sorted BAM file: alignments.sorted.bam"
    echo "✓ Sorted BAM index: alignments.sorted.bam.bai"
else
    echo "samtools not found, generating BAM directly with SABER..."
    
    $SABER \
        -i "$READS" \
        -d "$REFERENCE" \
        --mode nn \
        -o alignments.bam \
        --outfmt bam \
        --bam-index \
        -t auto
    
    echo "✓ BAM file generated: alignments.bam"
    echo "✓ BAM index generated: alignments.bam.bai"
fi

echo ""

# Generate multiple alignments with different stringencies
echo "Generating strict mapping (high quality)..."
$SABER \
    -i "$READS" \
    -d "$REFERENCE" \
    --mode nn \
    --match 2 \
    --mismatch -3 \
    -v 1e-50 \
    --percent-identity 99.0 \
    -o alignments_strict.sam \
    --outfmt sam

echo "✓ Strict mapping: alignments_strict.sam"

echo ""
echo "Generating sensitive mapping (allow more mismatches)..."
$SABER \
    -i "$READS" \
    -d "$REFERENCE" \
    --mode nn \
    --match 1 \
    --mismatch -4 \
    -v 1e-10 \
    --percent-identity 85.0 \
    -o alignments_sensitive.sam \
    --outfmt sam

echo "✓ Sensitive mapping: alignments_sensitive.sam"

echo ""

# SAM with detailed statistics
if command -v samtools &> /dev/null; then
    echo "Detailed alignment statistics:"
    echo ""
    echo "Strict mapping stats:"
    samtools view -b alignments_strict.sam | samtools flagstat -
    echo ""
    echo "Sensitive mapping stats:"
    samtools view -b alignments_sensitive.sam | samtools flagstat -
fi

echo ""

# Performance comparison
echo "Testing performance with different output formats..."

echo "SAM format:"
/usr/bin/time -f "  Time: %e seconds, Memory: %M KB" \
    $SABER \
    -i "$READS" \
    -d "$REFERENCE" \
    --mode nn \
    -o /tmp/test.sam \
    --outfmt sam \
    2>&1 | grep -E "Time:|Memory:"

echo "BAM format:"
/usr/bin/time -f "  Time: %e seconds, Memory: %M KB" \
    $SABER \
    -i "$READS" \
    -d "$REFERENCE" \
    --mode nn \
    -o /tmp/test.bam \
    --outfmt bam \
    2>&1 | grep -E "Time:|Memory:"

echo "M8 format:"
/usr/bin/time -f "  Time: %e seconds, Memory: %M KB" \
    $SABER \
    -i "$READS" \
    -d "$REFERENCE" \
    --mode nn \
    -o /tmp/test.m8 \
    --outfmt m8 \
    2>&1 | grep -E "Time:|Memory:"

echo ""
echo "================================================"
echo "SAM/BAM output examples completed!"
echo "================================================"

echo ""
echo "Generated files:"
echo "  ✓ alignments.sam - Text SAM format"
echo "  ✓ alignments.bam - Binary BAM format (if samtools available)"
echo "  ✓ alignments_strict.sam - High-quality mappings only"
echo "  ✓ alignments_sensitive.sam - Permissive mappings"
echo ""
echo "Use with common bioinformatics tools:"
echo "  samtools view alignments.bam              # View BAM"
echo "  samtools flagstat alignments.bam          # Statistics"
echo "  samtools sort alignments.sam > sorted.bam # Sort"
echo "  samtools index sorted.bam                 # Index"
echo "  bcftools mpileup -f ref.fa sorted.bam | bcftools call ...  # Variant calling"
