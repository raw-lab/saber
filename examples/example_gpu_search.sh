#!/bin/bash
# Example: GPU-Accelerated Search
# Demonstrates SABER with GPU acceleration for maximum performance
# on large-scale searches

set -e

SABER="./target/release/saber"
QUERIES="examples/data/large_query_set.fasta"
DATABASE="examples/data/large_database.fasta"

echo "================================================"
echo "SABER Example: GPU-Accelerated Search"
echo "================================================"

# Check GPU availability
echo "Checking GPU availability..."
if $SABER --gpu-info 2>/dev/null; then
    echo "✓ GPU detected and ready"
else
    echo "⚠ No GPU detected, falling back to CPU"
    CPU_ONLY=true
fi

echo ""

# GPU-accelerated search with auto thread optimization
if [ -z "$CPU_ONLY" ]; then
    echo "Running GPU-accelerated search..."
    $SABER \
        -i "$QUERIES" \
        -d "$DATABASE" \
        -m BLOSUM_62 \
        -o results_gpu.m8 \
        --outfmt m8 \
        --gpu \
        --gpu-device 0 \
        -t auto \
        --verbose

    echo "Results saved to: results_gpu.m8"
    echo ""

    # Multi-GPU search (if available)
    echo "Running multi-GPU search (if dual GPU system)..."
    $SABER \
        -i "$QUERIES" \
        -d "$DATABASE" \
        -m BLOSUM_62 \
        -o results_multi_gpu.m8 \
        --outfmt m8 \
        --gpu \
        --gpu-device 0,1 \
        -t auto 2>/dev/null || echo "Single GPU available, skipping multi-GPU test"

    echo ""

    # Comparison: GPU vs CPU
    echo "Running CPU-only search for comparison..."
    /usr/bin/time -v $SABER \
        -i "$QUERIES" \
        -d "$DATABASE" \
        -m BLOSUM_62 \
        -o results_cpu.m8 \
        --outfmt m8 \
        -t auto

    echo ""
    echo "GPU vs CPU performance comparison:"
    echo "Check timing information above for memory and execution time"

else
    echo "Running CPU-only search..."
    $SABER \
        -i "$QUERIES" \
        -d "$DATABASE" \
        -m BLOSUM_62 \
        -o results_cpu.m8 \
        --outfmt m8 \
        -t auto

    echo "Results saved to: results_cpu.m8"
fi

echo ""

# Large-scale batch processing with progress tracking
echo "Running large-scale batch search with progress..."
$SABER \
    -i "$QUERIES" \
    -d "$DATABASE" \
    -m BLOSUM_62 \
    -c 100000 \
    -o results_batch.m8 \
    --outfmt m8 \
    --gpu \
    -t auto \
    --verbose

echo ""
echo "Results saved to: results_batch.m8"

# Memory-efficient search with streaming
echo ""
echo "Running memory-efficient search (streaming mode)..."
$SABER \
    -i "$QUERIES" \
    -d "$DATABASE" \
    -m BLOSUM_62 \
    --streaming \
    -o results_streaming.m8 \
    --outfmt m8 \
    --gpu

echo "Results saved to: results_streaming.m8"

# Profile performance
echo ""
echo "Running with performance profiling..."
if command -v perf &> /dev/null; then
    perf stat -e cycles,instructions,cache-references,cache-misses,branches,branch-misses \
        $SABER \
        -i "$QUERIES" \
        -d "$DATABASE" \
        -m BLOSUM_62 \
        --gpu \
        -o results_profiled.m8 \
        --outfmt m8
else
    echo "perf not available, running without profiling"
    $SABER \
        -i "$QUERIES" \
        -d "$DATABASE" \
        -m BLOSUM_62 \
        --gpu \
        -o results_profiled.m8 \
        --outfmt m8
fi

echo ""
echo "================================================"
echo "GPU acceleration examples completed!"
echo "================================================"

# Summary
echo ""
echo "Summary of GPU features demonstrated:"
echo "  ✓ GPU detection and initialization"
echo "  ✓ Single GPU acceleration"
echo "  ✓ Multi-GPU support"
echo "  ✓ GPU vs CPU performance comparison"
echo "  ✓ Batch processing optimization"
echo "  ✓ Memory-efficient streaming mode"
echo "  ✓ Performance profiling"
