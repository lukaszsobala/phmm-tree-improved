#!/bin/bash

# Performance comparison script for PHMM-Tree parallelization
# Tests different thread counts and measures performance

echo "PHMM-Tree Parallelization Performance Test"
echo "=========================================="

# Check if test data exists
if [ ! -f "test_data.txt" ]; then
    echo "Creating sample test data..."
    # Create a simple distance matrix for testing
    cat > test_data.txt << EOF
5
species1  0.000 0.123 0.234 0.345 0.456
species2  0.123 0.000 0.345 0.456 0.567
species3  0.234 0.345 0.000 0.567 0.678
species4  0.345 0.456 0.567 0.000 0.789
species5  0.456 0.567 0.678 0.789 0.000
EOF
fi

# Function to run performance test
run_test() {
    local threads=$1
    echo "Testing with $threads thread(s)..."
    
    export OMP_NUM_THREADS=$threads
    
    # Time the execution (we'll simulate since we don't have full program)
    start_time=$(date +%s.%N)
    
    # Simulate processing (in real scenario, this would run the actual algorithm)
    echo "  Thread count: $threads"
    echo "  Simulated processing time varies with thread count"
    
    end_time=$(date +%s.%N)
    execution_time=$(echo "$end_time - $start_time" | bc -l)
    
    echo "  Execution time: ${execution_time} seconds"
    echo ""
}

# Test different thread configurations
echo "System information:"
echo "  CPU cores: $(nproc)"
echo "  OpenMP max threads: $(./test_parallel 2>&1 | grep "Maximum threads" | cut -d: -f2 | xargs)"
echo ""

# Run tests with different thread counts
for threads in 1 2 4 8 12; do
    if [ $threads -le $(nproc) ]; then
        run_test $threads
    fi
done

echo "Performance Notes:"
echo "- Speedup should be approximately linear with thread count for large datasets"
echo "- Small datasets may show minimal improvement due to overhead"
echo "- Optimal thread count is typically equal to the number of CPU cores"
echo "- Memory usage increases with thread count"

echo ""
echo "To run actual performance tests with real data:"
echo "1. Compile the full program: make"
echo "2. Use environment variable: OMP_NUM_THREADS=N ./phmm-tree"
echo "3. Compare execution times with different N values"
