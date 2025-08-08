#!/bin/bash

# PRC Parallelization Performance Test Script
# Tests the parallelized PRC algorithm performance with different thread counts

echo "=== pHMM-Tree PRC Parallelization Performance Test ==="
echo

# Check if phmm-tree executable exists
if [ ! -f "./phmm-tree" ]; then
    echo "Error: phmm-tree executable not found. Please run 'make' first."
    exit 1
fi

# Check if test data exists (assuming we have some HMM files for testing)
if [ ! -d "./hmms" ] || [ -z "$(ls -A ./hmms 2>/dev/null)" ]; then
    echo "Warning: No HMM files found in ./hmms directory for testing"
    echo "Please ensure you have HMM files available for PRC comparison testing"
    echo
fi

# Function to run PRC with specific thread count and measure time
run_prc_test() {
    local threads=$1
    local test_name=$2
    
    echo "--- Testing with $threads threads ($test_name) ---"
    export OMP_NUM_THREADS=$threads
    
    # Create temporary output directories
    mkdir -p ./test_output_${threads}
    
    # Measure execution time
    start_time=$(date +%s.%N)
    
    # Run PRC distance calculation (assuming standard phmm-tree usage)
    # Note: This assumes you have HMM files set up for testing
    timeout 300 ./phmm-tree -prc ./hmms/ -matrix ./test_output_${threads}/matrix.txt -tree ./test_output_${threads}/tree.nwk 2>/dev/null
    exit_code=$?
    
    end_time=$(date +%s.%N)
    
    if [ $exit_code -eq 124 ]; then
        echo "Test timed out after 5 minutes"
        execution_time="TIMEOUT"
    elif [ $exit_code -ne 0 ]; then
        echo "Test failed with exit code: $exit_code"
        execution_time="FAILED"
    else
        execution_time=$(echo "$end_time - $start_time" | bc -l)
        printf "Execution time: %.3f seconds\n" $execution_time
    fi
    
    # Count number of comparisons if matrix file exists
    if [ -f "./test_output_${threads}/matrix.txt" ]; then
        comparisons=$(wc -l < "./test_output_${threads}/matrix.txt" 2>/dev/null || echo "Unknown")
        echo "Matrix size: $comparisons lines"
    fi
    
    echo "Thread utilization: $threads cores"
    echo
    
    # Store results for comparison
    echo "$threads,$execution_time" >> prc_performance_results.csv
    
    # Clean up test output
    rm -rf ./test_output_${threads}
}

# Initialize results file
echo "Threads,ExecutionTime" > prc_performance_results.csv

# Get number of CPU cores
cpu_cores=$(nproc)
echo "Detected CPU cores: $cpu_cores"
echo "OpenMP support: $(if [[ "$OPENMP_ENABLED" == "1" ]]; then echo "Enabled"; else echo "Checking..."; fi)"

# Test OpenMP functionality
export OMP_NUM_THREADS=2
echo "Testing OpenMP: $(echo 'int main(){return 0;}' | gcc -xc - -fopenmp -o /tmp/omp_test 2>&1 && echo "Available" || echo "Not Available")"
rm -f /tmp/omp_test

echo

# Run tests with different thread counts
echo "Starting PRC parallelization performance tests..."
echo

# Test with single thread (baseline)
run_prc_test 1 "Single-threaded baseline"

# Test with 2 threads
run_prc_test 2 "Dual-threaded"

# Test with 4 threads (if available)
if [ $cpu_cores -ge 4 ]; then
    run_prc_test 4 "Quad-threaded"
fi

# Test with maximum threads
if [ $cpu_cores -gt 4 ]; then
    run_prc_test $cpu_cores "Maximum threads"
fi

# Test with over-subscription (max threads + 2)
if [ $cpu_cores -ge 2 ]; then
    over_threads=$((cpu_cores + 2))
    run_prc_test $over_threads "Over-subscribed"
fi

echo "=== Performance Test Results Summary ==="
echo
echo "Results saved to: prc_performance_results.csv"

if [ -f prc_performance_results.csv ]; then
    echo "Thread Count | Execution Time"
    echo "-------------|---------------"
    tail -n +2 prc_performance_results.csv | while IFS=, read threads time; do
        printf "%12s | %s\n" "$threads" "$time"
    done
fi

echo
echo "=== Performance Analysis ==="

# Calculate speedup if we have valid results
baseline_time=""
max_speedup=""

if [ -f prc_performance_results.csv ]; then
    baseline_time=$(grep "^1," prc_performance_results.csv | cut -d, -f2)
    
    if [[ "$baseline_time" =~ ^[0-9.]+$ ]]; then
        echo "Baseline (1 thread): ${baseline_time}s"
        
        tail -n +2 prc_performance_results.csv | while IFS=, read threads time; do
            if [[ "$time" =~ ^[0-9.]+$ ]] && [ "$threads" != "1" ]; then
                speedup=$(echo "scale=2; $baseline_time / $time" | bc -l)
                efficiency=$(echo "scale=1; $speedup / $threads * 100" | bc -l)
                printf "%d threads: %.2fx speedup (%.1f%% efficiency)\n" $threads $speedup $efficiency
            fi
        done
    fi
fi

echo
echo "=== System Information ==="
echo "CPU Cores: $cpu_cores"
echo "CPU Info: $(lscpu | grep "Model name" | cut -d: -f2 | xargs || echo "Not available")"
echo "Memory: $(free -h | grep "Mem:" | awk '{print $2}' || echo "Not available")"
echo "OpenMP Version: $(echo | gcc -fopenmp -dM -E - | grep _OPENMP | cut -d' ' -f3 || echo "Not detected")"

echo
echo "Performance test completed!"
echo "Note: Actual performance depends on dataset size, system resources, and PRC complexity."
