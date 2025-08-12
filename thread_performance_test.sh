#!/bin/bash

# Thread Performance Test Script for pHMM-Tree
# Tests various thread configurations with a specified dataset directory

# Check if directory argument is provided
if [ $# -ne 1 ]; then
    echo "Usage: $0 <hmm_directory>"
    echo "  hmm_directory: Path to directory containing HMM files for testing"
    echo ""
    echo "Example: $0 test_hmms_reduced/"
    exit 1
fi

TEST_DIR="$1"

echo "=========================================="
echo "pHMM-Tree Thread Performance Test"
echo "Testing with dataset: $TEST_DIR"
echo "=========================================="

# Check if test dataset exists
if [ ! -d "$TEST_DIR" ]; then
    echo "Error: Directory '$TEST_DIR' not found!"
    echo "Please provide a valid directory containing HMM files."
    exit 1
fi

# Check if directory contains HMM files
hmm_count=$(find "$TEST_DIR" -name "*.hmm" | wc -l)
if [ $hmm_count -eq 0 ]; then
    echo "Error: No .hmm files found in '$TEST_DIR'!"
    echo "Please provide a directory containing HMM files."
    exit 1
fi

echo "Found $hmm_count HMM files in $TEST_DIR"

# Create results directory
RESULTS_DIR="thread_test_results_$(date +%Y%m%d_%H%M%S)"
mkdir -p "$RESULTS_DIR"

echo "Results will be saved in: $RESULTS_DIR"
echo ""

# Initialize results table file
TABLE_FILE="$RESULTS_DIR/timing_results_table.txt"
echo "Thread Performance Results for Dataset: $TEST_DIR" > "$TABLE_FILE"
echo "Test Date: $(date)" >> "$TABLE_FILE"
echo "HMM Files: $hmm_count" >> "$TABLE_FILE"
echo "" >> "$TABLE_FILE"

# Table header
printf "%-15s %-15s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s\n" \
    "PRC_Threads" "Phylo_Threads" "PRC_Time" "Kitsch_Time" "Fitch_Time" "NJ_Time" "UPGMA_Time" "Total_Phylo" "Total_Time" "Status" >> "$TABLE_FILE"
echo "--------------------------------------------------------------------------------------------------------------------------------------------------------" >> "$TABLE_FILE"

# Function to run test with specific thread configuration
run_test() {
    local prc_threads=$1
    local phylo_threads=$2
    local test_name="prc${prc_threads}_phylo${phylo_threads}"
    
    echo "----------------------------------------"
    echo "Test: PRC threads=$prc_threads, Phylo threads=$phylo_threads"
    echo "Started at: $(date)"
    echo "----------------------------------------"
    
    # Clean up any previous test directories
    rm -rf prc_hmms_mode_*
    
    # Run the test and capture timing
    start_time=$(date +%s.%N)
    
    ./phmm-tree -prc -hmms -prc_threads $prc_threads -phylo_threads $phylo_threads "$TEST_DIR" \
        > "$RESULTS_DIR/${test_name}.log" 2>&1
    
    exit_code=$?
    end_time=$(date +%s.%N)
    
    # Calculate elapsed time
    elapsed=$(echo "$end_time - $start_time" | bc -l)
    
    if [ $exit_code -eq 0 ]; then
        echo "✓ Test completed successfully in ${elapsed} seconds"
        
        # Extract individual algorithm timings from the log
        prc_time=$(grep "PRC analysis completed in:" "$RESULTS_DIR/${test_name}.log" | sed -n 's/.*completed in: \([0-9.]*\) seconds.*/\1/p')
        kitsch_time=$(grep "Kitsch analysis completed in:" "$RESULTS_DIR/${test_name}.log" | sed -n 's/.*completed in: \([0-9.]*\) seconds.*/\1/p')
        fitch_time=$(grep "Fitch-Margoliash analysis completed in:" "$RESULTS_DIR/${test_name}.log" | sed -n 's/.*completed in: \([0-9.]*\) seconds.*/\1/p')
        nj_time=$(grep "Neighbor-Joining analysis completed in:" "$RESULTS_DIR/${test_name}.log" | sed -n 's/.*completed in: \([0-9.]*\) seconds.*/\1/p')
        upgma_time=$(grep "UPGMA analysis completed in:" "$RESULTS_DIR/${test_name}.log" | sed -n 's/.*completed in: \([0-9.]*\) seconds.*/\1/p')
        total_phylo_time=$(grep "Phylogenetic tree building completed in:" "$RESULTS_DIR/${test_name}.log" | sed -n 's/.*completed in: \([0-9.]*\) seconds.*/\1/p')
        
        # Use default values if extraction fails
        prc_time=${prc_time:-"N/A"}
        kitsch_time=${kitsch_time:-"N/A"}
        fitch_time=${fitch_time:-"N/A"}
        nj_time=${nj_time:-"N/A"}
        upgma_time=${upgma_time:-"N/A"}
        total_phylo_time=${total_phylo_time:-"N/A"}
        
        # Format elapsed time to 3 decimal places (handle long decimal numbers properly)
        elapsed_formatted=$(echo "$elapsed" | awk '{printf "%.3f", $1}')
        
        # Add row to results table
        printf "%-15s %-15s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s\n" \
            "$prc_threads" "$phylo_threads" "$prc_time" "$kitsch_time" "$fitch_time" \
            "$nj_time" "$upgma_time" "$total_phylo_time" "$elapsed_formatted" "SUCCESS" >> "$TABLE_FILE"
        
        # Create detailed report for this test
        echo "=== DETAILED TIMING REPORT for $test_name ===" > "$RESULTS_DIR/${test_name}_report.txt"
        echo "Thread Configuration: PRC=$prc_threads, Phylo=$phylo_threads" >> "$RESULTS_DIR/${test_name}_report.txt"
        echo "Total Execution Time: ${elapsed_formatted}s" >> "$RESULTS_DIR/${test_name}_report.txt"
        echo "" >> "$RESULTS_DIR/${test_name}_report.txt"
        echo "Individual Algorithm Timings:" >> "$RESULTS_DIR/${test_name}_report.txt"
        echo "  PRC Analysis: ${prc_time}s" >> "$RESULTS_DIR/${test_name}_report.txt"
        echo "  Kitsch Analysis: ${kitsch_time}s" >> "$RESULTS_DIR/${test_name}_report.txt"
        echo "  Fitch-Margoliash Analysis: ${fitch_time}s" >> "$RESULTS_DIR/${test_name}_report.txt"
        echo "  Neighbor-Joining Analysis: ${nj_time}s" >> "$RESULTS_DIR/${test_name}_report.txt"
        echo "  UPGMA Analysis: ${upgma_time}s" >> "$RESULTS_DIR/${test_name}_report.txt"
        echo "  Total Phylogenetic Time: ${total_phylo_time}s" >> "$RESULTS_DIR/${test_name}_report.txt"
        
        # Show brief summary immediately
        echo "  Timings: PRC=${prc_time}s, Kitsch=${kitsch_time}s, Fitch=${fitch_time}s, NJ=${nj_time}s, UPGMA=${upgma_time}s"
        
    else
        echo "✗ Test failed with exit code: $exit_code"
        printf "%-15s %-15s %-10s %-10s %-10s %-10s %-10s %-10s %-10s %-10s\n" \
            "$prc_threads" "$phylo_threads" "FAILED" "FAILED" "FAILED" "FAILED" "FAILED" "FAILED" "FAILED" "FAILED" >> "$TABLE_FILE"
    fi
    
    echo ""
}

# Start logging
echo "Thread Performance Test Started: $(date)" > "$RESULTS_DIR/summary.txt"
echo "Dataset: $TEST_DIR" >> "$RESULTS_DIR/summary.txt"
echo "HMM Files: $hmm_count" >> "$RESULTS_DIR/summary.txt"
echo "========================================" >> "$RESULTS_DIR/summary.txt"
echo "" >> "$RESULTS_DIR/summary.txt"

echo "Phase 1: Testing PRC thread variations (phylo threads = auto-detect)"
echo "===================================================================="

# Test PRC thread variations with phylo threads on auto-detect (0)
for prc_threads in 1 2 4 6 12; do
    run_test $prc_threads 0
done

echo ""
echo "Phase 2: Testing Phylogenetic thread variations (PRC threads = 12)"
echo "=================================================================="

# Test phylogenetic thread variations with PRC threads fixed at 12
for phylo_threads in 1 2 4 6 12; do
    run_test 12 $phylo_threads
done

echo ""
echo "=========================================="
echo "All tests completed!"
echo "=========================================="

# Finalize table
echo "" >> "$TABLE_FILE"
echo "Test completion time: $(date)" >> "$TABLE_FILE"

# Display the results table
echo ""
echo "=================== PERFORMANCE RESULTS TABLE ==================="
cat "$TABLE_FILE"
echo "=================================================================="
echo ""

# Display summary of findings
echo "PERFORMANCE ANALYSIS SUMMARY:"
echo "=============================="

# Show best performing configuration if we have successful tests
best_line=$(grep -v "FAILED" "$TABLE_FILE" | grep -v "PRC_Threads" | grep -v "^$" | grep -v "^Thread" | grep -v "^Test Date" | grep -v "^HMM Files" | grep -v "^---" | sort -k9 -n | head -1)

if [ -n "$best_line" ]; then
    echo "Best overall performance: $best_line"
    echo ""
fi

# Count successful and failed tests
successful_tests=$(grep -c "SUCCESS" "$TABLE_FILE")
failed_tests=$(grep -c "FAILED" "$TABLE_FILE")

echo "Tests completed: $((successful_tests + failed_tests))"
echo "Successful: $successful_tests"
echo "Failed: $failed_tests"
echo ""

echo "Detailed results are available in: $RESULTS_DIR"
echo "  - Timing table: $TABLE_FILE"
echo "  - Individual test reports: ${RESULTS_DIR}/*_report.txt"
echo "  - Full logs: ${RESULTS_DIR}/*.log"

echo ""
echo "To view the results table:"
echo "  cat $TABLE_FILE"
echo ""
echo "Test script completed successfully!"
