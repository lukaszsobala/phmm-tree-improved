# PRC Parallelization Documentation

## Overview
This document describes the parallelization of the PRC (Profile HMM Comparison) algorithm in phmm-tree to enable concurrent execution of multiple single-threaded PRC instances for faster all-vs-all distance matrix computation.

## Key Changes

### 1. prc_deal.cpp - Main Parallelization
- **Function**: `prc_each2()`
- **Strategy**: Converted nested sequential loops to OpenMP parallel processing of all HMM pairs
- **Implementation**: 
  - Pre-generates all HMM pairs to process
  - Uses `#pragma omp parallel for` with dynamic scheduling
  - Each PRC call runs as single-threaded with `OMP_NUM_THREADS=1` prefix

### 2. Thread-Safe Output Management
- **Progress Reporting**: Thread-safe progress output showing completion status
- **Format**: `[completed_count/total_pairs] hmm1_name vs hmm2_name (Thread thread_id)`
- **Synchronization**: `#pragma omp critical` sections for output and matrix updates

### 3. Matrix Update Safety  
- **Critical Section**: Matrix updates are protected with OpenMP critical sections
- **Function**: `matrix_get_each2_hmms_result_2()` called within critical region
- **Consistency**: Ensures thread-safe distance matrix population

## Technical Details

### Parallel Execution Flow
1. **Setup Phase**: 
   - Initialize HMM names and matrix structures
   - Generate list of all pairwise combinations
   
2. **Parallel Processing**:
   - Each thread processes different HMM pairs concurrently
   - Individual PRC calls run single-threaded to avoid conflicts
   - Results parsed and stored in thread-safe manner

3. **Synchronization**:
   - Critical sections protect shared data structures
   - Progress reporting synchronized to prevent output jumbling
   - Matrix updates coordinated across threads

### Command Construction
```cpp
// For hit-limited searches
str_cmd = "OMP_NUM_THREADS=1 prc -hits " + hit_num + " " + hmm1 + " " + hmm2;

// For standard searches  
str_cmd = "OMP_NUM_THREADS=1 prc " + hmm1 + " " + hmm2;
```

### Performance Benefits
- **Concurrent Execution**: Multiple PRC comparisons run simultaneously
- **CPU Utilization**: Full utilization of available CPU cores
- **Scalability**: Performance scales with number of available threads
- **Load Balancing**: Dynamic scheduling ensures optimal work distribution

### Output Format Changes
**Before**: 
```
PRC deal: 
hmm1_name  
. . . .
hmm2_name
. . . 
```

**After**:
```
PRC deal: Processing 15 pairwise comparisons...
[1/15] hmm1.hmm vs hmm2.hmm (Thread 0)
[2/15] hmm1.hmm vs hmm3.hmm (Thread 1)
[3/15] hmm2.hmm vs hmm3.hmm (Thread 0)
...
PRC deal completed: 15 comparisons processed.
```

## Usage Notes
- **Thread Detection**: Automatically uses available CPU cores
- **Thread Override**: Can be controlled with `OMP_NUM_THREADS` environment variable
- **PRC Requirements**: Each individual PRC call runs single-threaded to prevent conflicts
- **Output Clarity**: Thread-safe reporting prevents jumbled terminal output

## Limitations
- `prc_library()` function currently remains sequential (complex temporary file management)
- Memory usage increases with number of concurrent threads
- Requires sufficient system resources for concurrent PRC processes

## Future Enhancements
- Potential parallelization of `prc_library()` function
- Memory usage optimization for large datasets
- Enhanced error handling for concurrent execution scenarios
