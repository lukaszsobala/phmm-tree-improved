# PHMM-Tree Parallelization Guide

This version of PHMM-Tree has been enhanced with OpenMP parallelization for improved performance on multi-core systems.

## Supported Algorithms

### Fitch-Margoliash Method (fitch.c)
- Global rearrangements parallelized
- Distance calculations optimized
- Tree evaluation accelerated

### Kitsch Method (kitsch.c)  
- Fitch-Margoliash with contemporary tips
- Global optimization loops parallelized
- Statistical calculations accelerated

### Neighbor-Joining Method (neighbor.c)
- Matrix operations parallelized
- Distance calculations optimized
- Minimum finding operations accelerated

### UPGMA Method (upgma.c)
- Hierarchical clustering parallelized
- Matrix processing optimized
- Distance accumulation accelerated

## Features

### Automatic Thread Detection
- The program automatically detects the number of available CPU cores
- Sets the optimal number of threads without user intervention
- Falls back to sequential execution if OpenMP is not available

### Parallelized Components
1. **Global Rearrangements**: The most computationally intensive operations across all algorithms
2. **Distance Calculations**: Parallel computation of node and matrix distances
3. **Matrix Operations**: Parallel processing of distance matrices (Neighbor-Joining, UPGMA)
4. **Node Initialization**: Parallel setup of tree nodes and data structures
5. **Statistical Calculations**: Parallel computation of quality metrics and variance

### Performance Benefits
- Significant speedup for large datasets (>100 species)
- Optimal scaling up to the number of CPU cores
- Minimal overhead for small datasets (automatic threshold detection)

## Compilation

### Using Make
```bash
make                    # Compile with automatic OpenMP detection
make info              # Show compilation information
make test-parallel     # Test with different thread counts
make clean            # Clean build files
```

### Using CMake
```bash
mkdir build && cd build
cmake ..              # Configure with automatic OpenMP detection
make                  # Compile
```

### Manual Compilation
```bash
# With OpenMP support
gcc -O3 -fopenmp -DOPENMP_ENABLED *.c -o phmm-tree -lm -lstdc++

# Without OpenMP (sequential)
gcc -O3 *.c -o phmm-tree -lm -lstdc++
```

## Runtime Control

### Environment Variables
- `OMP_NUM_THREADS=N`: Set number of threads (overrides auto-detection)
- `OMP_SCHEDULE="type,chunk"`: Set scheduling policy

### Examples
```bash
# Use 4 threads
OMP_NUM_THREADS=4 ./phmm-tree

# Use dynamic scheduling
OMP_SCHEDULE="dynamic,1" ./phmm-tree

# Sequential execution
OMP_NUM_THREADS=1 ./phmm-tree
```

## Performance Guidelines

### When Parallelization Helps Most
- Large datasets (>100 species)
- Global rearrangement operations
- Complex tree topologies
- Multi-core systems with ample memory

### When to Use Sequential
- Small datasets (<50 species)
- Memory-constrained systems
- Single-core systems
- Debugging purposes

### Optimal Thread Count
- Generally equals the number of CPU cores
- For memory-intensive tasks, may be less than core count
- Hyperthreading may or may not help (dataset dependent)

## Technical Details

### Thread Safety
- Each thread operates on independent data structures
- Critical sections protect shared tree updates
- Reduction operations handle parallel accumulation

### Memory Usage
- Each thread may require additional memory for local trees
- Memory usage scales with thread count for global operations
- Automatic threshold prevents excessive memory usage

### Load Balancing
- Dynamic scheduling adapts to varying computation times
- Work-stealing helps with irregular workloads
- Chunk size automatically optimized for different operations

## Troubleshooting

### Performance Issues
```bash
# Check if OpenMP is working
export OMP_DISPLAY_ENV=TRUE
./phmm-tree

# Monitor thread usage
htop  # or top -H
```

### Compilation Issues
```bash
# Check compiler support
gcc -fopenmp --version

# Verify OpenMP availability
echo | gcc -fopenmp -E -dM - | grep -i openmp
```

### Runtime Issues
- Reduce thread count if memory errors occur
- Use `OMP_NUM_THREADS=1` for debugging
- Check system limits with `ulimit -a`

## Algorithm-Specific Notes

### Fitch-Margoliash Method
- Tree evaluation parallelized across interior nodes
- Distance matrix computations use SIMD-friendly loops
- Branch length optimization uses parallel iterations

### Performance Scaling
- Expected speedup: 2-8x on typical multi-core systems
- Best case: Linear scaling with core count
- Actual performance depends on dataset characteristics

## Backwards Compatibility
- Maintains full compatibility with original interface
- Produces identical results to sequential version
- Falls back gracefully when OpenMP unavailable
