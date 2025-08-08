/*
 * Test program to verify OpenMP parallelization in PHMM-Tree
 * Compile with: gcc -fopenmp -DOPENMP_ENABLED test_parallel.c -o test_parallel
 */

#include <stdio.h>
#include <stdlib.h>

#ifdef _OPENMP
#include <omp.h>
#endif

void test_thread_detection() {
    printf("=== Thread Detection Test ===\n");
    
#ifdef _OPENMP
    int max_threads = omp_get_max_threads();
    int num_procs = omp_get_num_procs();
    
    printf("OpenMP is available\n");
    printf("Maximum threads: %d\n", max_threads);
    printf("Number of processors: %d\n", num_procs);
    
    // Test parallel region
    printf("\nTesting parallel region:\n");
    #pragma omp parallel
    {
        int thread_id = omp_get_thread_num();
        int total_threads = omp_get_num_threads();
        
        #pragma omp critical
        {
            printf("Thread %d of %d running\n", thread_id, total_threads);
        }
    }
    
    // Test parallel for loop
    printf("\nTesting parallel for loop:\n");
    int i;
    #pragma omp parallel for schedule(static)
    for (i = 0; i < 10; i++) {
        #pragma omp critical
        {
            printf("Thread %d processing item %d\n", omp_get_thread_num(), i);
        }
    }
    
#else
    printf("OpenMP not available - sequential execution\n");
    printf("This is expected if compiled without -fopenmp\n");
#endif

    printf("\n=== Test Complete ===\n");
}

int main() {
    printf("PHMM-Tree Parallelization Test\n");
    printf("==============================\n");
    
    test_thread_detection();
    
    return 0;
}
