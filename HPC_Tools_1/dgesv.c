#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if defined(__APPLE__)
    #define USE_APPLE_DISPATCH
    #define USE_SIGNPOSTING
    #include <dispatch/dispatch.h>
    #include <pthread.h>
    #if defined(USE_SIGNPOSTING)
        #include <os/log.h>
        #include <os/signpost.h>
    #endif
#endif

#include "dgesv.h"

int *generate_int_vector(unsigned int length) {
  int *vector = (int *) malloc(sizeof(int) * length);
  return vector;
}

double *generate_double_vector(unsigned int length) {
  double *vector = (double *) malloc(sizeof(double) * length);
  return vector;
}

void swap(double *pa, double *pb) {
    double temp = *pa;
    *pa = *pb;
    *pb = temp;
}

void print_matrix(int n, double *a) {
    for (int i=0; i<n; i++) {
        for(int j=0; j<n; j++){
             printf("%lf     ", a[i*n + j]);
        }
        printf("\n");
    }
}

int my_dgesv(int n, double *a, double *b) {
    
    double *vv = generate_double_vector(n);
    
    double big, dum, temp;
    int imax;
    
    #if defined(USE_APPLE_DISPATCH)
    dispatch_queue_t concurrent_queue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0);
    int dispatch_stride = 16;
    dispatch_group_t coalesced_matrix_group = dispatch_group_create();
    #endif
    
    #if defined(USE_SIGNPOSTING)
    os_log_t log_handle = os_log_create("com.raulmonton.dgesv", OS_LOG_CATEGORY_POINTS_OF_INTEREST);
    os_signpost_id_t elimination_id = os_signpost_id_generate(log_handle);
    os_signpost_id_t inner_pivot_search = os_signpost_id_generate(log_handle);
    os_signpost_id_t inner_swapping_id = os_signpost_id_generate(log_handle);
    os_signpost_id_t inner_elimination_id = os_signpost_id_generate(log_handle);
    os_signpost_id_t matrix_coalesced_next_block_update_id = os_signpost_id_generate(log_handle);
    os_signpost_id_t matrix_coalesced_update_id = os_signpost_id_generate(log_handle);
    os_signpost_id_t triangular_elimination_id = os_signpost_id_generate(log_handle);
    #endif
    
    for (int i = 0; i<n; i++) {
        big = 0.0;
        for (int j = 0; j<n; j++) {
            if ((temp = fabs(a[i*n + j])) > big) {
                big = temp;
            }
        }
        
        // If the pivot element is zero, either the matrix is singular or we fucked up.
        if (big == 0.0) {
            printf("Singular matrix");
            exit(0);
        }
        
        vv[i] = 1.0 / big;
    }
    
    #if defined(USE_SIGNPOSTING)
    os_signpost_interval_begin(log_handle, elimination_id, "Elimination loop");
    #endif
    
    int block_size = 16;
    int block_end = n / block_size;
    if (n % block_size != 0) {
        block_end += 1;
    }
    for(int block = 0; block < block_end; block++) {
        
        int block_start_k = block * block_size;
        int block_end_k = (block + 1) * block_size;
        if (block_end_k > n) {
            block_end_k = n;
        }
        
        int next_block_end_k = block_start_k + (2 * block_size);
        if (next_block_end_k > n) {
            next_block_end_k = n;
        }
        int next_next_block_end_k = block_start_k + (3 * block_size);
        if (next_next_block_end_k > n) {
            next_next_block_end_k = n;
        }
        
        for(int k = block_start_k; k < block_end_k; k++) {
            
            // MARK: - Pivot search
            
            #if defined(USE_SIGNPOSTING)
            os_signpost_interval_begin(log_handle, inner_pivot_search, "Inner pivot search loop");
            #endif
            big = 0.0;
            // Search for pivot element only below the current row, on the current block or the next
            // one, but not any further.
            int end = k + block_size;
            if (end > n) {
                end = n;
            }
            for(int i=k; i<end; i++) {
                double temp = vv[i] * fabs(a[i*n + k]);
                if (temp > big) {
                    big = temp;
                    imax = i;
                }
            }
            #if defined(USE_SIGNPOSTING)
            os_signpost_interval_end(log_handle, inner_pivot_search, "Inner pivot search loop");
            #endif
            
            // MARK: - Row swapping
            // Check whether rows need to be swapped.
            if (k!=imax) {
                #if defined(USE_SIGNPOSTING)
                os_signpost_interval_begin(log_handle, inner_swapping_id, "Inner swapping loop");
                #endif
                
                // Swap rows...
                #pragma clang loop vectorize(enable) interleave(enable)
                for(int j=k; j<n; j++) {
                    swap(&a[imax*n + j], &a[k*n + j]);
                }
                #pragma clang loop vectorize(enable) interleave(enable)
                for(int j=0; j<n; j++) {
                    swap(&b[imax*n + j], &b[k*n + j]);
                }
                
                vv[imax] = vv[k];

                #if defined(USE_SIGNPOSTING)
                os_signpost_interval_end(log_handle, inner_swapping_id, "Inner swapping loop");
                #endif
            }
            
            // MARK: - Elimination
            // Update only the rows in this block and the next block
            #if defined(USE_SIGNPOSTING)
            os_signpost_interval_begin(log_handle, inner_elimination_id, "Inner elimination loop");
            #endif
            
            double pivot_value = a[k*n + k];
            
            for (int i = k+1; i < next_block_end_k; i++) {
                dum = a[i*n + k] / pivot_value;
                for (int j = k+1; j<n; j++) {
                    a[i*n + j] -= a[k*n + j] * dum;
                }
                for (int j = 0; j<n; j++) {
                    b[i*n + j] -= b[k*n + j] * dum;
                }
            }
            #if defined(USE_SIGNPOSTING)
            os_signpost_interval_end(log_handle, inner_elimination_id, "Inner elimination loop");
            #endif
        }
        
        // MARK: - Update the rest of the matrix
        
        if (block < (block_end - 2)) {
            
            // MARK: - Next block update
            // The previous big loop sees two blocks ahead. Here, we update the 3rd block.
            
            #if defined(USE_APPLE_DISPATCH)
            // This is essentially a memory fence, that waits until the coalesced matrix
            // update has finished updating the trailing blocks of the matrix.
            dispatch_group_wait(coalesced_matrix_group, DISPATCH_TIME_FOREVER);
            #endif
            
            #if defined(USE_SIGNPOSTING)
            os_signpost_interval_begin(log_handle, matrix_coalesced_next_block_update_id, "Next block coalesced update");
            #endif
            for(int k = block_start_k; k < block_end_k; k++) {
                double pivot_value = a[k*n + k];
                #pragma omp parallel for
                for (int i = next_block_end_k; i < next_next_block_end_k; i++) {
                    dum = a[i*n + k] / pivot_value;
                    for (int j = k+1; j<n; j++) {
                        a[i*n + j] -= a[k*n + j] * dum;
                    }
                    for (int j = 0; j<n; j++) {
                        b[i*n + j] -= b[k*n + j] * dum;
                    }
                }
            }
            #if defined(USE_SIGNPOSTING)
            os_signpost_interval_end(log_handle, matrix_coalesced_next_block_update_id, "Next block coalesced update");
            #endif
            
            // MARK: - Coalesced matrix update
            // Update the rest of the matrix. The elimination inner loop updates the current and
            // next blocks. The previous code segment updates the 3rd code block. Here, we update
            // from the 4th block until the end of the matrix.
            if (block < (block_end - 3)) {
                
                #if defined(USE_SIGNPOSTING)
                os_signpost_interval_begin(log_handle, matrix_coalesced_update_id, "Matrix coalesced update");
                #endif
                
                #if defined(USE_APPLE_DISPATCH)
                // On Apple platforms, execute concurrently in all available cores using Apple's
                // dispatch library.
                
                int number_of_threads = (n - next_next_block_end_k) / dispatch_stride;
                if ((n - next_next_block_end_k) % dispatch_stride != 0) {
                    number_of_threads += 1;
                }
                
                // The function dispatch_group_async schedules the code block inside asynchronously:
                // the code can continue executing at the call site, going to the next iteration of
                // the loop, while the code block inside executes in the background.
                // Then, if we need to touch any part of the matrix that has a data dependency on
                // this block, we'll call dispatch_group_wait to block execution until this block
                // has finished.
                dispatch_group_async(coalesced_matrix_group, concurrent_queue, ^{
                    dispatch_apply(number_of_threads, concurrent_queue, ^(size_t idx) {
                        for(int k = block_start_k; k < block_end_k; k++) {
                            double pivot_value = a[k*n + k];
                            
                            int start = next_next_block_end_k + idx * dispatch_stride;
                            int end = next_next_block_end_k + (idx + 1) * dispatch_stride;
                            if (end > n) {
                                end = n;
                            }
                            
                            for (int i = start; i < end; i++) {
                                double dum = a[i*n + k] / pivot_value;
                                for (int j = k+1; j<n; j++) {
                                    a[i*n + j] -= a[k*n + j] * dum;
                                }
                                for (int j = 0; j<n; j++) {
                                    b[i*n + j] -= b[k*n + j] * dum;
                                }
                            }
                        }
                    });
                    #if defined(USE_SIGNPOSTING)
                    os_signpost_interval_end(log_handle, matrix_coalesced_update_id, "Matrix coalesced update");
                    #endif
                });
                #else
                // On non-Apple platforms, use OpenMP to parallelize the loop instead.
                for(int k = block_start_k; k < block_end_k; k++) {
                    
                    double pivot_value = a[k*n + k];

                    #pragma omp parallel for
                    for (int i = next_next_block_end_k; i < n; i++) {
                        dum = a[i*n + k] / pivot_value;
                        for (int j = k+1; j<n; j++) {
                            a[i*n + j] -= a[k*n + j] * dum;
                        }
                        for (int j = 0; j<n; j++) {
                            b[i*n + j] -= b[k*n + j] * dum;
                        }
                    }
                }
                #endif
            }
        }
    }
    #if defined(USE_SIGNPOSTING)
    os_signpost_interval_end(log_handle, elimination_id, "Elimination loop");
    #endif
        
    // MARK: - Elimination up
    
    #if defined(USE_SIGNPOSTING)
    os_signpost_interval_begin(log_handle, triangular_elimination_id, "Triangular elimination loop");
    #endif
    #if defined(USE_APPLE_DISPATCH)
    // On Apple platforms, execute concurrently in all available cores using Apple's
    // dispatch library.
    for (int k = 0; k < n; k++) {
        double pivot_value = a[k*n + k];
        int number_of_threads = k / dispatch_stride;
        if (k % dispatch_stride != 0) {
            number_of_threads += 1;
        }
        dispatch_apply(number_of_threads, concurrent_queue, ^(size_t idx) {
            int start = idx * dispatch_stride;
            int end = (idx + 1) * dispatch_stride;
            if (end > k) {
                end = k;
            }
            for (int i = start; i < end; i++) {
                double dum = a[i*n + k] / pivot_value;
                for (int j = k; j<n; j++) {
                    a[i*n + j] -= a[k*n + j] * dum;
                }
                for (int j = 0; j<n; j++) {
                    b[i*n + j] -= b[k*n + j] * dum;
                }
            }
        });
    }
    
    int number_of_threads = n / dispatch_stride;
    if (n % dispatch_stride != 0) {
        number_of_threads += 1;
    }
    dispatch_apply(number_of_threads, concurrent_queue, ^(size_t idx) {
        int start = idx * dispatch_stride;
        int end = (idx + 1) * dispatch_stride;
        if (end > n) {
            end = n;
        }
        for (int i = start; i < end; i ++) {
            double diagonal = a[i*n + i];
            for (int j = 0; j < n; j++) {
                b[i*n + j] /= diagonal;
            }
        }
    });
    #else
    // On non-Apple platforms, use OpenMP to parallelize the loop instead.
    for (int k = 0; k < n; k++) {
        double pivot_value = a[k*n + k];
        #pragma omp parallel for
        for (int i = 0; i < k; i++) {
            dum = a[i*n + k] / pivot_value;
            for (int j = k; j<n; j++) {
                a[i*n + j] -= a[k*n + j] * dum;
            }
            for (int j = 0; j<n; j++) {
                b[i*n + j] -= b[k*n + j] * dum;
            }
        }
    }
    for (int i = 0; i < n; i ++) {
        double diagonal = a[i*n + i];
        for (int j = 0; j < n; j++) {
            b[i*n + j] /= diagonal;
        }
    }
    #endif
    #if defined(USE_SIGNPOSTING)
    os_signpost_interval_end(log_handle, triangular_elimination_id, "Triangular elimination loop");
    #endif

    return 0;
}
