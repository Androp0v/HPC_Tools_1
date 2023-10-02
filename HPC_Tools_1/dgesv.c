#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if defined(__APPLE__)
    #define USE_APPLE_DISPATCH
    #include <dispatch/dispatch.h>
    #include <os/log.h>
    #include <os/signpost.h>
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
    int dispatch_stride = 64;
    
    os_log_t log_handle = os_log_create("com.raulmonton.dgesv", OS_LOG_CATEGORY_POINTS_OF_INTEREST);
    os_signpost_id_t elimination_id = os_signpost_id_generate(log_handle);
    os_signpost_id_t inner_pivot_search = os_signpost_id_generate(log_handle);
    os_signpost_id_t inner_swapping_id = os_signpost_id_generate(log_handle);
    os_signpost_id_t inner_elimination_id = os_signpost_id_generate(log_handle);
    os_signpost_id_t matrix_coalesced_update_id = os_signpost_id_generate(log_handle);
    os_signpost_id_t backsubstitution_id = os_signpost_id_generate(log_handle);
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

    os_signpost_interval_begin(log_handle, elimination_id, "Elimination loop");
    
    int block_size = 16;
    for(int block=0; block<(n/block_size); block++) {
        
        int block_start_k = block * block_size;
        int block_end_k = (block + 1) * block_size;
        
        for(int k = block_start_k; k < block_end_k; k++) {
            
            // MARK: - Pivot search
            
            os_signpost_interval_begin(log_handle, inner_pivot_search, "Inner pivot search loop");
            big = 0.0;
            // Search for pivot element only below the current row, up to the end of the block.
            int end = k + 16;
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
            os_signpost_interval_end(log_handle, inner_pivot_search, "Inner pivot search loop");
            
            // Check whether rows need to be swapped.
            if (k != imax) {
                os_signpost_interval_begin(log_handle, inner_swapping_id, "Inner swapping loop");
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
                os_signpost_interval_end(log_handle, inner_swapping_id, "Inner swapping loop");
            }
            
            // MARK: - Elimination
            // Update only the rows in this block and the next block
            os_signpost_interval_begin(log_handle, inner_elimination_id, "Inner elimination loop");
            
            double pivot_value = a[k*n + k];
            
            int next_block_end_k = block_start_k + (2 * block_size);
            if (next_block_end_k > n) {
                next_block_end_k = n;
            }
            for (int i = k+1; i < next_block_end_k; i++) {
                dum = a[i*n + k] / pivot_value;
                for (int j = k+1; j<n; j++) {
                    a[i*n + j] -= a[k*n + j] * dum;
                }
                for (int j = 0; j<n; j++) {
                    b[i*n + j] -= b[k*n + j] * dum;
                }
            }
            os_signpost_interval_end(log_handle, inner_elimination_id, "Inner elimination loop");
        }
        
        // MARK: - Coalesced matrix update
        // Update the rest of the matrix
        os_signpost_interval_begin(log_handle, matrix_coalesced_update_id, "Matrix coalesced update");
        int next_block_end_k = block_start_k + (2 * block_size);
        if (next_block_end_k > n) {
            next_block_end_k = n;
        }
        
        #if defined(USE_APPLE_DISPATCH)
        // On Apple platforms, execute concurrently in all available cores using Apple's
        // dispatch library.
        int number_of_threadgroups = (n - next_block_end_k) / dispatch_stride;
        if ((n - next_block_end_k) % dispatch_stride != 0) {
            number_of_threadgroups += 1;
        }
        dispatch_apply(number_of_threadgroups, concurrent_queue, ^(size_t idx) {
            for(int k = block_start_k; k < block_end_k; k++) {
                
                double pivot_value = a[k*n + k];
                
                int start = next_block_end_k + idx * dispatch_stride;
                int end = next_block_end_k + (idx + 1) * dispatch_stride;
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
        #else
        // On non-Apple platforms, use OpenMP to parallelize the loop instead.
        for(int k = block_start_k; k < block_end_k; k++) {
            
            double pivot_value = a[k*n + k];
            
            int next_block_end_k = block_start_k + (2 * block_size);
            if (next_block_end_k > n) {
                next_block_end_k = n;
            }
            
            #pragma omp parallel for
            for (int i = next_block_end_k; i < n; i++) {
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
        os_signpost_interval_end(log_handle, matrix_coalesced_update_id, "Matrix coalesced update");
    }
    os_signpost_interval_end(log_handle, elimination_id, "Elimination loop");
    
    // MARK: - Backsubstitution
    os_signpost_interval_begin(log_handle, backsubstitution_id, "Backsubstitution loop");
    // On Apple platforms, execute concurrently in all available cores using Apple's
    // dispatch library.
    #if defined(USE_APPLE_DISPATCH)
    int number_of_threadgroups = n / dispatch_stride;
    if (n % dispatch_stride != 0) {
        number_of_threadgroups += 1;
    }
    dispatch_apply(number_of_threadgroups, concurrent_queue, ^(size_t idx) {
        int end = (idx + 1) * dispatch_stride;
        if (end > n) {
            end = n;
        }
        for (int i=n-1; i>=0; i--) {
            for (int k = idx * dispatch_stride; k < end; k++) {
                double sum = b[i*n + k];
                #pragma clang loop vectorize(enable) interleave(enable)
                for(int j=i+1; j<n; j++) {
                    sum -= a[i*n + j] * b[j*n + k];
                }
                b[i*n + k] = sum / a[i*n + i];
            }
        }
    });
    // On non-Apple platforms, use OpenMP to parallelize the loop instead.
    #else
    #pragma omp parallel for
    // On non-Apple platforms, use OpenMP to parallelize the loop instead.
    for (int i=n-1; i>=0; i--) {
        for (int k=0; k<n; k++) {
            double sum = b[i*n + k];
            for(int j=i+1; j<n; j++) {
                sum -= a[i*n + j] * b[j*n + k];
            }
            b[i*n + k] = sum / a[i*n + i];
        }
    }
    #endif
    os_signpost_interval_end(log_handle, backsubstitution_id, "Backsubstitution loop");
    
    return 0;
}
