#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#if defined(__APPLE__)
    #define USE_APPLE_DISPATCH
    #include <dispatch/dispatch.h>
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
    #endif
    
    for (int i = 0; i<n; i++) {
        big = 0.0;
        for (int j = 0; j<n; j++) {
            if ((temp = fabs(a[i*n + j])) > big) {
                big = temp;
            }
        }
        vv[i] = 1.0 / big;
    }

    for(int k=0; k<n; k++) {
        big = 0.0;
        // Search for pivot element only at or below the current row.
        for(int i=k; i<n; i++) {
            double temp = vv[i] * fabs(a[i*n + k]);
            if (temp > big) {
                big = temp;
                imax = i;
            }
        }
                
        // Check whether rows need to be swapped.
        if (k != imax) {
            // Swap rows...
            for(int j=0; j<n; j++) {
                swap(&a[imax*n + j], &a[k*n + j]);
                swap(&b[imax*n + j], &b[k*n + j]);
            }
            vv[imax] = vv[k];
        }

        // If the pivot element is zero, either the matrix is singular or we fucked up.
        if (a[k*n + k] == 0.0) {
            printf("Singular matrix");
            exit(0);
        }
        
        // MARK: - Elimination
        // On Apple platforms, execute concurrently in all available cores using Apple's
        // dispatch library.
        #if defined(USE_APPLE_DISPATCH)
        int number_of_threadgroups = ((n - (k+1)) / dispatch_stride) + 1;
        dispatch_apply(number_of_threadgroups, concurrent_queue, ^(size_t idx) {
            int start = k + 1 + idx * dispatch_stride;
            int end = k + 1 + (idx + 1) * dispatch_stride;
            if (end > n) {
                end = n;
            }
            for (int i = start; i<end; i++) {
                double dum = a[i*n + k] / a[k*n + k];
                #pragma clang loop vectorize(enable) interleave(enable)
                for (int j = k+1; j<n; j++) {
                    a[i*n + j] -= a[k*n + j] * dum;
                }
                #pragma clang loop vectorize(enable) interleave(enable)
                for (int j = 0; j<n; j++) {
                    b[i*n + j] -= b[k*n + j] * dum;
                }
            }
        });
        // On non-Apple platforms, use OpenMP to parallelize the loop instead.
        #else
        #pragma omp parallel for
        for (int i = k+1; i<n; i++) {
            temp = a[i*n + k] / a[k*n + k];
            for (int j = k+1; j<n; j++) {
                a[i*n + j] -= temp * a[k*n + j];
            }
            for (int j = 0; j<n; j++) {
                b[i*n + j] -= temp * b[k*n + j];
            }
        }
        #endif
    }
    
    // MARK: - Backsubstitution
    // On Apple platforms, execute concurrently in all available cores using Apple's
    // dispatch library.
    #if defined(USE_APPLE_DISPATCH)
    dispatch_apply((n / dispatch_stride) + 1, concurrent_queue, ^(size_t idx) {
        int end = (idx + 1) * dispatch_stride;
        if (end > n) {
            end = n;
        }
        for (int k = idx * dispatch_stride; k < end; k++) {
            for (int i=n-1; i>=0; i--) {
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
    for (int k=0; k<n; k++) {
        for (int i=n-1; i>=0; i--) {
            double sum = b[i*n + k];
            for(int j=i+1; j<n; j++) {
                sum -= a[i*n + j] * b[j*n + k];
            }
            b[i*n + k] = sum / a[i*n + i];
        }
    }
    #endif
        
    return 0;
}
