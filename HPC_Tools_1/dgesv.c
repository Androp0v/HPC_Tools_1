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

void swap(double *pa, double *pb) {
    double temp = *pa;
    *pa = *pb;
    *pb = temp;
}

int my_dgesv(int n, double *a, double *b) {
    
    int *ipiv = generate_int_vector(n);
    
    double big, dum, pivinv;
    int imax;
        
    #pragma clang loop vectorize(enable) interleave(enable)
    for(int j=0; j<n; j++) {
        ipiv[j] = 0;
    }
    
    #if defined(USE_APPLE_DISPATCH)
    dispatch_queue_t concurrentQueue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0);
    int dispatch_stride = 64;
    #endif
    
    for(int i=0; i<n; i++) {
        big = 0.0;
        for (int j=0; j<n; j++) {
            if (ipiv[j] == 0) {
                if (fabs(a[j*n + j]) >= big) {
                    big = fabs(a[j*n + j]);
                    imax = j;
                }
            }
        }
        ++(ipiv[imax]);

        if (a[imax*n + imax] == 0.0) {
            printf("Singular matrix");
            exit(0);
        }
        pivinv = 1.0 / a[imax*n + imax];
        a[imax*n + imax] = 1.0;
        
        #pragma clang loop vectorize(enable) interleave(enable)
        for(int l=0; l<n; l++) {
            a[imax*n + l] *= pivinv;
            b[imax*n + l] *= pivinv;
        }
                
        // On Apple platforms, execute concurrently in all available cores
        // using Apple's dispatch library.
        #if defined(USE_APPLE_DISPATCH)
        int number_of_threads = n / dispatch_stride;
        if (n % dispatch_stride != 0) {
            number_of_threads += 1;
        }
        dispatch_apply(number_of_threads, concurrentQueue, ^(size_t idx) {
            int max = (idx + 1) * dispatch_stride;
            if (max > n) {
                max = n;
            }
            for(int ll = idx * dispatch_stride; ll<max; ll++) {
                if (ll != imax) {
                    double dum = a[ll*n + imax];
                    a[ll*n + imax] = 0.0;
                    #pragma clang loop vectorize(enable) interleave(enable)
                    for (int l=0; l<n; l++) {
                        a[ll*n + l] -= a[imax*n + l] * dum;
                        b[ll*n + l] -= b[imax*n + l] * dum;
                    }
                }
            }
        });
        // On non-Apple platforms, use OpenMP to parallelize the loop instead.
        #else
        #pragma omp parallel for
        for(int ll=0; ll<n; ll++) {
            // Don't subtract from the pivot row...
            if (ll != imax) {
                double dum = a[ll*n + imax];
                a[ll*n + imax] = 0.0;
                #pragma clang loop vectorize(enable) interleave(enable)
                for (int l=0; l<n; l++) {
                    a[ll*n + l] -= a[imax*n + l] * dum;
                    b[ll*n + l] -= b[imax*n + l] * dum;
                }
            }
        }
        #endif
    }
        
    return 0;
}
