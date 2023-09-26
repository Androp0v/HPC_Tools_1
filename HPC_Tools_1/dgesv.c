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
    
    int *indxc = generate_int_vector(n);
    int *indxr = generate_int_vector(n);
    int *ipiv = generate_int_vector(n);
    
    double big, dum, pivinv;
    int icol, irow;
        
    #pragma clang loop vectorize(enable) interleave(enable)
    for(int j=0; j<n; j++) {
        ipiv[j] = 0;
    }
    
    #if defined(USE_APPLE_DISPATCH)
    dispatch_queue_t concurrentQueue = dispatch_get_global_queue(DISPATCH_QUEUE_PRIORITY_HIGH, 0);
    #endif
    
    for(int i=0; i<n; i++) {
        big = 0.0;
        for (int j=0; j<n; j++) {
            if (ipiv[j] != 1) {
                for(int k=0; k<n; k++) {
                    if (ipiv[k] == 0) {
                        if (fabs(a[j*n + k]) >= big) {
                            big = fabs(a[j*n + k]);
                            irow = j;
                            icol = k;
                        }
                    }
                }
            }
        }
        ++(ipiv[icol]);
        if (irow != icol) {
            for(int l=0; l<n; l++) {
                swap(&a[irow*n + l], &a[icol*n + l]);
                swap(&b[irow*n + l], &b[icol*n + l]);
            }
        }
        indxr[i] = irow;
        indxc[i] = icol;
        if (a[icol*n + icol] == 0.0) {
            printf("Singular matrix");
            exit(0);
        }
        pivinv = 1.0 / a[icol*n + icol];
        a[icol*n + icol] = 1.0;
        
        #pragma clang loop vectorize(enable) interleave(enable)
        for(int l=0; l<n; l++) {
            a[icol*n + l] *= pivinv;
            b[icol*n + l] *= pivinv;
        }
        
        // On Apple platforms, execute concurrently in all available cores
        // using Apple's dispatch library.
        #if defined(USE_APPLE_DISPATCH)
        dispatch_apply(n, concurrentQueue, ^(size_t ll) {
        #else
        for(int ll=0; ll<n; ll++) {
        #endif
            if (ll != icol) {
                double dum = a[ll*n + icol];
                a[ll*n + icol] = 0.0;
                #pragma clang loop vectorize(enable) interleave(enable)
                for (int l=0; l<n; l++) {
                    a[ll*n + l] -= a[icol*n + l] * dum;
                    b[ll*n + l] -= b[icol*n + l] * dum;
                }
            }
        #if defined(USE_APPLE_DISPATCH)
        });
        #else
        }
        #endif
    }

    for(int l=n-1; l>=0; l--) {
        if (indxr[l] != indxc[l]) {
            for(int k=0; k<n; k++) {
                swap(&a[k*n + indxr[l]], &a[k*n + indxc[l]]);
            }
        }
    }
        
    return 0;
}
