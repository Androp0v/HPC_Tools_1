#include <stdio.h>
#include <stdlib.h>
#include <math.h>

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
    
    int i, icol, irow, j, k, l, ll;
    double big, dum, pivinv;
    int *indxc = generate_int_vector(n);
    int *indxr = generate_int_vector(n);
    int *ipiv = generate_int_vector(n);
    
    for(j=0; j<n; j++) {
        ipiv[j] = 0;
    }
    
    for(i=0; i<n; i++) {
        big = 0.0;
        for (j=0; j<n; j++) {
            if (ipiv[j] != 1) {
                for(k=0; k<n; k++) {
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
            for(l=0; l<n; l++) {
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
        for(l=0; l<n; l++) {
            a[icol*n + l] *= pivinv;
            b[icol*n + l] *= pivinv;
        }
        for(ll=0; ll<n; ll++) {
            if (ll != icol) {
                dum = a[ll*n + icol];
                a[ll*n + icol] = 0.0;
                for (l=0; l<n; l++) {
                    a[ll*n + l] -= a[icol*n + l] * dum;
                    b[ll*n + l] -= b[icol*n + l] * dum;
                }
            }
        }
    }
    
    for(l=n-1; l>=0; l--) {
        if (indxr[l] != indxc[l]) {
            for(k=0; k<n; k++) {
                swap(&a[k*n + indxr[l]], &a[k*n + indxc[l]]);
            }
        }
    }
    
    return 0;
}
