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

int my_dgesv(int n, double a[n][n], double b[n][n]) {
    
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
                        if (fabs(a[j][k]) >= big) {
                            big = fabs(a[j][k]);
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
                swap(&a[irow][l], &a[icol][l]);
                swap(&b[irow][l], &b[icol][l]);
            }
        }
        indxr[i] = irow;
        indxc[i] = icol;
        if (a[icol][icol] == 0.0) {
            printf("Singular matrix");
            exit(0);
        }
        pivinv = 1.0 / a[icol][icol];
        a[icol][icol] = 1.0;
        for(l=0; l<n; l++) {
            a[icol][l] *= pivinv;
            b[icol][l] *= pivinv;
        }
        for(ll=0; ll<n; ll++) {
            if (ll != icol) {
                dum = a[ll][icol];
                a[ll][icol] = 0.0;
                for (l=0; l<n; l++) {
                    a[ll][l] -= a[icol][l] * dum;
                    b[ll][l] -= b[icol][l] * dum;
                }
            }
        }
    }
    
    for(l=n-1; l>=0; l--) {
        if (indxr[l] != indxc[l]) {
            for(k=0; k<n; k++) {
                swap(&a[k][indxr[l]], &a[k][indxc[l]]);
            }
        }
    }
    
    return 0;
}
