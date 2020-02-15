/* $Id: la.c,v 1.2 2005/04/25 13:15:54 elmar Exp $ */


#include <stdio.h>
#include <stdlib.h>
#include <g2c.h>
#include <clapack.h>
#include "define.h"

/* A: column-major order
 * i: row index
 * j: col index
 */
void printmatrix(doublecomplex *A, int n) {
        int i,j;
        for (i=0; i<n; i++) {
		printf("m: ");
                for (j=0; j<n; j++) {
                        printf("[%.5f %.5f]",A[j*n+i].r,A[j*n+i].i);
                }
                printf("\n");
        }
}

void printvector(doublecomplex *v, int n) {
        int i;
        for (i=0; i<n; i++) {
                printf("%.16f %.16f\n",v[i].r,v[i].i);
        }
}


/*
void matmul(doublecomplex *A, doublecomplex *B, doublecomplex *C, int N) {
        int i,j,k;
        for (i=0; i<N; i++) {
                for (j=0; j<N; j++) {
                        C[j*N+i].r = 0.0;
                        C[j*N+i].i = 0.0;
                        for (k=0; k<N; k++) {
                                C[j*N+i].r += A[k*N+i].r*B[j*N+k].r - A[k*N+i].i*B[j*N+k].i;
                                C[j*N+i].i += A[k*N+i].i*B[j*N+k].r + A[k*N+i].r*B[j*N+k].i;
                        }
                }
        }
}
*/

void eig(doublecomplex *A, doublecomplex *eigenvalue, doublecomplex *eigenvector, int n) {
        char    jobvl = 'N';    /* don't compute left eigenvectors */
        char    jobvr = 'V';    /* compute right eigenvectors */
        integer one = 1;
        integer lwork = 16*NMAX;
        integer info = 0;
        integer N;
        doublecomplex work[16*NMAX];
        doublereal rwork[2*NMAX];

        N = (integer)n;

        zgeev_(&jobvl, &jobvr, &N, A, &N, eigenvalue, (doublecomplex *)NULL, &one, eigenvector, &N, work, &lwork, rwork, &info);
        if (info) {
                fprintf(stderr,"error: zggev_ info: %d\n",(int)info);
                exit(1);
        }
}
