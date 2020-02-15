/* $Id: la.c,v 1.4 2005/09/05 14:22:10 elmar Exp $ */


#include <stdio.h>
#include <stdlib.h>
#include <g2c.h>
#include <clapack.h>
//#include <time.h>
#include "define.h"

/* A: column-major order
 * i: row index
 * j: col index
 */
void printmatrix(doublecomplex *A, int n) {
        int i,j;
        for (i=0; i<n; i++) {
                for (j=0; j<n; j++) {
                        printf("[%.5f %.5f]",A[j*n+i].r,A[j*n+i].i);
                }
                printf("\n");
        }
}
void fprintmatrix(doublecomplex *A, int n) {
        int i,j;
        for (i=0; i<n; i++) {
                for (j=0; j<n; j++) {
                        fprintf(stderr,"[%.5f %.5f]",A[j*n+i].r,A[j*n+i].i);
                }
                fprintf(stderr,"\n");
        }
}


void printvector(doublecomplex *v, int n) {
        int i;
        for (i=0; i<n; i++) {
                printf("%.16f %.16f\n",v[i].r,v[i].i);
        }
}


/* c = A*b 	c,b: complex[N], A: complex[NxN] lapack-style col-major
 * check before use!
//void matvec(doublecomplex *A, doublecomplex *b, doublecomplex *c, int N) {
//	int i,j;
//	for (i=0; i<N; i++) {
//		c[i].r = 0.0;
//		c[i].i = 0.0;
//		for (j=0; j<N; j++) {
//			c[i].r += A[i*N+j].r*b[i].r - A[i*N+j].i*b[i].i;
//			c[i].i += A[i*N+j].i*b[i].r + A[i*N+j].r*b[i].i;
//		}
//	}
//}
*/

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

/* blas version of matmul need blaswrap.h cblas..
void matmul_blas(doublecomplex *A, doublecomplex *B, doublecomplex *C, int N) {
        char no = 'n';
        doublecomplex alpha, beta;

        alpha.r = 1.0; alpha.i = 0.0;
        beta.r = 0.0; beta.i = 0.0;
        f2c_zgemm(&no,&no,&N,&N,&N,&alpha,A,&N,B,&N,&beta,C,&N);
}
*/

void eig(doublecomplex *A, doublecomplex *eigenvalue, doublecomplex *eigenvector, int n) {
        char    jobvl = 'N';    /* don't compute left eigenvectors */
        char    jobvr = 'V';    /* compute right eigenvectors */
        integer one = 1;
        integer lwork = 16*ESIZE;
        integer info = 0;
        integer N;
        doublecomplex work[16*ESIZE];
        doublereal rwork[2*ESIZE];

        N = (integer)n;

        zgeev_(&jobvl, &jobvr, &N, A, &N, eigenvalue, (doublecomplex *)NULL, &one, eigenvector, &N, work, &lwork, rwork, &info);
        if (info) {
                fprintf(stderr,"error: zggev_ info: %d\n",(int)info);
                exit(1);
        }
}


/*
void la_test(void) {
	int i,j,N;
	doublecomplex A[4*4];
	doublecomplex B[4*4];
	doublecomplex C[4*4];

	srand48((long)time(0));
	N = 4;
	for (i=0; i<N; i++) {
		for (j=0; j<N; j++) {
			A[j*N+i].r = drand48();
			A[j*N+i].i = drand48();
			B[j*N+i].r = drand48();
			B[j*N+i].i = drand48();
		}
	}
	printf("A=[\n");
	printmatrix(A,N);
	printf("]\nB=[\n");
	printmatrix(B,N);
	matmul(A,B,C,N);	
	printf("]\nC=[\n");
	printmatrix(C,N);
	printf("]\n");
}
*/
