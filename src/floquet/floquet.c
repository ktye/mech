/* $Id: floquet.c,v 1.10 2005/09/02 14:30:01 elmar Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <g2c.h>
#include "define.h"
#include "reinit.h"
#include "solve.h"
#include "la.h"
static char id[]="$Id: floquet.c,v 1.10 2005/09/02 14:30:01 elmar Exp $";

extern int NE;
extern int Gflag;
extern int gflag;
extern double	 ALPHA;

/* Y <- (A + c*B), c: real scalar, A,B,Y: complex matrix */
void matplusscal(doublecomplex *A, doublecomplex *B, double c, doublecomplex *Y) {
	int i,j;
	for (i=0; i<NE; i++) {
		for (j=0; j<NE; j++) {
			Y[j*NE+i].r = A[j*NE+i].r + c*B[j*NE+i].r;
			Y[j*NE+i].i = A[j*NE+i].i + c*B[j*NE+i].i;
		}
	}
}

void eye(doublecomplex *A) {
	int i,j;
	for (i=0; i<NE; i++) {
		for (j=0; j<NE; j++) {
			if (i==j) 	A[j*NE+i].r = 1.0;
			else		A[j*NE+i].r = 0.0;
			A[j*NE+i].i = 0.0;
		}
	}
}

void matscl(doublecomplex *A, double h) {
	int i,j;
	for (i=0; i<NE; i++) {
		for (j=0; j<NE; j++) {
			A[j*NE+i].r *= h;
			A[j*NE+i].i *= h;
		}
	}
}

/* runge kutta 4th order time-dependent matrix version, construct integral base */
/* dX/dt = A(t)*X,  X(0)=I, t:[0,2]
 * loop(t): 
 * K1 = h*A(t-h)  *   X(t-h)
 * K2 = h*A(t-h/2)* ( X(t-h) + 1/2*K1 )
 * K3 = h*A(t-h/2)* ( X(t-h) + 1/2*K2 )
 * K4 = h*A(t)    * ( X(t-h) +     K3 )
 * X(t) = X(t-h) + 1/6*K1 + 1/3*K2 + 1/3*K3 + 1/6*K4
 */
void floquet_rk(doublecomplex *X, diterator *Rey, iiterator *alphanum, diterator *t, diterator *k, iiterator *n, iiterator *N) {
	int i,j,T;
	double h, hscale;
	doublecomplex AA[ESIZE*ESIZE];
	doublecomplex A_[ESIZE*ESIZE];
	doublecomplex AB[ESIZE*ESIZE];
	doublecomplex *A0;
	doublecomplex *A1;
	doublecomplex K1[ESIZE*ESIZE];
	doublecomplex K2[ESIZE*ESIZE];
	doublecomplex K3[ESIZE*ESIZE];
	doublecomplex K4[ESIZE*ESIZE];
	doublecomplex X1[ESIZE*ESIZE];
	doublecomplex X2[ESIZE*ESIZE];
	doublecomplex X3[ESIZE*ESIZE];
	
	eye(X);
	reinit(Rey->val,alphanum->val,k->val,n->val,0.0,N->val);
	setmatrix(AA);
	A0 = AA;
	h = (t->max-t->min)/(double)(t->n-1);
	hscale = h * M_PI * Rey->val/ALPHA/ALPHA; /* Transform time variable */
	if (t->max - t->min != 2.0) { fprintf(stderr,"floquet_rk requires -t over period 2\n"); exit(1); }
	for (T=0; T<t->n; T++) {
		SETCNT(t,T);
		reinit(Rey->val,alphanum->val,k->val,n->val,t->val-h/2.0,N->val);
		setmatrix(A_);	
		if (T%2) {
			reinit(Rey->val,alphanum->val,k->val,n->val,t->val,N->val);
			setmatrix(AA);	
			A1 = AA;	
		} else {
			reinit(Rey->val,alphanum->val,k->val,n->val,t->val,N->val);
			setmatrix(AB);	
			A1 = AB;
		}
		matmul(A0,X,K1,NE); matscl(K1,hscale);
		matplusscal(X,K1,0.5,X1);
		matmul(A_,X1,K2,NE); matscl(K2,hscale);
		matplusscal(X,K2,0.5,X2);
		matmul(A_,X2,K3,NE); matscl(K3,hscale);
		matplusscal(X,K3,1.0,X3);
		matmul(A1,X3,K4,NE); matscl(K4,hscale);
		for (i=0; i<NE; i++) {
			for (j=0; j<NE; j++) {
				X[j*NE+i].r += K1[j*NE+i].r/6.0 + K2[j*NE+i].r/3.0 + K3[j*NE+i].r/3.0 + K4[j*NE+i].r/6.0;
				X[j*NE+i].i += K1[j*NE+i].i/6.0 + K2[j*NE+i].i/3.0 + K3[j*NE+i].i/3.0 + K4[j*NE+i].i/6.0;
			}
		}
	}
}

int cmp_gamma(const void *ap, const void *bp) {
	doublecomplex a,b;
	double x,y;
        a = *((doublecomplex *)ap);
        b = *((doublecomplex *)bp);
	x = hypot(a.r,a.i);	
	y = hypot(b.r,b.i);	
        if (x < y) return 1;
	else if (x > y) return -1;
	else return 0;
}

void floquet(diterator *Rey, iiterator *alphanum, diterator *t, diterator *k, iiterator *n, iiterator *N) {
	doublecomplex X[ESIZE*ESIZE];
	doublecomplex gamma[ESIZE];
	doublecomplex gamma_eigvec[ESIZE*ESIZE];

	N->val = N->max;
	floquet_rk(X, Rey, alphanum, t, k, n, N);
	if (Gflag) printmatrix(X,NE);
	eig(X, gamma, gamma_eigvec, NE);
	qsort(gamma, (size_t)NE, sizeof(doublecomplex), cmp_gamma);
	if (gflag) printvector(gamma,NE);
	printf("%.16f %.16f %.16f\n",gamma[0].r,gamma[0].i,hypot(gamma[0].r,gamma[0].i));
}
