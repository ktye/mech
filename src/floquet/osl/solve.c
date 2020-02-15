/* $Id: solve.c,v 1.7 2005/05/23 16:44:07 elmar Exp $ */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <g2c.h>

#include "define.h"
#include "base.h"
#include "baseflow.h"
#include "la.h"
#include "quad.h"

extern double	 R;
extern double	 T;
extern int	 N;
extern double	 K;
extern int	 vflag;
extern int	 mflag;

extern double	 BASEFLOW[YSIZE];
extern double	 BASEFLOW2[YSIZE];
extern double	 PHI[YSIZE][NMAX];
extern double	 PSI[YSIZE][NMAX];
extern double	 DPSI[YSIZE][NMAX];
extern double	 LAMBDA[NMAX];

void setmatrix(doublecomplex *A) {
	int	 i,j,Y;
	double	 a[YSIZE];

	for (j=0; j<N; j++) {
		for (i=0; i<N; i++) {
			for (Y=0; Y<YSIZE; Y++) {
				a[Y] = BASEFLOW[Y]*PHI[Y][i] - BASEFLOW2[Y]*PSI[Y][i];
				a[Y] *= PSI[Y][j];
			}
			if (i==j) {
				A[j*N+i].r = LAMBDA[i];
			} else {
				A[j*N+i].r = 0.0;
			}
			A[j*N+i].i = -K*quad(a);
		}
	}
	if (mflag)
		printmatrix(A,N);
}

void eig_fn(doublecomplex *eval, doublecomplex *evec) {
	int	 i, j, I;
	double	 sigma;
	double	 pr, pi, ur, ui, vr, vi, wr, wi;
	doublecomplex	*Ivec;

	I = 0;
	sigma = eval[0].r;
	for (i=1; i<N; i++) {
		if (eval[i].r > sigma) {
			sigma = eval[i].r;
			I = i;
		}
	}
	Ivec = (doublecomplex *)&(evec[I*N].r);
	for (i=0; i<YSIZE; i++) {
		pr = 0.0;
		pi = 0.0;
		ur = 0.0;
		ui = 0.0;
		vr = 0.0;
		vi = 0.0;
		wr = 0.0;
		wi = 0.0;
		for (j=0; j<N; j++) {
			pr += Ivec[j].r * PSI[i][j];
			pi += Ivec[j].i * PSI[i][j];
			ur += Ivec[j].r * DPSI[i][j];
			ui += Ivec[j].i * DPSI[i][j];
			vr -= Ivec[j].i * K*PSI[i][j];
			vi += Ivec[j].r * K*PSI[i][j];
			wr += Ivec[j].r * PHI[i][j];
			wi += Ivec[j].i * PHI[i][j];
		}
		printf("v: %g %g %g %g %g %g\n",SETY(i),hypot(hypot(ur,ui),hypot(vr,vi)),hypot(wr,wi),pr,pi,hypot(pr,pi));
	}

}

int cmp_eig(const void *ap, const void *bp) {
        doublecomplex a,b;
        a = *((doublecomplex *)ap);
        b = *((doublecomplex *)bp);
        if ((double)a.r < (double)b.r) return 1;
        else if ((double)a.r == (double)b.r) {
                if ((double)a.i < (double)b.i) return 1;
                else if ((double)a.i == (double)b.i) return 0;
                else return -1;
        } else return -1;
}

void solve(double *sigma, double *omega) {
	doublecomplex	 A[NMAX*NMAX];
	doublecomplex	 eigenvalue[NMAX];
	doublecomplex	 eigenvector[NMAX*NMAX];

	setmatrix(A);
//	printmatrix(A,N);
	eig(A, eigenvalue, eigenvector, N);
	if (vflag)
		eig_fn(eigenvalue, eigenvector);
	qsort(eigenvalue, (size_t)N, sizeof(doublecomplex), cmp_eig);
//	printvector(eigenvalue,N);
	*sigma = eigenvalue[0].r;
	*omega = eigenvalue[0].i;
}
