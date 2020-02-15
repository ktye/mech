/* $Id: solve.c,v 1.13 2005/07/28 10:03:51 elmar Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <g2c.h>
#include "define.h"
#include "quad.h"
#include "la.h"

static char id[]="$Id: solve.c,v 1.13 2005/07/28 10:03:51 elmar Exp $";



extern double	R;
extern double	K;
extern int 	Nphi;
extern int	NE;
extern int 	Mflag;
extern int	Lflag;
extern int	Sflag; 
extern int	Vflag; 
extern int	vflag; 

extern struct expvelfn {
        double u[RSIZE][ESIZE];
        double v[RSIZE][ESIZE];
        double w[RSIZE][ESIZE];
        double ev[ESIZE];
        int rsize;
        int esize;
} EXPVEL;

extern double SIGMA[ESIZE];
extern double OMEGA[ESIZE];

struct baseflowstruct {
        double     r[RSIZE];
        double    dr[RSIZE];
        double  real[RSIZE];
        double dreal[RSIZE];
        double  imag[RSIZE];
        double dimag[RSIZE];
} baseflow;


/*
 * A = lambda_i*delta_ij + <v_i, _convective terms_ >
 * 	lambda_i: eigenvalues of the expansion functions lambda = -1*(beta_i)/R:
 * 	v_i: 	ith expansion velocity
 * 	_convective terms_: [%iku*baseflow, %ikv*baseflow, %ikw*baseflow, + u*dbaseflow]
 * 		dbaseflow: derivative dbaseflow/dr
 * 		u: sum_k a_k*u_k, 	u_k: expansion velocity r-component
 * 		v: sum_k a_k*v_k, 	v_k: expansion velocity phi-component
 * 		w: sum_k a_k*w_k, 	w_k: expansion velocity z-component
 *
 * the inner product: <a, b> is defined as the integral of the scalar product 
 * of the conjugate complex of vector a with vector b, integrated over the domain r*dr:
 * 	<a,b> = integral cc(a) _dot_ b r dr; 	with r: [0,1], cc: conjugate complex
 *
 * since the expansion velocities are:
 * 	u = -%i^n*U
 * 	v = -%i^(n+1)*V
 * 	w = -%i^(n+1)*W
 *
 * the <v_i, _convective terms_ > leads to the matrix vector pair A*a with
 * A(i,j) += 
 * 	%i*k*integral( (U_i*U_j + V_i*V_j + W_i*W_j)*baseflow )*r*dr;
 * 	+ %i*integral( W_i*U_j*dbaseflow )*r*dr;
 */
void setmatrix(doublecomplex *A) {
	int i,j,r;
	double x,X[RSIZE];
	for (j=0; j<NE; j++) {
		for (i=0; i<NE; i++) {
			for (r=0; r<RSIZE; r++) {
				X[r] = EXPVEL.u[r][i]*EXPVEL.u[r][j]*baseflow.r[r];
				X[r]+= EXPVEL.v[r][i]*EXPVEL.v[r][j]*baseflow.r[r];
				X[r]+= EXPVEL.w[r][i]*EXPVEL.w[r][j]*baseflow.r[r];
				X[r]*= K;
				X[r]+= EXPVEL.w[r][i]*EXPVEL.u[r][j]*baseflow.dr[r];
				X[r]*= (double)r/(double)(RSIZE-1);
			}
			x = quad(X);
			A[j*NE+i].r = 0.0;
			A[j*NE+i].i = x;
			if (i==j) {
				A[j*NE+i].r -= (K*K + EXPVEL.ev[i]*EXPVEL.ev[i])/R;
			}
		}
		if (Lflag) {
			printf("lambda[%d]=%.16f\n",j,-(K*K + EXPVEL.ev[j]*EXPVEL.ev[j])/R);
		}
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

#define ISBAD(x) (isinf(x)||(isnan(x)))
int checkmatrix(doublecomplex *A, int N) {
	int i,j;
	for (i=0; i<N; i++) {
		for (j=0; j<N; j++) {
			if (ISBAD(A[j*N+i].r)) return 1;
			if (ISBAD(A[j*N+i].i)) return 1;
		}
	}
	return 0;
}


/*
 * print most unstable eigenvector
 */
void eig_vec(doublecomplex *eval, doublecomplex *evec) {
	int	 i, I;
	double	 sigma;
	doublecomplex *Ivec;

	I = 0;
	sigma = eval[0].r;
	for (i=1; i<NE; i++) {
		if (eval[i].r > sigma) {
			sigma = eval[i].r;
			I = i;		
		}
	}
	Ivec = (doublecomplex *)&(evec[I*NE].r);
	printf("v:");
	for (i=0; i<NE; i++) {
		printf(" %g",Ivec[i].r);
		printf(" %g",Ivec[i].i);
	}
	printf("\n");
}

/* 
 * print most unstable velocity eigenfunctions:
 * -> take unsorted eigenvalues and eigenvectors
 * -> look for the most unstable eigenvalue
 * -> take the corresponding eigenvector
 * -> multiply it with expansion velocities, respect complex notation
 *  	( u = -i^n U(r), v,w = i^{n+1}(V(r),W(r)) )
 */
void eig_fn(doublecomplex *eval, doublecomplex *evec) {
	int i, I, rc;
	double ur, vr, wr, ui, vi, wi;
	double sigma;
	double rf[RSIZE], f[RSIZE];
	doublecomplex *Ivec;

	I = 0;
	sigma = eval[0].r;
	for (i=1; i<NE; i++) {
		if (eval[i].r > sigma) {
			sigma = eval[i].r;
			I = i;		
		}
	}
	//printf("eig_fn I=%d sigma=%g\n",I,sigma);	
	//// Ith eigenvector
	//for (i=I*NE; i<I*NE+NE; i++) {
	//	printf("\t %g %g\n",evec[i].r,evec[i].i);
	//}
	Ivec = (doublecomplex *)&(evec[I*NE].r);
	for (rc=0; rc<RSIZE; rc++) {
		ur = 0.0; vr = 0.0; wr = 0.0;
		ui = 0.0; vi = 0.0; wi = 0.0;
		for (i=0; i<NE; i++) {
			switch (Nphi%4) {
			case 0:
				ur -= Ivec[i].r*EXPVEL.u[rc][i];
				ui -= Ivec[i].i*EXPVEL.u[rc][i];
				vr += Ivec[i].i*EXPVEL.v[rc][i];
				vi -= Ivec[i].r*EXPVEL.v[rc][i];
				wr += Ivec[i].i*EXPVEL.w[rc][i];
				wi -= Ivec[i].r*EXPVEL.w[rc][i];
				break;
			case 1:
				ur += Ivec[i].i*EXPVEL.u[rc][i];
				ui -= Ivec[i].r*EXPVEL.u[rc][i];
				vr += Ivec[i].r*EXPVEL.v[rc][i];
				vi += Ivec[i].i*EXPVEL.v[rc][i];
				wr += Ivec[i].r*EXPVEL.w[rc][i];
				wi += Ivec[i].i*EXPVEL.w[rc][i];
				break;
			case 2:
				ur += Ivec[i].r*EXPVEL.u[rc][i];
				ui += Ivec[i].i*EXPVEL.u[rc][i];
				vr -= Ivec[i].i*EXPVEL.v[rc][i];
				vi += Ivec[i].r*EXPVEL.v[rc][i];
				wr -= Ivec[i].i*EXPVEL.w[rc][i];
				wi += Ivec[i].r*EXPVEL.w[rc][i];
				break;
			case 3:
				ur -= Ivec[i].i*EXPVEL.u[rc][i];
				ui += Ivec[i].r*EXPVEL.u[rc][i];
				vr -= Ivec[i].r*EXPVEL.v[rc][i];
				vi -= Ivec[i].i*EXPVEL.v[rc][i];
				wr -= Ivec[i].r*EXPVEL.w[rc][i];
				wi -= Ivec[i].i*EXPVEL.w[rc][i];
				break;
			default: exit(1);
			}
		}
		printf("V: %g %g %g %g %g %g %g\n",(double)rc/(double)(RSIZE-1),ur,ui,vr,vi,wr,wi);
		f[rc] = ur*ur+ui*ui+vr*vr+vi*vi+wr*wr+wi*wi;
		rf[rc] = f[rc]*(double)rc/(double)(RSIZE-1);
	}
	printf("r: %.16f\n",quad(rf)/quad(f));
}

int solve(void) {
	int i;
	doublecomplex A[ESIZE*ESIZE];
	doublecomplex eigenvalue[ESIZE];
	doublecomplex eigenvector[ESIZE*ESIZE];

//	fprintf(stderr,"solve.");
//	fprintf(stderr,"setmatrix.");
	setmatrix(A);
//	fprintf(stderr,".setmatrix ");
	if (checkmatrix(A,NE)) {
		fprintf(stderr,"solve: ERROR\n");
		return 1;
	}
	if (Mflag) printmatrix(A,NE);
//	fprintf(stderr,"eig.");
	eig(A,eigenvalue,eigenvector,NE);
//	fprintf(stderr,".eig ");
	if (Vflag) eig_fn(eigenvalue, eigenvector);
	if (vflag) eig_vec(eigenvalue, eigenvector);
	qsort(eigenvalue, (size_t)NE, sizeof(doublecomplex), cmp_eig);
	for (i=0; i<NE; i++) {
		SIGMA[i] = (double)(eigenvalue[i].r);
		OMEGA[i] = (double)(eigenvalue[i].i);
	}
	if (Sflag) printvector(eigenvalue,NE);
//	fprintf(stderr,".solve\n");


	return 0;
}
