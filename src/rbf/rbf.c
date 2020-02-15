/* $Id: rbf.c,v 1.2 2006/02/18 18:18:05 elmar Exp $ */
/*
 * rbf - radial basis functions
 *
 * given:
 * 	N function values f_k at points x_k, x \in R^n
 *
 * interpolate value at point x:
 * 	f(x) = sum_{k=1}^N c_k phi(|x-x_k|)
 *
 * 	phi(d): kernel function (basis functions)
 * 	c_k: 	unknown coefficients
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "symsolve.h"

int	 DIM;
int	 NDIM;
double	 L;
double	*COEF;
double	*POINTS;


double kernel_exp(double xx, double ll, int derivative) {
	if (derivative)
		exit(1);
	else
		return exp(-xx/ll);
}
double kernel_poly2l(double xx, double ll, int derivative) {
	if (derivative)
		exit(1);
	else
		if (xx<4.0*ll)
			return 1.0-xx/ll/2.0+xx*xx/ll/ll/16.0;
		else
			return 0.0;
}
double (*kernel)() = kernel_exp;
/*
 * x:	array of input points (dimension (dim+1)*N)
 * dim:	space dimension
 * N:	Number of input points
 * l:	characteristic lengthscale
 */
void rbf_init(double *x, int dim, int N, double l, int kern) {
	int	 i, k, n, m, elements, Dim;
	double	 d;
	double 	*A;

	if (kern)
		kernel = kernel_exp;
	else
		kernel = kernel_poly2l;
	elements = (N*(N+1))/2;
	A = (double *)malloc(elements*sizeof(double));
	COEF = (double *)malloc(N*sizeof(double));

	Dim = dim + 1;
	DIM = dim;
	NDIM = N;
	POINTS = x;
	L = l*l;
	n = 0;
	for (i=0; i<N; i++) {
		for (k=0; k<=i; k++) {
			d = 0.0;
			for (m=0; m<dim; m++) {
				d += (x[Dim*i+m] - x[Dim*k+m])
				    *(x[Dim*i+m] - x[Dim*k+m]);
				fprintf(stderr,"i=%d k=%d m=%d d=%g\n",i,k,m,d);
			}
			A[n++] = kernel(d,L,0);	//exp(-d/L);
		}
	}
	for (i=0; i<N; i++) {
		COEF[i] = x[Dim*i+dim];
		fprintf(stderr,"RHS[%d] = %g\n", i, COEF[i]);
	}
	for (i=0; i<elements; i++) {
		fprintf(stderr,"A[%d] = %g\n", i, A[i]);
	}
	symsolve(A, COEF, N);
	for (i=0; i<N; i++) {
		fprintf(stderr,"COEF[%d] = %g\n", i, COEF[i]);
	}
}

double rbf(double *X) {
	int	 i, m;
	double	 f, d;

//	fprintf(stderr,"L = %g\n", L);
	for (i=0; i<NDIM; i++) {
//		fprintf(stderr,"COEF[%d] = %g\n", i, COEF[i]);
	}
	
	f = 0.0;
	for (i=0; i<NDIM; i++) {
		d = 0.0;
		for (m=0; m<DIM; m++)
			d += (POINTS[(DIM+1)*i+m] - X[m])
			    *(POINTS[(DIM+1)*i+m] - X[m]);
		f += COEF[i]*kernel(d,L,0);	//exp(-d/L);
	}
	return f;
}

double sph_approx(double *x, double *X, int dim, int N, double l) {
	int	 i,k;
	double	 f, d;	
	
	f = 0.0;
	for (i=0; i<N; i++) {
		d = 0.0;
		for (k=0; k<dim; k++) {
			d += (x[(dim+1)*i+k] - X[k])
			    *(x[(dim+1)*i+k] - X[k]);
		}
		f += x[(dim+1)*i+1] * exp(-d/l/l)/sqrt(M_PI);
	}
	return f;
}
