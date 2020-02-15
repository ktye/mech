/* $Id: base.c,v 1.5 2005/04/22 15:58:54 elmar Exp $ */
#include <math.h>
#include <stdio.h>

#include "define.h"
#include "quad.h"
extern double 	 K;
extern int 	 N;
extern double 	 R;
extern int	 S;

extern double	 PHI[YSIZE][NMAX];
extern double	 PSI[YSIZE][NMAX];
extern double	 DPSI[YSIZE][NMAX];
extern double	 LAMBDA[NMAX];

double	 BETA[NMAX];
double	 AMP[NMAX];


/*
 * base functions are solutions of the 
 * eigenvalue problem
 *
 * l phi = laplace phi
 * laplace phi = psi
 * rot psi = (u,v)
 * (u, v) = 0 @ y=1 & y=-1
 *
 * in 2d domain (x,y)
 *
 * ~exp(i k x), k: real
 * y [-1,1]
 *
 */

double disrel(double beta) {
	if (S)
		return beta*tan(beta) + K*tanh(K);
	else
		return K*tan(beta) - beta*tanh(K);
}

void disrel_init(void) {
	int	 n;
	double	 interval[2], x, y, yy;
	for (n=0; n<N; n++) {
		if (S) {
			interval[0] = ((double)n+0.5)*M_PI;
			interval[1] = (double)(n+1)*M_PI;
		} else {
			interval[0] = (double)(n+1)*M_PI;
			interval[1] = ((double)n+1.5)*M_PI;
		}
		yy = interval[0];
		y = interval[1];
		x = 0.0;
		while (y != yy) {
			x = 0.5*(interval[0]+interval[1]);
			yy = y;
			y = disrel(x);
			if (y>0.0)
				interval[1] = x;
			else
				interval[0] = x;
		}
		BETA[n] = x;
	}

	//
	/*
	for (n=0; n<YSIZE; n++) {
		y = (double)n*(((double)N+1.5)*M_PI)/(double)(YSIZE-1); 
		printf("%g %g\n",y,y*tan(y));
	}
	printf("\n:1:2:yellow\n");
	for (n=0; n<YSIZE; n++) {
		y = (double)n*(((double)N+1.5)*M_PI)/(double)(YSIZE-1); 
		printf("%g %g\n",y,-K*tanh(K));
	}
	for (n=0; n<YSIZE; n++) {
		y = (double)n*(((double)N+1.5)*M_PI)/(double)(YSIZE-1); 
		printf("%g %g\n",y,disrel(y));
	}
	printf("\n:1:2:red\n");
	for (n=0; n<N; n++) {
		printf("%g %g\n%g %g\n\n",BETA[n],-1.0,BETA[n],1.0);
		fprintf(stderr,"BETA[%d] = %g\n",n,BETA[n]);
	}
	exit(0);
	*/
	//

}

double phifn(double y, int n) {
	double	 b, l;
	b = BETA[n];
	l = - b*b - K*K;
	if (S)
		return l*cosh(K)/cos(b)*cos(b*y);
	else
		return l*sinh(K)/sin(b)*sin(b*y);
}

double psifn(double y, int n) {
	double	 b;
	b = BETA[n];
	if (S)
		return cosh(K)*cos(b*y)/cos(b)-cosh(K*y);
	else
		return sin(b*y)*sinh(K)/sin(b)-sinh(K*y);
}

double dpsifn(double y, int n) {
	double	 b, l;
	b = BETA[n];
	l = - b*b - K*K;
	if (S)
		return -b*cosh(K)*sin(b*y)/cos(b)-K*sinh(K*y);
	else
		return b*cos(b*y)*sinh(K)/sin(b)-K*cosh(K*y);
}

double phi(double y, int n) {
	return phifn(y,n)/AMP[n];
}

double psi(double y, int n) {
	return psifn(y,n)/AMP[n];
}

double dpsi(double y, int n) {
	return dpsifn(y,n)/AMP[n];
}

double lambda(int n) {
	return -(BETA[n]*BETA[n]+K*K)/R;
}

void base_init(void) {
	double	 a[YSIZE], y;
	int	 n, i;
	disrel_init();
	for (n=0; n<N; n++) {
		for (i=0; i<YSIZE; i++) {
			a[i] = phifn((double)i/(double)(YSIZE-1),n);
			a[i] *= -psifn((double)i/(double)(YSIZE-1),n);
		}
		AMP[n] = sqrt(quad(a));
	}

	for (n=0; n<N; n++) {
		for (i=0; i<YSIZE; i++) {
			y = SETY(i);
			PHI[i][n] = phi(y,n);
			PSI[i][n] = psi(y,n);
			DPSI[i][n] = dpsi(y,n);
		}
		LAMBDA[n] = lambda(n);
	}
	
	/* orthogonality <psi_n, psi_m> = - integral( phi_n, psi_m ) = I;
	for (i=0; i<N; i++) {
		for (j=0; j<N; j++) {
			for (n=0; n<YSIZE; n++) {
				a[n] = -PHI[n][i]*PSI[n][j];
			}
			printf("%g ",quad(a));
		}
		printf("\n");
	}
	*/
}


