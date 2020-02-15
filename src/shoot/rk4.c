/* $Id: rk4.c,v 1.4 2004/09/10 09:37:56 elmar Exp $ */

/*
 * Runge Kutta 4th order, with maxima storage
 */

#include <math.h>
#include "rk4.h"

#define MAXORDER 12
extern int flag_print;
extern unsigned int Nsteps;
extern unsigned int Lsteps;


void ode_printvec(double *y, int n, double x) {
	int i;
	printf("%g ",x);
	for (i=0; i<n; i++)
		printf("%g ",y[i]);
	printf("%g\n");

}


void rk4(double *y, int ny, double (*f)(), double *ymax) {
	double 	t1[MAXORDER], t2[MAXORDER], t3[MAXORDER],
		k1[MAXORDER], k2[MAXORDER], k3[MAXORDER],k4[MAXORDER];
	int 	i,j;
	double 	d;

	for (i=0; i<ny; i++) ymax[i] = 0.0;
	d = (double)(Nsteps-1);
	i = 0;
	if (flag_print)	ode_printvec(y,ny, SSTEP(i) );
	for (i=0; i<Nsteps-1; i++) {
		for (j=0; j<ny; j++) t1[j] = y[j]+0.5*(k1[j]=f(y, j, ny, 2*i)/d);
		for (j=0; j<ny; j++) t2[j] = y[j]+0.5*(k2[j]=f(t1, j, ny, 2*i+1)/d);
		for (j=0; j<ny; j++) t3[j] = y[j]+    (k3[j]=f(t2, j, ny, 2*i+1)/d);
		for (j=0; j<ny; j++) k4[j] =                 f(t3, j, ny, 2*i+2)/d;

		for (j=0; j<ny; j++) y[j] += (k1[j]+2.0*k2[j]+2.0*k3[j]+k4[j])/6.0;
		if (flag_print)	ode_printvec(y,ny,SSTEP(2*i+2));
		for (j=0; j<ny; j++) if (fabs(y[j])>ymax[j]) ymax[j] = fabs(y[j]);
	}
}
