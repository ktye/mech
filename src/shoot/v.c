/* $Id: v.c,v 1.10 2004/09/10 09:37:56 elmar Exp $ */

/*
 * v - differential equation n=0, phi-mode
 */

#include <stdlib.h>
#include <math.h>
#include "baseflow.h"
#include "rk4.h"
extern int flag_print;
extern int debug;

extern double OMEGAR;
extern double OMEGAI;
extern double K;
extern double R;
extern void (*ode)();

/*
 * ode: V'' + V'/r + (a-1/r^2)V = 0
 * with a=iR(omega-k*baseflow)-k^2
 * shifted: V(r) = r*v(s), s=r^2
 * --> v'' + 2v'/s + av/4s = 0
 *
 * y[0] + i*y[1]: v
 * y[2] + i*y[3]: v'
 *
 * initial values: 
 * 	v(0) = 1, v'(0) = -a[0]/8
 * regularity:
 * 	v''(0) = a[0]^2/96 - a[1]/12
 */
double v(double *y, int i, int ny, int step) {
	double a0r, a0i, a1i, ar, ai, W, dW, W0, W1;
	double s;

	W0 = BASEFLOW[0];
	W1 = BASEFLOW1[0];
	a0r = -K*K/(R*R) -OMEGAI;
	a0i = OMEGAR - K*W0;
	a1i = -K*W1;
	
	W = BASEFLOW[step];
	dW= BASEFLOW1[step];
	ar = a0r;
	ai = OMEGAR - K*W;

	s = SSTEP(step);	
	if (0==step) {
		switch (i) {
		case 0: return -a0r/8.0;
		case 1: return -a0i/8.0;
		case 2: return (a0r*a0r-a0i*a0i)/96.0;
		case 3: return a0r*a0i/48.0-a1i/12.0;
		default: exit(1);
		}
	} else {
		switch (i) {
		case 0: return y[2];
		case 1: return y[3];
		case 2: return (-2.0*y[2] + (ai*y[1]-ar*y[0])/4.0)/s;
		case 3: return (-2.0*y[3] - (ai*y[0]+ar*y[1])/4.0)/s;
		default: exit(1);
		}
	}
}

double shootv(void) {
	int i;
	double y[4], ymax[4];
	double a0r, a0i, W0;


	W0 = BASEFLOW[0];
	a0r = -K*K/(R*R) -OMEGAI;
	a0i = OMEGAR - K*W0;

	flag_print = 0; if (debug==1)	flag_print = 1;
	y[0] = 1.0;
	y[1] = 0.0;
	y[2] = -a0r/8.0;
	y[3] = -a0i/8.0;
	ode(y, 4, v, ymax);
	return hypot(y[0],y[1]);
}
