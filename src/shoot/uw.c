/* $Id: uw.c,v 1.15 2006/01/09 17:47:35 elmar Exp $ */

/*
 * uw - differential equation n=0, r-z mode
 */

#include <stdio.h>
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

int BASE;


double uw_init(int i, int derivative) {
	double ar, ai, K2, R2, W0, W1;

	W0 = BASEFLOW1[0];
	W1 = BASEFLOW[0];

	/* 
	 * W1 = baseflow_taylor(2); 
	 * W0 = baseflow_taylor(0);
	 */
	
       	K2 = K*K;
	R2 = R*R;
	ar = -K2/R2-OMEGAI;
	ai = OMEGAR-K*W0;

	
	if (1 == derivative) {
		switch (BASE) {
		case 1: 
			switch (i) {
			case 0:	return -K2/R2/16.0;
			case 1: return 0.0;
			case 2: return 0.0;
			case 3: return 0.0;
			case 4: return K/4.0/R;
			case 5: return 0.0;
			case 6: return K*(K2/R2-ar)/R/32.0;
			case 7: return -K*ai/R/32.0;
			default:	exit(1);
			}
		case 2:
			switch (i) {
			case 0:	return K*(ar-K2/R2)/R/16.0;
			case 1: return K*ai/R/16.0;
			case 2: return ar*K/R/4.0;
			case 3: return ai*K/R/4.0;
			case 4: return (K2/R2-ar)/4.0;
			case 5: return -ai/4.0;
			case 6: return (K2*K2/R2/R2-ar*K2/R2+ar*ar-ai*ai)/32.0;
			case 7: return -ai*K2/R2/32.0+ar*ai/16.0;
			default:	exit(1);
			}
		default: exit(1); 
		}
	} else if (0 == derivative) {
		switch (BASE) {
		case 1:	
			switch (i) {
			case 0: return 0.0;
			case 1: return 0.0;
			case 2: return 1.0;
			case 3: return 0.0;
			case 4: return 0.0;
			case 5: return 0.0;
			case 6: return K/R/4.0;
			case 7: return 0.0;
			default: exit(1);
			}
		case 2:
			switch (i) {
			case 0: return -K/R/2.0;
			case 1: return 0.0;
			case 2: return 0.0;
			case 3: return 0.0;
			case 4: return 1.0;
			case 5: return 0.0;
			case 6: return (K2/R2-ar)/R/4.0; /* ??? warum /R? */
			case 7: return -ai/4.0;
			default: exit(1);
			}
		default:	exit(1);
		}
	} else {
		exit(1);
	}

}

/*
 * shifted version: s = r^2
 * U -> r*u(s)
 * f -> f(s)
 * g -> 2*r*g(s)
 * W -> w(s)
 */
double uw(double *y, int i, int ny, int step) {		//ny=8
	double W, dW, K2, R2;
	double x, ar, ai;

	x = SSTEP(step);
	W = BASEFLOW[step];
	dW= BASEFLOW1[step];

	/*
	 * W  = baseflow_s(x,0);
	 * dW = baseflow_s(x,1);
	 */

	K2 = K*K;
	R2 = R*R;
	ar = -K2/R2-OMEGAI;
	ai = OMEGAR-K*W;


	if (0==step) {
		return uw_init(i, 1);
	}
	switch (i) {
	case 0:	return -(y[4]*K/R/2.0+y[0])/x;
	case 1: return -(y[5]*K/R/2.0+y[1])/x;
	case 2: return (y[1]*ai-y[0]*ar)/2.0;
	case 3: return (-y[1]*ar-y[0]*ai)/2.0;
	case 4: return y[6];
	case 5: return y[7];
	case 6: return -y[1]*R*dW/2.0 +( y[2]*K/R/4.0 +(y[4]*(K2/R2-ar)+y[5]*ai)/4.0 -y[6] )/x;
	case 7: return y[0]*R*dW/2.0 +(y[3]*K/R/4.0 +(y[5]*(K2/R2-ar) -y[4]*ai)/4.0 -y[7])/x;
	default:	exit(1);
	}
}

double shootuw(void) {
	int i;
	double y[8], ymax[8];
	double det[8];	//  { U_1.real, U_1.imag,  U_2.real, U_2.imag, W_1.real, W_1.imag, W_2.real, W_2.imag }
	double dreal, dimag;

	BASE=1;
	flag_print = 0; if (debug==1)	flag_print = 1;
	for (i=0; i<8; i++)	y[i] = uw_init(i, 0);
	ode(y, 8, uw, ymax);;
	det[0] = y[0];	// U1.real
	det[1] = y[1];	// U1.imag
	det[4] = y[4];  // W1.real
	det[5] = y[5];	// W1.imag

	BASE=2;
	flag_print = 0; if (debug==2)	flag_print = 1;
	for (i=0; i<8; i++)	y[i] = uw_init(i, 0);
	ode(y, 8, uw, ymax);
	det[2] = y[0]; 	// U2.real
	det[3] = y[1];	// U2.imag
	det[6] = y[4];	// W2.real
	det[7] = y[5];	// W2.imag

	// determinant
	dreal = -det[1]*det[7] + det[0]*det[6] + det[3]*det[5] - det[2]*det[4];
	dimag =  det[0]*det[7] + det[1]*det[6] - det[2]*det[5] - det[3]*det[4];
	return hypot(dreal, dimag);
}
