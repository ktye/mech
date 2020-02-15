/* $Id: uvw1.c,v 1.7 2004/09/10 09:37:56 elmar Exp $ */

/*
 * uvw1 - differential equation n = 1
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

double dBASE[3];	/* [u(0), w(0), f(0)] */

double uvw1init(int i, int derivative) {
	double gamma0, gamma1, K2, R2, a0r, a0i, a1r, a1i; 

	gamma0 = BASEFLOW[0];
	gamma1 = BASEFLOW1[0];
	/*
	 * gamma0  = baseflow_taylor(0);
	 * gamma1 = baseflow_taylor(2);
	 */

	K2 = K*K;
	R2 = R*R;
	a0r = -K2/R2-OMEGAI;
	a0i = OMEGAR-K*gamma0;
	a1r = 0.0;
	a1i = -K*gamma1;

	if (0 == derivative) {
		switch(i) {
		case 0: return dBASE[0];
		case 1: return 0.0;
		case 2: return dBASE[2];
		case 3: return 0.0;
		case 4: return dBASE[0];
		case 5: return 0.0;
		case 6: return -a0r*dBASE[0]-dBASE[2];
		case 7: return -a0i*dBASE[0]-dBASE[2];
		case 8: return dBASE[1];
		case 9: return 0.0;
		case 10: return 2.0*dBASE[1];
		case 11: return 0.0;
		default: exit(1);
		}
	} else if (1 == derivative) {
		switch (i) {
		case 0: return - dBASE[1]*K/R/3.0 - dBASE[0]*a0r/6.0 - dBASE[2];
		case 1: return - dBASE[0]*a0i/6.0;
		case 2: return (a0r*K2/R2/144.0-7.0*a1r/24.0+(a0r*a0r-a0i*a0i)/36.0)*dBASE[0] + (7.0*a0r*K/72.0/R + K2*K/72.0/R/R2)*dBASE[1] + (K2/R2/36.0 + a0r/36.0)*dBASE[2];
		case 3: return (gamma1*K/24.0+a0i*K2/R2/144.0-7.0*a1i/24.0+a0i*a0r/18.0)*dBASE[0] + (7.0*a0i*K/72.0/R)*dBASE[1] + (a0i/36.0)*dBASE[2];
		case 4: return -a0r*dBASE[0]/2.0 - dBASE[2]/2.0;
		case 5: return -a0i*dBASE[0]/2.0;
		case 6: return (-a0r*K2/144.0/R2 - 5.0*a1r/24.0 +5.0*(a0r*a0r-a0i*a0i)/36.0)*dBASE[0] + (-a0r*K/72.0/R-K*K2/72.0/R/R2)*dBASE[1] + (5.0*a0r/36.0-K2/36.0/R2)*dBASE[2];
		case 7: return (-a0i*K2/144.0/R2 - gamma1*K/24.0 - 5.0*a1i/24.0 + 5.0*a0i*a0r/18.0)*dBASE[0] + (-a0i*K/72.0/R)*dBASE[1] + (5.0*a0i/36.0)*dBASE[2];   
		case 8: return (a0r*K/24.0/R)*dBASE[0] + (K2/12.0/R2-a0r/8.0)*dBASE[1] + (K/6.0/R)*dBASE[2];
		case 9: return (gamma1*R/4.0+a0i*K/24.0/R)*dBASE[0] + (-a0i/8.0)*dBASE[1];
		case 10: return (a0r*K/6.0/R)*dBASE[0] + (K2/3.0/R2-a0r/2.0)*dBASE[1] + (2.0*K/3.0/R)*dBASE[2];
		case 11: return (gamma1*R+a0i*K/6/R)*dBASE[0] + (-a0i/2.0)*dBASE[1];
		default: exit(1);
		}
	} else {
		exit(1);
	}
}

/*
 * u:      y[0]+ %i*y[1]$
 * f:      y[2]+ %i*y[3]$
 * v:      y[4]+ %i*y[5]$
 * h:      y[6]+ %i*y[7]$
 * w:      y[8]+ %i*y[9]$
 * g:      y[10]+ %i*y[11]$
 */
double uvw1(double *y, int i, int ny, int step) {
	double s, K2, R2, s2, ar, ai, W, dW;

	if (0==step) {
		return uvw1init(i, 1);
	}

	s = SSTEP(step);

	W = BASEFLOW[step];
	dW= BASEFLOW1[step];
	/*
	 * W  = baseflow_s(s,0);
	 * dW = baseflow_s(s,1);
	 */
	K2 = K*K;
	R2 = R*R;
	ar = -K2/R2-OMEGAI;
	ai = OMEGAR-K*W;
	s2 = s*s;

	switch (i) {
	case 0: return (y[4]-y[0])/s/2.0-y[8]*K/R/2.0;
	case 1: return (y[5]-y[1])/s/2.0-y[9]*K/R/2.0;
	case 2: return y[8]*K/(R*s)/2.0+(y[1]*ai-y[0]*ar)/s/2.0-y[2]/s/2.0+3.0*(y[0]-y[4])/(2.0*s2);
	case 3: return y[9]*K/(R*s)/2.0+(-y[1]*ar-y[0]*ai)/s/2.0-y[3]/s/2.0+3.0*(y[1]-y[5])/(2.0*s2);
	case 4: return y[6]/2.0;
	case 5: return y[7]/2.0;
	case 6: return -y[8]*K/(R*s)/2.0+(y[5]*ai-y[4]*ar)/s/2.0-y[6]/s-y[2]/s/2.0+3.0*(y[4]-y[0])/(2.0*s2);
	case 7: return -y[9]*K/(R*s)/2.0+(-y[5]*ar-y[4]*ai)/s/2.0-y[7]/s-y[3]/s/2.0+3.0*(y[5]-y[1])/(2.0*s2);
	case 8: return y[10]/s/2.0-y[8]/s;
	case 9: return y[11]/s/2.0-y[9]/s;
	case 10: return (y[0]-y[4])*K/(R*s)/2.0-y[1]*dW*R+y[2]*K/R/2.0+(y[8]*(K2/R2-ar)+y[9]*ai)/2.0;
	case 11: return (y[1]-y[5])*K/(R*s)/2.0+y[0]*dW*R+y[3]*K/R/2.0+(y[9]*(K2/R2-ar)-y[8]*ai)/2.0;

	}
}

/*
* u:      y[0]+ %i*y[1]$
* f:      y[2]+ %i*y[3]$
* v:      y[4]+ %i*y[5]$
* h:      y[6]+ %i*y[7]$
* w:      y[8]+ %i*y[9]$
* g:      y[10]+ %i*y[11]$
 */
double shootuvw1(void) {
	int i;
	int base;
	double scale;
	double y[12], ymax[12];
	double dreal, dimag;
	double ur[3],ui[3],vr[3],vi[3],wr[3],wi[3];

	for (base=1; base<=3; base++) {
		for (i=0; i<3; i++) dBASE[i] = 0.0;
		dBASE[base-1]=1.0;
		flag_print = 0; 
		for (i=0; i<12; i++) y[i] = uvw1init(i,0);
		ode(y, 12, uvw1, ymax);

		scale = (hypot(ymax[0],ymax[1]) + hypot(ymax[4],ymax[5]) + hypot(ymax[8],ymax[9]))/3.0;
		if (scale>0.0)
			dBASE[base-1]=1.0/scale;
		if (debug==base)	flag_print = 1;
		for (i=0; i<12; i++) y[i] = uvw1init(i,0);
		ode(y, 12, uvw1, ymax);
		ur[base-1] = y[0];
		ui[base-1] = y[1];
		vr[base-1] = y[4];
		vi[base-1] = y[5];
		wr[base-1] = y[8];
		wi[base-1] = y[9];
	}

	dreal = ur[0]*(vr[1]*wr[2]-vi[1]*wi[2]-wr[1]*vr[2]+wi[1]*vi[2])
		-ui[0]*(vi[1]*wr[2]+vr[1]*wi[2]-wi[1]*vr[2]-wr[1]*vi[2])
		-vr[0]*(ur[1]*wr[2]-ui[1]*wi[2]-wr[1]*ur[2]+wi[1]*ui[2])
		+vi[0]*(ui[1]*wr[2]+ur[1]*wi[2]-wi[1]*ur[2]-wr[1]*ui[2])
		+wr[0]*(ur[1]*vr[2]-ui[1]*vi[2]-vr[1]*ur[2]+vi[1]*ui[2])
		-wi[0]*(ui[1]*vr[2]+ur[1]*vi[2]-vi[1]*ur[2]-vr[1]*ui[2]);
	dimag = ui[0]*(vr[1]*wr[2]-vi[1]*wi[2]-wr[1]*vr[2]+wi[1]*vi[2])
		+ur[0]*(vi[1]*wr[2]+vr[1]*wi[2]-wi[1]*vr[2]-wr[1]*vi[2])
		-vi[0]*(ur[1]*wr[2]-ui[1]*wi[2]-wr[1]*ur[2]+wi[1]*ui[2])
		-vr[0]*(ui[1]*wr[2]+ur[1]*wi[2]-wi[1]*ur[2]-wr[1]*ui[2])
		+wi[0]*(ur[1]*vr[2]-ui[1]*vi[2]-vr[1]*ur[2]+vi[1]*ui[2])
		+wr[0]*(ui[1]*vr[2]+ur[1]*vi[2]-vi[1]*ur[2]-vr[1]*ui[2]);
	return hypot(dreal, dimag);
}
