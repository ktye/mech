/* $Id: baseflow.c,v 1.13 2004/09/30 08:32:16 elmar Exp $ */

/*
 * baseflow - laminar velocity profile
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "readfile.h"

double *BASEFLOW;	/* global array with baseflow as a function of s=r^2 */
double *BASEFLOW1;	/* derivative of baseflow with respect to s */
double *BASESIGNAL;	/* baseflow and derivatives in complex form */

extern int HPFLOW;
extern unsigned int Nsteps;
extern unsigned int Lsteps;
extern double ALPHA;
extern int ALPHANUM;
extern double TIME;

extern char *basefile;

void initbaseflow() {
	int i;
	double alpha;
	double t;
	unsigned int N;

	alpha = 0.0;
	t = 0.0;
	if (HPFLOW) {
		BASEFLOW = (double *)malloc((size_t)Lsteps*sizeof(double));
		BASEFLOW1= (double *)malloc((size_t)Lsteps*sizeof(double));
		if ((!BASEFLOW)||(!BASEFLOW1)) {
			fprintf(stderr,"error: malloc baseflow [%u]\n",Lsteps);
			exit(1);
		}
		for (i=0; i<Lsteps; i++) {
			BASEFLOW[i] =  1.0 - (double)i/(double)(Lsteps-1);
			BASEFLOW1[i]= -1.0;
		}
		return;
	} else {
		BASESIGNAL = readfile(basefile, &alpha, &N);
		Lsteps = N;
		Nsteps = Lsteps/2 + 1;
		ALPHA = alpha;

		BASEFLOW = (double *)malloc((size_t)Lsteps*sizeof(double));
		BASEFLOW1= (double *)malloc((size_t)Lsteps*sizeof(double));
		for (i=0; i<Lsteps; i++) {
			BASEFLOW[i] = BASESIGNAL[2*i]*cos(M_PI*TIME) + BASESIGNAL[2*i+1]*sin(M_PI*TIME);
			BASEFLOW1[i] = BASESIGNAL[2*i+2*Lsteps]*cos(M_PI*TIME) + BASESIGNAL[2*i+1+2*Lsteps]*sin(M_PI*TIME);
		}
	}
	return;
}

void reinitbaseflow(void) {
	int i;
	for (i=0; i<Lsteps; i++) {
		BASEFLOW[i] = BASESIGNAL[2*i]*cos(M_PI*TIME) + BASESIGNAL[2*i+1]*sin(M_PI*TIME);
		BASEFLOW1[i] = BASESIGNAL[2*i+2*Lsteps]*cos(M_PI*TIME) + BASESIGNAL[2*i+1+2*Lsteps]*sin(M_PI*TIME);
	}
}

double baseflow(double x, int i) {
	perror("attic baseflow");
	exit(1);
	return 0.0;
}

double baseflow_s(double x, int i) {
	perror("attic baseflow_s");
	exit(1);
	return 0.0;
}

double baseflow_taylor(int i) {
	perror("attic baseflow_taylor");
	exit(1);
	return 0.0;
}

void printbaseflow(void) {
	int i;
	for (i=0; i<Nsteps; i++)
		printf("%g %g %g\n",(double)i/(double)(Nsteps-1),BASEFLOW[2*i],BASEFLOW1[2*i]);
}
