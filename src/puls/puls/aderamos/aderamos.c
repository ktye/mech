/* $Id: aderamos.c,v 1.1 2003/12/08 17:02:34 elmar Exp $ */

/*
 * aderamos - pulsatile pipe flow
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "amoswrap.h"
#include "complex.h"

#define VECDIM 250

double scale (double a, double b, double A, double B, double x) {
	return (x*(A-B)/(a-b)+(a*B-b*A)/(a-b));
}

double w(double r, double t, double alpha, double *R, int NR, double *phase) {
	double Jr_r, Jr_i;
	double J_r, J_i;
	double JJ_r, JJ_i;
	double Wr, Wi;
	double exp_r, exp_i;
	double wr, wi;
	double W,abs,p;
	double N;
	int n,k;

	W = 0.0;
	*phase = 0.0;
	for (k=1; k<NR; k+=2) {
		n = (1+k)/2;
		N=(double)n;

		J0(sqrt(N)*alpha*M_SQRT1_2*r,sqrt(N)*alpha*M_SQRT1_2 *r,&Jr_r,&Jr_i);
		J0(sqrt(N)*alpha*M_SQRT1_2,sqrt(N)*alpha*M_SQRT1_2 ,&J_r,&J_i);
		cdiv(Jr_r,Jr_i,J_r,J_i,&JJ_r,&JJ_i);
		JJ_r -= 1.0;
		
		cmult(-R[k+1]/(N*pow(alpha,4.0)),R[k]/(N*pow(alpha,4.0)),JJ_r,JJ_i,&Wr,&Wi);

		exp_r = cos(N*t);
		exp_i = -sin(N*t);

		wr = 0.0;
		wi = 0.0;
		cmult(Wr,Wi,exp_r,exp_i,&wr,&wi);
		W += wr;
		cxyrp(wr, wi, &abs, &p);
		*phase += p;
	}

	return W+R[0]*(1.0-r*r)/(4.0*alpha*alpha);
}
void usage(char **argv) {
	fprintf(stderr,"usage: %s alpha R TIME\n",argv[0]);
	exit(1);
}
int main(int args, char **argv) {
	double R[3];	
	double t,alpha;
	double x,y,xmin,xmax,phase;
	int i,steps;
	char *var;
	R[0] = 0.0; 	// Zero mean flow
	R[2] = 0.0; 	// No imag part

	if (args!=4)
		usage(argv);
	alpha = atof(argv[1]);
	R[1] = atof(argv[2]);
	t = atof(argv[3]);
	xmin=0.0;
	xmax=1.0;
	steps=100;

	if (var = getenv("steps"))	steps=atoi(var);
	if (var = getenv("xmin")) xmin=atof(var);
	if (var = getenv("xmax")) xmax=atof(var);

	for (i=0; i<steps; i++) {
		x=scale(0,steps-1,xmin,xmax,i);
		y=w(x,t,alpha,R,3,&phase);
		printf("%.20f %.20f %.20f\n",t,x,y);
	}

	return 0;
}
