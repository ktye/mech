/* $Id: quasisteady.c,v 1.12 2006/01/18 15:52:45 elmar Exp $ */

#include <stdio.h>
#include <math.h>
#include <g2c.h>
//#include "hooke.h"
#include "define.h"
#include "reinit.h"
#include "solve.h"


static char id[]="$Id: quasisteady.c,v 1.12 2006/01/18 15:52:45 elmar Exp $";

#define ISBAD(x) (isinf(z)||isnan(z))	

extern int	oflag;
extern int	Tflag;
extern double R;
extern double ALPHA;
extern double K;
extern double TIME;
extern int Nphi;
extern int NE;

extern double SIGMA[ESIZE];
extern double OMEGA[ESIZE];

#ifdef DEBUG
extern struct statstruct {
        unsigned long i;
        unsigned long jobs;
        time_t starttime;
        diterator *t, *R, *k;
        iiterator *n, *N, *alpha;
	double Nerr;
	int Nadapt;
} stat;
#endif /* DEBUG */

int increase(iiterator *N, int inc) {
	if (N->val > ESIZE-1)
		return 1;
	if (N->max - inc >= N->val) {
		N->val += inc;
		return 0;
	} else if (N->max > N->val) {
		N->val = N->max;
		return 0;
	} else 	return 1;
}
int decrease(iiterator *N, int dec) {
	if (N->min + dec <= N->val) {
		N->val -= dec;
		return 0;
	} else if (N->min < N->val) {
		N->val = N->min;
		return 0;
	} else 	return 1;
}

double relmax(double a, double b, double c, double d) {
	double A,B;
	if (a!=0.0)  	A = fabs((b-a)/a);
	else		A = fabs(b);
	if (c!=0.0)	B = fabs((d-c)/c);
	else		B = fabs(d);
	return ( (A) > (B) ? (A):(B));
}

void quasisteady(diterator *Rey, iiterator *alphanum, diterator *t, diterator *k, iiterator *n, iiterator *N, tuple *result) {
	int i;
	int unstable_times = 0;

	result->R = Rey->val;
	result->n = n->val;
	result->N = N->val;
	result->a = alphanum->val;
	result->sigma = -10E17;

	for (i = 0; i<t->n; i++) {
		SETCNT(t,i);
		reinit(Rey->val,alphanum->val,k->val,n->val,t->val,N->val);
		if (solve()) {
			fprintf(stderr,"floating point exception for R=%g A=%g n=%d K=%.16f t=%.16f\n", R,ALPHA,Nphi,K, t->val);
		} else if (!i) {
			result->sigma = SIGMA[0];
			result->omega = OMEGA[0];
			result->t = TIME;
			result->k = K;
		} else if (SIGMA[0] > result->sigma) {
			result->sigma = SIGMA[0];
			result->omega = OMEGA[0];
			result->t = TIME;
			result->k = K;
		}
		if (SIGMA[0] > 0.0)
			++unstable_times;
//		fprintf(stderr,"%.16f %.16f %.16f %.16f\n", SIGMA[0], OMEGA[0], TIME, K);
		if (Tflag)
			printf("T: %g %g %d %.16f %.16f %.16f %.16f\n",R,ALPHA,Nphi,K,SIGMA[0],OMEGA[0],TIME);
#ifdef DEBUG
		++stat.i;
#endif /* DEBUG */
	}
	/* omega has opposite sign! */
	if (!oflag) printf("%g %g %d %.16f %.16f %.16f %.16f %g\n", result->R,ALPHA,result->n,result->k,result->sigma, result->omega, result->t, 100.0*(double)unstable_times/(double)t->n);
}

/*
double opti_fn(double x[2], int n) {
	double k,t;
	t = x[0];
	k = x[1];
	reinit(R, (int)(ALPHA/5.0), k, Nphi, t, NE);
	if (solve()) return 1000.0;
	else return -SIGMA[0];
}

void quasisteady_opti(tuple *result) {
	int iters, itermax = 100;
	double start[2], end[2];
	double epsilon = 1E-5, rho = 0.5;

	start[0] = result->t;
	start[1] = result->k;
	reinit(result->R,result->a,result->k,result->n,result->t,result->N);
	iters = hooke(2, start, end, rho, epsilon, itermax, opti_fn);
	result->t = end[0];
	result->k = end[1];
	reinit(result->R,result->a,result->k,result->n,result->t,result->N);
	solve();
	result->sigma = SIGMA[0];
	result->omega = OMEGA[0];
	fprintf(stderr,"iters=%d t: %g --> %g\tk:%g --> %g\n",iters,start[0],end[0],start[1],end[1]);
	printf("%g %g %d %.16f %.16f %.16f %.16f\n", result->R,ALPHA,result->n,result->k,result->sigma, result->omega, result->t);
}
*/
