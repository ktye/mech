/* $Id: loop.c,v 1.6 2005/04/29 14:35:07 elmar Exp $ */
#include <stdio.h>
#include <stdlib.h>
#include "define.h"
#include "base.h"
#include "baseflow.h"
#include "solve.h"

extern double	 R;
extern double	 K;
extern double	 T;
extern double	 L;
extern int	 tflag, kflag, rflag;


void setmax(double *max, double r, double k, double t, double sigma, double omega) {
	max[0] = r;
	max[1] = k;
	max[2] = t;
	max[3] = sigma;
	max[4] = omega;
}

void maxcp(double *src, double *dst) {
	dst[0] = src[0];
	dst[1] = src[1];
	dst[2] = src[2];
	dst[3] = src[3];
	dst[4] = src[4];
}

void printmax(double *x, double dt) {
	printf("%g %.16f %.16f %.16f %.16f %g\n", x[0]/L, x[1]/L, x[2], x[3]/L, x[4]/L, dt);
}

/*
void loop(diterator *rey, diterator *kay, diterator *time) {
	int	 iR, iK, iT;
	double	 sigma, omega;
	for (iR=0; iR<rey->n; iR++) {
		SETCNT(rey, iR);
		R = rey->val;
		base_init();
		for (iK=0; iK<kay->n; iK++) {
			SETCNT(kay, iK);
			K = kay->val;
			base_init();
			for (iT=0; iT<time->n; iT++) {
				SETCNT(time, iT);
				T = time->val;
				baseflow_init();
				solve(&sigma, &omega);
				printf("%g %.16f %.16f %.16f %.16f\n", R, K, T, sigma, omega);
			}
		}
	}
}
*/

void loop(diterator *rey, diterator *kay, diterator *time) {
	int	 iR, iK, iT, dt, dtmax = 0;
	double	 sigma, omega;
	double	 tmax[5], kmax[5], rmax[5];
	for (iR=0; iR<rey->n; iR++) {
		SETCNT(rey, iR);
		R = rey->val;
		base_init();
		for (iK=0; iK<kay->n; iK++) {
			SETCNT(kay, iK);
			K = kay->val;
			base_init();
			dt = 0;
			for (iT=0; iT<time->n; iT++) {
				SETCNT(time, iT);
				T = time->val;
				baseflow_init();
				solve(&sigma, &omega);
				if (sigma > 0.0)
					dt++;
				if (!iT) {
					setmax(tmax, R, K, T, sigma, omega);
				} else {
					if (sigma > tmax[3])
						setmax(tmax, R, K, T, sigma, omega);
				}
				if (!tflag)
					printf("%g %.16f %.16f %.16f %.16f\n", R/L, K/L, T, sigma/L, omega/L);
			}
			if (!iK) {
				maxcp(tmax, kmax);
				dtmax = dt;
			} else {
				if (tmax[3] > kmax[3]) {
					maxcp(tmax,kmax);
					dtmax = dt;
				}
			}
			if (!kflag&&tflag)
				printmax(tmax,(double)dtmax/(double)time->n);

		}
		if (!iR) {
			maxcp(kmax, rmax);
		} else {
			if (kmax[3] > rmax[3])
				maxcp(kmax,rmax);
		}
		if (!rflag&&kflag)
			printmax(kmax,(double)dtmax/(double)time->n);
	}
	if (rflag)
		printmax(rmax,0.0);
}
