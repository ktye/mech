/* $Id: reinit.c,v 1.6 2005/05/26 08:52:27 elmar Exp $ */

#include <stdlib.h>
#include <stdio.h>
#include "expvel.h"
#include "readbase.h"

static char id[]="$Id: reinit.c,v 1.6 2005/05/26 08:52:27 elmar Exp $";

extern double 	R;
extern double 	ALPHA;
extern char	*FILENAME;
extern double 	K;
extern double 	TIME;
extern int 	Nphi;
extern int 	NE;
extern int	Ninit;

/* Reynolds Number */
void reinit_R(double r) {
	R = r;
}

/* Womersley number alpha, filename with baseflow */
void reinit_alpha(int alpha) {
	if ((alpha < 0) || (alpha>999)) {
		fprintf(stderr,"error in reinit_alpha: alpha out of [0,99]\n");
		exit(1);
	}
	readbase(alpha);
	ALPHA = (double)alpha/5.0;
	readbase_settime();
}

/* axial wavenumber k */
/* angular wavenumber n */
/* number of expansion functions N */
void reinit_nk(double k, int n, int N) {
	K = k;
	Nphi = n;
	NE = N;
	expvel_init(Nphi, K, NE);
}

/* time t */
void reinit_t(double t) {
	TIME = t;
	readbase_settime();
}

/* reinit */
void reinit(double Reynolds, int alphanum, double k, int n, double t, int N) {
	//fprintf(stderr,"reinit:\n\tR=%g [%g]\n\talpha=%g [%g]\n\tk=%g [%g]\n\tn=%d [%d]\n\tt=%g [%g]\n\tN=%d [%d]\n",Reynolds,R,(double)alphanum/5.0,ALPHA,k,K,n,Nphi,t,TIME,N,NE);
	if (R!=Reynolds) 
		reinit_R(Reynolds);
	if (ALPHA!=(double)(alphanum)/5.0) {
		reinit_alpha(alphanum);
	}
	if ((K!=k)||(Nphi!=n)) 
		reinit_nk(k, n, N);
	if (N!=NE) {
		if (N>Ninit)
			reinit_nk(k, n, N);
		else
			NE = N;
	}
	if (TIME!=t) 
		reinit_t(t);
}
