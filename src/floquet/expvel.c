/* $Id: expvel.c,v 1.9 2004/12/01 10:16:21 elmar Exp $ */

/*
 * expvel - expansion velocity functions
 */

static char id[]="$Id: expvel.c,v 1.9 2004/12/01 10:16:21 elmar Exp $";
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "In.h"
#include "define.h"
#include "quad.h"

extern struct expvelfn { double   u[RSIZE][ESIZE]; double   v[RSIZE][ESIZE]; double   w[RSIZE][ESIZE]; double ev[ESIZE]; int rsize; int esize; } EXPVEL;
extern int 	Mflag;
extern int 	Eflag;

double *readzerojn(int n, int maxindex) {
	int i;
	double *zerojn;
	FILE *fp;
	char filename[1024];
	snprintf(filename,sizeof(filename),ZEROJNFILE,n);

	fp = fopen(filename,"r");
	if (!fp) {
		perror("readzerojn: fopen");
		fprintf(stderr,"error in readzerojn: cannot open file %s\n",filename);
		exit(1);
	}
	zerojn = (double *)malloc((size_t)(maxindex)*sizeof(double));
	if (!zerojn) {
		perror("readzerojn: malloc");
		exit(1);
	}
	for (i=0; i<maxindex; i++) {
		if (fscanf(fp,"%lf",&zerojn[i]) != 1) {
			perror("readzerojn: scanf");
			exit(1);
		}
	}
	fclose(fp);
	return zerojn;
}


double expvel_disrel(int N, double k, double beta) {
	double n;	
	n = (double)N;

	if (!N)		
		return k*iv(0.0,k)*j1(beta) - beta*iv(1.0,k)*j0(beta);
	if (k==0.0) 
		return beta*jn(N,beta)*jn(N+1,beta);
	return (  2.0*k* iv(n,    k)/iv(n-1.0,k) * jn(N-1,beta) * jn(N+1,beta)
		- beta * iv(n+1.0,k)/iv(n-1.0,k) * jn(N-1,beta) * jn(N,  beta)
		+ beta *                           jn(N+1,beta) * jn(N,  beta));
}

void solve_disrel0(int n, int maxindex) {
	int i;
	double *zerojn, *zerojnplus;

	if (n) {
		zerojn = readzerojn(n, maxindex);
		zerojnplus = readzerojn(n+1, maxindex);
		for (i=0; i<maxindex-1; i++) {
			EXPVEL.ev[i] = zerojn[i/2];
			i++;
			EXPVEL.ev[i] = zerojnplus[i/2];
		}
		if (i<maxindex) {
			EXPVEL.ev[i] = zerojn[i/2];
		}
		free(zerojn);
		free(zerojnplus);
	} else { // n == 0 --> beta_i = j2,i 
		zerojn = readzerojn(2, maxindex);
		for (i=0; i<maxindex; i++)
			EXPVEL.ev[i] = zerojn[i];
	}
}

int disrel_fzero(int n, double k, int i, double *X, double dx, double *initerror) {
	double x0, x1, y0, y1;
	double x, y, xold;
	int iters = 0;
	double EPS = 10E-12;

	x0 = *X - dx;
	x1 = *X + dx;
	y0 = expvel_disrel(n, k, x0);
	y1 = expvel_disrel(n, k, x1);
	if (((y0>0.0)&&(y1>0.0))||((y0<0.0)&&(y1<0.0))) {
		// fprintf(stderr,"disrel_fzero: wrong intervall start=[%g %g] end=[%g %g]\n",x0,y0,x1,y1);
		return 0;
	}
	x = (x0+x1)/2.0;
	xold = -1.0;
	while (x1-x0 > EPS) {
		xold = x;
		x = (x0+x1)/2.0;
		y = expvel_disrel(n,k,x);
		if (((y0>0.0)&&(y>0.0))||((y0<0.0)&&(y<0.0)))
			x0 = x;
		else
			x1 = x;
		iters++;
	}
	*initerror = fabs(*X-x);
	*X = x;
	return iters;
}

void solve_disrel(int n, double K, int maxindex) {
	int i;
	double dx, k, dk, initerror=0;
	int iters;
	dx = 0.1;
	solve_disrel0(n, maxindex);
	for (i=0; i<maxindex; i++) {
		iters = 0;
		k = 0.0;
		dk = 0.1;
		//fprintf(stderr,"n=%d k=%g i=%d\tev: %.16f --> ",n,k,i,EXPVELEV[i]);
		while (k<=K) {
			iters = 0; 
			iters = disrel_fzero(n, k, i, &(EXPVEL.ev[i]), dx, &initerror);
			while (!iters) {
				//fprintf(stderr,"k=%g dk=%g\n",k,dk);
				k -= dk;
				dk /= 2.0;
			}
			if (initerror/dx<0.1) {
				dk *= 2.0;
			} else if (initerror/dx>0.8) {
				dk /= 2.0;
			}
			if (k==K) break;
			k += dk; if (k>K) k = K;
		}
		//fprintf(stderr,"%.16f iters=%d dk=%g k=%g\n",EXPVELEV[i],iters,dk,k);
				/*
			if (disrel_fzero(n, k, i, &EXPVELEV[i], dx)) {
				wrong_code;
			}
			*/
	}
}

/*
void expansion_pressure(int index, int n, double k, double r, double *u, double *v, double *w) { 
	double Inplus_k, Inplus_kr;
	double Inminus_k,Inminus_kr;
	double In_kr;
	Inplus_k   = iv((double)(n+1), k);	
	Inplus_kr  = iv((double)(n+1), k*r);
	Inminus_k  = iv((double)(n-1), k);
	Inminus_kr = iv((double)(n-1), k*r);
	In_kr      = iv((double)(n),   k*r);	        
	*u = k*(-Inminus_kr - Inplus_kr)/2.0;
	*v = k*( Inminus_kr - Inplus_kr)/2.0;
	*w = k*In_kr;
}
*/

void expvel(int index, int n, double k, double r, double *u, double *v, double *w) {
	double Inplus_k,  Jnplus_b,    Inplus_kr,    Jnplus_br;
	double Inminus_k, Jnminus_b,   Inminus_kr,   Jnminus_br;
	/*
	double  dInplus_kr,  dJnplus_br,  dInminus_kr,  dJnminus_br;
	double ddInplus_kr, ddJnplus_br, ddInminus_kr, ddJnminus_br;
	*/
	double   In_kr,  Jn_br;
	/*
	double  dIn_kr, dJn_br;
	double ddIn_kr,ddJn_br;
	*/
	double beta;

	beta = EXPVEL.ev[index];

	Inplus_k   = iv((double)(n+1), k);	
	Jnplus_b   = jn(         n+1,  beta);
	Inplus_kr  = iv((double)(n+1), k*r);
	Jnplus_br  = jn(         n+1,  beta*r);
	Inminus_k  = iv((double)(n-1), k);
	Jnminus_b  = jn(         n-1,  beta);
	Inminus_kr = iv((double)(n-1), k*r);
	Jnminus_br = jn(         n-1,  beta*r);
	In_kr      = iv((double)(n),   k*r);	        
	Jn_br      = jn(n,             beta*r);	


	if (n==0) {
		*u = - k*Inplus_kr + k*Inplus_k/Jnplus_b*Jnplus_br;
		*v = 0.0;
		*w = k*In_kr - beta*Inplus_k/Jnplus_b*Jn_br;
		return;
	}
	if (k==0.0) {
		fprintf(stderr,"special case k->0 not implemented\n");
		exit(1);
	}

	/*
		if (beta == 0.0) { 	// n>0 beta=0  -->  [[ beta must be != 0 ]] (I guess)
		*u = k*(pow(r,(double)(n+1))*Inplus_k + pow(r,(double)(n-1))*Inminus_k - Inminus_kr - Inplus_kr)/2.0;
		*v = k*(pow(r,(double)(n+1))*Inplus_k - pow(r,(double)(n-1))*Inminus_k + Inminus_kr - Inplus_kr)/2.0;
		*w = k*In_kr - (double)(n+1)*pow(r,(double)n)*Inplus_k - (double)(n<1)*Inminus_k;
		} else {		
	*/
		// n>0 beta>0
	//if (derivative==0) {
	*u = k*( Inplus_k*Jnplus_br/Jnplus_b + Inminus_k*Jnminus_br/Jnminus_b - Inminus_kr - Inplus_kr)/2.0;
	*v = k*( Inplus_k*Jnplus_br/Jnplus_b - Inminus_k*Jnminus_br/Jnminus_b + Inminus_kr - Inplus_kr)/2.0;
	*w = k*In_kr - beta*(Inplus_k/Jnplus_b - Inminus_k/Jnminus_b)*Jn_br/2.0;
	/*
	} else if (derivative==1) {
		dInplus_kr   = k*(In_kr + iv((double)(n+2), k*r))/2.0;
		dInminus_kr  = k*(iv((double)(n-2),k*r) + In_kr)/2.0;
		dIn_kr       = k*(Inminus_kr + Inplus_kr)/2.0;
		dJnplus_br   = beta*(Jn_br - jn(n+2,beta*r))/2.0;
		dJnminus_br  = beta*(jn(n-2,beta*r) - Jn_br)/2.0;
		dJn_br       = beta*(Jnminus_br-Jnplus_br)/2.0;
		*u = k*( Inplus_k*dJnplus_br/Jnplus_b + Inminus_k*dJnminus_br/Jnminus_b - dInminus_kr - dInplus_kr)/2.0;
		*v = k*( Inplus_k*dJnplus_br/Jnplus_b - Inminus_k*dJnminus_br/Jnminus_b + dInminus_kr - dInplus_kr)/2.0;
		*w = k*dIn_kr - beta*(Inplus_k/Jnplus_b - Inminus_k/Jnminus_b)*dJn_br/2.0;
	} else if (derivative==2) {
		ddInplus_kr   = k*k*(Inminus_kr+2.0*Inplus_kr+iv((double)(n+3),k*r))/4.0;
		ddInminus_kr  = k*k*(iv((double)(n-3),k*r),+2.0*Inminus_kr+Inplus_kr)/4.0;
		ddIn_kr       = k*k*(iv((double)(n-2),k*r)+2.0*In_kr+iv((double)(n+2),k*r))/4.0;
		ddJnplus_br   = beta*beta*(Jnminus_br-2.0*Jnplus_br+jn(n+3,beta*r))/4.0;
		ddJnminus_br  = beta*beta*(jn(n-3,beta*r)-2.0*Jnminus_br+Jnplus_br)/4.0;
		ddJn_br       = beta*beta*(jn(n-2,beta*r)-2.0*Jn_br+jn(n+2,beta*r))/4.0;
		*u = k*( Inplus_k*ddJnplus_br/Jnplus_b + Inminus_k*ddJnminus_br/Jnminus_b - ddInminus_kr - ddInplus_kr)/2.0;
		*v = k*( Inplus_k*ddJnplus_br/Jnplus_b - Inminus_k*ddJnminus_br/Jnminus_b + ddInminus_kr - ddInplus_kr)/2.0;
		*w = k*ddIn_kr - beta*(Inplus_k/Jnplus_b - Inminus_k/Jnminus_b)*ddJn_br/2.0;
	} else {
		fprintf(stderr,"error in expvel: unknown derivative\n");
		exit(1);
	}
	*/

	/* multiply vel->r by factor (-i^n), multiply vel->phi and vel->z by factor (i^(n+1))
	switch (n%4) {
	case 0:
		vel->r.r *= -1.0;		vel->r.i = 0.0;
		vel->phi.i = vel->phi.r;	vel->phi.r = 0.0;
		vel->z.i = vel->z.r;		vel->z.r = 0.0;
		break;
	case 1:
		vel->r.i = - vel->r.r;		vel->r.r = 0.0;
		vel->phi.r *= -1.0;		vel->phi.i = 0.0;
		vel->z.r *= -1.0;		vel->z.i = 0.0;
		break;
	case 2:
		vel->r.i = 0.0;
		vel->phi.i = -vel->phi.r;	vel->phi.r = 0.0;
		vel->z.i = -vel->z.r;		vel->z.r = 0.0;
		break;
	case 3:
		vel->r.i = vel->r.r;		vel->r.r = 0.0;
		vel->phi.i = 0.0;
		vel->z.i = 0.0;
		break;
	default:
		fprintf(stderr,"internal error\n");
		exit(1);
	}
	*/
}

void expvel_init(int n, double k, int N) {
	int i, j;
	double r, s, t;
	double x[RSIZE];
//	EXPVEL.ev = malloc((size_t)(N)*sizeof(double));
	if (N>ESIZE) {
		fprintf(stderr,"expvel_init: error N > ESIZE; recompile\n");
		exit(1);
	}
	solve_disrel(n,k,N);
	if (Eflag) {
		for (i=0; i<N; i++) {
			printf("EXPVELEV[%d] = %.16f\n",i,EXPVEL.ev[i]);
		}
	}

	for (j=0; j<N; j++) {
		for (i=0; i<RSIZE; i++) {
			r = (double)i/(double)(RSIZE-1);
			expvel(j, n, k, r, &(EXPVEL.u[i][j]), &(EXPVEL.v[i][j]), &(EXPVEL.w[i][j]));
//			expansion_pressure(j, n, k, r, &(Pr[i][j]), &(Pphi[i][j]), &(Pz[i][j]));
			s = EXPVEL.u[i][j]; s *= s;
			t = EXPVEL.v[i][j]; s += t*t;
			t = EXPVEL.w[i][j]; s += t*t;
			x[i] = r*s;
		}
		s = quad(x);
		for (i=0; i<RSIZE; i++) {
			EXPVEL.u[i][j] /= sqrt(s);
			EXPVEL.v[i][j] /= sqrt(s);
			EXPVEL.w[i][j] /= sqrt(s);
			r = (double)i/(double)(RSIZE-1);
		}
	}
}

void expvel_ortho(N) {
	int i,j,r;
	double x, X[RSIZE], max = 0.0;
	for (i=0; i<N; i++) {
		for (j=0; j<N; j++) {
			for (r=0; r<RSIZE; r++) {
				X[r] = EXPVEL.u[r][i]*EXPVEL.u[r][j];
				X[r]+= EXPVEL.v[r][i]*EXPVEL.v[r][j];
				X[r]+= EXPVEL.w[r][i]*EXPVEL.w[r][j];
				X[r]*= (double)r/(double)(RSIZE-1);
			}
			x = quad(X);
			if (i!=j) if (fabs(x) > max) max = fabs(x);
			if (i==j) if (fabs(1.0-x)>max) max = fabs(1.0-x);
			printf("%.5f ",x);
		}
		printf("\n");
	}
	fprintf(stderr,"max_ij { I - <v_i, v_j> } = %g\n",max);
}
