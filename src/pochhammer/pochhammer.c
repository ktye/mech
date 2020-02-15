/* $Id: pochhammer.c,v 1.4 2004/03/08 19:57:25 elmar Exp $ */

/*
 * pochhammer - pochhammer dispersion relation shooting algorithm
 */

static char id[]="$Id: pochhammer.c,v 1.4 2004/03/08 19:57:25 elmar Exp $";
int flag_print = 1;

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "ode.h"
#include "cubic.h"

/* global parameter */
double K = 0.0;		/* wave number squared k^2 */
//double OMEGA=0.0;	/* frequency	*/
double A=0.0;		/* OMEGA^2/A^2	*/
double B=0.0;		/* OMEGA^2/B^2	*/
double C=0.0;		/* a^2/b^2	*/
double D=0.0;		/* b^2/a^2	*/

/*
 * RHS 	Pochhammer n=0 	Longitudinal waves
 * y[0]: A
 * y[1]: U (\bar{U})
 * y[2]: W
 * y[3]: Z
 */
double schwinger(double *y, int j, int ny, double x) {
	if (j==0)	return y[1];
	if (j==1)	return -y[0];
}

/* Pochhammer frequency relation
 * omega2	omega*omega
 * kk		k^2
 * aa		a^2
 * bb		b^2
 */
/* FEHER:	Komplexe Berechnung ist hier nötig!! Sonst sind alle c<1 verloren! */
double dispersion(double omega2, double kk, double aa, double bb, int n) {
	double a,b,c;
	double d,e,f;
	double g,h,i;
	double alpha, beta, aalpha, bbeta, cc,k;
	if (omega2/aa-kk < 0.0) {
		fprintf(stderr,"error: alpha^2 < 0: %g, w^2=%g\t,aa=%g\t,kk=%g\n",omega2/aa-kk,omega2,aa,kk);
		return 10.0;
		//exit(1);
	}
	if (omega2/bb-kk < 0.0) {
		fprintf(stderr,"error: beta^2 < 0: %g, w^2=%g\t,bb=%g\t,kk=%g\n",omega2/bb-kk,omega2,bb,kk);
		return 10.0;
		//exit(1);
	}
	alpha = sqrt(omega2/aa-kk);
	aalpha = omega2/aa-kk;
	beta = sqrt(omega2/bb-kk);
	bbeta = omega2/bb-kk;
	k = sqrt(kk);
	cc = aa/bb;

//	printf("\nalpha=%f, beta=%f\naalpha=%f, bbeta=%f\nk=%f,cc=%f\n",alpha,beta,aalpha,bbeta,k,cc);

	if (n==0) {
		a = (cc*(kk + aalpha) - 2.0*kk)*j0(alpha) - 2.0*alpha*j1(alpha);
		b = 2.0*beta*j0(beta) - 2.0*k*j1(beta);
		c = 2.0*alpha*k*j1(alpha);
		d = (kk-bbeta)*j1(beta);
		return a*d-b*c;
	} else if (n==1) {
		a = ((2.0-cc)*kk+(2.0-aalpha)*cc)*j1(alpha)-2.0*alpha*j0(alpha);
		b = 2.0*j1(alpha);
		c = 2.0*alpha*k*j0(alpha)-2.0*k*j1(alpha);
		d = 2.0*cc*j1(beta)+2.0*(1.0-cc)*beta*j0(beta);
		e = (2.0-bbeta)*j1(beta) - 2.0*beta*j0(beta);
		f = -k*j1(beta);
		g = 2.0*beta*k*(1.0-cc)*j1(beta)+k*(4.0-2.0*cc)*j0(beta);
		h = -k*beta*j1(beta)-2.0*k*j0(beta);
		i = (bbeta-kk)*j0(beta) - beta * j1(beta);
		return a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g);
	} else {
		fprintf(stderr,"error: n>1 not implemented\n");
		exit(1);
	}
}
double pochhammer(double *y, int j, int ny, double r) { return 0.0; }
double pochhammer_uw(double *y, int j, int ny, double r) {
	double rr;
	rr = r*r;
	if (r==0.0) {

		/* debug */
		fprintf(stderr,"error: r=0\n");
		exit(1);

		switch (j) {
		case 0:
			return 0.0;	/* U_2 = 0 */
		case 1:
			return 0.0;
		case 2:
			return 0.0;
		case 3:
			return 1.0; 	/* bin zu faul zum rechnen (aber nicht sicher) */
		}
	}
	switch (j) {
	case 0:
		return (r*(C*K-B)+1.0/r)*y[1] + r*sqrt(K)*(D-1.0)*y[3];
		break;
	case 1:
		return y[0]/r;	
		break;
	case 2:
		return y[3];
		break;
	case 3:
		return -y[3]/r + (B+A*K)*y[2] + sqrt(K)*(A-1.0)*(y[0]+y[1])/r;
		break;
	default:
		fprintf(stderr,"ERROR: your wrong here.\n");
		exit(1);
	}
}

/* error at boundary */
double rb(double *y, int n) {
	double rb1, rb2;
	if (n==0) {
		rb1 = A*y[0] + (A - 2.0)*(y[1] + sqrt(K)*y[2]);
		rb2 = sqrt(K)*y[1] - y[3];
		return (fabs(rb1) + fabs(rb2));
	} else {
		return 0.0;
	}
}

/* initial values */
void initialise(double *y, int n) {
	if (n==0) {
		y[0] = 0.0;
		y[1] = 0.0;
		y[2] = 1.0;
		y[3] = 0.0;
	} else {
	}
}
void testr0() {
	int i;
	double y[4],x,r;
	initialise(y,0);
	for (r=10.0; r>0.00000001; r/=2.0) {
		printf("\n%g: ",r);
		for (i=0; i<4; i++) {
			x = pochhammer_uw(y, i, 4, r) ;
			printf("%g ",x);
		}
	}
	exit(0);
}

double max(double x, double y) {
	if (x>y) return x;
	else return y;
}

double aim(double x, double range) {
	return x + range*(drand48()-0.5);
}



double quadroot(double p, double q, int i) {
	double x;
	x = p*p/4.0 -q;
	if (x<0) return 0.0;
	if (i) return -p/2 + sqrt(x);
	else return -p/2 - sqrt(x);
}
double dispers(double kk, double aa, double bb, int i) {
	double x,f,a,b,c,d,co[4],root[3];
	double kappa = 0.9;
	switch(i) {
	case 0:
		// rayleigh love
		a = (3.0*aa-4.0*bb)*bb/(aa-bb);
		b = (aa-2.0*bb)/(4.0*(aa-bb));
		f = a*kk/(1.0+b*kk);
		return f;
	case 1:
		// rayleigh love pochhammer
		a = (3.0*aa-4.0*bb)*bb/(aa-bb);
		b = (2.0*aa-4.0*bb)*bb/(4.0*(aa-bb));
		c = (6.0*aa*aa-3.0*aa*bb-4.0*bb*bb)/(8.0*aa*(aa-bb));
		f = a*kk/(1.0-b*kk*kk+c*kk); 
		return f;
	case 2:
		// ebb
		f = (2.0*aa*bb-4.0*bb*bb)*kk*kk/(2.0*aa-4.0*bb);
		return f;
	case 3:
		// timo1
		a = (((aa*bb-bb*bb)*kappa)+3.0*aa*bb-4.0*bb)/(aa-bb);
		b = (3.0*aa*bb*bb-4.0*bb*bb)*kappa/(aa-bb);
		c = 4.0*bb*kappa;
		f = quadroot(a*kk-c, b*kk*kk, 0);
		return f;
	case 4:
		// timo2
		a = (((aa*bb-bb*bb)*kappa)+3.0*aa*bb-4.0*bb)/(aa-bb);
		b = (3.0*aa*bb*bb-4.0*bb*bb)*kappa/(aa-bb);
		c = 4.0*bb*kappa;
		f = quadroot(a*kk-c, b*kk*kk, 1);
		return f;
	case 5:
		a = (3.0*aa*bb+bb*bb)/(3.0*aa+3.0*bb);
		b = aa*bb*bb/(3.0*aa+3.0*bb);
		c = (79.0*aa*bb*bb*bb*bb-159.0*aa*aa*bb*bb*bb+54.0*aa*aa*aa*bb*bb)/(36.0*bb*bb*bb-144.0*aa*bb*bb-108.0*aa*aa*bb+72.0*aa*aa*aa);
		c = 0.0;
		f = quadroot(-a*kk, b*kk*kk+c*kk*kk*kk, 0);
		// timo1 pc
		return f;
	case 6:
		a = (3.0*aa*bb+bb*bb)/(3.0*aa+3.0*bb);
		b = aa*bb*bb/(3.0*aa+3.0*bb);
		c = (79.0*aa*bb*bb*bb*bb-159.0*aa*aa*bb*bb*bb+54.0*aa*aa*aa*bb*bb)/(36.0*bb*bb*bb-144.0*aa*bb*bb-108.0*aa*aa*bb+72.0*aa*aa*aa);
		c = 0.0;
		f = quadroot(-a*kk, b*kk*kk+c*kk*kk*kk, 1);
		// timo2 pc
		return f;
	case 7:
		//rayleigh surface waves (1st root, should all be the same)
		co[0] = 1.0;
		co[1] = -8.0*bb*kk;
		co[2] = -16.0*bb*bb*bb*kk*kk/aa + 24.0*bb*bb*kk*kk;
		co[3] = 16.0*bb*bb*bb*bb*kk*kk*kk/aa - 16.0*bb*bb*bb*kk*kk*kk;
		solve3(co,root);
		return root[0];
	case 8:
		// rayleigh surface short wave limit
		return kk*4.0*(1.0-bb*sqrt(bb)/sqrt(aa));
	}
}
void usage(char **argv) {
	fprintf(stderr,"usage: %s kmax N aa bb mode\n", argv[0]);
	exit(1);
}
int main (int args, char **argv) {
	/*
	int i,j,n,N;
	double c,k,dk, kk, ww, w0, aa, bb, omegamin;
	double f,try;
	double range,error,best;
	int mag,maxmag;
	*/
	int i, j, N, mode;
	double k, kk, kmax, dk, w2, c, aa, bb;
//	double y[6];
	if (args != 6) usage(argv);
	kmax=atof(argv[1]);
	N = atoi(argv[2]);
	aa = atof(argv[3]);
	bb = atof(argv[4]);
	mode = atoi(argv[5]);

	

	for (i=1; i<N; i++) {
		k = (double)i*kmax/(double)(N-1);
		printf("%g\t",k);
		kk = k*k;
		for (j=0; j<9; j++) {
			w2 = dispers(kk, aa, bb, j);
			if (mode)	printf("%g\t",sqrt(w2));	
			else		printf("%g\t",sqrt(w2)/k);	
		}
		printf("\n");
	}

	/*
	w0 = atof(argv[1]);
	range = atof(argv[2]);
	dk = atof(argv[3]);
	aa = atof(argv[4]);
	bb = atof(argv[5]);
	n = atoi(argv[6]);

	error=0.001;
	maxmag=10000;
	
	N = atoi(argv[7]);
	srand48(time(0));
	*/
	

	/* plot asymptotik */


	/* bild mode
	c=0.0;
	for (i=0; i<N; i++) {
		k=0.0;
		for (j=0; j<N; j++) {
			kk = k*k;
			ww = c*c*kk;
			if ((c*c<aa)||(c*c<bb)) f=10.0;
			else f = fabs(dispersion(ww,kk,aa,bb,n));
			printf("%g\n",f);
			k+=dk;
		}
		c += range;
	}
	return 0;
	bild mode end */


	/* shooting mode
	kk=0.0;
//	w=14.57*bb; 	// j^2_1,1 = 14.57
	ww=w0;
	for (i=0; i<N; i++) {
		mag=0;
		best=100000;
//		printf("%g %g (%g)\n",kk,w,dispersion(w,kk,aa,bb,n));
		f=1.0;
		while (f>error) {
			mag++;
			try = aim(ww,range);
			if ((try<kk*aa)||(try<kk*bb)) f=10.0;
			else f = fabs(dispersion(try,kk,aa,bb,n));
			if (f<best)	best=f;
			if (mag>maxmag) if (!(mag%100000)) fprintf(stderr,"error: mag=%d kk=%g\ttry=%g\tf=%g\tbest=%g\n",mag,kk,try,f,best);
			if (mag>1000*maxmag) return 1;
		}
		ww = try;
		printf("%g %g (%d) (%g)\n",kk,ww,mag,f);
		kk+=dk;	
	}
	*/

	return 0;
}

/* fuer die odes
 	A = omega*omega/aa;
	B = omega*omega/bb;
	C = aa/bb;
	D = bb/aa;
	oomega = omega*omega;
initialise(y,n);
heun(y, 4, pochhammer_uw, 0.0001, 1.0/(double)N, N);
flag_print = 0;
fprintf(stderr,"heun\t--> %g\n",rb(y,n));
*/
