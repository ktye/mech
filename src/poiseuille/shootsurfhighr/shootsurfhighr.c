/* $Id: shootsurfhighr.c,v 1.1 2004/04/21 13:53:03 elmar Exp $ */

/*
 * shootsurfhighr - shooting algorithm high Reynolds limit
 */

static char id[]="$Id: shootsurfhighr.c,v 1.1 2004/04/21 13:53:03 elmar Exp $";
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

struct doublecomplex {
	double r;
	double i;
};
typedef struct doublecomplex cdouble;

double scale (double a, double b, double A, double B, double x) {
	return (x*(A-B)/(a-b)+(a*B-b*A)/(a-b));
}

/* f'' + 1/r*f' + (a+b*r^2-1/r^2)*f = 0			*/
/* f(0) = 0, f'(0) = 1					*/
/* solve initial value problem, return end value @ r=1 	*/

/* transformed ode f = r*g(s), s=r^2 	*/
/* 4sg'' + 8g' +(a+bs)g = 0 		*/
/* g(0) = 1, g'(0) = -a/12		*/
double shoot2(double kz, double omegar, double omegai, int N) {
	int k;
	double K, h; 
	double r, k_n;
	cdouble g, g_, g__; 
	cdouble a;
	double bi;
	h = 1.0/(double)N;

	a.r = -omegai;
	a.i = omegar-kz;
	bi = kz;

	g__.r = 1.0;
	g__.i = 0.0;
	g_.r  = 1.0-a.r*h/8.0;
	g_.i  = -a.i*h/8.0;
	for (k=1; k<N; k++) {
		K = (double)k;
		g.r = (1.0-K)/(1.0+K)*g__.r + g_.r*(8.0*K-a.r*h)/(4.0*K+4.0) + g_.i*(bi*h*K/(double)N+a.i*h)/(4.0*K+4.0);
		g.i = (1.0-K)/(1.0+K)*g__.i + g_.i*(8.0*K-a.r*h)/(4.0*K+4.0) - (bi*h*K/(double)N+a.i*h)/(4.0*K+4.0)*g_.r;
		g__.r = g_.r;
		g__.i = g_.i;
		g_.r = g.r;
		g_.i = g.i;
	}
	return sqrt(g.r*g.r+g.i*g.i);
}


void usage(char **s) {
	fprintf(stderr,"usage %s: omega.r_min omega.r_max omega.i_min omega.i_max kmax NR NI\n",s[0]);
	exit(1);
}

int main (int args, char **argv) {
	int i,j,lc;
	int N; 				/* discretisation number */
	int NR, NI;
	double k;
	double omegar, omegai, wrmin, wrmax, wimin, wimax;
	double val, minval,kmax;


	double x,y;

	if (args!=8)	usage(argv);
	wrmin = atof(argv[1]);
	wrmax = atof(argv[2]);
	wimin = atof(argv[3]);
	wimax = atof(argv[4]);
	kmax = atof(argv[5]);
	NR = atoi(argv[6]);
	NI = atoi(argv[7]);


	N = 5096;

	fprintf(stderr,"k[0 %f]\n",kmax);
	fprintf(stderr,"omegar=[%f %f]x%d\n",wrmin,wrmax,NR);
	fprintf(stderr,"omegai=[%f %f]x%d\n",wimin,wimax,NI);

	for (i=0; i<NI; i++) {
		omegai = scale(0,NI,wimin,wimax,i);
		for (j=0; j<NR; j++) {
			omegar = scale(0,NR,wrmin,wrmax,j);
			minval = shoot2(0.0,omegar,omegai,N);
			for (k=0.0; k<kmax; k+=kmax/1024.0) {
				val  = shoot2(k,omegar,omegai,N);	
				if (val<minval) minval = val; 
			}
			printf("%g\n",minval);
		}
	}

	
	return 0;
}
