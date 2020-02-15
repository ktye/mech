/* $Id: shoothighr.c,v 1.1 2004/02/17 18:45:36 elmar Exp $ */

/*
 * shoothighr - shooting algorithm for high Reynolds limit
 */

static char id[]="$Id: shoothighr.c,v 1.1 2004/02/17 18:45:36 elmar Exp $";
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define MAXLOOP		4096
struct doublecomplex {
	double r;
	double i;
};
typedef struct doublecomplex cdouble;

double scale (double a, double b, double A, double B, double x) {
	return (x*(A-B)/(a-b)+(a*B-b*A)/(a-b));
}

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

/* get random value pair (*X, *Y) near (x,y) within range rho */
void aim(double x, double y, double *X, double *Y, double rho) {
	double alpha,s;
	s = drand48();		/* s [0,1] */
	s = rho*(2.0*s-1.0);	/* s [-rho, rho] */
	*X = x+s;
	s = drand48();
	s = rho*(2.0*s-1.0);
	*Y = y+s;
}

void usage(char **s) {
	fprintf(stderr,"usage %s: kmax \n",s[0]);
	exit(1);
}

int main (int args, char **argv) {
	int i,j,lc;
	int N; 				/* discretisation number */
	double omegar, omegai, kz, k, dk;
	double bestr, besti, bestval, val;
	double range, shots;
	double maxerr;


	double x,y;

	if ((args!=2))	usage(argv);
	kz = atof(argv[1]);

	omegar = 0.0;
	omegai = -3.83171*3.83171;

	srand48(time(0));
	N = 2048;
	range = 0.02;
	shots = 20;
	maxerr = 0.0001;
	dk = 0.1;

	bestr = omegar;
	besti = omegai;
	k = 0.0;
	bestval = shoot2(k,bestr,besti,N);
	while (k<kz) {
		lc = 0;
		while(bestval>maxerr) {	
			aim(bestr,besti,&omegar,&omegai,range);
			for (j=0; j<shots; j++) {
				val  = shoot2(k,omegar,omegai,N);	
				if (val<bestval) {
					bestval = val;
					bestr = omegar;
					besti = omegai;
				}
			}
			lc++;
			if (lc>MAXLOOP) {
				fprintf(stderr,"break @ k=%f\n",k);
				return 1;
			}
		}
		fprintf(stderr,"k=%f omega=[%f %f] val=%f\n",k,bestr,besti,bestval);
		printf("%f %f\n",bestr,besti);
		k += dk;
		bestval=maxerr+1.0;
	}
	return 0;
}
