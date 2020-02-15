/* $Id: shootk.c,v 1.2 2004/02/20 19:52:38 elmar Exp $ */

/*
 * shootk - shooting algorithm konst k lines
 *
 * usage %s: k omega.r omega.i
 * (omega.r, omega.i) = values for R->infty @ konst k
 *
 * output is scaled by Reynolds number [R*omegar, R*omegai]
 */

static char id[]="$Id: shootk.c,v 1.2 2004/02/20 19:52:38 elmar Exp $";
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

/* f'' + 1/r*f' + (a+b*r^2-1/r^2)*f = 0			*/
/* f(0) = 0, f'(0) = 1					*/
/* solve initial value problem, return end value @ r=1 	*/

/* transformed ode f = r*g(s), s=r^2 	*/
/* 4sg'' + 8g' +(a+bs)g = 0 		*/
/* g(0) = 1, g'(0) = -a/12		*/
double shoot2(double R, double kz, double omegar, double omegai, int N) {
	int k;
	double K, h; 
	double r, k_n;
	cdouble g, g_, g__; 
	cdouble a;
	double bi;
	h = 1.0/(double)N;

	a.r = -kz*kz-omegai*R;
	a.i = omegar*R-kz*R;
	bi = kz*R;

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
	/* alpha=4.6/(rho*rho); */
	*X = x+s;
	s = drand48();
	s = rho*(2.0*s-1.0);
	*Y = y+s;
}

/* get Reynolds number [infty, 0] */
double getReynolds(int i, int N) {
	double R0 = 1.0;
	double Rinfty = 1000.0;
	double rate = 10.0;
	return R0 + (Rinfty-R0)*exp(-rate*(double)i/(double)N);
}

void usage(char **s) {
	fprintf(stderr,"usage %s: k omega.r omega.i\n",s[0]);
	exit(1);
}

int main (int args, char **argv) {
	int i,j,l,lc;
	int N,NR; 				/* discretisation number */
	double omegar, omegai, k, Rk, R;
	double bestr, besti, bestval, val;
	double range, shots;
	double maxerr;


	double x,y;

	if ((args!=4))	usage(argv);

	srand48(time(0));
	N = 2048;
	NR = 2048;
	shots = 20;
	maxerr = 0.001;
	bestval = maxerr + 1.0;

	Rk = atof(argv[1]);
	omegar=atof(argv[2])/getReynolds(0,NR);
	omegai=atof(argv[3])/getReynolds(0,NR);
	fprintf(stderr,"Rinfty: omega=(%f,%f) @ Rk=%f\n",omegar,omegai,Rk);
	
	bestr = omegar;
	besti = omegai;
	for (l=0; l<NR; l++) {
		R = getReynolds(l,NR);
		lc = 0;
		while(bestval>maxerr) {	
			range = 0.02/sqrt(R);
			aim(bestr,besti,&omegar,&omegai,range);
			for (j=0; j<shots; j++) {
				k = Rk/R;
				val  = shoot2(R,k,omegar,omegai,N);	
				if (val<bestval) {
					bestval = val;
					bestr = omegar;
					besti = omegai;
				}
			}
			lc++;
			if (lc>MAXLOOP) {
				fprintf(stderr,"break @ k=%f R=%f, omega=[%f, %f] range=%f\n",k,R,omegar,omegai,range);
				return 1;
			}
		}
		fprintf(stderr,"k=%f R=%f omega=[%f %f] val=%f\n",k,R,bestr,besti,bestval);
		printf("%f %f\n",bestr*R,besti*R);
		bestval=maxerr+1.0;
	}
	return 0;
}
