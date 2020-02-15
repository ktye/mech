/* $Id: shootxx.cc,v 1.1 2004/04/21 13:53:04 elmar Exp $ */

/*
 * shootxx.cc - c++ version of shoot
 */

static char id[]="$Id: shootxx.cc,v 1.1 2004/04/21 13:53:04 elmar Exp $";
#include <iostream>
#include <complex>
#include <math.h>
using namespace std;


#define MAXLOOP		4096
double scale (double a, double b, double A, double B, double x) {
	return (x*(A-B)/(a-b)+(a*B-b*A)/(a-b));
}


/* transformed ode f = r*g(s), s=r^2 	*/
/* 4sg'' + 8g' +(a+bs)g = 0 		*/
/* g(0) = 1, g'(0) = -a/12		*/
double shoot2(double R, double kz, complex <double> omega, int N) {
	int k;
	double K, h; 
	double r, k_n;
	complex<double> g; 
	complex<double> a(-kz*kz-R*omega.imag(),R*omega.real()-R*kz);
	complex<double> b(0,kz*R);
	h = 1.0/(double)N;

	complex<double> g__(1.0,0.0);
	complex<double> g_(1.0,0.0);
	g_ -= a*h/8.0;

	for (k=1; k<N; k++) {
		K = (double)k;
		g = (1.0-K)/(1.0+K)*g__ + (2.0*K/(1.0+K)-h*(a+b*K/(double)N)/(4.0*K+4.0))*g_;
		g__ = g_;
		g_ = g;

	}
	return abs(g);
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

void usage(char **s) {
	fprintf(stderr,"usage %s: R k [omega.r omega.i]\n",s[0]);
	exit(1);
}

int main (int args, char **argv) {
	int i,j,lc;
	int N; 				/* discretisation number */
	double kz, k, dk, R;
	double bestval, val;
	complex <double> best;
	double range, shots;
	double maxerr;
	double omegar,omegai;


	double x,y;

	if ((args!=3))	usage(argv);
	R = atof(argv[1]);
	kz = atof(argv[2]);
	complex <double> omega(0.0, -3.83171*3.83171/R);


	srand48(time(0));
	N = 2048;
	range = 0.02/R;
	shots = 20;
	maxerr = 0.001;
	dk = 0.1/sqrt(R);

	if (args==5) {
		val  = shoot2(R,kz,omega,N);	
		printf("%f\n",val);
		return 0;
	}

	best = omega;
	k = 0.0;
	bestval = shoot2(R,k,best,N);
	while (k<kz) {
		lc = 0;
		while(bestval>maxerr) {	
			aim(best.real(),best.imag(),&omegar, &omegai,range);
			omega = complex<double>(omegar,omegai);
			for (j=0; j<shots; j++) {
				val  = shoot2(R,k,omega,N);	
				if (val<bestval) {
					bestval = val;
					best = omega;
				}
			}
			lc++;
			if (lc>MAXLOOP) {
				fprintf(stderr,"break @ k=%f\n",k);
				return 1;
			}
		}
		fprintf(stderr,"k=%f omega=[%f %f] val=%f\n",k,best.real(),best.imag(),bestval);
		printf("%f %f\n",best.real(),best.imag());
		k += dk;
		bestval=maxerr+1.0;
	}
	return 0;
}
