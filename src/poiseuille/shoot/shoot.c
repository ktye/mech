/* $Id: shoot.c,v 1.4 2004/02/20 19:52:36 elmar Exp $ */

/*
 * shoot - shooting algorithm
 */

static char id[]="$Id: shoot.c,v 1.4 2004/02/20 19:52:36 elmar Exp $";
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
double shoot(double R, double kz, double omegar, double omegai, int N) {
	int k;
	double K, h; 
	double r, k_n;
	cdouble f, f_, f__; 
	cdouble a;
	double bi;
	h = 1.0/(double)N;

	a.r = -kz*kz-omegai*R;
	a.i = omegar*R-kz*R;
	bi = kz*R;

	f__.r = 0.0;
	f__.i = 0.0;
	f_.r  = h;
	f_.i  = 0.0;
	//printf("%f %f %f\n",0.0,0.0,0.0);
	//printf("%f %f %f\n",h,h,0.0);
	for (k=1; k<N; k++) {
		K = (double)k;
		k_n = K/(double)N;
		f.r = -f__.r*(2.0*K-1.0)/(2.0*K+1.0) + f_.r*(1.0+2.0*K*K-a.r*k_n*k_n)/(K*K+K/2.0);
		f.r += f_.i*(a.i*k_n*k_n+bi*pow(k_n,4))/(K*K+K/2.0);
		f.i = -f__.i*(2.0*K-1.0)/(2.0*K+1.0) + f_.i*(1.0+2.0*K*K-a.r*k_n*k_n)/(K*K+K/2.0);
		f.i += f_.r*(-a.i*k_n*k_n-bi*pow(k_n,4))/(K*K+K/2.0);
		r = (K+1.0)*h;
	//	printf("%f %f %f\n",r,f.r,f.i);
		f__.r = f_.r;
		f__.i = f_.i;
		f_.r = f.r;
		f_.i = f.i;
	}
	//fprintf(stderr,"%f\n",sqrt(f.r*f.r+f.i*f.i));
	return sqrt(f.r*f.r+f.i*f.i);
}

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
//		r = sqrt((K+1.0)*h);
//		printf("%f %f %f\n",r,r*g.r,r*g.i);
		g__.r = g_.r;
		g__.i = g_.i;
		g_.r = g.r;
		g_.i = g.i;
	}
	//fprintf(stderr,"%f\n",sqrt(f.r*f.r+f.i*f.i));
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

/* get initial values for omega */
void startomega(double R, double kz, double *omegar, double *omegai) {
	*omegar = 0.0;
	*omegai = -3.83171*3.83171/R;
}

void usage(char **s) {
	fprintf(stderr,"usage %s: R k [omega.r omega.i]\n",s[0]);
	exit(1);
}

int main (int args, char **argv) {
	int i,j,lc;
	int N; 				/* discretisation number */
	double omegar, omegai, kz, k, dk, R;
	double bestr, besti, bestval, val;
	double range, shots;
	double maxerr;


	double x,y;

	if ((args!=3)&&(args!=5))	usage(argv);
	R = atof(argv[1]);
	kz = atof(argv[2]);
	if (args==5) {
		omegar=atof(argv[3]);
		omegai=atof(argv[4]);
	} else {
		startomega(R, kz, &omegar, &omegai);
	}


	srand48(time(0));
	N = 2048;
	range = 0.02/R;
	shots = 20;
	maxerr = 0.001;
	dk = 0.1/R;

	if (args==5) {
		val  = shoot2(R,kz,omegar,omegai,N);	
		printf("%f\n",val);
		return 0;
	}

	bestr = omegar;
	besti = omegai;
	k = 0.0;
	bestval = shoot2(R,k,bestr,besti,N);
	while (k<kz) {
		lc = 0;
		while(bestval>maxerr) {	
			aim(bestr,besti,&omegar,&omegai,range);
			for (j=0; j<shots; j++) {
				val  = shoot2(R,k,omegar,omegai,N);	
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
		printf("%f %f %f %f\n",R, k, bestr,besti);
		k += dk;
		bestval=maxerr+1.0;
	}
	return 0;
}
