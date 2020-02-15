#include <stdio.h>
#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <ginac/ginac.h>
using namespace std;
using namespace GiNaC;

//#include "poly.h"

ex scale(ex x0, ex x1, ex y0, ex y1, ex scalex) {
	return scalex*(y0-y1)/(x0-x1) + (x0*y1-x1*y0)/(x0-x1);
}

void usage(char **argv) {
	cerr << "usage: " << argv[0] << " alpha R TIME N(tayler-expansion)" << endl;
	exit(1);
}
int main(int args, char **argv)
{
	symbol x("x"), y("y");
	unsigned int N = 5;	/* order of expansion system */
	int steps=100;
	double xmin=0.0,xmax=1.0;
	
	char *var;
	if (var = getenv("steps"))      steps=atoi(var);
	if (var = getenv("xmin")) xmin=atof(var);
	if (var = getenv("xmax")) xmax=atof(var);

	if (args != 5)
		usage(argv);

	ex a = 9;	/* Womersley number */
	ex R = 10;	/* Reynolds number */
	ex n = 1;	/* order of excitation ~exp(i*n*t) */
	ex t;	/* time */


	a = atof(argv[1]);
	R = atof(argv[2]);
	t = atof(argv[3]);
	N = atoi(argv[4]);

	matrix A(N,N);
	for (int i=1; i<=N; i++) {
		for (int k=1; k<=N; k++) {
			if (i==k)
				A(i-1,k-1) = 4*i*i - I*n*a*a;
			if (i==k+1)
				A(i-1,k-1) = I*n*a*a;
			if (i+1==k)
				A(i-1,k-1) = -4*i*i;
		}
	}
	matrix b(N,1);
//	for (int i=0; i<N; i++)
//		b(i,0) = R/(a*a);
	b(0,0) = R/(a*a);


	ex l = symbolic_matrix(N,1,"l");
	l = A.solve(ex_to<matrix>(l),b);


	ex w = 0;
	for (int i=0; i<N; i++) {
		w = w + l[i]*pow(x,2*i);
	}
	w = w*(x*x-1);

	cerr << w << endl;
	/*
	ex coeffi;
	numeric co;
	for (int i=0; i<=N+2; i+=2) {
		co = ex_to<numeric>(expand(w).coeff(x,i));
		cout << i << " " << real(co)*real(co)+imag(co)*imag(co) << endl;
	}
	*/

	ex X;
	ex W;
	for (int i=0; i<steps; i++) {
		X = evalf(scale(0,steps-1,xmin,xmax,i));
		W = w*exp(-I*n*t);
		//cout << t << " " << X << " " << real(ex_to<numeric>(w.subs(x==X).evalf())) << endl;
		//cout << t << " " << pqert << " " << real(ex_to<numeric>(w.subs(x==pqert))) << endl;
		cout << t << " " << X << " " << real(ex_to<numeric>(W.subs(x==X))) << endl;
		//cout << t << " " << X << " " << w.subs(x==X) << endl;
	}


	/*
	lst wlst(0);
	int timesteps=10;
	for (int i=0; i<timesteps; i++) {
		t = scale(0,timesteps-1,0,Pi,i);
		wlst.append(w*exp(I*n*t));
	}
	*/
}
