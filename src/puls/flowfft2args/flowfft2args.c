/* $Id: flowfft2args.c,v 1.2 2003/10/02 18:13:44 elmar Exp $ */

/*
 * flowfft2args - convert flow fourier components to ./puls arguments (alpha R[n])
 *
 * argv[1] - Womersley parameter alpha
 * argv[2] - mean component
 * argv[3] - mean.imag part (discarded) should be zero
 * argv[n] - real and imag parts of oscillating components
 *
 * out:
 * alpha
 * R0
 * R1.r
 * R1.i
 * [Rn.r
 *  Rn.i]
 */

/*
 * Complex:
 * R[n] = i*n*sqrt(i*n) * alpha^5 * C[n]/(1-J_1(sqrt(i*n)*alpha)/J_0(sqrt(i*n)*alpha))
 */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "amoswrap.h"
#include "complex.h"


void R(double alpha, int n, double Creal, double Cimag, double *real, double *imag) {
	double J0r=0.0, J0i=0.0, J1r=0.0, J1i=0.0;
	double JJr=0.0, JJi=0.0;
	double aCr=0.0, aCi=0.0;
	double A=0.0, B=0.0;
	double N;


	N=(double)n;

	J1(sqrt(N)*M_SQRT1_2*alpha,sqrt(N)*M_SQRT1_2*alpha,&J1r,&J1i);
	J0(sqrt(N)*M_SQRT1_2*alpha,sqrt(N)*M_SQRT1_2*alpha,&J0r,&J0i);
	cdiv(J1r,J1i,J0r,J0i,&JJr,&JJi);

	cmult(M_SQRT1_2*sqrt(N)/(N*N*pow(alpha,5)),M_SQRT1_2*sqrt(N)/(N*N*pow(alpha,5)),JJr,JJi,&A,&B);
	B -= 0.5/(N*pow(alpha,4));


	*real = (A*Creal - B*Cimag)/(A*A+B*B)/(2.0*M_PI);
	*imag = - (B*Creal + A*Cimag)/(A*A+B*B)/(2.0*M_PI);
}

void usage(char **argv) {
	printf("usage %s alpha meanflow 0 [Cn.real Cn.imag] ...\n",argv[0]);
	exit(1);
}
int main (int args, char **argv) {
	int i;
	double alpha;
	double R0;
	double Cr, Ci;
	double R1=0.0, R2=0.0;

	if ((args-4)%2 != 0)	usage(argv);

	for (i=1; i<args; ++i) {
		if (i==1) {
			alpha=atof(argv[i]);
			printf("%f\n",alpha);
		} else if (i==2) {
			R0=(8.0/M_PI)*alpha*alpha*atof(argv[i]);
			printf("%f\n",R0);
		}
	}
	for (i=4; i<args; i+=2) {
		Cr=atof(argv[i]);
		Ci=atof(argv[i+1]);
		R(alpha, (i-2)/2, Cr, Ci, &R1, &R2);
		printf("%f\n%f\n",R1,R2);
	}

	return 0;
}
