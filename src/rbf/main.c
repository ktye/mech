/* $Id: main.c,v 1.2 2006/02/18 18:18:05 elmar Exp $ */

/*
 * rbf - radial basis functions
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>

#include "rbf.h"
char	*progname;

void usage(void) {
	fprintf(stderr,"usage: %s -d dimensions [-l lengthscale] x0 [y0 ...] f0 [...]\n",progname);
	exit(1);
}

int main(int args, char **argv) {
	int	 i;
	extern char	*optarg;
	extern int	 optind;
	int	 c, d = 1;
	double	 l = 1.0;
	double	*x;
	double  *X, f;
	int	 kflag = 0;

	progname = argv[0];
	while ((c = getopt(args, argv, "d:l:k")) != -1) {
		switch (c) {
		case 'd':
			d = atoi(optarg);
			break;
		case 'l':
			l = atof(optarg);
			break;
		case 'k':
			kflag = 1;
			break;
		default:
			usage();
		}
	}
	args -= optind;
	argv += optind;

	if (d<1) { fprintf(stderr,"error: d<1\n"); return 1; }
	if (args<1) { fprintf(stderr,"error: no points given\n"); return 1; }
	if (args%(d+1)) { fprintf(stderr,"error: wrong number of points\n"); return 1; }

	x = (double *)malloc(args*sizeof(double));
	for (i=0; i<args; i++)
		x[i] = atof(argv[i]);

	printf(":1:2:red\n");
	for (i=0; i<args; i+=2)
		printf("%g 0\n%g %g\n\n",x[i],x[i],x[i+1]);
	printf(":1:2:green\n");

	rbf_init(x, d, args/(d+1), l, kflag);

	X = (double *)malloc(d*sizeof(double));
	for (;;) {
		for (i=0; i<d; i++)
			if (scanf("%lf",&X[i]) != 1)
				return 0;
		f = rbf(X);
		//f = sph_approx(x, X, d, args/(d+1), l);
		for (i=0; i<d; i++)
			printf("%g ", X[i]);
		printf("%g\n",f);
	}



	return 0;
}
