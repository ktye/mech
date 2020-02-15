/* $Id: sym.c,v 1.1 2006/02/17 18:29:48 elmar Exp $ */

/*
 * sym - test symsolve
 */


#include <stdio.h>
#include <stdlib.h>
#include "symsolve.h"
int nomain(int args, char **argv) {
	int	 i, n, elements;
	double	*A, *x;

	args--;
	argv++;

	printf("args=%d\n",args);
	for (n=1; n<args; n++)
		if (args == (n*(n+1))/2 + n)
			break;
	if (n==args)
		return 1;

	printf("n=%d\n",n);


	elements = (n*(n+1))/2;
	A = (double *)malloc(elements*sizeof(double));
	x = (double *)malloc(n*sizeof(double));
	
	printf("A = ");
	for (i=0; i<elements; i++) {
		A[i] = atof(argv[i]);
		printf("%g ",A[i]);
	}
	printf("\nrhs = ");
	for (i=elements; i<args; i++) {
		x[i-elements] = atof(argv[i]);
		printf("%g ",x[i-elements]);
	}
	printf("\n");

	symsolve(A, x, n);
	printf("x = ");
	for (i=0; i<n; i++)
		printf("%g ", x[i]);
	printf("\n");

	
	return 0;
}
