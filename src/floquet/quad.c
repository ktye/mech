/* $Id: quad.c,v 1.2 2004/11/24 15:01:34 elmar Exp $ */
#include <stdio.h>
#include "define.h"

static char id[]="$Id: quad.c,v 1.2 2004/11/24 15:01:34 elmar Exp $";


/* simpson_integration on interval [0,1] with RSIZE steps */
/* array version, RSIZE must be odd! */
double quad(double *f) {
	int c, k;
	double sum, h;

	/*
	if (!(N%2)) {
		fprintf(stderr,"error quad_simpson: N must be odd\n");
		exit(1);
	}
	*/
	c = 2;
	h = 1.0/(double)(RSIZE-1);
	sum = f[0];
	for (k=1; k<RSIZE-1; k++) {
		c = 6-c;
		sum += (double)c*f[k];
	}
	return (sum + f[RSIZE-1])*h/3.0;
}
