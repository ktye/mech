/* $Id: quad.c,v 1.2 2005/04/18 17:11:51 elmar Exp $ */
#include <stdio.h>
#include "define.h"

/* simpson_integration on interval [-1,1] with YSIZE steps */
/* array version, YSIZE must be odd! */
double quad(double *f) {
	int c, k;
	double sum, h;

	c = 2;
	h = 2.0/(double)(YSIZE-1);
	sum = f[0];
	for (k=1; k<YSIZE-1; k++) {
		c = 6-c;
		sum += (double)c*f[k];
	}
	return (sum + f[YSIZE-1])*h/3.0;
}
