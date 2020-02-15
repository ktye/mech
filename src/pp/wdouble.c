/* $Id: wdouble.c,v 1.2 2004/06/10 17:08:09 elmar Exp $ */

/*
 * wdouble - write double
 */

static char id[]="$Id: wdouble.c,v 1.2 2004/06/10 17:08:09 elmar Exp $";
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

void usage(char **argv) {
	fprintf(stderr,"usage: %s \n", argv[0]);
	exit(1);
}

int main (int args, char **argv) {
	int i,j,MAX;
	double x,y,f;
	MAX=100;


	x=0.0;
	for (j=0; j<MAX; j++) {
		x = (double)j*2.0*M_PI/(double)MAX;
		for (i=0; i<MAX; i++) { 
			y = (double)i*4.0*M_PI/(double)MAX;
			f = sin(x*y);
			fwrite(&f,sizeof(double),1,stdout);
		}
	}
	return 0;
}


