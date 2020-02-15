/* $Id: main.c,v 1.1 2004/10/28 16:56:20 elmar Exp $ */

/*
 * colbin - column data to binary translation
 *
 * little endian version
 *
 * converts one column of 2d-tabled data to binary stream
 */

static char id[]="$Id: main.c,v 1.1 2004/10/28 16:56:20 elmar Exp $";
#include <stdlib.h>
#include <stdio.h>

char *progname;
void usage(void) {
	fprintf(stderr,"usage: %s xdim x0 xe xlabel ydim y0 ye ylabel zlabel col maxcols\n",progname);
	exit(1);
}
int main (int args, char **argv) {
	int i,j,c;
	int xdim, ydim, col, cols;
	char *xlabel, *ylabel, *zlabel;
	double x0, xe, y0, ye, z;

	progname = argv[0];
	if (args != 12) usage();
	xdim = atoi(argv[1]);
	ydim = atoi(argv[5]);
	xlabel=argv[4];
	ylabel=argv[8];
	zlabel=argv[9];
	x0=atof(argv[2]);
	xe=atof(argv[3]);
	y0=atof(argv[6]);
	ye=atof(argv[7]);
	col = atoi(argv[10]);
	cols = atoi(argv[11]);
	
	printf("little-endian double 2d %d %d\n",xdim,ydim);
	printf("%s %g %g %s %g %g %s\n\n\n\n\n\n\n\n.\n",xlabel,x0,xe,ylabel,y0,ye,zlabel);
	for (i=0; i<xdim; i++) {
		for (j=0; j<ydim; j++) {
			for (c=0; c<cols; c++) {
				scanf("%lg",&z);
				if (c==col) {
					fwrite(&z, sizeof(double), 1, stdout);
				}
			}
		}
	}


	return 0;
}
