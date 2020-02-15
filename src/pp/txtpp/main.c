/* $Id: main.c,v 1.1 2004/10/01 08:55:28 elmar Exp $ */

/*
 * main - ascii table to pp-data converter
 */

static char id[]="$Id: main.c,v 1.1 2004/10/01 08:55:28 elmar Exp $";
#include <stdio.h>
#include <stdlib.h>



int main (int args, char **argv) {
	int i,j;
	int x = 1;
	int xdim,ydim,rows,row;
	double xmin,ymin,xmax,ymax;
	char *xstr,*ystr;
	double z;

	if (args != 11) {
		fprintf(stderr,"usage: %s xdim ydim xstr xmin xmax ystr ymin ymax rows row\n",argv[0]);
		exit(0);
	}
	xdim=atoi(argv[1]);
	ydim=atoi(argv[2]);
	xstr=argv[3];
	xmin=atof(argv[4]);
	xmax=atof(argv[5]);
	ystr=argv[6];
	ymin=atof(argv[7]);
	ymax=atof(argv[8]);
	rows=atoi(argv[9]);
	row=atoi(argv[10]);

	if (*(char *)&x == 1) printf("little-endian ");
	else	printf("big-endian ");
	printf("double 2d %d %d\n",ydim,xdim);
	printf("%s %g %g %s %g %g row%d\n",xstr,xmin,xmax,ystr,ymin,ymax,row);
	printf("\n\n\n\n\n\n\n.\n");

	z = 0.0;
	for (j=0; j<xdim*ydim; j++) {
		for (i=0; i<rows; i++) {
			scanf("%lg",&z);
			if (i==row-1) {
				fwrite(&z, sizeof(double), 1, stdout);
			}
		}
	}

	return 0;
}
