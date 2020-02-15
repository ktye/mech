/* $Id: membrane.c,v 1.1 2004/04/21 13:53:05 elmar Exp $ */

/*
 * membrane - 1,2 mode circular mem.
 */

static char id[]="$Id: membrane.c,v 1.1 2004/04/21 13:53:05 elmar Exp $";
#include <stdio.h>
#include <math.h>

#define ROWS	600
#define COLS	600
#define	FRAMES	20

int mendian(void) {
	int x = 1;
	if(*(char *)&x == 1) return 0; // little-endian
	else return 1; // big-endian
}

write_header(FILE *file, double time) {
	if (mendian()) 
		fprintf(file,"big-endian double 2d %d %d\n",ROWS,COLS);
	else 
		fprintf(file,"little-endian double 2d %d %d\n",ROWS,COLS);
	fprintf(file,"x 0 1 y 0 1 z\n\n\n\n\n\n\n\n\n");
}

int main (int args, char **argv) {
	int i,j,f;
	double x,y,z,phi;
	char filename[10];
	FILE *fp;

	for (f=0; f<FRAMES; f++) {
		snprintf(filename,sizeof(filename),"out.%.3d",f);
		fp = fopen(filename,"w");
		if (fp == NULL) {
			perror(NULL);
			exit(1);
		}
		write_header(fp, 0.0);

		for (j=0; j<ROWS; j++) {
			y = (double)j/(double)(ROWS-1);
			for (i=0; i<COLS; i++) {
				x = (double)i/(double)(COLS-1);
				phi = atan2(y-0.5,x-0.5);
				// kleine Drehung mit der Periode
				z = j1(2.0*7.01558*sqrt((x-0.5)*(x-0.5)+(y-0.5)*(y-0.5)))*cos(2.0*M_PI*f/(double)(FRAMES-1))*sin(phi-2.0*M_PI*f/(double)(FRAMES-1));
				if (((y-0.5)*(y-0.5) + (x-0.5)*(x-0.5)) > 0.25)
					z = 0.0;
				fwrite(&z, sizeof(double), 1, fp);
			}
		}
		fclose(fp);
	}

	return 0;
}
