/* $Id: mkpic.c,v 1.1 2004/04/21 13:53:05 elmar Exp $ */

/*
 * membrane - 1,2 mode circular mem.
 */

static char id[]="$Id: mkpic.c,v 1.1 2004/04/21 13:53:05 elmar Exp $";
#include <stdio.h>
#include <math.h>

#define	FRAMES	2
int ROWS, COLS;

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

	if (args != 2) exit(1);
	ROWS = atoi(argv[1]);
	COLS = ROWS;

	for (f=0; f<FRAMES-1; f++) {
		snprintf(filename,sizeof(filename),"out-%.3d",f);
		fp = fopen(filename,"w");
		if (fp == NULL) {
			perror(NULL);
			exit(1);
		}
		write_header(fp, 0.0);

		z = 0.0;
		for (j=0; j<ROWS; j++) {
			y = (double)j/(double)(ROWS-1);
			for (i=0; i<COLS; i++) {
				z += 1.0;
				fwrite(&z, sizeof(double), 1, fp);
			}
		}
		fclose(fp);
	}

	return 0;
}
