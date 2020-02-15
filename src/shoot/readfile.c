/* $Id: readfile.c,v 1.5 2004/09/30 08:32:16 elmar Exp $ */

/*
 * readbase - read baseflow binary file
 *
 */

#include <stdio.h>
#include <stdlib.h>
#include <libgen.h>
#include <math.h>
#include "endian.h"

extern int ALPHANUM;
double *readfile(const char *path, double *alpha, unsigned int *N) {
	int i;
	FILE *fp;
	char dchar[3];
	double x;
	unsigned int nX;
	int log2steps, endian;
	int fileok;
	char Endian;
	int swapendian = 0;
	double *w;
	char *file;

	file = basename(path);
	fileok = 0;
	if (strlen(file) == 12) {
		if (!strncmp("base-",file,5)) {
			dchar[0] = file[5];
			dchar[1] = file[6];
			dchar[2] = '\0';
			ALPHANUM = atoi(dchar);
			dchar[0] = file[8];
			dchar[1] = file[9];
			dchar[2] = '\0';
			log2steps = atoi(dchar);
			Endian = file[11];
			if (Endian == 'L') endian = 0;
			else if (Endian == 'B') endian = 1;
			else { fprintf(stderr,"endian error\n"); exit(1); }
			fileok = 1;
		}
	}
	if (!fileok) { fprintf(stderr, "filename cannot be scanned\n"); exit(1); }

	if ((log2steps<1)||(log2steps>15)) { fprintf(stderr,"error: log2steps not int [1,15]\n"); exit(1); }
	nX = (1<<log2steps)+1;


	if (endian != mendian()) swapendian = 1;

	w = (double *)malloc((size_t)(4*nX)*sizeof(double));
	if (!w) { perror("w memory allocation"); exit(1); }

	fp = fopen(path,"r");
	if (!fp) { perror("main fopen"); exit(1); }
	if (fread(alpha, sizeof(double), 1, fp) != 1) { perror("main fread"); exit(1); }
	if (swapendian) dswap(alpha);
	for (i=0; i<4*nX; i++) {
		if (fread(&w[i], sizeof(double), 1, fp) != 1) { perror("main fread"); exit(1); }
		if (swapendian) dswap(&w[i]);
	}
	if (fread(&x, sizeof(double), 1, fp) == 1) { fprintf(stderr, "error: file too long\n"); exit(1); }
	fclose(fp);

	*N = nX;
	return w;
}
