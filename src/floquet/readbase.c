/* $Id: readbase.c,v 1.6 2005/05/26 08:52:27 elmar Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "define.h"
static char id[]="$Id: readbase.c,v 1.6 2005/05/26 08:52:27 elmar Exp $";

extern struct baseflowstruct {
        double     r[RSIZE];
        double    dr[RSIZE];
        double  real[RSIZE];
        double dreal[RSIZE];
        double  imag[RSIZE];
        double dimag[RSIZE];
} baseflow;

extern double 	TIME;

extern int	Bflag;
extern char 	*BASEPREFIX;

void readbase(int alphanum) {
	FILE *fp;
	int i;
	double r,wr,wi,w1r,w1i;
	char n[4];
	char filename[1024];

	snprintf((char *)&n,4,"%.3d",alphanum);
	(void)strncpy(filename, BASEPREFIX, sizeof(filename)-1);
	filename[sizeof(filename)-1] = '\0';
	(void)strncat(filename, (char *)&n, sizeof(filename)-1-strlen(filename));
	/*
	n[0] = filename[strlen(filename)-3];
	n[1] = filename[strlen(filename)-2];
	n[2] = filename[strlen(filename)-1];
	n[3] = '\0';
	*/

	fp = fopen(filename,"r");
	if (!fp) {
		fprintf(stderr,"error in readbase nsff: %s\n",filename);
		exit(1);
	}
	for (i=0; i<RSIZE; i++) {
		if (fscanf(fp, "%lf %lf %lf %lf %lf",&r,&wr,&wi,&w1r,&w1i) !=5) {
			fprintf(stderr,"readbase: error at line %d in %s\n",i,filename);
			exit(1);
		}
		baseflow.real[i]  = wr;
		baseflow.imag[i]  = wi;
		baseflow.dreal[i] = w1r;
		baseflow.dimag[i] = w1i;
	}
	if (fscanf(fp,"%lf",&r) > 0) {
		fprintf(stderr,"readbase: file too long: %s\n",filename);
		exit(1);
	}
	fclose(fp);
}

void readbase_settime(void) {
	int i;
	for (i=0; i<RSIZE; i++) {
		baseflow.r[i]  = cos(M_PI*TIME)*baseflow.real[i]  + sin(M_PI*TIME)*baseflow.imag[i];
		baseflow.dr[i] = cos(M_PI*TIME)*baseflow.dreal[i] + sin(M_PI*TIME)*baseflow.dimag[i];
		if (Bflag) {
			printf("baseflow[%d] %.16f %.16f\n",i,baseflow.r[i],baseflow.dr[i]);
		}
	}
}
