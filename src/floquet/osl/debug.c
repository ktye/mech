/* $Id: debug.c,v 1.1 2005/04/20 15:07:29 elmar Exp $ */
#include <stdio.h>
#include <string.h>

#include "define.h"
#include "base.h"
#include "baseflow.h"

extern double	 BASEFLOW[YSIZE], BASEFLOW2[YSIZE];
extern double	 PHI[YSIZE][NMAX], PSI[YSIZE][NMAX];
extern double	 LAMBDA[NMAX];
extern int	 N;

void debug(char *c) {
	int	 i, n;

	baseflow_init();
	base_init();
	if (!strcmp(c,"baseflow")) {
		for (i=0; i<YSIZE; i++) {
			printf("%g %g %g\n",SETY(i), BASEFLOW[i], BASEFLOW2[i]);
		}
	} else if (!strcmp(c,"phi")) {
		for (n=0; n<N; n++) {
			for (i=0; i<YSIZE; i++) {
				printf("%g %g\n",SETY(i), PHI[i][n]);
			}
			printf("\n");
		}
	} else if (!strcmp(c,"psi")) {
		for (n=0; n<N; n++) {
			for (i=0; i<YSIZE; i++) {
				printf("%g %g\n",SETY(i), PSI[i][n]);
			}
			printf("\n");
		}
	} else {
		printf("debug: unkown command\n");
	}
	exit(0);
}
