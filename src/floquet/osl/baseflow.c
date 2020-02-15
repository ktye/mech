/* $Id: baseflow.c,v 1.7 2005/05/25 15:27:06 elmar Exp $ */
#include <math.h>

#include "define.h"

extern int	 bflag;
extern double	 L;
extern double	 T;
extern double 	 BASEFLOW[YSIZE];
extern double 	 BASEFLOW2[YSIZE];
double baseflow(double Y, int d) {
	double	 eta;
	//y = L*(Y + 1.0);
	eta = L*(Y+1.0)/M_SQRT2;
	if (d==0) {
		if (bflag == 1)		/* Oscillating stokes layer (OSL) */
			return -0.5*exp(-eta)*cos(M_PI*T-eta);
		else if (bflag == 2)	/* Poiseuille flow */
			return 1.0 - Y*Y;
		else			/* pressure driven OSL */
			return 0.5*(cos(M_PI*T) - exp(-eta)*cos(M_PI*T-eta));
	} else if (d==2) {
		if (bflag == 2)		/* Poiseuille flow */
			return -2.0;
		else			/* OSL */
			return 0.5*L*L*exp(-eta)*sin(M_PI*T-eta);
	} else
		exit(1);
	return 0.0;
}

void baseflow_init(void) {
	int	 i;
	double	 y;
	for (i=0; i<YSIZE; i++) {
		y = SETY(i);
		BASEFLOW[i] = baseflow(y, 0);
		BASEFLOW2[i] = baseflow(y, 2);
	}
}
