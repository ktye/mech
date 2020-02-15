/* $Id: hooke.c,v 1.4 2004/10/29 15:01:24 elmar Exp $ */

/* Algorithm adapted from:
 * http://www.netlib.org/opt/hooke.c
 * see URL for description of parameters (rho) */


#include <stdio.h>
#include <math.h>

#define      VARS	2

extern double shoot(void);
extern double OMEGAR;
extern double OMEGAI;
extern double K;

/*
 *	double f(double x[VARS], int n) {
 *		double y;
 *		OMEGAR=x[0];
 *		OMEGAI=x[1];
 *		return shoot();
 *	}
 */

double best_nearby(double delta[VARS], double point[VARS], double prevbest, int nvars, double (*f)()) {
	   double	   z[VARS];
	   double	   minf, ftmp;
	   int		   i;
	   minf = prevbest;
	   for (i = 0; i < nvars; i++)
		   z[i] = point[i];
	   for (i = 0; i < nvars; i++) {
		   z[i] = point[i] + delta[i];
		   ftmp = f(z, nvars);
		   if (ftmp < minf)
			   minf = ftmp;
		   else {
			   delta[i] = 0.0 - delta[i];
			   z[i] = point[i] + delta[i];
			   ftmp = f(z, nvars);
			   if (ftmp < minf)
				   minf = ftmp;
			   else
				   z[i] = point[i];
		   }
	   }
	   for (i = 0; i < nvars; i++)
		   point[i] = z[i];
	   return (minf);
}

int hooke(int nvars, double startpt[VARS], double endpt[VARS], double rho, double epsilon, int itermax, double (*f)()) {
	   double	   delta[VARS];
	   double	   newf, fbefore, steplength, tmp;
	   double	   xbefore[VARS], newx[VARS];
	   int		   i, keep;
	   int		   iters, iadj;
	   for (i = 0; i < nvars; i++) {
		   newx[i] = xbefore[i] = startpt[i];
		   delta[i] = fabs(startpt[i] * rho);
		   if (delta[i] == 0.0)
			   delta[i] = rho;
	   }
	   iadj = 0;
	   steplength = rho;
	   iters = 0;
	   fbefore = f(newx, nvars);
	   newf = fbefore;
	   while ((iters < itermax) && (steplength > epsilon)) {
		   iters++;
		   iadj++;
		   for (i = 0; i < nvars; i++) {
			   newx[i] = xbefore[i];
		   }
		   newf = best_nearby(delta, newx, fbefore, nvars, f);
		   keep = 1;
		   while ((newf < fbefore) && (keep == 1)) {
			   iadj = 0;
			   for (i = 0; i < nvars; i++) {
				   if (newx[i] <= xbefore[i])
					   delta[i] = 0.0 - fabs(delta[i]);
				   else
					   delta[i] = fabs(delta[i]);
				   tmp = xbefore[i];
				   xbefore[i] = newx[i];
				   newx[i] = newx[i] + newx[i] - tmp;
			   }
			   fbefore = newf;
			   newf = best_nearby(delta, newx, fbefore, nvars, f);
			   if (newf >= fbefore)
				   break;
			   keep = 0;
			   for (i = 0; i < nvars; i++) {
				   keep = 1;
				   if (fabs(newx[i] - xbefore[i]) >
				       (0.5 * fabs(delta[i])))
					   break;
				   else
					   keep = 0;
			   }
		   }
		   if ((steplength >= epsilon) && (newf >= fbefore)) {
			   steplength = steplength * rho;
			   for (i = 0; i < nvars; i++) {
				   delta[i] *= rho;
			   }
		   }
	   }
	   for (i = 0; i < nvars; i++)
		   endpt[i] = xbefore[i];
	   return (iters);
}
