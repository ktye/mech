/* $Id: rk4.h,v 1.3 2004/09/10 09:37:56 elmar Exp $ */

extern unsigned int Lsteps;
#define SSTEP(i)	((double)(i)/(double)(Lsteps-1))

void ode_printvec(double *y, int n, double x);

/* integration from 0 to 1 in Nsteps steps */
void rk4(double *y, int ny, double (*f)(), double *ymax);
