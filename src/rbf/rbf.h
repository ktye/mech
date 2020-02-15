/* $Id: rbf.h,v 1.2 2006/02/18 18:18:05 elmar Exp $ */
void rbf_init(double *x, int dim, int N, double l, int kern);
double rbf(double *x);
double sph_approx(double *x, double *X, int dim, int N, double l);
