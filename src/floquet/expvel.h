/* $Id: expvel.h,v 1.2 2004/11/26 19:04:59 elmar Exp $ */
void expvel_init(int n, double k, int N);
double expvel_disrel(int n, double k, double beta);
void expvel(int index, int n, double k, double r, double *u, double *v, double *w);
void expvel_ortho(int N);
