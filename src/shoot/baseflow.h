/* $Id: baseflow.h,v 1.6 2004/09/30 08:32:16 elmar Exp $ */

extern double *BASEFLOW;
extern double *BASEFLOW1;

void initbaseflow(void);
void reinitbaseflow(void);
void printbaseflow(void);
double baseflow(double, int);
double baseflow_s(double, int);
double baseflow_taylor(int);
