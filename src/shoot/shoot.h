/* $Id: shoot.h,v 1.6 2004/11/04 12:26:04 elmar Exp $ */


double shoot(void);
void single(double wrmin, double wrmax, double wimin, double wimax, int Nw);
void inter(double wrmin, double wrmax, double wimin, double wimax, int Nw, int Nk);
void plot(double rho, int Nk);
void plot_all(double rho, int Nk);
void tmaxima(double t[3], double k[9],double omegar[9],double omegai[9]);
