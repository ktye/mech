/*
 * Operations on complex numbers
 */

#include <math.h>
/* multiplication */
void 
cmult(double xr, double xi, double yr, double yi, double *zr, double *zi) {
	*zr = xr*yr - xi*yi;
	*zi = xi*yr + xr*yi;
}

/* division */
void 
cdiv(double xr, double xi, double yr, double yi, double *zr, double *zi) {
	*zr = (xr*yr + xi*yi)/(yr*yr+yi*yi);
	*zi = (xi*yr - xr*yi)/(yr*yr+yi*yi);
}

/* Cartesian to Polar */
void 
cxyrp(double xr, double xi, double *r, double *p) {
	*r = sqrt(xi*xi + xr*xr);
	*p = atan2(xi,xr);
}

/* Polar to Cartesian */
void 
crpxy(double r, double p, double *xr, double *xi) {
	*xr = r * cos(p);
	*xi = r * sin(p);
}
