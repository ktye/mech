/* $Id: define.h,v 1.7 2005/01/21 10:37:18 elmar Exp $ */

/* filenames for data files with zeros of Bessel functions Jn */
#define ZEROJNFILE	"/home/elmar/dat/zerojn/zerojn-%.3d-100"

/* filename prefix for baseflow data */
#define DEFBASEFILE	"/home/elmar/dat/baseflow/baseflow-"

/* Steps of radial integration; length of expansion velocities */
#define RSIZE 1025 /* RSIZE must be odd, due to quad.c integration algorithm */

/* Maximal number of expansion functions */
#define ESIZE 128

/* Number of Time steps for Floquet time integration
#define TSIZE 1025
*/

/* set value in loop */
#define SETCNT(x,i) if (x->n==1) x->val = x->min; else x->val = (x->min + (double)i*(x->max - x->min)/(double)(x->n-1))

typedef struct {
	double  val;
	double  min;
	double  max;
	int     n;
	int     ndef;
} diterator;

typedef struct {
	int  val;
	int  min;
	int  max;
	int  n;
	int  ndef;
} iiterator;

typedef struct {
	double 	R;
	int	a;
	int 	n;
	double 	k;
	double	t;
	int	N;
	double 	sigma;
	double 	omega;
} tuple;
