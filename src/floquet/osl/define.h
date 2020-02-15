/* $Id: define.h,v 1.2 2005/04/18 17:11:51 elmar Exp $ */
#define NMAX 100
#define YSIZE 1025		/* YSIZE must be ODD */

#define SETY(y)	(-1.0 + 2.0*(double)(y)/(double)(YSIZE-1))
#define SETCNT(x,i) if (x->n==1) x->val = x->min; else x->val = (x->min + (double)i*(x->max - x->min)/(double)(x->n-1))

typedef struct {
        double  val;
        double  min;
        double  max;
        int     n;
        int     ndef;
} diterator;

