/* $Id: loop.c,v 1.9 2004/12/03 15:57:08 elmar Exp $ */

#include <stdio.h>
#include "define.h"
#include "floquet.h"
#include "quasisteady.h"

extern int 	Qflag;
extern int 	Fflag;
extern int	oflag;

extern double	ALPHA;

#ifdef DEBUG
extern struct statstruct {
	unsigned long i;
	unsigned long jobs;
	time_t starttime;
	diterator *t, *R, *k;
	iiterator *n, *N, *alpha;
	double Nerr;
	int Nadapt;
} stat;
#endif /* DEBUG */

unsigned int numberof(iiterator *x) {
	if (x->max-x->min) {
		return (unsigned long)((x->max-x->min)/x->n + 1);
	} else return (unsigned long)1;
}

void printtuple(tuple *x,const char *s) {
	fprintf(stderr,"tuple %s> R=%g a=%d n=%d N=%d k=%.16f t=%.16f sigma=%.16f omega=%.16f\n", s, x->R, x->a, x->n, x->N, x->k, x->t, x->sigma, x->omega);
}

void tuplecp(tuple *src, tuple *dest) {
	dest->R = src->R;
	dest->a = src->a;
	dest->n = src->n;
	dest->k = src->k;
	dest->t = src->t;
	dest->N = src->N;
	dest->sigma = src->sigma;
	dest->omega = src->omega;
}

void ditercp(diterator *src, diterator *dest) {
	dest->val = src->val;
	dest->min = src->min;
	dest->max = src->max;
	dest->n   = src->n;
	dest->ndef= src->ndef;
}

/* hooke version (does not work)
void loop(diterator *Rey, iiterator *alpha, diterator *t, diterator *k, iiterator *n, iiterator *N) {
	int iR, ik, o;
	tuple new,max;
	N->val = N->min;
	for (alpha->val=alpha->min; alpha->val<=alpha->max; alpha->val+=alpha->n) {
		for (n->val=n->min; n->val<=n->max; n->val += n->n) {
			for (iR=0; iR<Rey->n; iR++) {
				SETCNT(Rey,iR);
				for (ik=0; ik<k->n; ik++) {
					SETCNT(k,ik);
					if (Qflag) quasisteady(Rey,alpha,t,k,n,N,&new);
					if (Fflag) floquet(Rey,alpha,t,k,n,N);
					if (oflag) {
						if (!ik) {
							tuplecp(&new,&max);
						} else if (new.sigma > max.sigma) {
							tuplecp(&new,&max);
						}
					}
				}
				if (oflag) {
					printtuple(&max,"start");
					max.N = N->max;
					quasisteady_opti(&max);
					printtuple(&max,"end");
				}
			}
		}
	}
}
*/

void loop(diterator *Rey, iiterator *alpha, diterator *t, diterator *k, iiterator *n, iiterator *N) {
	int iR, ik, o;
	double dk,dt;
	diterator origtime;
	diterator origk;
	tuple max, new;

#ifdef DEBUG
	stat.jobs  = (unsigned long)1;
	stat.jobs *= numberof(alpha);
	stat.jobs *= numberof(n);
	stat.jobs *= (unsigned long)(Rey->n);
	stat.jobs *= (unsigned long)(k->n);
	stat.jobs *= (unsigned long)(t->n);
	stat.i = (unsigned long)0;
#endif /* DEBUG */
	N->val = N->min;
	for (alpha->val=alpha->min; alpha->val<=alpha->max; alpha->val+=alpha->n) {
		for (n->val=n->min; n->val<=n->max; n->val += n->n) {
			for (iR=0; iR<Rey->n; iR++) {
				SETCNT(Rey,iR);
				dk = 1.0;
				if (oflag) {
					ditercp(t,&origtime);
					ditercp(k,&origk);
				}
				for (o=0; o<=oflag; o++) {
//					fprintf(stderr,"==> o=%d\n",o);
					if (oflag) {
						N->val = N->min + (o*(N->max-N->min))/oflag;
//						fprintf(stderr,"N=%d\n",N->val);
						if (o) {
							dk = (k->max - k->min)/10.0;
							k->min = max.k - dk; if (k->min <= 0.0) k->min = 0.0000001;
							k->max = max.k + dk;
							k->n = 11;
							dt = (t->max - t->min)/10.0;
							t->min = max.t - dt;
							t->max = max.t + dt;
							t->n = 11;
						}
//						printf("dk=%g dt=%g N=%d\n",dk,dt,N->val);
					}
					for (ik=0; ik<k->n; ik++) {
						SETCNT(k,ik);
						if (Qflag) {
//							fprintf(stderr,"Q: R=%g a=%d t=[%g %g] k=%g n=%d N=%d\n",Rey->val,alpha->val,t->min,t->max,k->val,n->val,N->val);
							quasisteady(Rey,alpha,t,k,n,N,&new);
							if (oflag) {
								if (!ik) {
									tuplecp(&new,&max);
								} else if (new.sigma > max.sigma) {
									tuplecp(&new,&max);
								}
							}
						}
						if (Fflag) floquet(Rey,alpha,t,k,n,N);
					}
					if (!oflag) dk = 0.0;
//					printtuple(&max,">>>");
				}
				if (oflag) {
					printf("%g %g %d %.16f %.16f %.16f %.16f\n",max.R,ALPHA,max.n,max.k,max.sigma,max.omega,max.t);
					fflush(stdout);
					ditercp(&origtime,t);
					ditercp(&origk,k);
				}
			}
		}
	}
}
