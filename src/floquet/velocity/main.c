/* $Id: main.c,v 1.4 2005/09/12 13:52:39 elmar Exp $ */

/*
 * velocity - disturbance velocity
 *
 * stdin:
 * 	time q0.r q0.i q1.r q1.i ... qN.r qN.i
 */

static char id[]="$Id: main.c,v 1.4 2005/09/12 13:52:39 elmar Exp $";
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "expvel.h"
#include "define.h"
#include "quad.h"

#define OPTIONS "R:k:n:N:z:eh"

struct expvelfn {
	double u[RSIZE][ESIZE];
	double v[RSIZE][ESIZE];
	double w[RSIZE][ESIZE];
	double ev[ESIZE];
	int rsize;
	int esize;
} EXPVEL;

double 	 R = 0.0;
char	*BASEPREFIX;
double   K = 0.0;
double 	 TIME = -1.0;
double	 Z = 0.0;
int	 Nphi = -1;
int	 NE = 0;

int	 Bflag = 0;
int	 Eflag = 0;
int	 eflag = 0;

char 	*progname;

void usage(void) {
	printf("usage: %s OPTIONS\n",progname);
	printf("\tOPTIONS: %s\n",OPTIONS);
	exit(1);
}

double retrans(double *qr, double *qi, double *u, double *v, double *w, int rc) {
        int 	 i;
        double	 ur, vr, wr, ui, vi, wi;

	ur = 0.0; vr = 0.0; wr = 0.0;
	ui = 0.0; vi = 0.0; wi = 0.0;
	for (i=0; i<NE; i++) {
		switch (Nphi%4) {
		case 0:
			ur -= qr[i]*EXPVEL.u[rc][i];
			ui -= qi[i]*EXPVEL.u[rc][i];
			vr += qi[i]*EXPVEL.v[rc][i];
			vi -= qr[i]*EXPVEL.v[rc][i];
			wr += qi[i]*EXPVEL.w[rc][i];
			wi -= qr[i]*EXPVEL.w[rc][i];
			break;
		case 1:
			ur += qi[i]*EXPVEL.u[rc][i];
			ui -= qr[i]*EXPVEL.u[rc][i];
			vr += qr[i]*EXPVEL.v[rc][i];
			vi += qi[i]*EXPVEL.v[rc][i];
			wr += qr[i]*EXPVEL.w[rc][i];
			wi += qi[i]*EXPVEL.w[rc][i];
			break;
		case 2:
			ur += qr[i]*EXPVEL.u[rc][i];
			ui += qi[i]*EXPVEL.u[rc][i];
			vr -= qi[i]*EXPVEL.v[rc][i];
			vi += qr[i]*EXPVEL.v[rc][i];
			wr -= qi[i]*EXPVEL.w[rc][i];
			wi += qr[i]*EXPVEL.w[rc][i];
			break;
		case 3:
			ur -= qi[i]*EXPVEL.u[rc][i];
			ui += qr[i]*EXPVEL.u[rc][i];
			vr -= qr[i]*EXPVEL.v[rc][i];
			vi -= qi[i]*EXPVEL.v[rc][i];
			wr -= qr[i]*EXPVEL.w[rc][i];
			wi -= qi[i]*EXPVEL.w[rc][i];
			break;
		default: exit(1);
		}
	}
	*u = ur*cos(Z*M_PI) - ui*sin(Z*M_PI);
	*v = vr*cos(Z*M_PI) - vi*sin(Z*M_PI);
	*w = wr*cos(Z*M_PI) - wi*sin(Z*M_PI);
	if (eflag)
		return hypot(ur,ui) + hypot(vr,vi) + hypot(wr,wi);
	else
		return 0.0;
                //printf("%g %g %g %g\n",(double)rc/(double)(RSIZE-1),off+u,off+v,off+w);
}

void plot(double *qr, double *qi) {
	int rc;
	double	 u, v, w;
	double	 off;

	off = 10.0*Z;
	printf("\n:1:2:green\n");
        for (rc=0; rc<RSIZE; rc+=32) {
		retrans(qr, qi, &u, &v, &w, rc);
                printf("%g %g\n",u+off,(double)rc/(double)(RSIZE-1));
	}
	printf("\n:1:2:yellow\n");
        for (rc=0; rc<RSIZE; rc+=32) {
		retrans(qr, qi, &u, &v, &w, rc);
                printf("%g %g\n",v+off,(double)rc/(double)(RSIZE-1));
	}	
	printf("\n:1:2:red\n");
        for (rc=0; rc<RSIZE; rc+=32) {
		retrans(qr, qi, &u, &v, &w, rc);
                printf("%g %g\n",w+off,(double)rc/(double)(RSIZE-1));
	}
	printf("\n");
}

void energy(double *qr, double *qi) {
	int	 rc;
	double	 E[RSIZE];
	double	 e, u, v, w;
        for (rc=0; rc<RSIZE; rc++)
		E[rc] = (double)rc/(double)(RSIZE-1)*retrans(qr, qi, &u, &v, &w, rc);
	e = quad(E);
	printf("%g %g\n",TIME,e);
}

/* split arguments of the form  a:b[:c] */
void dsplitarg(diterator *X, const char *str) {
        char s[1024], *p, *tokens[5], *last;
        int i;

        i = 0;
        (void)strncpy(s,str,sizeof(s)-1);
        s[sizeof(s)-1] = '\0';
        for ((p = strtok_r(s, ":", &last)); p;
                (p = strtok_r(NULL, ":,", &last))) {
                        if (i < 5 - 1)
                                tokens[i++] = p;
        }
        tokens[i] = NULL;
        if (1==i) {
                X->min = atof(str);
                X->max = atof(str);
                X->n = 1;
        } else if (2==i) {
                X->min = atof(tokens[0]);
                X->max = atof(tokens[1]);
                X->n = X->ndef;
        } else if (3==i) {
                X->min = atof(tokens[0]);
                X->max = atof(tokens[1]);
                X->n = atoi(tokens[2]);
        } else {
                fprintf(stderr,"error in dsplitarg: cannot split %s\n",str);
                exit(1);
        }
        X->val = X->min;
        if ((X->max < X->min)||(X->n<1)) {
                fprintf(stderr,"error: wrong argument %s\n",str);
                exit(1);
        }
}

int main(int args, char **argv) {
	int	 c;
	extern char 	*optarg;
	char 	defaultbase[]=DEFBASEFILE;
	double	 qr[ESIZE], qi[ESIZE];
	double	 time;
	diterator	 z;
	z.val = 0.0;	z.min = 0.0; 	z.max = 0.0;	z.n = 1;	z.ndef = 1;

	progname = argv[0];
	BASEPREFIX = defaultbase;
	while ((c=getopt(args, argv, OPTIONS)) != -1) {
		switch(c) {
		case 'R':
			R = atof(optarg);
			break;
		case 'k':
			K = atof(optarg);
			break;
		case 'n':
			Nphi = atoi(optarg);
			break;
		case 'N':
			NE = atoi(optarg);
			break;
		case 'z':
			dsplitarg(&z,optarg);
			break;
		case 'e':
			eflag = 1;
			break;
		case 'h':
			usage();
		default:
			fprintf(stderr,"getopt default:\n");
			usage();
		}
	}
	if (R<=0.0)  { fprintf(stderr,"main error: option -R missing [R>0.0]\n"); exit(1); }
	if (K<=0.0)    { fprintf(stderr,"main error: option -k missing [K>0.0]\n"); exit(1); }
	if (Nphi<0)       { fprintf(stderr,"main error: option -n missing [n>=0]\n"); exit(1); }
	if (NE<1)       { fprintf(stderr,"main error: option -N missing [N>0]\n"); exit(1); }


	expvel_init(Nphi, K, NE);
	while (scanf("%lf",&time) != EOF) {
		for (c=0; c<NE; c++) {
			scanf("%lf",&(qr[c]));
			scanf("%lf",&(qi[c]));
		}
		if (eflag) {
			TIME=time;
			energy(qr, qi);
		} else {
			for (c=0; c<z.n; c++) {
				if (z.n==1) z.val = z.min; else z.val = (z.min + (double)c*(z.max - z.min)/(double)(z.n-1));
				Z = z.val;
				plot(qr, qi);
				fflush(stdout);
			}
			printf(":frame\n");
		}
	}


	return 0;
}
