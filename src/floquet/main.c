/* $Id: main.c,v 1.21 2006/01/09 09:58:22 elmar Exp $ */
/* main.c | elmar@lighthill OpenBSD i386 | Thu Nov 11 14:24:25 CET 2004 */

/*
 * floquet - floquet and quasi--steady galerkin method for pipe flow stability analysis
 */

static char id[]="$Id: main.c,v 1.21 2006/01/09 09:58:22 elmar Exp $";
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <time.h>
#include <string.h>
#include <signal.h>
#include "define.h"
#include "reinit.h"
#include "floquet.h"
#include "expvel.h"
#include "loop.h"
#ifdef DEBUG
#include "status.h"
#endif /* DEBUG */

#define OPTIONS "QFTo:GgMBELSVvAOf:a:R:k:n:N:t:h"

struct expvelfn {
	double u[RSIZE][ESIZE];
	double v[RSIZE][ESIZE];
	double w[RSIZE][ESIZE];
	double ev[ESIZE];
	int rsize;
	int esize;
} EXPVEL;

struct baseflowstruct {
	double     r[RSIZE];
	double    dr[RSIZE];
	double  real[RSIZE];
	double dreal[RSIZE];
	double  imag[RSIZE];
	double dimag[RSIZE];
} baseflow;

#ifdef DEBUG
struct statstruct {
	unsigned long i;
	unsigned long jobs;
	time_t starttime;
	diterator *t, *R, *k;
	iiterator *n, *N, *alpha;
	double Nerr;
	int Nadapt;
} stat;
#endif /* DEBUG */

double SIGMA[ESIZE]; /* Eigenvalues of quasi-steady time evolution      */
double OMEGA[ESIZE]; /* Imag part;  ~= exp [ (SIGMA + %I*OMEGA)*t ]	*/

double 	R = 0.0;
double 	ALPHA = -1.0;
char	*BASEPREFIX;
double  K = 0.0;
double 	TIME = -1.0;
int	Nphi = -1;
int	NE = 0;
int	Ninit = 0;

int	Qflag;
int	Fflag;
int 	oflag;
int	Gflag;
int	gflag;
int	Mflag;
int	Bflag;
int	Eflag;
int	Lflag;
int	Sflag;
int	Tflag;
int	Vflag;
int	vflag;
int 	Aflag;
int 	Oflag;

char *progname;

void usage(void) {
	printf("usage: %s [-M] -f basefile -R Reynolds -k wavnumber -n phimode -N system_order [-t time]\n",progname);
	printf("\tOPTIONS: %s\n",OPTIONS);
	exit(1);
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

void isplitarg(iiterator *X, const char *str) {
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
		X->min = atoi(str);
		X->max = atoi(str);
		X->n = 1;
	} else if (2==i) {
                X->min = atoi(tokens[0]);
                X->max = atoi(tokens[1]);
		X->n = X->ndef;
        } else if (3==i) {
                X->min = atoi(tokens[0]);
                X->max = atoi(tokens[1]);
                X->n = atoi(tokens[2]);
        } else {
		fprintf(stderr,"error in isplitarg: cannot split %s\n",str);
		exit(1);
	}
	X->val = X->min;
	if ((X->max < X->min)||(X->n<1)) {
		fprintf(stderr,"error: wrong argument %s\n",str);
		exit(1);
	}
}

int main(int args, char **argv) {
	int c;
	extern char *optarg;
	char defaultbase[]=DEFBASEFILE;
	diterator k, Rey, t;
	iiterator n, N, alpha;

	printf("# "); for (c=0; c<args; c++) printf("%s ",argv[c]); printf("\n");

	progname = argv[0];
#ifdef DEBUG
	(void)signal(SIGHUP,debug);
	(void)signal(SIGTERM,terminate);
	(void)signal(SIGINT,terminate);
	(void)signal(SIGUSR1,status);
	stat.starttime = time(0);
	stat.t = &t;
	stat.k = &k;
	stat.R = &Rey;
	stat.n = &n;
	stat.N = &N;
	stat.alpha = &alpha;
#endif /* DEBUG */
	Rey.val = -1.0;	Rey.ndef = 2; 
	k.val = 0.0; 	k.ndef = 10;
	t.val = 0.0;	t.min = 0.0; 	t.max = 1.0; 	t.n = 128; 	t.ndef = 128;
	alpha.val=0; 	alpha.min = 0; 	alpha.max=99;	alpha.n = 1;	alpha.ndef = 1;
	n.val = -1;	n.ndef = 1;
	N.val = 0;	N.ndef = 10;
	Qflag = 0;
	Fflag = 0;
	Gflag = 0;
	gflag = 0;
	Mflag = 0;
	Bflag = 0;
	Eflag = 0;
	Lflag = 0;
	Sflag = 0;
	Tflag = 0;
	Vflag = 0;
	vflag = 0;
	Aflag = 0;
	Oflag = 0;
	oflag = 0;
	BASEPREFIX = defaultbase;
	while ((c=getopt(args, argv, OPTIONS)) != -1) {
		switch(c) {
		case 'Q': /* quasi-steady */
			Qflag = 1;
			break;
		case 'F': /* floquet */
			Fflag = 1;
			break;
		case 'T': /* output at all times */
			Tflag = 1;
			break;
		case 'o': /* optimise */
			oflag = atoi(optarg);
			break;
		case 'G': /* print floquet matrix */
			Gflag = 1;
			break;
		case 'g': /* print floquet eigenvalues */
			gflag = 1;
			break;
		case 'M': /* print matrix */
			Mflag = 1;
			break;
		case 'B': /* print baseflow */
			Bflag = 1;
			break;
		case 'E': /* print expansion eigenvalues EXPVEL.ev[] */
			Eflag = 1;
			break;
		case 'L': /* print expansion lambda (eigenvalues) lambda = -((EXPVEL.ev[])^2+K^2)/R */
			Lflag = 1;
			break;
		case 'S': /* print quasi-steady eigenvalues SIGMA, OMEGA */
			Sflag = 1;
			break;
		case 'V': /* print most unstable eigenfunction */
			Vflag = 1;
			break;
		case 'v': /* print most unstable eigenvector */
			vflag = 1;
			break;
		case 'A': /* print expansion ansatz(eigen)function N */
			Aflag = 1;
			break;
		case 'O': /* check orthogonality of eigenfunctions */
			Oflag = 1;
			break;
		case 'f':
			BASEPREFIX = optarg;
			break;
		case 'a':
			isplitarg(&alpha,optarg);
			break;
		case 'R':
			dsplitarg(&Rey,optarg);
			break;
		case 'k':
			dsplitarg(&k,optarg);
			break;
		case 'n':
			isplitarg(&n,optarg);
			break;
		case 'N':
			isplitarg(&N,optarg);
			break;
		case 't':
			dsplitarg(&t,optarg);
			break;
		case 'h':
			usage();
		default:
			fprintf(stderr,"getopt default:\n");
			usage();
		}
	}
	if (Rey.val<0.0)  { fprintf(stderr,"main error: option -R missing [R>0.0]\n"); exit(1); }
	if (k.val<=0.0)    { fprintf(stderr,"main error: option -k missing [K>0.0]\n"); exit(1); }
	if (alpha.val<0)   { fprintf(stderr,"main error: option -f missing\n"); exit(1); }
	if (n.val<0)       { fprintf(stderr,"main error: option -n missing [n>=0]\n"); exit(1); }
	if (N.val<1)       { fprintf(stderr,"main error: option -N missing [N>0]\n"); exit(1); }
	if (oflag<0)	{ fprintf(stderr,"main error: oflag < 0\n"); exit(1); }

	/*
	for (i=0; i<k.n; i++) {
		k.val = (k.n==1)||(k.min + (double)i*(k.max - k.min)/(double)(k.n-1));
		printf("i=%d k=%g\n",i,k.val);
	}
	return 0;
	*/

	reinit(Rey.val,alpha.val,k.val,n.val,t.val,N.val);
	if (Aflag) { for (c=0; c<RSIZE; c++) printf("%.16f %.16f %.16f %.16f\n",(double)c/(double)(RSIZE-1),EXPVEL.u[c][N.val-1],EXPVEL.v[c][N.val-1],EXPVEL.w[c][N.val-1]); return 0; }
	if (Bflag) for (c=0; c<RSIZE; c++) printf("%.16f %.16f %.16f\n",(double)c/(double)(RSIZE-1),baseflow.r[c],baseflow.dr[c]);
	if (Oflag) { expvel_ortho(N.val); exit(0); } 
	if (Fflag) { t.min=0.0; t.max=2.0; }
	loop(&Rey,&alpha,&t,&k,&n,&N);

	return 0;
}

// ./o.floquet -M -R3000 -k0.01 -n2 -N5 -t0.0 -f ~/dat/baseflow/baseflow-60
