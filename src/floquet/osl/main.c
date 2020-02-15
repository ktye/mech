/* $Id: main.c,v 1.7 2005/05/23 16:44:07 elmar Exp $ */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <string.h>

#include "define.h"
#include "debug.h"
#include "loop.h"

#define	OPTIONS	"K:L:N:R:S:T:tkrd:vmb:"
char	*progname;

double	 K;
double	 L;
int	 N;
double	 R;
int	 S;
double	 T;
int	 tflag, kflag, rflag, vflag, mflag, bflag;

double	 BASEFLOW[YSIZE];
double	 BASEFLOW2[YSIZE];
double	 PHI[YSIZE][NMAX];
double	 PSI[YSIZE][NMAX];
double	 DPSI[YSIZE][NMAX];
double	 LAMBDA[NMAX];

void usage(void) {
	fprintf(stderr,"usage: %s %s\n",progname,OPTIONS);
	exit(1);
}

void outdom(const char *var, const char *domain) {
	fprintf(stderr,"error: option %s out of domain %s\n",var,domain);
	exit(1);
}

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
	extern char	*optarg;
	extern int	 optind;
	int	 c;
	char	*debc;
	diterator	 time, kay, rey;

	progname = argv[0];
	L = 10.0;
	N = 0;
	S = 0;
	tflag = 0;
	kflag = 0;
	rflag = 0;
	vflag = 0;
	bflag = 0;
	debc = (char *)NULL;
	time.val = 0.0;	time.min = 0.0; time.max = 1.0; time.n = 128; time.ndef = 128;
	kay.val = 0.0;	kay.min = 0.0; kay.max = 1.0; kay.n = 128; kay.ndef = 128;
	rey.val = 0.0;	rey.min = 0.0; rey.max = 1.0; rey.n = 128; rey.ndef = 128;
	while ((c = getopt(args, argv, OPTIONS)) != -1) {
		switch (c) {
		case 'K':	/* axial wavenumber in boundary layer thickness lengthscale */
			dsplitarg(&kay, optarg);
			break;
		case 'L':	/* infinity */
			L = atof(optarg);
			break;
		case 'N':	/* system order */
			N = atoi(optarg);
			break;
		case 'R':	/* Reynolds number */
			dsplitarg(&rey, optarg);
			break;
		case 'S':	/* S=1: symmetric, S=0 antimetric */
			S = atoi(optarg);
			break;
		case 'T':	/* time [0,1) */
			dsplitarg(&time, optarg);
			break;
		case 't':	/* maximise over T-loop */
			tflag = 1;
			break;
		case 'k':	/* maximise over T and K-loop */
			kflag = 1;
			tflag = 1;
			break;
		case 'r':	/* maximise over T,K and R-loop */
			rflag = 1;
			kflag = 1;
			tflag = 1;
			break;
		case 'd':
			debc = optarg;
			break;
		case 'v':
			vflag = 1;
			break;
		case 'm':
			mflag = 1;
			break;
		case 'b':	/* type of baseflow, 1:OSL, 2: Poiseuille, else: pressure driven OSL */
			bflag = atoi(optarg);
			break;
		default:
			usage();
		}
	}
	args -= optind;
	argv += optind;

	if (kay.val <= 0.0) 	outdom("K",">0 (double)");
	if (L < 1.0)	outdom("L",">>1 (double)");
	if ((N <= 0)||(N>NMAX)) 	outdom("N","[0,NMAX] (int)");
	if (rey.val<0.0)	outdom("R",">0 (double)");
	if ((S != 0)&&(S!=1))	outdom("S","[0|1]");
//	if (time.val<0.0)	outdom("T","(0,1)");
	
	/* scaling to fit [-1,1] channel width */
	rey.min *= L;
	rey.max *= L;
	rey.val *= L;
	kay.min *= L;
	kay.max *= L;
	kay.val *= L;

	R = rey.val;
	K = kay.val;
	T = time.val;

	if (debc)
		debug(debc);
	loop(&rey, &kay, &time);
	
	return 0;
}
