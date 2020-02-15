/* $Id: status.c,v 1.2 2004/12/02 09:13:35 elmar Exp $ */

#ifdef DEBUG
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>
#include <time.h>
#include "define.h"
#include <g2c.h>
#include "solve.h"
#include "expvel.h"
#include "readbase.h"

extern double	R;
extern double	K;
extern double	TIME;
extern double	ALPHA;
extern int	NE;
extern int	Nphi;

extern struct statstruct {
	unsigned long i;
	unsigned long jobs;
	time_t starttime;
	diterator *t, *R, *k;
	iiterator *n, *N, *alpha;
	double Nerr;
	int Nadapt;
} stat;


extern struct expvelfn { double   u[RSIZE][ESIZE]; double   v[RSIZE][ESIZE]; double   w[RSIZE][ESIZE]; double ev[ESIZE]; int rsize; int esize; } EXPVEL;

extern char 	*progname;
extern int 	Bflag;
extern int 	Mflag;
extern int 	Sflag;
extern int 	NE;
extern double 	R;
extern double	K;

void status(int sig) {
	time_t now, end;
	long computed;
	double rate;

	now = time(0);
	computed = (long)now - (long)stat.starttime;
	rate = (double)stat.jobs/(double)stat.i;
	end = (long)(rate*(double)((long)now-(long)stat.starttime))+(long)stat.starttime;
	fprintf(stderr,"+---------------------------------------\n");
	fprintf(stderr,"| SIGINFO %s PID %d\n",progname,(int)getpid());
	fprintf(stderr,"+---------------------------------------\n");
	fprintf(stderr,"| start: %s",ctime(&stat.starttime));
	fprintf(stderr,"| now  : %s",ctime(&now));
	fprintf(stderr,"| end  : %s",ctime(&end));
	fprintf(stderr,"+---------------------------------------\n");
	fprintf(stderr,"| R: %g \t[%g %g]x(%d)\n",stat.R->val,stat.R->min,stat.R->max,stat.R->n);
	fprintf(stderr,"| k: %g \t[%g %g]x(%d)\n",stat.k->val,stat.k->min,stat.k->max,stat.k->n);
	fprintf(stderr,"| t: %g \t[%g %g]x(%d)\n",stat.t->val,stat.t->min,stat.t->max,stat.t->n);
	fprintf(stderr,"| alpha: %d \t[%d %d](++%d)\n",stat.alpha->val,stat.alpha->min,stat.alpha->max,stat.alpha->n);
	fprintf(stderr,"| n: %d \t[%d %d](++%d)\n",stat.n->val,stat.n->min,stat.n->max,stat.n->n);
	fprintf(stderr,"| N: %d \t[%d %d](++%d)\n",stat.N->val,stat.N->min,stat.N->max,stat.N->n);
	fprintf(stderr,"+---------------------------------------\n");
	fprintf(stderr,"| %lu of %lu (%3.2f %%)\n",stat.i,stat.jobs,100.0*(double)stat.i/(double)stat.jobs);
	if (stat.Nadapt) fprintf(stderr,"+--status Nadapt, Nerr=%.5g---------\n",stat.Nerr);
	else fprintf(stderr,"+---------------------------------------\n");
}

void terminate(int sig) {
	fprintf(stderr,"signal: %d\n",sig);
	exit(sig);
}

void debug(int sig) {
	int c,i,quit;
	status(0);
	fprintf(stderr,"signal: %d --> DEBUG MODE\nenter h for help\n: ",sig);
	quit = 0;
	while ((c=getchar())!=EOF && !quit) {
		switch (c) {
		case 'h':
			fprintf(stderr,"'h'\thelp\n's'\tstatus\n'M'\tmatrix\n'B'\tbaseflow\n'E'\texpansion eigenvalues\n'L'\texpansion lambda\n'S'\tquasi-steady eigenvalues (SIGMA + i OMEGA)\n'O'\torthogonality of eigenfunctions\n'q'\tquit\n");
			fprintf(stderr,": ");
			break;
		case 's':
			status(0);
			fprintf(stderr,": ");
			break;
		case 'M':
			Mflag = 1;
			solve();
			Mflag = 0;
			fprintf(stderr,": ");
			break;
		case 'B':
			Bflag = 1;
			readbase_settime();
			Bflag = 0;
			fprintf(stderr,": ");
			break;
		case 'E':
			for (i=0; i<NE; i++)
				fprintf(stderr,"EXPVELEV[%d] = %.16f\n",i,EXPVEL.ev[i]);
			fprintf(stderr,": ");
			break;
		case 'L':
			for (i=0; i<NE; i++)
				fprintf(stderr,"lambda[%d]=%.16f\n",i,-(K*K + EXPVEL.ev[i]*EXPVEL.ev[i])/R);
			fprintf(stderr,": ");
			break;
		case 'S':
			Sflag = 1;
			solve();
			Sflag = 0;
			fprintf(stderr,": ");
			break;
		case 'O':
			expvel_ortho(NE);
			fprintf(stderr,": ");
			break;
		case 'q':
			quit = 1;
			break;
		default:
			if (c!='\n') {
				fprintf(stderr,"unknown character '%c', press 'h' for help\n",c);
				fprintf(stderr,": ");
			}
		}
	}
	fprintf(stderr,"quit DEBUG, continue.\n");
}
#endif /* DEBUG */
