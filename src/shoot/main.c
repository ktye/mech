/* $Id: main.c,v 1.16 2004/11/09 15:12:23 elmar Exp $ */

/*
 * shoot - pulsatile pipe flow stability analysis (shooting algorithm)
 */

static char id[]="$Id: main.c,v 1.16 2004/11/09 15:12:23 elmar Exp $";
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include "shoot.h"
#include "rk4.h"
#include "endian.h"
#include "baseflow.h"

int	flag_print = 0;
int	debug = 0;
int	HPFLOW = 0;

int 	Nphi = -9;
unsigned int 	Nsteps = 1025; 
unsigned int	Lsteps = 2049;
int	Iterations;
double 	OMEGAR = 0.0;
double 	OMEGAI = 0.0;
double 	K = 0.0;
double 	R = 0.0;
char 	*basefile;

double	ALPHA = -1.0;
int	ALPHANUM = -1;
double	TIME = -1.0;
int	TIMENUM = -1;
int	RNUM = -1;

void 	(*ode)() = rk4;
const char 	*progname;

void header(int rows, int cols) {
	if (mendian() == 0)	printf("little-endian ");
	else if (mendian() == 1)	printf("big-endian ");
	else {
		fprintf(stderr,"endian error\n");
		exit(1);
	}
	printf("double 2d %d %d\n",rows,cols);
}

int splitomega(double *wr, double *wi, double *wrmin, double *wrmax, double *wimin, double *wimax, const char *str, int *I) {
	char s[1024], *p, *tokens[5], *last;
	int i,j;
	i = *I;

	(void)strncpy(s,str,sizeof(s)-1);
	s[sizeof(s)-1] = '\0';
	for ((p = strtok_r(s, ":,", &last)); p; 
		(p = strtok_r(NULL, ":,", &last))) {
			if (i < 5 - 1)
				tokens[i++] = p;
	}
	tokens[i] = NULL;
	*I = i;
	if (2==i) {
		*wr = atof(tokens[0]);
		*wi = atof(tokens[1]);
		*wrmin = 0.0; *wimin = 0.0; *wrmax = 0.0; *wrmax = 0.0;
		return 0;
	} else if (4==i) {
		*wrmin = atof(tokens[0]);
		*wimin = atof(tokens[1]);
		*wrmax = atof(tokens[2]);
		*wimax = atof(tokens[3]);
		*wr = 0.0; *wi = 0.0;
		return 0;
	}
	return 1;
}  

void usage(void) {
	fprintf(stderr,"%s\nusage: %s [-S|-I|-P] OPTIONS\n",id,progname);
	fprintf(stderr,"usage: %s OPTIONS\n",progname);
	fprintf(stderr,"OPTIONS\n");
	fprintf(stderr,"	Without -S|-I|-P the ode solution of one base (option -b) is printed to \n");
	fprintf(stderr,"	stdout and the determinant to stderr\n");
	fprintf(stderr,"	-n mode\n");
	fprintf(stderr,"		n=-1 angular symmetric mode\n");
	fprintf(stderr,"		n= 0 radial-axial symmetric mode\n");
	fprintf(stderr,"		n= 1 first nonsymmetric (antisymmetric) mode\n\n");
	fprintf(stderr,"	-w omega_r,omega_i\n");
	fprintf(stderr,"		real and imaginary part of frequency\n\n");
	fprintf(stderr,"	-b base\n");
	fprintf(stderr,"		integral base; for n=0: [1,2]; for n=1: [1,2,3]\n\n");
	fprintf(stderr,"	-N log2(steps-1)\n");
	fprintf(stderr,"		[1-15] default: 10, unused if -f\n\n");
	fprintf(stderr,"	-R Reynol's number\n\n");
	fprintf(stderr,"	-f baseflowfile | -H\n");
	fprintf(stderr,"		baseflow input file or -H stationary Hagen-Poiseuille flow\n\n");
	fprintf(stderr,"usage: %s -S -n mode -w omegarange -W number_of_omega -k wavenumber -R reynolds -N steps [-F file | -H]\n",progname);
	fprintf(stderr,"	Option -S will produce a binary output of the determinant for a given \n");
	fprintf(stderr,"	range of the complex omega plane.\n");
	fprintf(stderr,"	-n mode\n\n");
	fprintf(stderr,"	-w min_omega_r,min_omega_i:max_omega_r,max_omega_i\n");
	fprintf(stderr,"		defines the complex omega plane\n\n");
	fprintf(stderr,"	-W number of omega steps (integer)\n\n");
	fprintf(stderr,"	-k wavenumber	\n\n");
	fprintf(stderr,"	-R Reynold's number\n\n");
	fprintf(stderr,"	-N log2(steps-1), unused if -f\n\n");
	fprintf(stderr,"	-F baseflowfile | -H\n\n");
	fprintf(stderr,"usage: %s -I -n mode -w omegarange -W No_omega -k wavenumber -K No_k -R reynolds -N steps [-f file | -H]\n",progname);
	fprintf(stderr,"	Option -I will produce a binary output of the determinant for a given \n");
	fprintf(stderr,"	range of the complex omega plane spanned over a interval of the \n");
	fprintf(stderr,"	wavenumbers using the minimal determinant. Options are the same as for\n");
	fprintf(stderr,"	-S with the additional\n");
	fprintf(stderr,"	-K number of steps from k=0 to the wavenumber (option -k)\n\n");
	fprintf(stderr,"usage: %s -P -n mode [-k wavenumber] [-K steps] -R reynolds [-N steps] [-r rho] [-f baseflowfile | -H] [-t time]\n",progname);
	fprintf(stderr,"	Option -P will search the complex omega plane for a minimal determinant\n");
	fprintf(stderr,"	starting at known starting values for k=0. It increases k and prints\n");
	fprintf(stderr,"	\"k omega_r omega_i\" which lead to the minimal determinant.\n");
	fprintf(stderr,"	Options are\n");
	fprintf(stderr,"	-n mode, default all modes\n\n");
	fprintf(stderr,"	-k maximal wavenumber\n\n");
	fprintf(stderr,"	-K initial K steps\n\n");
	fprintf(stderr,"	-R Reynolds number\n\n");
	fprintf(stderr,"	-N log2(steps-1), unused if -f\n\n");
	fprintf(stderr,"	-r rho Hooke-Jeeves parameter (minima search) [0,1]\n");
	fprintf(stderr,"	-f baseflowfile | -H\n");
	fprintf(stderr,"	-t time [0,1] (half periode)\n\n");

	exit(1);
}

int main(int args, char **argv) {
	int i, mode;
	double rho, wr, wi, wrmin, wrmax, wimin, wimax;
	extern char *optarg;
	extern int optint;
	int c, NoW, Nw, Nk, Ngiven, wgiven, log2steps, nphigiven, tgiven;
	char _basefile[1024];
	double kary[9], tary[3], wrary[9], wiary[9];

	progname = argv[0];
	mode = -1;
	HPFLOW = 0;
	debug = 0;
	Nphi = -9;
	nphigiven = 0;
	K = -1.0;
	Ngiven = 0;
	log2steps = 0;
	wgiven = 0;
	Nsteps = 1025;
	Lsteps = 2049;
	wr = 0.0; wi = 0.0; wrmin = 0.0; wrmax = 0.0; wimin = 0.0; wimax = 0.0;
	basefile = (char *)NULL;
	R = 0.0; RNUM=(int)R;
	NoW = 0;
	Nw = 0;
	Nk = 0;
	rho = -1.0;
	TIME = -1.0;
	tgiven = 0;
	while ((c = getopt(args, argv, "vSIPMHdhb:n:w:k:N:f:R:W:K:r:t:")) != -1 ) {
		switch (c) {
		case 'v':
			printf("%s\n",id);
			return 0;
		case 'S': 
			mode = 0; 
			break;
		case 'I': 
			mode = 1; 
			break;
		case 'P': 
			mode = 2; 
			break;
		case 'M':
			mode = 3;
			break;
		case 'H': 
			HPFLOW = 1; 
			break;
		case 'd': 
			debug = 1; 
			break;
		case 'h': 
			usage(); 
			break;
		case 'b': 
			debug = atoi(optarg); 
			break;
		case 'n': 
			nphigiven = 1;
			Nphi = atoi(optarg); 
			break;
		case 'w': 
			if (splitomega(&wr,&wi,&wrmin,&wrmax,&wimin,&wimax,optarg,&NoW)) {
				fprintf(stderr,"wrong argument to -w\n"); 
			}
			wgiven = 1; 
			break;
		case 'k': 
			K = atof(optarg); 
			break;
		case 'N': 
			Ngiven = 1; 
			log2steps = atoi(optarg); 
			if ((log2steps<0)||(log2steps>15)) { fprintf(stderr, "error: -N argument not in [1,15] \n"); exit(1); }
			Nsteps = (1<<log2steps)+1;
			Lsteps = (2<<log2steps)+1;
			break;
		case 'f': 
			strncpy(_basefile,optarg,sizeof(_basefile)); 
			basefile=_basefile;
			break;
		case 'R': 
			R = atof(optarg); 
			RNUM = (int)R;
			break;
		case 'W': 
			Nw = atoi(optarg); 
			break;
		case 'K': 
			Nk = atoi(optarg); 
			break;
		case 'r': 
			rho = atof(optarg); 
			break;
		case 't':
			TIME = atof(optarg);
			tgiven = 1;
			break;
		default:
			fprintf(stderr,"getopt default:\n");
			usage();
		}
	}
	args -= optind;
	argv += optind;

	if ((basefile == NULL)&&(HPFLOW == 0)) { fprintf(stderr,"missing option: -H|-f baseflow\n"); exit(1); }
	else if ((basefile != NULL)&&(HPFLOW == 1)) { fprintf(stderr,"concurring options: -H, -f baseflow\n"); exit(1); }
	if ((Nphi < -1)||(Nphi >1)||(!nphigiven)) { fprintf(stderr,"-n out of [-1,0,1]\n"); exit(1); }
	if ((K==-1.0)&&(mode!=2)&&(mode!=3)) { fprintf(stderr,"missing option -k\n"); exit(1); }
	if ((Ngiven)&&(basefile!=NULL)) {fprintf(stderr,"concurring optins -N, -f\n"); exit(1); }
	if (R == 0.0) { fprintf(stderr,"missing option -R != 0.0 R=%g\n",R); exit(1); }
	if (!wgiven&&(mode!=2)&&(mode!=3)) { fprintf(stderr,"missing option -w\n"); exit(1); }
	if ((tgiven)&&((TIME>2.0)||(TIME<0.0))) {fprintf(stderr, "missing option -t time in [0,1]\n"); exit(1); }

	initbaseflow();
	switch(mode) {
	case (-1): 
		if (debug == 0) debug = 1;
		if ((wrmin!=0.0)||(wrmax!=0.0)||(wimin!=0.0)||(wimax!=0.0)) { fprintf(stderr,"wrong argument to -w\n"); exit (1); }
		OMEGAR = wr;
		OMEGAI = wi;
fprintf(stderr,"DEGUB: shoot n=%d w=(%g,%g) K=%g R=%g N=%u f=%s alpha=%g time=%g HP=%d\n",Nphi,OMEGAR,OMEGAI,K,R,Nsteps,basefile,ALPHA,TIME,HPFLOW);
		fprintf(stderr,"%g\n", shoot()); return 0;
	case ( 0): 
		if (((wr!=0.0)||(wi!=0.0))||((wrmin==0.0)&&(wrmax==0.0))) { fprintf(stderr,"wrong argument to -w\n"); exit(1); }
		if ((wrmin>=wrmax)||(wimin>=wimax)) { fprintf(stderr,"error: (wrmin>=wrmax)||(wimin>=wimax)\n"); exit(1); }
		if (Nw<1) { fprintf(stderr,"wrong argument to -W\n"); exit(1); }
fprintf(stderr,"DEGUB: single n=%d w=(%g,%g)x(%g,%g) K=%g Nw=%d R=%g N=%u f=%s alpha=%g time=%g HP=%d\n",Nphi,wrmin,wimin,wrmax,wimax,K,Nw,R,Nsteps,basefile,ALPHA,TIME,HPFLOW);
		single(wrmin, wrmax, wimin, wimax, Nw); return 0;
	case ( 1): 
		if (((wr!=0.0)||(wi!=0.0))||((wrmin==0.0)&&(wrmax==0.0))) { fprintf(stderr,"wrong argument to -w\n"); exit(1); }
		if ((wrmin>=wrmax)||(wimin>=wimax)) { fprintf(stderr,"error: (wrmin>=wrmax)||(wimin>=wimax)\n"); exit(1); }
		if (Nw<1) { fprintf(stderr,"wrong argument to -W\n"); exit(1); }
		if (Nk<1) { fprintf(stderr,"wrong argument to -K\n"); exit(1); }
fprintf(stderr,"DEGUB: inter n=%d w=(%g,%g)x(%g,%g) K=%g Nw=%d Nk=%d R=%g N=%u f=%s alpha=%g time=%g HP=%d\n",Nphi,wrmin,wimin,wrmax,wimax,K,Nw,Nk,R,Nsteps,basefile,ALPHA,TIME,HPFLOW); return 0;
		inter(wrmin, wrmax, wimin, wimax, Nw, Nk);  return 0;
	case ( 2): 
		if (tgiven) {
			TIMENUM=(int)(100.0*TIME);
			if (debug) {
				printbaseflow();
				return 0;
			}
			plot(rho, Nk);
		} else {
			for (i=0; i<100; i++) {
				TIMENUM=i;
				TIME=(double)i/100.0;
				if (TIME>0.5) {
					TIME = 0.25*(1.0+sqrt(16.0*TIME-7.0));
				} else {
					TIME = 3.0/4.0 - sqrt(9.0/16.0 - TIME);
				}
				reinitbaseflow();
				plot(rho, Nk);
			}
		}
		printf("#status terminated\n");
		return 0;
	case ( 3):
		scanf("%lg",&tary[0]);
		scanf("%lg",&kary[0]);
		scanf("%lg",&kary[1]);
		scanf("%lg",&kary[2]);
		scanf("%lg",&wrary[0]);
		scanf("%lg",&wrary[1]);
		scanf("%lg",&wrary[2]);
		scanf("%lg",&wiary[0]);
		scanf("%lg",&wiary[1]);
		scanf("%lg",&wiary[2]);
		scanf("%lg",&tary[1]);
		scanf("%lg",&kary[3]);
		scanf("%lg",&kary[4]);
		scanf("%lg",&kary[5]);
		scanf("%lg",&wrary[3]);
		scanf("%lg",&wrary[4]);
		scanf("%lg",&wrary[5]);
		scanf("%lg",&wiary[3]);
		scanf("%lg",&wiary[4]);
		scanf("%lg",&wiary[5]);
		scanf("%lg",&tary[2]);
		scanf("%lg",&kary[6]);
		scanf("%lg",&kary[7]);
		scanf("%lg",&kary[8]);
		scanf("%lg",&wrary[6]);
		scanf("%lg",&wrary[7]);
		scanf("%lg",&wrary[8]);
		scanf("%lg",&wiary[6]);
		scanf("%lg",&wiary[7]);
		scanf("%lg",&wiary[8]);
		tmaxima(tary,kary,wrary,wiary);
		return 0;
	default:
		perror("internal error: switch mode");
		return 1;
	}
}
