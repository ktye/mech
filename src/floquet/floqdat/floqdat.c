/* $Id: floqdat.c,v 1.7 2006/01/18 15:52:46 elmar Exp $ */

/*
 * floqdat - print data from floquet calculations
 *
 * 	derived from 
 * 		src/smooth
 * 		src/intergrid
 */

/*
 * BUG: smoothlevel 1=2 3=4 5=6 ...
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>

#define	 Cols	 8
#define	 Xdim	75
#define	 Ydim	99

char	*progname;

char	*filename;
char def0[] = "/home/elmar/dat/floquet/job-05/dat.0";
char def1[] = "/home/elmar/dat/floquet/job-05/dat.1";

double	 mem[Xdim][Ydim][Cols];
double	 mem0[Xdim][Ydim][Cols];

int	 flipmem = 0;
int	 cflag;

/* 
 * COLS:
 *	0 R
 *	1 alpha
 *	2 k
 *	3 sigma
 *	4 omega
 *	5 tsin
 *	6 dt
 *	7 rho
 */

const char	*var[]  = { "k","kappa","sigma","omega","Sigma","Omega","t","tsin","dt","epsilon","rho","rhod","sigmad","omegad","sigmaw"};
int	 vars, varn;
/*
 * [0]  k	wavenumber 	[radius]
 * [1]  kappa	wavenumber 	[boundary layer thickness]
 * [2]  sigma	amplification 	[transport time]
 * [3]  omega	oscillation	[transport time] (wrong sign)
 * [4]  Sigma	amplification	[diffusive time]
 * [5]  Omega	oscillation	[diffusive time]
 * [6]  t	time		[cos(pi*t)]
 * [7]  tsin	time		[sin(pi*t)]
 * [8]  dt	unstable times
 * [9]  epsilon	1/quasi-steady
 * [10]	rho	radius of act.	[radius]
 * [11] rhod	1-radius of act	[boundary layer thickness]
 * [12] sigmad	amplification	[boundary layer time]
 * [13] omegad	oscillation	[boundary layer time]
 * [14] sigmaw	amplification	[oscillation-time]
 */

int	 smoothlevel;
int 	 pflag;

double getval(double *col) {
	double	 z;
	switch(varn) {
		case 0: 
			z = col[2];
			if ((cflag)&&(col[3]<0.0))	z = 0.0;
			break;
		case 1:
			z = col[2]/col[1];
			if ((cflag)&&(col[3]<0.0))	z = 0.0;
			break;
		case 2:
			z = col[3];
			if ((cflag)&&(col[3]<0.0))	z = 0.0;
			break;
		case 3:
			z = col[4];
			if ((cflag)&&(col[3]<0.0))	z = 0.0;
			break;
		case 4:
			z = col[0]*col[3];
			if ((cflag)&&(col[3]<0.0))	z = 0.0;
			break;
		case 5:
			z = fabs(col[0]*col[4]);
			if ((cflag)&&(col[3]<0.0))	z = 0.0;
			break;
		case 6:
			if ((z = col[5] + 0.5) > 1.0)
				z -= 1.0;
			if ((cflag)&&(col[3]<0.0))	z = 0.0;
			break;
		case 7:
			z = col[5];
			if ((cflag)&&(col[3]<0.0))	z = 0.0;
			break;
		case 8:
			z = col[6];
			if ((cflag)&&(col[3]<0.0))	z = 0.0;
			break;
		case 9:
			z = col[0]*hypot(col[3],col[4])/(col[1]*col[1]);
			if ((cflag)&&(col[3]<0.0))	z = 0.0;
			break;
		case 10:
			z = col[7];
			if ((cflag)&&(col[3]<0.0))	z = 0.0;
			break;
		case 11:
			z = (1.0-col[7])*col[1];
			if ((cflag)&&(col[3]<0.0))	z = 0.0;
			break;
		case 12:
			z = col[3]/col[1];
			if ((cflag)&&(col[3]<0.0))	z = 0.0;
			break;
		case 13:
			z = col[4]/col[1];
			if ((cflag)&&(col[3]<0.0))	z = 0.0;
			break;
		case 14:
			z = col[3]*col[0]/col[1]/col[1];
			if ((cflag&&(col[3]<0.0)))	z = 0.0;
			break;
		default:
			exit(1);
			break;
	}
	return z;
}

void print(void) {
	int	 i,j,c;
	double	 z;
	for (j=0; j<Ydim; j++) {
		for (i=0; i<Xdim; i++) {
			if (flipmem) for (c=0; c<Cols; c++) mem[i][j][c] = mem0[i][j][c];
			z = getval(mem[i][j]);
			printf("%g %g %g\n",mem[i][j][1],mem[i][j][0],z);
		}
		printf("\n");
	}
}

void intergrid(double R, double alpha) {
	int	 i,j,I=0,J=0;
	double	 col[Cols], xi, et;

	col[0] = R;
	col[1] = alpha;
	if ((R<mem[0][0][0])||(R>mem[Xdim-1][0][0])||(alpha<mem[0][0][1])||(alpha>mem[0][Ydim-1][1]))
		return ;
	for (i=1; i<Xdim; i++)
		if ((R>=mem[i-1][0][0])&&(R<=mem[i][0][0]))
			I = i;
	for (j=1; j<Ydim; j++)
		if ((alpha>=mem[0][j-1][1])&&(alpha<=mem[0][j][1]))
			J = j;
	xi  = (R-mem[I-1][0][0])/(mem[I][0][0]-mem[I-1][0][0]);
	et  = (alpha-mem[0][J-1][1])/(mem[0][J][1]-mem[0][J-1][1]);
	for (i=2; i<Cols; i++)
		col[i] = xi*et*mem[I][J][i]
			+xi*(1.0-et)*mem[I][J-1][i]
			+(1.0-xi)*et*mem[I-1][J][i]
			+(1.0-xi)*(1.0-et)*mem[I-1][J-1][i];
	printf("%g %g %g\n",alpha,R,getval(col));
}

double smooth(int i, int j, int c, int flip) {
        int a,b,x,y;
        double z;

        z = 0.0;
        for (a=i-1; a<=i+1; a++) {
                x = a;
                if (a<0) x=0;
                if (a>=Xdim) x = Xdim-1;
                for (b=j-1; b<=j+1; b++) {
                        y = b;
                        if (b<0) y=0;
                        if (b>=Ydim) y = Ydim-1;
			if (flip)
				z += mem0[x][y][c];
			else
				z += mem[x][y][c];
                }
        }
        return z/9.0;
}


void usage(void) {
	fprintf(stderr,"usage: %s [-c] [-f filename | -n0|1 ] -v variable [-s smoothlevel] [R1 Alpha1 [R2 Alpha2 ...]]\n",progname);
	fprintf(stderr,"usage: %s -f filename -v variable [-s smoothlevel] -p   (read from stdin)\n",progname);
	fprintf(stderr,"usage: %s -v?	print values\n",progname);
	exit(1);
}

int main(int args, char **argv) {
	int	 i,j,k;
	extern char	*optarg;
	extern int	 optind;
	int	 c;
	int	 n;
	double	 z1, z2;
	FILE	*fp;


	vars = (int)(sizeof(var)/sizeof(void *));
	progname = argv[0];
	filename = (char *)NULL;
	smoothlevel = 0;
	varn = -1;
	pflag = 0;
	cflag = 0;
	n = -1;
	while ((c = getopt(args, argv, "cn:f:v:s:p")) != -1) {
		switch (c) {
		case 'c':	/* clip */
			cflag = 1;
			break;
		case 'n':	/* angular wavenumber */
			n = atoi(optarg);
			if (n == 0)
				filename = def0;
			else if (n == 1)
				filename = def1;
			else {
				fprintf(stderr,"wrong option to -n: 0|1\n");
				exit(1);
			}
			break;
		case 'f':	/* input file */
			filename = optarg;
			break;
		case 's':	/* smoothlevel */
			smoothlevel = atoi(optarg);
			flipmem = smoothlevel%2;
			break;
		case 'v':	/* name of variable to print */
			for (i=0; i<vars; i++) {
				if (!strcmp(var[i],optarg))
					varn = i;
			}
			if (varn<0) {
				fprintf(stderr,"unknown variable -v %s\nOptions are\n",optarg);
				for (i=0; i<vars; i++) {
					fprintf(stderr,"\t%s\n",var[i]);
				}
				exit(1);
			}
			break;
		case 'p':	/* pipe mode, (R,alpha) from stdin */
			pflag = 1;
			break;
		default:
			usage();
		}
	}
	args -= optind;
	argv += optind;

	if (varn < 0) { fprintf(stderr,"error: -v missing\n"); usage(); }

	/* read file to memory */
	if (!filename) { fprintf(stderr,"error: -n | -f missing\n"); usage(); }
	if ((fp = fopen(filename,"r")) != NULL) {
		for (j=0; j<Ydim; j++)	
			for (i=0; i<Xdim; i++)
				for (k=0; k<Cols; k++) {
					if (fscanf(fp,"%lg",&(mem[i][j][k])) != 1) {
						fprintf(stderr,"error: fscanf\n");
						exit(1);
					}
				}
	} else {
		fprintf(stderr,"error: fopen %s\n",filename);
		exit(1);
	}
	fclose(fp);

	if (smoothlevel < 0) { fprintf(stderr,"error: smoothlevel is negative\n"); exit(1); }
	if (smoothlevel) 
		for (j=0; j<Ydim; j++)
			for (i=0; i<Xdim; i++)
				for (c=0; c<Cols; c++)
					mem0[i][j][c] = mem[i][j][c];
	while (smoothlevel--)
		for (j=0; j<Ydim; j++)
			for (i=0; i<Xdim; i++) {
				if (smoothlevel%2) {
					mem[i][j][0] = mem0[i][j][0];
					mem[i][j][1] = mem0[i][j][1];
				} else {
					mem0[i][j][0] = mem[i][j][0];
					mem0[i][j][1] = mem[i][j][1];
				}
				for (c=2; c<Cols; c++) {
					if (smoothlevel%2)
						mem[i][j][c] = smooth(i,j,c,1);
					else
						mem0[i][j][c] = smooth(i,j,c,0);
				}
			}

	if ((!args)&&(!pflag))
		print();

	for (i=0; i<args-1; i+=2)
		intergrid(atof(argv[i]),atof(argv[i+1]));

	if (pflag)
		while (scanf("%lf %lf",&z1,&z2) == 2)
			intergrid(z1,z2);	

	return 0;
}
