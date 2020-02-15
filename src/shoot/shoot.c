/* $Id: shoot.c,v 1.40 2004/11/04 12:26:04 elmar Exp $ */

/*
 * shoot - pulsatile pipe flow stability analysis (shooting algorithm)
 */

static char id[]="$Id: shoot.c,v 1.40 2004/11/04 12:26:04 elmar Exp $";
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "v.h"
#include "uw.h"
#include "uvw1.h"
#include "rk4.h"
#include "hooke.h"

extern int	flag_print;
extern int	debug;
extern int 	Nphi;
extern int 	Iterations;
extern double 	OMEGAR;
extern double 	OMEGAI;
extern double 	K;
extern double 	R;
extern double	ALPHA, TIME;
extern int	RNUM, TIMENUM, ALPHANUM;
extern void (*ode)();

double scale (double a, double b, double A, double B, double x) {
	return (x*(A-B)/(a-b)+(a*B-b*A)/(a-b));
}

double shoot(void) {
	double x;
	Iterations++;
	switch (Nphi) {
	case (-1): return shootv();
	case (0): return shootuw();
	case (1): return shootuvw1();
	default: exit(1);
	}
}

void single(double wrmin, double wrmax, double wimin, double wimax, int Nw) {
	int 	i,j;
	double 	f;
	for (i=0; i<Nw; i++) {
		OMEGAI = scale(0,Nw-1,wimin,wimax,i);
		for (j=0; j<Nw; j++) {
			OMEGAR = scale(0,Nw-1,wrmin,wrmax,j);
			f = 0.0;
			f = shoot();
			if (!debug)	fwrite(&f,sizeof(double),1,stdout);
		}
	}
}

void inter(double wrmin, double wrmax, double wimin, double wimax, int Nw, int Nk) {
	int 	i, j, k, kmax;
	double 	fmin, f;
	kmax = K;
	for (i=0; i<Nw; i++) {
		OMEGAI = scale(0,Nw-1,wimin,wimax,i);
		for (j=0; j<Nw; j++) {
			OMEGAR = scale(0,Nw-1,wrmin,wrmax,j);
			fmin = 999999999.0;
			for (k=0; k<Nk; k++) {
				K = scale(0, Nw-1, 0.0, kmax, k);       // k in [0,kmax]
				f = shoot();
				if (f<fmin)     fmin = f;
			}
			//printf("%g\n",fmin);
			if (!debug)	fwrite(&fmin,sizeof(double),1,stdout);
		}
	}
}

double shootfunction(double x[2], int n) {
	double y;
	OMEGAR=x[0];
	OMEGAI=x[1];
	return shoot();
}

void plot(double rho, int Nk) {
	int i, iters, itermax;
	double startpt[2];
	double endpt[2];
	double epsilon;
	double dk, kmax, bigk;
	double error; 
	double lastomegar, lastomegai, lastk, llastomegar, llastomegai, llastk, delta, omegarguess, omegaiguess;
	double quit = 0;
	double tol = 1.5;		/* reduce stepsize if OMEGA-OMEGAold > tol */
	double dkmin = 1E-4;		/* break if stepsize < dkmin */
	double w1, k1, w2, k2, wrhat, wihat, khat, shifttime; /* w1,k1 (start of instability), w2,k2 (end of instability), w[ri]hat, khat (maxima), shifttime (shifted +0.5, periodic continuation) */
	double bottomband, topband, sideband;
	int region = 0; 		/* 0 start (stable), 1 (unstable), 2 (reentrance of stability) */
	int first, instability;
	int breakstatus = 0;	/* 0 nobreak, 1 stepsize, 2 bottom, 3 top, 4 side, 5 kmax, 6 ok*/
	static char *breakstr[] = {"nobreak","stepsize","bottom","top","side","kmax","reentrance"};
	
	OMEGAR = 0.0;
	if (Nphi == -1) OMEGAI = -14.682;
	else if (Nphi == 0) OMEGAI = -26.37;
	else if (Nphi == 1) OMEGAI= -14.682;
	else { fprintf(stderr,"starting value unknown\n"); exit(1); }
	if (rho == -1.0) rho = 0.5;
	if ((rho>1.0)||(rho<0.0)) { fprintf(stderr,"error: -r rho must be in [0,1]\n"); exit(1); }

	if (!Nk)	Nk = 50;
	epsilon = 1E-6;
	itermax = 3000;
	if ( K>=0.0 ) kmax = K;
	else {
		K = 100.0;
		kmax = 100.0*R;
	}
	bottomband = -150.0;
	topband = 300.0;
	sideband = 10.0*R;
	dk = K/(double)Nk;
	K = 0.0;
	printf("#shootplot: n=%d R=%d (%g) alpha=%d (%g) t=%d (%g) rho=%g Nk=%d kmax=%g\n",Nphi,RNUM,R,ALPHANUM,ALPHA,TIMENUM,TIME,rho,Nk,kmax,R);
	printf("%g %g %g\n",0.0,OMEGAR,OMEGAI);
	wrhat = 0.0;
	wihat = OMEGAI;
	khat = 0.0;
	quit = 0;
	llastk = 0.0;
	lastk = 0.0;
	llastomegar = OMEGAR;
	llastomegai = OMEGAI;
	lastomegar = OMEGAR;
	lastomegai = OMEGAI;
	omegarguess = OMEGAR;
	omegaiguess = OMEGAI;
	first = 1;
	instability = 0;
	while (1) {
		error = 2.0 * tol;
		bigk = K;
		dk *= 2.0;
		while (error > tol) {
			dk /= 2.0;
			if (dk < dkmin) {
				breakstatus = 1;
				break;
			}
			if (!first) {
				delta = dk/(lastk-llastk);
				omegarguess = (1.0 + delta)*lastomegar - delta*llastomegar;
				omegaiguess = (1.0 + delta)*lastomegai - delta*llastomegai;
				// fprintf(stderr,"guess: dk=%g llk=%g lk=%g delta=%g llastw=(%g,%g) lastw=(%g,%g) -> (%g,%g)\n",dk,llastk,lastk,delta,llastomegar,llastomegai,lastomegar,lastomegai, omegarguess, omegaiguess);
			}
			K = bigk + dk;
			startpt[0] = omegarguess;
			startpt[1] = omegaiguess;
			Iterations=0;
			iters = hooke(2, startpt, endpt, rho, epsilon, itermax, shootfunction);	
			error = hypot(omegarguess-endpt[0],omegaiguess-endpt[1]);
			// fprintf(stderr,"stepsize control: dk=%g, error=%g, K=%g w=(%g,%g) wstart=(%g,%g)\n",dk,error,K,endpt[0],endpt[1],startpt[0],startpt[1]);
		}
		if (error<tol/3.0) dk *= 2.0;
		if (breakstatus)	break;
		if (first) first = 0; 
		llastomegar = lastomegar;
		llastomegai = lastomegai;
		llastk = lastk;
		lastomegar = OMEGAR;
		lastomegai = OMEGAI;
		lastk = K;
		OMEGAR = endpt[0];
		OMEGAI = endpt[1];
		if (OMEGAI >= 0.0) instability = 1;
		if (OMEGAI > wihat) {
			wrhat = OMEGAR;
			wihat = OMEGAI;
			khat = K;
		}
		if ((!region)&&(OMEGAI >= 0.0)) {
			w1 = llastomegar + llastomegai*(llastomegar-OMEGAR)/(OMEGAI-llastomegai);
			k1 = llastk + llastomegai*(llastk-K)/(OMEGAI-llastomegai);
			region = 1;
		} else if ((region==1)&&(OMEGAI<=0.0)) {
			w2 = llastomegar + llastomegai*(llastomegar-OMEGAR)/(OMEGAI-llastomegai);
			k2 = llastk + llastomegai*(llastk-K)/(OMEGAI-llastomegai);
			region = 2;
			breakstatus = 6;
		}
		printf("%g %g %g\n",K,OMEGAR,OMEGAI);
		if (OMEGAI <= bottomband) { breakstatus = 2; break; } 
		if (OMEGAI >= topband) { breakstatus = 3; break; } 
		if (fabs(OMEGAR) >= sideband) { breakstatus = 4; break; } 
		if ((K>kmax)&&(region!=1)) { breakstatus = 5; break; }
		if (breakstatus) break;
	}
	if ((shifttime = TIME+0.5) > 1.0) shifttime -= 1.0;
	printf("#status instability=%d break=%s region=%d khat=%g wrhat=%g wihat=%g k1=%g w1=%g k2=%g w2=%g R=%g (%d) T=%g (%d) alpha=%g (%d) shifttime=%g\n",instability,breakstr[breakstatus],region,khat,wrhat,wihat,k1,w1,k2,w2,R,RNUM,TIME,TIMENUM,ALPHA,ALPHANUM,shifttime);
}



struct tkrw {
	double k;
	double wr;
	double wi;
};

int kmaxima(struct tkrw *start, struct tkrw *middle, struct tkrw *end, int depth) {
	int i, iters, itermax = 3000;
	struct tkrw left, right, lleft, center, rright;
	double startpt[2], endpt[2], rho = 0.5, epsilon = 1E-6, errorr, errori;

	if (middle->k == 0.0) {
		start->k = 0.0;
		middle->k = 0.0;
		end->k = 0.0;
		start->wr = 0.0;
		middle->wr = 0.0;
		end->wr = 0.0;
		start->wi = -26.374387;
		middle->wi = -26.374387;
		end->wi = -26.374387;
		return 0;
	}
	lleft.k = start->k;
	lleft.wr = start->wr;
	lleft.wi = start->wi;
	center.k = middle->k;
	center.wr = middle->wr;
	center.wi = middle->wi;
	rright.k = end->k;
	rright.wr = end->wr;
	rright.wi = end->wi;
	for (i=0; i<depth+5; i++) {
		left.k = (lleft.k + center.k)/2.0;
		left.wr = (lleft.wr + center.wr)/2.0;
		left.wi = (lleft.wi + center.wi)/2.0;
		right.k = (rright.k + center.k)/2.0;
		right.wr = (rright.wr + center.wr)/2.0;
		right.wi = (rright.wi + center.wi)/2.0;

		K = left.k;
		startpt[0] = left.wr;
		startpt[1] = left.wi;
	// fprintf(stderr,"left: t=%g k=%g w=(%g %g) -> ",TIME,K,startpt[0],startpt[1]);
		iters = hooke(2, startpt, endpt, rho, epsilon, itermax, shootfunction);	
		/*
		errorr = (startpt[0] - endpt[0])/(startpt[0]+endpt[0]);
		errori = (startpt[1] - endpt[1])/(startpt[1]+endpt[1]);
		*/
		errorr = 0.0;
		errori = 0.0;
		if (hypot(errorr,errori)>0.2) { 
			fprintf (stderr,"\n too far away (%g %g)->(%g %g)\n",startpt[0],startpt[1],endpt[0],endpt[1]); 
			left.wr = startpt[0];
			left.wi = startpt[1];
		} else {
			left.wr = endpt[0];
			left.wi = endpt[1];
		}
	// fprintf(stderr,"(%g %g) i=%d (%d)\n",endpt[0],endpt[1],i,depth);

		K = right.k;
		startpt[0] = right.wr;
		startpt[1] = right.wi;
	// fprintf(stderr,"right: t=%g k=%g w=(%g %g) -> ",TIME,K,startpt[0],startpt[1]);
		iters = hooke(2, startpt, endpt, rho, epsilon, itermax, shootfunction);	
		/*
		errorr = (startpt[0] - endpt[0])/(startpt[0]+endpt[0]);
		errori = (startpt[1] - endpt[1])/(startpt[1]+endpt[1]);
		*/
		errorr = 0.0;
		errori = 0.0;
		if (hypot(errorr,errori)>0.2) { 
			fprintf (stderr,"\n too far away (%g %g)->(%g %g)\n",startpt[0],startpt[1],endpt[0],endpt[1]); 
			right.wr = startpt[0];
			right.wi = startpt[1];
		} else {
			right.wr = endpt[0];
			right.wi = endpt[1];
		}
	// fprintf(stderr,"(%g %g)\n",endpt[0],endpt[1]);

		if ((center.wi > left.wi)&&(center.wi > right.wi)) {
			lleft.k = left.k;
			lleft.wr = left.wr;
			lleft.wi = left.wi;
			rright.k = right.k;
			rright.wr = right.wr;
			rright.wi = right.wi;
		} else if (left.wi > center.wi) {
			rright.k = center.k;
			rright.wr = center.wr;
			rright.wi = center.wi;
			center.k = left.k;
			center.wr = left.wr;
			center.wi = left.wi;
		} else if (right.wi > center.wi) {
			lleft.k = center.k;
			lleft.wr = center.wr;
			lleft.wi = center.wi;
			center.k = right.k;
			center.wr = right.wr;
			center.wi = right.wi;
		} else {
			fprintf(stderr,"WARNING: kmaxima (else)\n");
		}
		fprintf(stderr,"%d [%d], band=%g\n",i,depth,fabs(lleft.k-rright.k));
		if ( (fabs(lleft.k - rright.k)) < 1E-4 ) break;
		
	}
	if ( (fabs(start->k-center.k) < 0.1)||(fabs(end->k-center.k) < 0.1) ) {
		fprintf(stderr,"WARNING kmaxima (near edge)\n");
	}
	start->k = lleft.k;
	start->wr = lleft.wr;
	start->wi = lleft.wi;
	middle->k = center.k;
	middle->wr = center.wr;
	middle->wi = center.wi;
	end->k = rright.k;
	end->wr = rright.wr;
	end->wi = rright.wi;
	fprintf(stderr,"==> %g %g (%g %g)\n",TIME,middle->k,middle->wr,middle->wi);
	return 0;
}

struct timekrw {
	double time;
	struct tkrw left, center, right;
};

void timekrwcpy(struct timekrw *dest, struct timekrw *src) {
	dest->time = src->time;
	dest->left.k = src->left.k;
	dest->left.wr = src->left.wr;
	dest->left.wi = src->left.wi;
	dest->center.k = src->center.k;
	dest->center.wr = src->center.wr;
	dest->center.wi = src->center.wi;
	dest->right.k = src->right.k;
	dest->right.wr = src->right.wr;
	dest->right.wi = src->right.wi;
}

void tmaxima(double t[3], double k[9], double wr[9], double wi[9]) {
	int depth = 0;
	double dt, shifttime, epsilon_1;
	struct timekrw prev, prior, middle, later, past; 

	prev.time = t[0];
	prev.left.k = k[0]; 	prev.center.k = k[1]; 	prev.right.k = k[2];
	prev.left.wr = wr[0]; 	prev.center.wr = wr[1]; prev.right.wr = wr[2];
	prev.left.wi = wi[0]; 	prev.center.wi = wi[1]; prev.right.wi = wi[2];
	TIME = prev.time;
	kmaxima(&prev.left, &prev.center, &prev.right, depth);

	middle.time = t[1];
	middle.left.k = k[3]; 	middle.center.k = k[4]; 	middle.right.k = k[5];
	middle.left.wr = wr[3]; middle.center.wr = wr[4]; 	middle.right.wr = wr[5];
	middle.left.wi = wi[3]; middle.center.wi = wi[4]; 	middle.right.wi = wi[5];
	TIME = middle.time;
	kmaxima(&middle.left, &middle.center, &middle.right, depth);

	past.time = t[2];
	past.left.k = k[6]; 	past.center.k = k[7]; 	past.right.k = k[8];
	past.left.wr = wr[6]; 	past.center.wr = wr[7]; past.right.wr = wr[8];
	past.left.wi = wi[6]; 	past.center.wi = wi[7]; past.right.wi = wi[8];
	TIME = past.time;
	kmaxima(&past.left, &past.center, &past.right, depth);

	dt = 1.0;
	while ((dt > 1E-4)||(depth<5)) {
		prior.time = (prev.time + middle.time)/2.0;
		prior.left.k = (prev.left.k + middle.left.k)/2.0;
		prior.left.wr = (prev.left.wr + middle.left.wr)/2.0;
		prior.left.wi = (prev.left.wi + middle.left.wi)/2.0;
		prior.center.k = (prev.center.k + middle.center.k)/2.0;
		prior.center.wr = (prev.center.wr + middle.center.wr)/2.0;
		prior.center.wi = (prev.center.wi + middle.center.wi)/2.0;
		prior.right.k = (prev.right.k + middle.right.k)/2.0;
		prior.right.wr = (prev.right.wr + middle.right.wr)/2.0;
		prior.right.wi = (prev.right.wi + middle.right.wi)/2.0;

		later.time = (past.time + middle.time)/2.0;
		later.left.k = (past.left.k + middle.left.k)/2.0;
		later.left.wr = (past.left.wr + middle.left.wr)/2.0;
		later.left.wi = (past.left.wi + middle.left.wi)/2.0;
		later.center.k = (past.center.k + middle.center.k)/2.0;
		later.center.wr = (past.center.wr + middle.center.wr)/2.0;
		later.center.wi = (past.center.wi + middle.center.wi)/2.0;
		later.right.k = (past.right.k + middle.right.k)/2.0;
		later.right.wr = (past.right.wr + middle.right.wr)/2.0;
		later.right.wi = (past.right.wi + middle.right.wi)/2.0;

		TIME = prior.time;
		reinitbaseflow();
		kmaxima(&prior.left,&prior.center,&prior.right, depth);

		TIME = later.time;
		reinitbaseflow();
		kmaxima(&later.left,&later.center,&later.right, depth);
	
		if ((middle.center.wi > prior.center.wi)&&(middle.center.wi > later.center.wi)) {
			timekrwcpy(&prev,&prior);
			timekrwcpy(&past,&later);
		} else if (prior.center.wi > middle.center.wi) {
			timekrwcpy(&past,&middle);
			timekrwcpy(&middle,&prior);
		} else if (later.center.wi > middle.center.wi) {
			timekrwcpy(&prev,&middle);
			timekrwcpy(&middle,&later);
		} else {
			fprintf(stderr,"WARNING tmaxima (else)\n");
		}
		dt = fabs(prev.time - past.time);
		fprintf(stderr,"===> dt=%g depth=%d [%g %g (%g %g)]<===\n",dt,depth,middle.time,middle.center.k,middle.center.wr,middle.center.wi);
		depth++;
	}

	TIME = middle.time;
	K = middle.center.k;
	OMEGAR = middle.center.wr;
	OMEGAI = middle.center.wi;
	fprintf(stderr,"t=%g k=%g w=(%g %g)\n",TIME,K,OMEGAR,OMEGAI);
	if ((shifttime = TIME+0.5) > 1.0) shifttime -= 1.0;
	epsilon_1 = hypot(OMEGAR,OMEGAI)/(ALPHA*ALPHA);
	printf("%g %g %g %g %g %g %g %g %g\n",R,ALPHA,TIME,shifttime,K/R,K,OMEGAR,OMEGAI,epsilon_1);
}
