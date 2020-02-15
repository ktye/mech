/* $Id: puls.c,v 1.5 2003/12/11 18:07:55 elmar Exp $ */

/*
 * main - pulsatile pipe flow
 */

#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include <X11/Xlib.h>
#include "plot.h"
#include "amoswrap.h"
#include "complex.h"
#include "timer.h"

#define VECDIM 250
/*
int WIDTH=1280;
int HEIGHT=960;
*/

// void J0(double zr, double zi, double *yr, double *yi);


/* 
 * Dimensionslose Parameter
 *
 * alpha: 	Womersleyzahl
 * R0:		Reynoldszahl des stationären Anteils
 * R[n]:	Reynoldszahl der harmonischen Anteile
 * 
 * R[n] = G[n]*L^3 / (rho*nu^2)
 *
 * Skala:
 * Lange (Radius): 	L
 * Geschwindigkeit: 	U=L*omega
 * Zeit:		T=1/omega
 *
 * Druckgradient:	[rho*U^2/L]
 * Dichte:		rho
 * Zähigkeit:		nu
 *
 */
/* Strömungsprofil w=w(r,t) */

double scale (double a, double b, double A, double B, double x) {
	return (x*(A-B)/(a-b)+(a*B-b*A)/(a-b));
}

double w(double r, double t, double alpha, double *R, int NR, double *phase) {
	double Jr_r, Jr_i;
	double J_r, J_i;
	double JJ_r, JJ_i;
	double Wr, Wi;
	double exp_r, exp_i;
	double wr, wi;
	double W,abs,p;
	double N;
	int n,k;

	W = 0.0;
	*phase = 0.0;
	for (k=1; k<NR; k+=2) {
		n = (1+k)/2;
		N=(double)n;

		J0(sqrt(N)*alpha*M_SQRT1_2*r,sqrt(N)*alpha*M_SQRT1_2 *r,&Jr_r,&Jr_i);
		J0(sqrt(N)*alpha*M_SQRT1_2,sqrt(N)*alpha*M_SQRT1_2 ,&J_r,&J_i);
		cdiv(Jr_r,Jr_i,J_r,J_i,&JJ_r,&JJ_i);
		JJ_r -= 1.0;
		
		cmult(-R[k+1]/(N*pow(alpha,4.0)),R[k]/(N*pow(alpha,4.0)),JJ_r,JJ_i,&Wr,&Wi);

		exp_r = cos(N*t);
		exp_i = -sin(N*t);

		wr = 0.0;
		wi = 0.0;
		cmult(Wr,Wi,exp_r,exp_i,&wr,&wi);
		W += wr;
		cxyrp(wr, wi, &abs, &p);
		*phase += p;
	}

	return W+R[0]*(1.0-r*r)/(4.0*alpha*alpha);
}

double dwdr(double r, double t) {
	return 0;
}

double d2wdr2(double r, double t, double alpha, double *R, int NR) {
	double Jr_r, Jr_i;
	double J_r, J_i;
	double JJ_r, JJ_i;
	double J1r_r, J1r_i;
	double J1J0_r, J1J0_i;
	double Wr, Wi;
	double exp_r, exp_i;
	double wr, wi;
	double W,abs,p;
	double N;
	double Cr, Ci, CRr, CRi, CRJ1J0r, CRJ1J0i, RJJr,RJJi;
	int n,k;

	W = 0.0;
	for (k=1; k<NR; k+=2) {
		n = (1+k)/2;
		N = (double)n;
		J0(sqrt(N)*alpha*M_SQRT1_2*r,sqrt(N)*alpha*M_SQRT1_2 *r,&Jr_r,&Jr_i);

		J0(sqrt(N)*alpha*M_SQRT1_2,sqrt(N)*alpha*M_SQRT1_2 ,&J_r,&J_i);

		cdiv(Jr_r,Jr_i,J_r,J_i,&JJ_r,&JJ_i);

		J1(sqrt(N)*alpha*M_SQRT1_2*r,sqrt(N)*alpha*M_SQRT1_2 *r,&J1r_r,&J1r_i);

		cdiv(J1r_r,J1r_i,J_r,J_i,&J1J0_r,&J1J0_i);

		Cr = -M_SQRT1_2/(r*sqrt(N)*pow(alpha,3));
		Ci = M_SQRT1_2/(r*sqrt(N)*pow(alpha,3));

//		if (k%2) {
			cmult(Cr,Ci,R[k],R[k+1],&CRr,&CRi);
			cmult(CRr,CRi,J1J0_r,J1J0_i,&CRJ1J0r,&CRJ1J0i);
			cmult(R[k]/(alpha*alpha),R[k+1]/(alpha*alpha),JJ_r,JJ_i,&RJJr,&RJJi);
			Wr = CRJ1J0r + RJJr;
			Wi = CRJ1J0i + RJJi;
//		}

		exp_r = cos(N*t);
		exp_i = -sin(N*t);

		wr = 0.0;
		wi = 0.0;
		cmult(Wr,Wi,exp_r,exp_i,&wr,&wi);
		W += wr;
	}

	return W-R[0]/(2.0*alpha*alpha);
}

int inflexion(double y, int init) {
	int sign;
	static int last_sign = 0;
	int inflex = 0;

	if (y>0)	sign = 1;
	else		sign = 0;
	
	if (last_sign != sign)	inflex = 1;
	last_sign = sign;
	if (init) return 0;
	return inflex;
}

double flow(double t, double alpha, double *R, int N) {
	double r,dr = 0.01;
	double f=0.0;
	double phase=0.0;

	for (r=0.0; r<1.0; r+=dr) {
		f += r*dr*w(r, t, alpha, R, N, &phase);
	}
	return 2.0*M_PI*f;
}

double flow_fourier(double t, double alpha, double *R, int N) {
	int i,nn;
	double n;
	double U;
	double A=0.0,B=0.0,Cr=0.0,Ci=0.0;
	double J1r=0.0,J1i=0.0,J0r=0.0,J0i=0.0,JJr=0.0,JJi=0.0;

	U = R[0]/(16.0*alpha*alpha);
	for (i=1; i<N; i++) {
		nn = (i+1)/2;
		n = (double)nn;

		J1(sqrt(n)*M_SQRT1_2*alpha,sqrt(n)*M_SQRT1_2*alpha,&J1r,&J1i);
		J0(sqrt(n)*M_SQRT1_2*alpha,sqrt(n)*M_SQRT1_2*alpha,&J0r,&J0i);
		cdiv(J1r,J1i,J0r,J0i,&JJr,&JJi);

		cmult(M_SQRT1_2*sqrt(n)/(n*n*pow(alpha,5.0)),M_SQRT1_2*sqrt(n)/(n*n*pow(alpha,5.0)),JJr,JJi,&A,&B);
		B -= 0.5/(n*pow(alpha,4.0));
		//printf("n=%d, A=%f, B=%f\n",nn,A,B);



		if (i%2) {
			Cr = A*R[i] - B*R[i+1];
			U += Cr*cos(n*t); 
		} else	{
			Ci = -B*R[i-1]-A*R[i];
			U -= Ci*sin(n*t);
		}
	}
	return 2.0*M_PI*U;
}

double pressure_gradient(double t, double *R, int N) {
	int i,omega;
	double p,pi;

	p = R[0];
	for (i=1; i<N; i++) {
		omega = (1+i)/2;
		if (i%2)	p+=R[i]*cos(omega*t);
		else		p-=R[i]*sin(omega*t);
	}
	return p;
}


double infl(double alpha, double *R, int N) {
        double y;
        int *X;
        double dr=0.01;
        double dt;
        int k;
        int nmax;
        int T;
        int DT=100;
        int i,j,I,J;
        int infl;
        int inflexions=0;
        int init;
        double t,r;

        k=N;
        nmax = 1+k/2;
        I = nmax*DT;

        dt = 2.0*M_PI/((double)I-1.0);
        X = (int *)malloc((size_t)(I) * sizeof(int));

        J = 100;
        dr = 1.0/((double)J-1.0);


        t=0.0;
        for (i=0; i<I; i++) {
                t += dt;
                X[i]=0;
                r=0.0;
                for (j=0; j<J; j++) {
                        r+=dr;
                        if (r==0.0)     init = 1;
                        else            init = 0;
                        y=d2wdr2(r, t, alpha, R, N);
                        if (inflexion(y,init)) {
                                inflexions++;
                                break;
                        }
                }
        }
        return (double)inflexions/(double)I;
}

/* print velocity profile(r,t_i), pressure_gradient(t), mean_flow(t) 
 * to stdout, stderr, stderr. r: [-1,1], t: [0,2*pi] */
void profile(double alpha, double *R, int N) {
	int i,j,inflx,linemode;
	int Nx=512;
	int Nt=100;
	double x,y,t;
	double phase=0;
	char *var;

	if (var = getenv("steps"))	Nx=atoi(var);
	if (var = getenv("timesteps"))	Nt=atoi(var);

	for (j=0; j<Nt; j++) {
		for (i=0; i<Nx; i++) {
			x=scale(0,Nx-1,-1.0,1.0,i);
			t=scale(0,Nt-1,0,2.0*M_PI,j);
			y=w(x,t,alpha,R,N,&phase);
			printf("w:%f %f %f\n",t,x,y);
		}
		printf("w:\n");
	}
	printf("\n");
	for (i=0; i<Nx; i++) {
		t=scale(0,Nx-1,0,2.0*M_PI,i);
		y=pressure_gradient(t,R,N);
		printf("g:%f %f\n",t,y);
	}
	printf("\n");
	for (i=0; i<Nx; i++) {
		t=scale(0,Nx-1,0,2.0*M_PI,i);
		y=flow_fourier(t,alpha,R,N);
		printf("u:%f %f\n",t,y);
	}
	printf("\n");
	for (i=0; i<Nx; i++) {
		t = (double)i*2.0*M_PI/(Nx-1.0);
		inflx=0;
		for (x=0.0; x<=1.0001; x+=0.02) {
			y=d2wdr2(x,t,alpha,R,N);
			if (x==0.0)	linemode = 1;
			else 		linemode = 0;	
			if (inflexion(y,linemode)) 	inflx=1;
		}
		if (inflx) 	printf("i:%f 1\n",t);
		else		printf("i:%f 0\n",t);
	}
}


int main(int args, char **argv) {
	char *fgcolor = "green";
	char *bgcolor = "black";
	char *xlabel = "X";
	char *ylabel = "Y";
	char *title = "-- (puls) -- ";
	double xmin, xmax, ymin, ymax, time;
	int width, height;
	double x, y;
	int i,n,N;
	int linemode = 2;
	int linestyle = 1;
	int linewidth = 2;
	char linecolor[100];
	char buf[1024];
	char *var;
	double alpha;
	double t;
	double *R;
	double p,pmax,pi,wmax,w2max,f,ff,flowmax;
	double flowvec[VECDIM];
	double ffvec[VECDIM];
	int inflexionvec[VECDIM], inflx;
	double phase=0.0;
	
	strncpy (linecolor,"green",99);
	linecolor[99]='\0';

	width=0;
	height=0;

	if ((args<3) || ((args>3)&&!(args%2))) {
		fprintf(stderr,"usage: %s alpha R0 [R1.r R1.i] [R2.r R2.i] [...] [Rn.r Rn.i]]\n",argv[0]);
		exit(1);
	}
	else {
		alpha = atof(argv[1]);
		R = (double *)malloc((size_t)(args-2) * sizeof(double));
		N = args-2;
		for (i=2;i<args;i++)
			R[i-2] = atof(argv[i]);
	}

	/* only calculate inflexion ratio */
	if (var = getenv("inflexion")) {
		for (i=1; i<args; i++) printf("%s ",argv[i]);
		printf ("\t%f\n",infl(alpha, R, N));
		return 0;
	}

	/* print profile */
	if (var = getenv("profile")) {
		profile(alpha, R, N);
		return 0;
	}



	xmin=-1.0;
	xmax=1.0;
	ymin=-1.0;
	ymax=1.0;
	time=1.0;
//	if (var = getenv("xmin")) xmin=atof(var);
//	if (var = getenv("xmax")) xmax=atof(var);
	if (var = getenv("ymin")) ymin=atof(var);
	if (var = getenv("ymax")) ymax=atof(var);
	if (var = getenv("time")) time=atof(var);
	
	plotinit (xmin, xmax, ymin, ymax, fgcolor, bgcolor, xlabel, 
			ylabel, argv[0], width, height);

	pmax=0.0;
	pi = 0.0;
	for (i=0; i<N; i++) {
		if (i<2) pmax += fabs(R[i]);
		else {
			if (!(i%2)) pi=R[i]*R[i];
			else pmax += sqrt(pi + R[i]*R[i]);
		}	
	}
	flowmax=0.0;
	for (i=0; i<VECDIM; i++) {
		t = (double)i*2.0*M_PI/(VECDIM-1.0);
		f=flow(t,alpha,R,N);
		ff=flow_fourier(t,alpha,R,N);
		flowvec[i] = f;
		ffvec[i] = ff;
		if (fabs(f)>flowmax)	flowmax=fabs(f);
	}
	wmax=0.0;
	for (t=0.0; t<2.0*M_PI; t+=0.05) {
		for (x=xmin; x<=xmax+0.0001; x+=0.02) {
			y=w(x,t,alpha,R,N,&phase);
			if (fabs(y)>wmax)	wmax=fabs(y);
		}
	}
	wmax *= 1.1;

	w2max=0.0;
	for (t=0.0; t<2.0*M_PI; t+=0.05) {
		for (x=xmin; x<=xmax+0.0001; x+=0.02) {
			y=d2wdr2(x,t,alpha,R,N);
			if (fabs(y)>w2max)	w2max=fabs(y);
		}
	}
	for (i=0; i<VECDIM; i++) {
		t = (double)i*2.0*M_PI/(VECDIM-1.0);
		inflx=0;
		for (x=xmin; x<=xmax+0.0001; x+=0.02) {
			y=d2wdr2(x,t,alpha,R,N);
			if (x==xmin)	linemode = 1;
			else 		linemode = 0;	
			if (inflexion(y,linemode)) 	inflx=1;
		}
		if (inflx) 	inflexionvec[i] = 1;
		else		inflexionvec[i] = 0;
	}
	
	t=0.0;
	while (1) {
		t += time*elapsed_time();
		if (t>2.0*M_PI)	t-=2.0*M_PI;
	//	plotbackground (xmin, xmax, ymin, ymax, fgcolor, bgcolor, xlabel, ylabel, argv[0], width, height);
		plotbackground (xmin, xmax, -wmax, wmax, fgcolor, bgcolor, xlabel, ylabel, argv[0], width, height);
		for (x=xmin; x<=xmax+0.0001; x+=0.02) {
			y=w(x,t,alpha,R,N,&phase);
			if (x==xmin)	linemode = 1;
			else 		linemode = 0;	
		//	plot(x,y,xmin,xmax,ymin,ymax,linemode,linestyle,linewidth,"green");
			plot(x,y,xmin,xmax,-wmax,wmax,linemode,linestyle,linewidth,"green");
		}

		/* Print 2nd derivative of profile ********
		for (x=xmin; x<=xmax+0.0001; x+=0.02) {
			y=d2wdr2(x,t,alpha,R,N);
			if (x==xmin)	linemode = 1;
			else 		linemode = 0;	
			if (inflexion(y,linemode))	plotinflexion(x,xmin,xmax);
			plot(x,y*wmax/w2max,xmin,xmax,-wmax,wmax,linemode,linestyle,linewidth,"yellow");
		}
		*/

		for (x=0.0; x<=2.0*M_PI; x+=0.02) {
			p = pressure_gradient(x,R,N);
			if (x==0.0) plotpressure(x, p, 0.0, 2.0*M_PI, -pmax*1.1, pmax*1.1, 1);
			else	plotpressure(x, p, 0.0, 2.0*M_PI, -pmax*1.1, pmax*1.1, 0); 
		}
		for (i=0; i<VECDIM; i++) {
			x = (double)i*2.0*M_PI/(VECDIM-1.0);
			p = flowvec[i];
			if (x==0.0) plotflow(x, p, 0.0, 2.0*M_PI, -flowmax*1.1, flowmax*1.1, 1);
			else	plotflow(x, p, 0.0, 2.0*M_PI, -flowmax*1.1, flowmax*1.1, 0); 
		}	

		/* flowfourier */
		for (i=0; i<VECDIM; i++) {
			x = (double)i*2.0*M_PI/(VECDIM-1.0);
			p = ffvec[i];
			if (x==0.0) plotfflow(x, p, 0.0, 2.0*M_PI, -flowmax*1.1, flowmax*1.1, 1);
			else	plotfflow(x, p, 0.0, 2.0*M_PI, -flowmax*1.1, flowmax*1.1, 0); 
		}
		plot_inflexionline(inflexionvec,VECDIM);
		plot_time(t,0,2.0*M_PI);
		plotdisplay();
		usleep((useconds_t)10000);
		clearplot();
	}

	plotexit();

	return 0;
	// Window mode
	if ((width==0)&&(height==0)) sleep(2);
	return 0;
}
