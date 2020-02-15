/*
 * ode: 	y'[i] = f(y,x)
 * ode backend
 * 	euler 	- euler explicit
 * 	heun	- modified euler
 * 	rk2	- runge kutta 2nd order
 * 	rk4	- runge kutta 4th order
 * ARGUMENTS:
 * 	*y	array (initialised with initial values)
 * 	ny	*y dim
 *	f	function pointer to f(double *y, int i, int ny, double x)
 *	x0	start of x
 *	h	stepsize, 	N*h = 1 	
 * 	N	timesteps
 */

extern int flag_print;

/* Maximal number of equations */
#define	MAXORDER	12

void printvec(double *y, int n, double x) {
	int i;
	printf("%g ",x);
	for (i=0; i<n; i++)
		printf("%g ",y[i]);
	printf("\n");

}

void euler(double *y, int ny, double (*f)(), double x0, double h, int N) {
	int i,j;
	double x;
	double y_[MAXORDER];
	x = x0;
	if (flag_print)	printvec(y,ny,x);
	for (i=1; i<=N; i++) {
		for (j=0; j<ny; j++) y_[j] = y[j];
		for (j=0; j<ny; j++) {
			y[j] = y_[j] + h*f(y_, j, ny, x);
		}
		if (flag_print)	printvec(y,ny,x+h);
		x += h;
	}
}

void heun(double *y, int ny, double (*f)(), double x0, double h, int N) {
	int i,j;
	double x;
	double y_[MAXORDER];
	double yhf[MAXORDER];
	x = x0;
	if (flag_print)	printvec(y,ny,x);
	for (i=1; i<N; i++) {
		for (j=0; j<ny; j++) y_[j] = y[j];
		for (j=0; j<ny; j++) yhf[j] = y_[j] + h*f(y_, j, ny, x);
		for (j=0; j<ny; j++) {
			y[j] = y_[j] + 0.5*h*( f(y_, j, ny, x) + f(yhf, j, ny, x+h));
		}
		if (flag_print)	printvec(y,ny,x+h);
		x += h;
	}
}

void rk2(double *y, int ny, double (*f)(), double x0, double h, int N) {
	int i,j;
	double x;
	double k1[MAXORDER];
	double k2[MAXORDER];
	double y_[MAXORDER];
	x = x0;
	if (flag_print)	printvec(y,ny,x);
	for (i=1; i<N; i++) {
		for (j=0; j<ny; j++) y_[j] = y[j];
		for (j=0; j<ny; j++) k1[j] = y_[j] + 0.5*h*f(y_, j, ny, x);
		for (j=0; j<ny; j++) k2[j] = h*f(k1, j, ny, x+h/2.0);
		for (j=0; j<ny; j++) y[j] = y_[j] + k2[j];
		if (flag_print)	printvec(y,ny,x+h);
		x += h;
	}
}

void rk4(double *y, int ny, double (*f)(), double x0, double h, int N) {
	int i,j;
	double x;
	double k1[MAXORDER];
	double k2[MAXORDER];
	double k3[MAXORDER];
	double k4[MAXORDER];
	double yk1[MAXORDER];
	double yk2[MAXORDER];
	double yk3[MAXORDER];
	double y_[MAXORDER];
	x = x0;
	if (flag_print)	printvec(y,ny,x);
	for (i=1; i<N; i++) {
		for (j=0; j<ny; j++) y_[j] = y[j];
		for (j=0; j<ny; j++) k1[j] = h*f(y_, j, ny, x);
		for (j=0; j<ny; j++) yk1[j] = y_[j]+0.5*k1[j];
		for (j=0; j<ny; j++) k2[j] = h*f(yk1, j, ny, x+h/2.0);
		for (j=0; j<ny; j++) yk1[j] = y_[j]+0.5*k2[j];
		for (j=0; j<ny; j++) k3[j] = h*f(yk2, j, ny, x+h/2.0);
		for (j=0; j<ny; j++) yk3[j] = y_[j]+k3[j];
		for (j=0; j<ny; j++) k4[j] = h*f(yk3, j, ny, x+h);
		for (j=0; j<ny; j++) y[j] = y_[j] + k1[j]/6.0 + k2[j]/3.0 + k3[j]/3.0 + k4[j]/6.0;
		if (flag_print)	printvec(y,ny,x+h);
		x += h;
	}
}

