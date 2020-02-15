/*
 * amoswrap
 */
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <g2c.h>


/*
 * *_exp: are exponentially scaled versions (see amos sources for details; KODE=2)
 */

/* 
 * Bessel functions of first kind Jn, Jn_exp, J0, J1
 */
void Jn(double v, double zr, double zi, double *yr, double *yi) {
	double yreal, yimag;
	int nz, ierr;
	int n = 1;
	int kode = 1;

	zbesj_(&zr, &zi, &v, &kode, &n, &yreal, &yimag, &nz, &ierr);
	*yr = yreal;
	*yi = yimag;
}
void Jn_exp(double v, double zr, double zi, double *yr, double *yi) {
	double yreal, yimag;
	int nz, ierr;
	int n = 1;
	int kode = 2;

	zbesj_(&zr, &zi, &v, &kode, &n, &yreal, &yimag, &nz, &ierr);
	*yr = yreal;
	*yi = yimag;
}
void J0(double zr, double zi, double *yr, double *yi) {
	double yreal, yimag;
	int nz, ierr;
	int n = 1;
	int kode = 1;
	double v = 0.0;

	zbesj_(&zr, &zi, &v, &kode, &n, &yreal, &yimag, &nz, &ierr);
	*yr = yreal;
	*yi = yimag;
}
void J1(double zr, double zi, double *yr, double *yi) {
	double yreal, yimag;
	int nz, ierr;
	int n = 1;
	int kode = 1;
	double v = 1.0;

	zbesj_(&zr, &zi, &v, &kode, &n, &yreal, &yimag, &nz, &ierr);
	*yr = yreal;
	*yi = yimag;
}

/* 
 * Bessel functions of second kind (Weber functions, Neumann functions) Yn, Yn_exp
 */
void Yn(double v, double zr, double zi, double *yr, double *yi) {
	double yreal, yimag;
	double wr, wi;
	int nz, ierr;
	int n = 1;
	int kode = 1;

	zbesy_(&zr, &zi, &v, &kode, &n, &yreal, &yimag, &nz, &wr, &wi,&ierr);
	*yr = yreal;
	*yi = yimag;
}
void Yn_exp(double v, double zr, double zi, double *yr, double *yi) {
	double yreal, yimag;
	double wr, wi;
	int nz, ierr;
	int n = 1;
	int kode = 2;

	zbesy_(&zr, &zi, &v, &kode, &n, &yreal, &yimag, &nz, &wr, &wi, &ierr);
	*yr = yreal;
	*yi = yimag;
}

/*
 * Bessel functions of third kind (Hankel functions)
 * Hn1, Hn2, Hn1_exp, Hn2_exp
 * 
 *
 */
void Hn1(double v, double zr, double zi, double *yr, double *yi) {
	double yreal, yimag;
	int nz, ierr;
	int n = 1;
	int m = 1;
	int kode = 1;

	zbesh_(&zr, &zi, &v, &kode, &m, &n, &yreal, &yimag, &nz, &ierr);
	*yr = yreal;
	*yi = yimag;
}
void Hn2(double v, double zr, double zi, double *yr, double *yi) {
	double yreal, yimag;
	int nz, ierr;
	int n = 1;
	int m = 2;
	int kode = 1;

	zbesh_(&zr, &zi, &v, &kode, &m, &n, &yreal, &yimag, &nz, &ierr);
	*yr = yreal;
	*yi = yimag;
}
void Hn1_exp(double v, double zr, double zi, double *yr, double *yi) {
	double yreal, yimag;
	int nz, ierr;
	int n = 1;
	int m = 1;
	int kode = 2;

	zbesh_(&zr, &zi, &v, &kode, &m, &n, &yreal, &yimag, &nz, &ierr);
	*yr = yreal;
	*yi = yimag;
}
void Hn2_exp(double v, double zr, double zi, double *yr, double *yi) {
	double yreal, yimag;
	int nz, ierr;
	int n = 1;
	int m = 2;
	int kode = 2;

	zbesh_(&zr, &zi, &v, &kode, &m, &n, &yreal, &yimag, &nz, &ierr);
	*yr = yreal;
	*yi = yimag;
}
/*
 * modified Bessel functions In, Kn, In_exp, Kn_exp
 */
void In(double v, double zr, double zi, double *yr, double *yi) {
	double yreal, yimag;
	int nz, ierr;
	int n = 1;
	int kode = 1;

	zbesi_(&zr, &zi, &v, &kode, &n, &yreal, &yimag, &nz, &ierr);
	*yr = yreal;
	*yi = yimag;
}
void In_exp(double v, double zr, double zi, double *yr, double *yi) {
	double yreal, yimag;
	int nz, ierr;
	int n = 1;
	int kode = 2;

	zbesi_(&zr, &zi, &v, &kode, &n, &yreal, &yimag, &nz, &ierr);
	*yr = yreal;
	*yi = yimag;
}
void Kn(double v, double zr, double zi, double *yr, double *yi) {
	double yreal, yimag;
	int nz, ierr;
	int n = 1;
	int kode = 1;

	zbesk_(&zr, &zi, &v, &kode, &n, &yreal, &yimag, &nz, &ierr);
	*yr = yreal;
	*yi = yimag;
}
void Kn_exp(double v, double zr, double zi, double *yr, double *yi) {
	double yreal, yimag;
	int nz, ierr;
	int n = 1;
	int kode = 2;

	zbesk_(&zr, &zi, &v, &kode, &n, &yreal, &yimag, &nz, &ierr);
	*yr = yreal;
	*yi = yimag;
}
/*
 * Airy Functions Ai, Ai1, Bi, Bi1, Ai_exp, Ai1_exp, Bi_exp, Bi1_exp
 *
 * ODE:		w'' - z*w = 0
 * Linear independent solutions:	Ai(z),Bi(z)
 *
 * Ai1 := Ai'
 * Bi1 := Bi'
 *
 *
 */
void Ai(double zr, double zi, double *yr, double *yi) {
	double yreal, yimag;
	int nz, ierr;
	int kode = 1;
	int id = 0;
	zairy_(&zr, &zi, &id, &kode, &yreal, &yimag, &nz, &ierr);
	*yr = yreal;
	*yr = yimag;
}
void Ai1(double zr, double zi, double *yr, double *yi) {
	double yreal, yimag;
	int nz, ierr;
	int kode = 1;
	int id = 1;
	zairy_(&zr, &zi, &id, &kode, &yreal, &yimag, &nz, &ierr);
	*yr = yreal;
	*yr = yimag;
}
void Bi(double zr, double zi, double *yr, double *yi) {
	double yreal, yimag;
	int nz, ierr;
	int n = 1;
	int kode = 1;
	int id = 0;
	zbiry_(&zr, &zi, &id, &kode, &yreal, &yimag, &nz, &ierr);
	*yr = yreal;
	*yr = yimag;
}
void Bi1(double zr, double zi, double *yr, double *yi) {
	double yreal, yimag;
	int nz, ierr;
	int n = 1;
	int kode = 1;
	int id = 1;
	zbiry_(&zr, &zi, &id, &kode, &yreal, &yimag, &nz, &ierr);
	*yr = yreal;
	*yr = yimag;
}
void Ai_exp(double zr, double zi, double *yr, double *yi) {
	double yreal, yimag;
	int nz, ierr;
	int n = 1;
	int kode = 2;
	int id = 0;
	zairy_(&zr, &zi, &id, &kode, &yreal, &yimag, &nz, &ierr);
	*yr = yreal;
	*yr = yimag;
}
void Ai1_exp(double zr, double zi, double *yr, double *yi) {
	double yreal, yimag;
	int nz, ierr;
	int n = 1;
	int kode = 2;
	int id = 1;
	zairy_(&zr, &zi, &id, &kode, &yreal, &yimag, &nz, &ierr);
	*yr = yreal;
	*yr = yimag;
}
void Bi_exp(double zr, double zi, double *yr, double *yi) {
	double yreal, yimag;
	int nz, ierr;
	int n = 1;
	int kode = 2;
	int id = 0;
	zbiry_(&zr, &zi, &id, &kode, &yreal, &yimag, &nz, &ierr);
	*yr = yreal;
	*yr = yimag;
}
void Bi1_exp(double zr, double zi, double *yr, double *yi) {
	double yreal, yimag;
	int nz, ierr;
	int n = 1;
	int kode = 2;
	int id = 1;
	zbiry_(&zr, &zi, &id, &kode, &yreal, &yimag, &nz, &ierr);
	*yr = yreal;
	*yr = yimag;
}

/*
 * logarithm Log
 */
void Log(double zr, double zi, double *yr, double *yi) {
	double yreal, yimag;
	int ierr;
	azlog_(&zr, &zi, &yreal, &yimag, &ierr);
	*yr = yreal;
	*yi = yimag;
}
