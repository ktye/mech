/* $Id: symsolve.c,v 1.1 2006/02/17 18:29:48 elmar Exp $ */
/*
 * symsolve - solve real symmetric system
 *
 * A*x = b
 * A: symmetric NxN
 *
 * overwrites A and x
 */

#include <stdlib.h>
#include <g2c.h>
#include <clapack.h>
void symsolve(double *A, double *x, int N) {
	char	 uplo = 'U';
	integer  nrhs = 1;
	integer	*ipiv;
	integer  ldb;
	integer  info = 0;

	ipiv = (integer *)malloc(N*sizeof(integer));
	ldb = (integer)N;
	dspsv_(&uplo, (integer *)&N, &nrhs, (doublereal *)A, ipiv, (doublereal *)x, &ldb, &info);
}
