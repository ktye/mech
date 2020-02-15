/* $Id: la.h,v 1.1 2005/04/18 17:11:51 elmar Exp $ */
void printmatrix(doublecomplex *A, int n);
void printvector(doublecomplex *v, int n);
void matvec(doublecomplex *A, doublecomplex *b, doublecomplex *c, int N);
void matmul(doublecomplex *A, doublecomplex *B, doublecomplex *C, int N);
void eig(doublecomplex *A, doublecomplex *eigenvalue, doublecomplex *eigenvector, int n);
