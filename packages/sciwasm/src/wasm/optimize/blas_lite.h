/**
 * Minimal BLAS/LAPACK declarations for L-BFGS-B.
 *
 * Only the routines actually called by lbfgsb.c are declared here.
 * These use Fortran-style interfaces (pass-by-pointer).
 */

#ifndef BLAS_LITE_H
#define BLAS_LITE_H

/* BLAS Level 1 */
void   daxpy_(int* n, double* alpha, double* x, int* incx, double* y, int* incy);
void   dscal_(int* n, double* alpha, double* x, int* incx);
void   dcopy_(int* n, double* x, int* incx, double* y, int* incy);
double dnrm2_(int* n, double* x, int* incx);
double ddot_(int* n, double* x, int* incx, double* y, int* incy);

/* LAPACK */
void   dpotrf_(char* uplo, int* n, double* a, int* lda, int* info);
void   dtrtrs_(char* uplo, char* trans, char* diag, int* n, int* nrhs,
               double* a, int* lda, double* b, int* ldb, int* info);

#endif /* BLAS_LITE_H */
