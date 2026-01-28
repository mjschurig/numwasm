/**
 * Minimal BLAS/LAPACK implementations for L-BFGS-B.
 *
 * These are simple reference implementations of the routines needed
 * by lbfgsb.c. They use Fortran-style interfaces (pass-by-pointer)
 * to match the calling convention expected by the L-BFGS-B code.
 */

#include "blas_lite.h"
#include <math.h>

/* ============================== */
/* BLAS Level 1                   */
/* ============================== */

/**
 * daxpy: y := alpha*x + y
 */
void daxpy_(int* n, double* alpha, double* x, int* incx, double* y, int* incy) {
    int nn = *n;
    double a = *alpha;
    int ix = *incx;
    int iy = *incy;

    if (nn <= 0 || a == 0.0) return;

    if (ix == 1 && iy == 1) {
        for (int i = 0; i < nn; i++) {
            y[i] += a * x[i];
        }
    } else {
        int jx = ix > 0 ? 0 : -(nn - 1) * ix;
        int jy = iy > 0 ? 0 : -(nn - 1) * iy;
        for (int i = 0; i < nn; i++) {
            y[jy] += a * x[jx];
            jx += ix;
            jy += iy;
        }
    }
}

/**
 * dscal: x := alpha*x
 */
void dscal_(int* n, double* alpha, double* x, int* incx) {
    int nn = *n;
    double a = *alpha;
    int ix = *incx;

    if (nn <= 0 || ix <= 0) return;

    if (ix == 1) {
        for (int i = 0; i < nn; i++) {
            x[i] *= a;
        }
    } else {
        for (int i = 0; i < nn * ix; i += ix) {
            x[i] *= a;
        }
    }
}

/**
 * dcopy: y := x
 */
void dcopy_(int* n, double* x, int* incx, double* y, int* incy) {
    int nn = *n;
    int ix = *incx;
    int iy = *incy;

    if (nn <= 0) return;

    if (ix == 1 && iy == 1) {
        for (int i = 0; i < nn; i++) {
            y[i] = x[i];
        }
    } else {
        int jx = ix > 0 ? 0 : -(nn - 1) * ix;
        int jy = iy > 0 ? 0 : -(nn - 1) * iy;
        for (int i = 0; i < nn; i++) {
            y[jy] = x[jx];
            jx += ix;
            jy += iy;
        }
    }
}

/**
 * dnrm2: ||x||_2
 */
double dnrm2_(int* n, double* x, int* incx) {
    int nn = *n;
    int ix = *incx;

    if (nn < 1 || ix < 1) return 0.0;
    if (nn == 1) return fabs(x[0]);

    double scale = 0.0;
    double ssq = 1.0;

    int j = 0;
    for (int i = 0; i < nn; i++) {
        if (x[j] != 0.0) {
            double absxi = fabs(x[j]);
            if (scale < absxi) {
                ssq = 1.0 + ssq * (scale / absxi) * (scale / absxi);
                scale = absxi;
            } else {
                ssq += (absxi / scale) * (absxi / scale);
            }
        }
        j += ix;
    }
    return scale * sqrt(ssq);
}

/**
 * ddot: x^T * y
 */
double ddot_(int* n, double* x, int* incx, double* y, int* incy) {
    int nn = *n;
    int ix = *incx;
    int iy = *incy;

    if (nn <= 0) return 0.0;

    double s = 0.0;
    if (ix == 1 && iy == 1) {
        for (int i = 0; i < nn; i++) {
            s += x[i] * y[i];
        }
    } else {
        int jx = ix > 0 ? 0 : -(nn - 1) * ix;
        int jy = iy > 0 ? 0 : -(nn - 1) * iy;
        for (int i = 0; i < nn; i++) {
            s += x[jx] * y[jy];
            jx += ix;
            jy += iy;
        }
    }
    return s;
}


/* ============================== */
/* LAPACK                         */
/* ============================== */

/**
 * dpotrf: Cholesky factorization of a real symmetric positive definite matrix.
 *
 * Computes the Cholesky factorization of a real symmetric positive
 * definite matrix A: A = L * L^T  (lower triangular).
 *
 * Only the lower triangle is used/modified when uplo = 'L'.
 *
 * @param uplo  'U' or 'L' — only 'L' is used by lbfgsb.c
 * @param n     Order of the matrix
 * @param a     On entry: the symmetric matrix. On exit: the factor L.
 *              Column-major, stride lda.
 * @param lda   Leading dimension of a
 * @param info  0 = success, j > 0 = not positive definite (leading minor j)
 */
void dpotrf_(char* uplo, int* n, double* a, int* lda, int* info) {
    int nn = *n;
    int ld = *lda;
    *info = 0;

    if (nn <= 0) return;

    /* Lower triangular Cholesky: A = L * L^T */
    if (*uplo == 'L' || *uplo == 'l') {
        for (int j = 0; j < nn; j++) {
            double sum = a[j + ld * j];
            for (int k = 0; k < j; k++) {
                sum -= a[j + ld * k] * a[j + ld * k];
            }
            if (sum <= 0.0) {
                *info = j + 1;
                return;
            }
            a[j + ld * j] = sqrt(sum);
            double ajj = a[j + ld * j];

            for (int i = j + 1; i < nn; i++) {
                double s = a[i + ld * j];
                for (int k = 0; k < j; k++) {
                    s -= a[i + ld * k] * a[j + ld * k];
                }
                a[i + ld * j] = s / ajj;
            }
        }
    } else {
        /* Upper triangular Cholesky: A = U^T * U */
        for (int j = 0; j < nn; j++) {
            double sum = a[j + ld * j];
            for (int k = 0; k < j; k++) {
                sum -= a[k + ld * j] * a[k + ld * j];
            }
            if (sum <= 0.0) {
                *info = j + 1;
                return;
            }
            a[j + ld * j] = sqrt(sum);
            double ajj = a[j + ld * j];

            for (int i = j + 1; i < nn; i++) {
                double s = a[j + ld * i];
                for (int k = 0; k < j; k++) {
                    s -= a[k + ld * j] * a[k + ld * i];
                }
                a[j + ld * i] = s / ajj;
            }
        }
    }
}

/**
 * dtrtrs: Solve a triangular system of equations.
 *
 * Solves A * X = B  or  A^T * X = B,
 * where A is a triangular matrix and B is a general matrix.
 *
 * @param uplo   'U' for upper triangular, 'L' for lower triangular
 * @param trans  'N' for no transpose, 'T' for transpose
 * @param diag   'N' for non-unit diagonal, 'U' for unit diagonal
 * @param n      Order of the matrix A
 * @param nrhs   Number of right-hand sides (columns of B)
 * @param a      The triangular matrix (column-major, stride lda)
 * @param lda    Leading dimension of a
 * @param b      On entry: right-hand side. On exit: solution. (stride ldb)
 * @param ldb    Leading dimension of b
 * @param info   0 = success, -i = i-th argument error, j > 0 = singular
 */
void dtrtrs_(char* uplo, char* trans, char* diag, int* n, int* nrhs,
             double* a, int* lda, double* b, int* ldb, int* info) {
    int nn = *n;
    int nr = *nrhs;
    int ld_a = *lda;
    int ld_b = *ldb;
    int lower = (*uplo == 'L' || *uplo == 'l');
    int notrans = (*trans == 'N' || *trans == 'n');
    int nounit = (*diag == 'N' || *diag == 'n');

    *info = 0;

    /* Check for singularity */
    if (nounit) {
        for (int j = 0; j < nn; j++) {
            if (a[j + ld_a * j] == 0.0) {
                *info = j + 1;
                return;
            }
        }
    }

    /* Solve for each right-hand side */
    for (int k = 0; k < nr; k++) {
        double* bk = &b[ld_b * k];

        if (lower && notrans) {
            /* L * x = b : forward substitution */
            for (int i = 0; i < nn; i++) {
                double s = bk[i];
                for (int j = 0; j < i; j++) {
                    s -= a[i + ld_a * j] * bk[j];
                }
                bk[i] = nounit ? s / a[i + ld_a * i] : s;
            }
        } else if (lower && !notrans) {
            /* L^T * x = b : backward substitution */
            for (int i = nn - 1; i >= 0; i--) {
                double s = bk[i];
                for (int j = i + 1; j < nn; j++) {
                    s -= a[j + ld_a * i] * bk[j];
                }
                bk[i] = nounit ? s / a[i + ld_a * i] : s;
            }
        } else if (!lower && notrans) {
            /* U * x = b : backward substitution */
            for (int i = nn - 1; i >= 0; i--) {
                double s = bk[i];
                for (int j = i + 1; j < nn; j++) {
                    s -= a[i + ld_a * j] * bk[j];
                }
                bk[i] = nounit ? s / a[i + ld_a * i] : s;
            }
        } else {
            /* U^T * x = b : forward substitution */
            for (int i = 0; i < nn; i++) {
                double s = bk[i];
                for (int j = 0; j < i; j++) {
                    s -= a[j + ld_a * i] * bk[j];
                }
                bk[i] = nounit ? s / a[i + ld_a * i] : s;
            }
        }
    }
}
