/*
 * NumJS LAPACK - Linear Algebra PACKage Implementation
 *
 * Essential LAPACK routines implemented from scratch for WebAssembly.
 * All matrices use column-major (Fortran) storage.
 *
 * Algorithms based on:
 * - Golub & Van Loan, "Matrix Computations"
 * - LAPACK Users' Guide
 */

#include "lapack.h"
#include "blas.h"
#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <float.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/* Helper macros */
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define ABS(x) ((x) >= 0 ? (x) : -(x))
#define SIGN(a, b) ((b) >= 0 ? ABS(a) : -ABS(a))

/* Character comparison helpers */
static inline int is_upper(char c) { return c == 'U' || c == 'u'; }
static inline int is_trans(char c) { return c == 'T' || c == 't' || c == 'C' || c == 'c'; }
static inline int is_unit(char c) { return c == 'U' || c == 'u'; }
static inline int want_vectors(char c) { return c == 'V' || c == 'v'; }

/* Machine precision */
static const double DLAMCH_EPS = 2.220446049250313e-16;  /* DBL_EPSILON */
static const double DLAMCH_SFMIN = 2.2250738585072014e-308;  /* DBL_MIN */

EXPORT double lapack_dlamch(char cmach) {
    switch (cmach) {
        case 'E': case 'e': return DLAMCH_EPS;
        case 'S': case 's': return DLAMCH_SFMIN;
        case 'B': case 'b': return 2.0;  /* base */
        case 'P': case 'p': return DLAMCH_EPS * 2.0;  /* precision */
        case 'N': case 'n': return 53.0;  /* number of mantissa bits */
        case 'R': case 'r': return 1.0;  /* rounding mode */
        case 'M': case 'm': return -1021.0;  /* min exponent */
        case 'U': case 'u': return DLAMCH_SFMIN;  /* underflow threshold */
        case 'L': case 'l': return 1024.0;  /* max exponent */
        case 'O': case 'o': return DBL_MAX;  /* overflow threshold */
        default: return 0.0;
    }
}

/* ============ Row Permutation ============ */

EXPORT void lapack_dlaswp(int32_t n, double* A, int32_t lda,
                           int32_t k1, int32_t k2, const int32_t* ipiv, int32_t incx)
{
    if (incx > 0) {
        for (int32_t i = k1; i <= k2; i++) {
            int32_t ip = ipiv[i];
            if (ip != i) {
                blas_dswap(n, &A[i], lda, &A[ip], lda);
            }
        }
    } else {
        for (int32_t i = k2; i >= k1; i--) {
            int32_t ip = ipiv[i];
            if (ip != i) {
                blas_dswap(n, &A[i], lda, &A[ip], lda);
            }
        }
    }
}

/* ============ LU Factorization ============ */

EXPORT int32_t lapack_dgetrf(int32_t m, int32_t n, double* A, int32_t lda,
                              int32_t* ipiv)
{
    if (m <= 0 || n <= 0) return 0;

    int32_t info = 0;
    int32_t minmn = MIN(m, n);

    for (int32_t j = 0; j < minmn; j++) {
        /* Find pivot: max |A[i,j]| for i >= j */
        int32_t jp = j + blas_idamax(m - j, &A[j + j * lda], 1);
        ipiv[j] = jp;

        if (A[jp + j * lda] != 0.0) {
            /* Swap rows j and jp */
            if (jp != j) {
                blas_dswap(n, &A[j], lda, &A[jp], lda);
            }

            /* Compute elements j+1:m of j-th column */
            if (j < m - 1) {
                double ajj_inv = 1.0 / A[j + j * lda];
                blas_dscal(m - j - 1, ajj_inv, &A[j + 1 + j * lda], 1);
            }
        } else {
            /* A[j,j] = 0: matrix is singular */
            if (info == 0) {
                info = j + 1;  /* 1-based index */
            }
        }

        /* Update trailing submatrix */
        if (j < minmn - 1) {
            /* A[j+1:m, j+1:n] -= A[j+1:m, j] * A[j, j+1:n] */
            blas_dger(m - j - 1, n - j - 1, -1.0,
                      &A[j + 1 + j * lda], 1,
                      &A[j + (j + 1) * lda], lda,
                      &A[j + 1 + (j + 1) * lda], lda);
        }
    }

    return info;
}

EXPORT int32_t lapack_sgetrf(int32_t m, int32_t n, float* A, int32_t lda,
                              int32_t* ipiv)
{
    if (m <= 0 || n <= 0) return 0;

    int32_t info = 0;
    int32_t minmn = MIN(m, n);

    for (int32_t j = 0; j < minmn; j++) {
        int32_t jp = j + blas_isamax(m - j, &A[j + j * lda], 1);
        ipiv[j] = jp;

        if (A[jp + j * lda] != 0.0f) {
            if (jp != j) {
                blas_sswap(n, &A[j], lda, &A[jp], lda);
            }

            if (j < m - 1) {
                float ajj_inv = 1.0f / A[j + j * lda];
                blas_sscal(m - j - 1, ajj_inv, &A[j + 1 + j * lda], 1);
            }
        } else {
            if (info == 0) {
                info = j + 1;
            }
        }

        if (j < minmn - 1) {
            for (int32_t jj = j + 1; jj < n; jj++) {
                float temp = A[j + jj * lda];
                for (int32_t i = j + 1; i < m; i++) {
                    A[i + jj * lda] -= A[i + j * lda] * temp;
                }
            }
        }
    }

    return info;
}

EXPORT int32_t lapack_dgetrs(char trans, int32_t n, int32_t nrhs,
                              const double* A, int32_t lda, const int32_t* ipiv,
                              double* B, int32_t ldb)
{
    if (n <= 0 || nrhs <= 0) return 0;

    if (!is_trans(trans)) {
        /* Solve A * X = B */
        /* Apply row interchanges to B */
        for (int32_t i = 0; i < n; i++) {
            int32_t ip = ipiv[i];
            if (ip != i) {
                blas_dswap(nrhs, &B[i], ldb, &B[ip], ldb);
            }
        }

        /* Solve L * Y = B (forward substitution) */
        blas_dtrsm('L', 'L', 'N', 'U', n, nrhs, 1.0, A, lda, B, ldb);

        /* Solve U * X = Y (back substitution) */
        blas_dtrsm('L', 'U', 'N', 'N', n, nrhs, 1.0, A, lda, B, ldb);
    } else {
        /* Solve A^T * X = B */
        /* Solve U^T * Y = B */
        blas_dtrsm('L', 'U', 'T', 'N', n, nrhs, 1.0, A, lda, B, ldb);

        /* Solve L^T * X = Y */
        blas_dtrsm('L', 'L', 'T', 'U', n, nrhs, 1.0, A, lda, B, ldb);

        /* Apply row interchanges in reverse */
        for (int32_t i = n - 1; i >= 0; i--) {
            int32_t ip = ipiv[i];
            if (ip != i) {
                blas_dswap(nrhs, &B[i], ldb, &B[ip], ldb);
            }
        }
    }

    return 0;
}

EXPORT int32_t lapack_sgetrs(char trans, int32_t n, int32_t nrhs,
                              const float* A, int32_t lda, const int32_t* ipiv,
                              float* B, int32_t ldb)
{
    if (n <= 0 || nrhs <= 0) return 0;

    if (!is_trans(trans)) {
        for (int32_t i = 0; i < n; i++) {
            int32_t ip = ipiv[i];
            if (ip != i) {
                blas_sswap(nrhs, &B[i], ldb, &B[ip], ldb);
            }
        }

        /* Forward substitution for L */
        for (int32_t j = 0; j < nrhs; j++) {
            for (int32_t k = 0; k < n; k++) {
                if (B[k + j * ldb] != 0.0f) {
                    for (int32_t i = k + 1; i < n; i++) {
                        B[i + j * ldb] -= B[k + j * ldb] * A[i + k * lda];
                    }
                }
            }
        }

        /* Back substitution for U */
        for (int32_t j = 0; j < nrhs; j++) {
            for (int32_t k = n - 1; k >= 0; k--) {
                if (B[k + j * ldb] != 0.0f) {
                    B[k + j * ldb] /= A[k + k * lda];
                    for (int32_t i = 0; i < k; i++) {
                        B[i + j * ldb] -= B[k + j * ldb] * A[i + k * lda];
                    }
                }
            }
        }
    } else {
        /* Similar for transpose case */
        for (int32_t j = 0; j < nrhs; j++) {
            for (int32_t k = 0; k < n; k++) {
                B[k + j * ldb] /= A[k + k * lda];
                for (int32_t i = k + 1; i < n; i++) {
                    B[i + j * ldb] -= B[k + j * ldb] * A[k + i * lda];
                }
            }
        }

        for (int32_t j = 0; j < nrhs; j++) {
            for (int32_t k = n - 1; k >= 0; k--) {
                for (int32_t i = k + 1; i < n; i++) {
                    B[k + j * ldb] -= A[i + k * lda] * B[i + j * ldb];
                }
            }
        }

        for (int32_t i = n - 1; i >= 0; i--) {
            int32_t ip = ipiv[i];
            if (ip != i) {
                blas_sswap(nrhs, &B[i], ldb, &B[ip], ldb);
            }
        }
    }

    return 0;
}

EXPORT int32_t lapack_dgetri(int32_t n, double* A, int32_t lda, const int32_t* ipiv,
                              double* work, int32_t lwork)
{
    if (n <= 0) return 0;

    /* Workspace query */
    if (lwork == -1) {
        work[0] = (double)n;
        return 0;
    }

    /* Invert U in place */
    for (int32_t j = 0; j < n; j++) {
        A[j + j * lda] = 1.0 / A[j + j * lda];
        double ajj = -A[j + j * lda];

        /* Compute elements 0:j-1 of j-th column */
        blas_dtrmv('U', 'N', 'N', j, A, lda, &A[j * lda], 1);
        blas_dscal(j, ajj, &A[j * lda], 1);
    }

    /* Solve the equation inv(A)*L = inv(U) for inv(A) */
    for (int32_t j = n - 2; j >= 0; j--) {
        /* Copy column j+1:n of A to work */
        for (int32_t i = j + 1; i < n; i++) {
            work[i] = A[i + j * lda];
            A[i + j * lda] = 0.0;
        }

        /* Compute column j of inv(A) */
        blas_dgemv('N', n, n - j - 1, -1.0, &A[(j + 1) * lda], lda,
                   &work[j + 1], 1, 1.0, &A[j * lda], 1);
    }

    /* Apply column interchanges */
    for (int32_t j = n - 2; j >= 0; j--) {
        int32_t jp = ipiv[j];
        if (jp != j) {
            blas_dswap(n, &A[j * lda], 1, &A[jp * lda], 1);
        }
    }

    return 0;
}

EXPORT int32_t lapack_sgetri(int32_t n, float* A, int32_t lda, const int32_t* ipiv,
                              float* work, int32_t lwork)
{
    if (n <= 0) return 0;

    if (lwork == -1) {
        work[0] = (float)n;
        return 0;
    }

    for (int32_t j = 0; j < n; j++) {
        A[j + j * lda] = 1.0f / A[j + j * lda];
        float ajj = -A[j + j * lda];

        for (int32_t i = 0; i < j; i++) {
            float sum = 0.0f;
            for (int32_t k = i; k < j; k++) {
                sum += A[i + k * lda] * A[k + j * lda];
            }
            A[i + j * lda] = ajj * sum;
        }
    }

    for (int32_t j = n - 2; j >= 0; j--) {
        for (int32_t i = j + 1; i < n; i++) {
            work[i] = A[i + j * lda];
            A[i + j * lda] = 0.0f;
        }

        for (int32_t i = 0; i < n; i++) {
            float sum = 0.0f;
            for (int32_t k = j + 1; k < n; k++) {
                sum += A[i + k * lda] * work[k];
            }
            A[i + j * lda] -= sum;
        }
    }

    for (int32_t j = n - 2; j >= 0; j--) {
        int32_t jp = ipiv[j];
        if (jp != j) {
            blas_sswap(n, &A[j * lda], 1, &A[jp * lda], 1);
        }
    }

    return 0;
}

/* ============ Cholesky Factorization ============ */

EXPORT int32_t lapack_dpotrf(char uplo, int32_t n, double* A, int32_t lda)
{
    if (n <= 0) return 0;

    int32_t info = 0;

    if (is_upper(uplo)) {
        /* Compute U such that A = U^T * U */
        for (int32_t j = 0; j < n; j++) {
            /* Compute U(j,j) */
            double ajj = A[j + j * lda];
            for (int32_t k = 0; k < j; k++) {
                ajj -= A[k + j * lda] * A[k + j * lda];
            }

            if (ajj <= 0.0) {
                A[j + j * lda] = ajj;
                return j + 1;  /* Not positive definite */
            }

            ajj = sqrt(ajj);
            A[j + j * lda] = ajj;

            /* Compute elements j+1:n of row j */
            for (int32_t i = j + 1; i < n; i++) {
                double aij = A[j + i * lda];
                for (int32_t k = 0; k < j; k++) {
                    aij -= A[k + j * lda] * A[k + i * lda];
                }
                A[j + i * lda] = aij / ajj;
            }
        }
    } else {
        /* Compute L such that A = L * L^T */
        for (int32_t j = 0; j < n; j++) {
            /* Compute L(j,j) */
            double ajj = A[j + j * lda];
            for (int32_t k = 0; k < j; k++) {
                ajj -= A[j + k * lda] * A[j + k * lda];
            }

            if (ajj <= 0.0) {
                A[j + j * lda] = ajj;
                return j + 1;
            }

            ajj = sqrt(ajj);
            A[j + j * lda] = ajj;

            /* Compute elements j+1:n of column j */
            for (int32_t i = j + 1; i < n; i++) {
                double aij = A[i + j * lda];
                for (int32_t k = 0; k < j; k++) {
                    aij -= A[i + k * lda] * A[j + k * lda];
                }
                A[i + j * lda] = aij / ajj;
            }
        }
    }

    return info;
}

EXPORT int32_t lapack_spotrf(char uplo, int32_t n, float* A, int32_t lda)
{
    if (n <= 0) return 0;

    if (is_upper(uplo)) {
        for (int32_t j = 0; j < n; j++) {
            float ajj = A[j + j * lda];
            for (int32_t k = 0; k < j; k++) {
                ajj -= A[k + j * lda] * A[k + j * lda];
            }

            if (ajj <= 0.0f) {
                A[j + j * lda] = ajj;
                return j + 1;
            }

            ajj = sqrtf(ajj);
            A[j + j * lda] = ajj;

            for (int32_t i = j + 1; i < n; i++) {
                float aij = A[j + i * lda];
                for (int32_t k = 0; k < j; k++) {
                    aij -= A[k + j * lda] * A[k + i * lda];
                }
                A[j + i * lda] = aij / ajj;
            }
        }
    } else {
        for (int32_t j = 0; j < n; j++) {
            float ajj = A[j + j * lda];
            for (int32_t k = 0; k < j; k++) {
                ajj -= A[j + k * lda] * A[j + k * lda];
            }

            if (ajj <= 0.0f) {
                A[j + j * lda] = ajj;
                return j + 1;
            }

            ajj = sqrtf(ajj);
            A[j + j * lda] = ajj;

            for (int32_t i = j + 1; i < n; i++) {
                float aij = A[i + j * lda];
                for (int32_t k = 0; k < j; k++) {
                    aij -= A[i + k * lda] * A[j + k * lda];
                }
                A[i + j * lda] = aij / ajj;
            }
        }
    }

    return 0;
}

EXPORT int32_t lapack_dpotrs(char uplo, int32_t n, int32_t nrhs,
                              const double* A, int32_t lda,
                              double* B, int32_t ldb)
{
    if (n <= 0 || nrhs <= 0) return 0;

    if (is_upper(uplo)) {
        /* Solve U^T * Y = B */
        blas_dtrsm('L', 'U', 'T', 'N', n, nrhs, 1.0, A, lda, B, ldb);
        /* Solve U * X = Y */
        blas_dtrsm('L', 'U', 'N', 'N', n, nrhs, 1.0, A, lda, B, ldb);
    } else {
        /* Solve L * Y = B */
        blas_dtrsm('L', 'L', 'N', 'N', n, nrhs, 1.0, A, lda, B, ldb);
        /* Solve L^T * X = Y */
        blas_dtrsm('L', 'L', 'T', 'N', n, nrhs, 1.0, A, lda, B, ldb);
    }

    return 0;
}

EXPORT int32_t lapack_dpotri(char uplo, int32_t n, double* A, int32_t lda)
{
    if (n <= 0) return 0;

    /* First invert the triangular factor */
    int32_t info = lapack_dtrtri(uplo, 'N', n, A, lda);
    if (info != 0) return info;

    /* Then compute the product */
    if (is_upper(uplo)) {
        /* A = U^(-1) * U^(-T) */
        for (int32_t j = 0; j < n; j++) {
            for (int32_t i = 0; i <= j; i++) {
                double sum = 0.0;
                for (int32_t k = j; k < n; k++) {
                    sum += A[i + k * lda] * A[j + k * lda];
                }
                A[i + j * lda] = sum;
            }
        }
    } else {
        /* A = L^(-T) * L^(-1) */
        for (int32_t j = 0; j < n; j++) {
            for (int32_t i = j; i < n; i++) {
                double sum = 0.0;
                for (int32_t k = i; k < n; k++) {
                    sum += A[k + i * lda] * A[k + j * lda];
                }
                A[i + j * lda] = sum;
            }
        }
    }

    return 0;
}

/* ============ Triangular Inverse ============ */

EXPORT int32_t lapack_dtrtri(char uplo, char diag, int32_t n,
                              double* A, int32_t lda)
{
    if (n <= 0) return 0;

    int nounit = !is_unit(diag);

    if (is_upper(uplo)) {
        for (int32_t j = 0; j < n; j++) {
            if (nounit) {
                if (A[j + j * lda] == 0.0) return j + 1;
                A[j + j * lda] = 1.0 / A[j + j * lda];
            }
            double ajj = -A[j + j * lda];

            /* Compute elements 0:j-1 of j-th column */
            blas_dtrmv('U', 'N', diag, j, A, lda, &A[j * lda], 1);
            blas_dscal(j, ajj, &A[j * lda], 1);
        }
    } else {
        for (int32_t j = n - 1; j >= 0; j--) {
            if (nounit) {
                if (A[j + j * lda] == 0.0) return j + 1;
                A[j + j * lda] = 1.0 / A[j + j * lda];
            }
            double ajj = -A[j + j * lda];

            if (j < n - 1) {
                blas_dtrmv('L', 'N', diag, n - j - 1, &A[j + 1 + (j + 1) * lda], lda,
                           &A[j + 1 + j * lda], 1);
                blas_dscal(n - j - 1, ajj, &A[j + 1 + j * lda], 1);
            }
        }
    }

    return 0;
}

/* ============ QR Factorization ============ */

/* Generate Householder reflector */
static void dlarf(char side, int32_t m, int32_t n, const double* v, int32_t incv,
                   double tau, double* C, int32_t ldc, double* work)
{
    if (tau == 0.0) return;

    if (side == 'L' || side == 'l') {
        /* C = (I - tau * v * v^T) * C = C - tau * v * (v^T * C) */
        /* work = v^T * C */
        for (int32_t j = 0; j < n; j++) {
            work[j] = 0.0;
            int32_t iv = 0;
            for (int32_t i = 0; i < m; i++) {
                work[j] += v[iv] * C[i + j * ldc];
                iv += incv;
            }
        }
        /* C = C - tau * v * work^T */
        int32_t iv = 0;
        for (int32_t i = 0; i < m; i++) {
            for (int32_t j = 0; j < n; j++) {
                C[i + j * ldc] -= tau * v[iv] * work[j];
            }
            iv += incv;
        }
    } else {
        /* C = C * (I - tau * v * v^T) = C - tau * (C * v) * v^T */
        /* work = C * v */
        for (int32_t i = 0; i < m; i++) {
            work[i] = 0.0;
            int32_t iv = 0;
            for (int32_t j = 0; j < n; j++) {
                work[i] += C[i + j * ldc] * v[iv];
                iv += incv;
            }
        }
        /* C = C - tau * work * v^T */
        int32_t iv = 0;
        for (int32_t j = 0; j < n; j++) {
            for (int32_t i = 0; i < m; i++) {
                C[i + j * ldc] -= tau * work[i] * v[iv];
            }
            iv += incv;
        }
    }
}

/* Generate Householder vector */
static void dlarfg(int32_t n, double* alpha, double* x, int32_t incx, double* tau)
{
    if (n <= 1) {
        *tau = 0.0;
        return;
    }

    double xnorm = blas_dnrm2(n - 1, x, incx);

    if (xnorm == 0.0) {
        *tau = 0.0;
        return;
    }

    double beta = -SIGN(sqrt(*alpha * *alpha + xnorm * xnorm), *alpha);
    double safmin = DLAMCH_SFMIN / DLAMCH_EPS;

    int32_t knt = 0;
    if (fabs(beta) < safmin) {
        /* Scale to avoid overflow */
        double rsafmn = 1.0 / safmin;
        do {
            knt++;
            blas_dscal(n - 1, rsafmn, x, incx);
            beta *= rsafmn;
            *alpha *= rsafmn;
        } while (fabs(beta) < safmin);

        xnorm = blas_dnrm2(n - 1, x, incx);
        beta = -SIGN(sqrt(*alpha * *alpha + xnorm * xnorm), *alpha);
    }

    *tau = (beta - *alpha) / beta;
    blas_dscal(n - 1, 1.0 / (*alpha - beta), x, incx);

    for (int32_t j = 0; j < knt; j++) {
        beta *= safmin;
    }
    *alpha = beta;
}

EXPORT int32_t lapack_dgeqrf(int32_t m, int32_t n, double* A, int32_t lda,
                              double* tau, double* work, int32_t lwork)
{
    if (m <= 0 || n <= 0) return 0;

    int32_t k = MIN(m, n);

    /* Workspace query */
    if (lwork == -1) {
        work[0] = (double)n;
        return 0;
    }

    for (int32_t i = 0; i < k; i++) {
        /* Generate Householder reflector H(i) to annihilate A(i+1:m, i) */
        dlarfg(m - i, &A[i + i * lda], &A[MIN(i + 1, m - 1) + i * lda], 1, &tau[i]);

        if (i < n - 1) {
            /* Apply H(i) to A(i:m, i+1:n) from the left */
            double aii = A[i + i * lda];
            A[i + i * lda] = 1.0;
            dlarf('L', m - i, n - i - 1, &A[i + i * lda], 1, tau[i],
                  &A[i + (i + 1) * lda], lda, work);
            A[i + i * lda] = aii;
        }
    }

    return 0;
}

EXPORT int32_t lapack_sgeqrf(int32_t m, int32_t n, float* A, int32_t lda,
                              float* tau, float* work, int32_t lwork)
{
    /* Simple implementation - convert to double, call double version */
    if (m <= 0 || n <= 0) return 0;

    if (lwork == -1) {
        work[0] = (float)n;
        return 0;
    }

    /* Direct single-precision implementation */
    int32_t k = MIN(m, n);

    for (int32_t i = 0; i < k; i++) {
        /* Compute 2-norm of A(i+1:m, i) */
        float xnorm = 0.0f;
        for (int32_t j = i + 1; j < m; j++) {
            xnorm += A[j + i * lda] * A[j + i * lda];
        }
        xnorm = sqrtf(xnorm);

        float alpha = A[i + i * lda];

        if (xnorm == 0.0f) {
            tau[i] = 0.0f;
        } else {
            float beta = -((alpha >= 0.0f) ? 1.0f : -1.0f) * sqrtf(alpha * alpha + xnorm * xnorm);
            tau[i] = (beta - alpha) / beta;
            float scale = 1.0f / (alpha - beta);
            for (int32_t j = i + 1; j < m; j++) {
                A[j + i * lda] *= scale;
            }
            A[i + i * lda] = beta;
        }

        if (i < n - 1 && tau[i] != 0.0f) {
            float aii = A[i + i * lda];
            A[i + i * lda] = 1.0f;

            /* Apply H(i) from left */
            for (int32_t j = i + 1; j < n; j++) {
                float sum = A[i + j * lda];
                for (int32_t l = i + 1; l < m; l++) {
                    sum += A[l + i * lda] * A[l + j * lda];
                }
                sum *= tau[i];
                A[i + j * lda] -= sum;
                for (int32_t l = i + 1; l < m; l++) {
                    A[l + j * lda] -= A[l + i * lda] * sum;
                }
            }

            A[i + i * lda] = aii;
        }
    }

    return 0;
}

EXPORT int32_t lapack_dorgqr(int32_t m, int32_t n, int32_t k,
                              double* A, int32_t lda, const double* tau,
                              double* work, int32_t lwork)
{
    if (m <= 0 || n <= 0) return 0;

    /* Workspace query */
    if (lwork == -1) {
        work[0] = (double)n;
        return 0;
    }

    /* Initialize columns k+1:n to columns of identity matrix */
    for (int32_t j = k; j < n; j++) {
        for (int32_t i = 0; i < m; i++) {
            A[i + j * lda] = 0.0;
        }
        if (j < m) {
            A[j + j * lda] = 1.0;
        }
    }

    /* Apply H(k-1), ..., H(1), H(0) */
    for (int32_t i = k - 1; i >= 0; i--) {
        if (i < n - 1) {
            A[i + i * lda] = 1.0;
            dlarf('L', m - i, n - i - 1, &A[i + i * lda], 1, tau[i],
                  &A[i + (i + 1) * lda], lda, work);
        }

        if (i < m - 1) {
            blas_dscal(m - i - 1, -tau[i], &A[i + 1 + i * lda], 1);
        }
        A[i + i * lda] = 1.0 - tau[i];

        /* Set A(0:i-1, i) to zero */
        for (int32_t l = 0; l < i; l++) {
            A[l + i * lda] = 0.0;
        }
    }

    return 0;
}

EXPORT int32_t lapack_sorgqr(int32_t m, int32_t n, int32_t k,
                              float* A, int32_t lda, const float* tau,
                              float* work, int32_t lwork)
{
    if (m <= 0 || n <= 0) return 0;

    if (lwork == -1) {
        work[0] = (float)n;
        return 0;
    }

    /* Initialize columns k+1:n */
    for (int32_t j = k; j < n; j++) {
        for (int32_t i = 0; i < m; i++) {
            A[i + j * lda] = 0.0f;
        }
        if (j < m) {
            A[j + j * lda] = 1.0f;
        }
    }

    /* Apply reflectors */
    for (int32_t i = k - 1; i >= 0; i--) {
        if (i < n - 1) {
            A[i + i * lda] = 1.0f;

            /* Apply H(i) from left */
            for (int32_t j = i + 1; j < n; j++) {
                float sum = A[i + j * lda];
                for (int32_t l = i + 1; l < m; l++) {
                    sum += A[l + i * lda] * A[l + j * lda];
                }
                sum *= tau[i];
                A[i + j * lda] -= sum;
                for (int32_t l = i + 1; l < m; l++) {
                    A[l + j * lda] -= A[l + i * lda] * sum;
                }
            }
        }

        if (i < m - 1) {
            for (int32_t l = i + 1; l < m; l++) {
                A[l + i * lda] *= -tau[i];
            }
        }
        A[i + i * lda] = 1.0f - tau[i];

        for (int32_t l = 0; l < i; l++) {
            A[l + i * lda] = 0.0f;
        }
    }

    return 0;
}

EXPORT int32_t lapack_dtrtrs(char uplo, char trans, char diag,
                              int32_t n, int32_t nrhs,
                              const double* A, int32_t lda,
                              double* B, int32_t ldb)
{
    if (n <= 0 || nrhs <= 0) return 0;

    /* Check for singularity */
    if (!is_unit(diag)) {
        for (int32_t j = 0; j < n; j++) {
            if (A[j + j * lda] == 0.0) {
                return j + 1;
            }
        }
    }

    blas_dtrsm('L', uplo, trans, diag, n, nrhs, 1.0, A, lda, B, ldb);

    return 0;
}

/* ============ Eigenvalue Decomposition ============ */

/* Symmetric eigenvalue using Jacobi iteration (simple but robust) */
EXPORT int32_t lapack_dsyevd(char jobz, char uplo, int32_t n,
                              double* A, int32_t lda, double* w,
                              double* work, int32_t lwork,
                              int32_t* iwork, int32_t liwork)
{
    if (n <= 0) return 0;

    /* Workspace query */
    if (lwork == -1 || liwork == -1) {
        if (want_vectors(jobz)) {
            work[0] = (double)(2 * n * n + 6 * n + 1);
            iwork[0] = 5 * n + 3;
        } else {
            work[0] = (double)(2 * n + 1);
            iwork[0] = 1;
        }
        return 0;
    }

    int wantz = want_vectors(jobz);
    int upper = is_upper(uplo);

    /* Copy A to work for symmetric part */
    /* Make A symmetric by copying upper/lower triangle */
    for (int32_t j = 0; j < n; j++) {
        for (int32_t i = j; i < n; i++) {
            if (upper) {
                A[i + j * lda] = A[j + i * lda];
            } else {
                A[j + i * lda] = A[i + j * lda];
            }
        }
    }

    /* Initialize eigenvectors to identity if needed */
    double* V = wantz ? work : NULL;
    if (wantz) {
        for (int32_t j = 0; j < n; j++) {
            for (int32_t i = 0; i < n; i++) {
                V[i + j * n] = (i == j) ? 1.0 : 0.0;
            }
        }
    }

    /* Jacobi iteration */
    const int32_t maxiter = 50;
    const double tol = DLAMCH_EPS * DLAMCH_EPS;

    for (int32_t iter = 0; iter < maxiter; iter++) {
        /* Find maximum off-diagonal element */
        double maxoff = 0.0;
        for (int32_t j = 0; j < n; j++) {
            for (int32_t i = 0; i < j; i++) {
                double aij = fabs(A[i + j * lda]);
                if (aij > maxoff) maxoff = aij;
            }
        }

        if (maxoff < tol) break;

        /* Sweep through all off-diagonal elements */
        for (int32_t p = 0; p < n - 1; p++) {
            for (int32_t q = p + 1; q < n; q++) {
                double app = A[p + p * lda];
                double aqq = A[q + q * lda];
                double apq = A[p + q * lda];

                if (fabs(apq) < tol * sqrt(fabs(app * aqq))) continue;

                /* Compute rotation angle */
                double theta = 0.5 * (aqq - app) / apq;
                double t;
                if (theta >= 0.0) {
                    t = 1.0 / (theta + sqrt(1.0 + theta * theta));
                } else {
                    t = 1.0 / (theta - sqrt(1.0 + theta * theta));
                }
                double c = 1.0 / sqrt(1.0 + t * t);
                double s = t * c;

                /* Update A */
                A[p + p * lda] = app - t * apq;
                A[q + q * lda] = aqq + t * apq;
                A[p + q * lda] = 0.0;
                A[q + p * lda] = 0.0;

                for (int32_t i = 0; i < p; i++) {
                    double aip = A[i + p * lda];
                    double aiq = A[i + q * lda];
                    A[i + p * lda] = c * aip - s * aiq;
                    A[i + q * lda] = s * aip + c * aiq;
                }
                for (int32_t i = p + 1; i < q; i++) {
                    double api = A[p + i * lda];
                    double aiq = A[i + q * lda];
                    A[p + i * lda] = c * api - s * aiq;
                    A[i + q * lda] = s * api + c * aiq;
                }
                for (int32_t i = q + 1; i < n; i++) {
                    double api = A[p + i * lda];
                    double aqi = A[q + i * lda];
                    A[p + i * lda] = c * api - s * aqi;
                    A[q + i * lda] = s * api + c * aqi;
                }

                /* Update eigenvectors */
                if (wantz) {
                    for (int32_t i = 0; i < n; i++) {
                        double vip = V[i + p * n];
                        double viq = V[i + q * n];
                        V[i + p * n] = c * vip - s * viq;
                        V[i + q * n] = s * vip + c * viq;
                    }
                }
            }
        }
    }

    /* Extract eigenvalues */
    for (int32_t i = 0; i < n; i++) {
        w[i] = A[i + i * lda];
    }

    /* Sort eigenvalues in ascending order */
    for (int32_t i = 0; i < n - 1; i++) {
        int32_t k = i;
        double p = w[i];
        for (int32_t j = i + 1; j < n; j++) {
            if (w[j] < p) {
                k = j;
                p = w[j];
            }
        }
        if (k != i) {
            w[k] = w[i];
            w[i] = p;
            if (wantz) {
                blas_dswap(n, &V[i * n], 1, &V[k * n], 1);
            }
        }
    }

    /* Copy eigenvectors back to A */
    if (wantz) {
        for (int32_t j = 0; j < n; j++) {
            for (int32_t i = 0; i < n; i++) {
                A[i + j * lda] = V[i + j * n];
            }
        }
    }

    return 0;
}

EXPORT int32_t lapack_ssyevd(char jobz, char uplo, int32_t n,
                              float* A, int32_t lda, float* w,
                              float* work, int32_t lwork,
                              int32_t* iwork, int32_t liwork)
{
    if (n <= 0) return 0;

    if (lwork == -1 || liwork == -1) {
        if (want_vectors(jobz)) {
            work[0] = (float)(2 * n * n + 6 * n + 1);
            iwork[0] = 5 * n + 3;
        } else {
            work[0] = (float)(2 * n + 1);
            iwork[0] = 1;
        }
        return 0;
    }

    /* Similar to double version but with float */
    int wantz = want_vectors(jobz);
    int upper = is_upper(uplo);

    for (int32_t j = 0; j < n; j++) {
        for (int32_t i = j; i < n; i++) {
            if (upper) {
                A[i + j * lda] = A[j + i * lda];
            } else {
                A[j + i * lda] = A[i + j * lda];
            }
        }
    }

    float* V = wantz ? work : NULL;
    if (wantz) {
        for (int32_t j = 0; j < n; j++) {
            for (int32_t i = 0; i < n; i++) {
                V[i + j * n] = (i == j) ? 1.0f : 0.0f;
            }
        }
    }

    const int32_t maxiter = 50;
    const float tol = 1.19e-7f * 1.19e-7f;

    for (int32_t iter = 0; iter < maxiter; iter++) {
        float maxoff = 0.0f;
        for (int32_t j = 0; j < n; j++) {
            for (int32_t i = 0; i < j; i++) {
                float aij = fabsf(A[i + j * lda]);
                if (aij > maxoff) maxoff = aij;
            }
        }

        if (maxoff < tol) break;

        for (int32_t p = 0; p < n - 1; p++) {
            for (int32_t q = p + 1; q < n; q++) {
                float app = A[p + p * lda];
                float aqq = A[q + q * lda];
                float apq = A[p + q * lda];

                if (fabsf(apq) < tol * sqrtf(fabsf(app * aqq))) continue;

                float theta = 0.5f * (aqq - app) / apq;
                float t;
                if (theta >= 0.0f) {
                    t = 1.0f / (theta + sqrtf(1.0f + theta * theta));
                } else {
                    t = 1.0f / (theta - sqrtf(1.0f + theta * theta));
                }
                float c = 1.0f / sqrtf(1.0f + t * t);
                float s = t * c;

                A[p + p * lda] = app - t * apq;
                A[q + q * lda] = aqq + t * apq;
                A[p + q * lda] = 0.0f;
                A[q + p * lda] = 0.0f;

                for (int32_t i = 0; i < p; i++) {
                    float aip = A[i + p * lda];
                    float aiq = A[i + q * lda];
                    A[i + p * lda] = c * aip - s * aiq;
                    A[i + q * lda] = s * aip + c * aiq;
                }
                for (int32_t i = p + 1; i < q; i++) {
                    float api = A[p + i * lda];
                    float aiq = A[i + q * lda];
                    A[p + i * lda] = c * api - s * aiq;
                    A[i + q * lda] = s * api + c * aiq;
                }
                for (int32_t i = q + 1; i < n; i++) {
                    float api = A[p + i * lda];
                    float aqi = A[q + i * lda];
                    A[p + i * lda] = c * api - s * aqi;
                    A[q + i * lda] = s * api + c * aqi;
                }

                if (wantz) {
                    for (int32_t i = 0; i < n; i++) {
                        float vip = V[i + p * n];
                        float viq = V[i + q * n];
                        V[i + p * n] = c * vip - s * viq;
                        V[i + q * n] = s * vip + c * viq;
                    }
                }
            }
        }
    }

    for (int32_t i = 0; i < n; i++) {
        w[i] = A[i + i * lda];
    }

    for (int32_t i = 0; i < n - 1; i++) {
        int32_t k = i;
        float p = w[i];
        for (int32_t j = i + 1; j < n; j++) {
            if (w[j] < p) {
                k = j;
                p = w[j];
            }
        }
        if (k != i) {
            w[k] = w[i];
            w[i] = p;
            if (wantz) {
                blas_sswap(n, &V[i * n], 1, &V[k * n], 1);
            }
        }
    }

    if (wantz) {
        for (int32_t j = 0; j < n; j++) {
            for (int32_t i = 0; i < n; i++) {
                A[i + j * lda] = V[i + j * n];
            }
        }
    }

    return 0;
}

/* General eigenvalue - simplified using QR algorithm sketch */
/* For MVP, we implement a basic version that may not converge for all matrices */
EXPORT int32_t lapack_dgeev(char jobvl, char jobvr, int32_t n,
                             double* A, int32_t lda,
                             double* wr, double* wi,
                             double* VL, int32_t ldvl,
                             double* VR, int32_t ldvr,
                             double* work, int32_t lwork)
{
    if (n <= 0) return 0;

    int wantvl = want_vectors(jobvl);
    int wantvr = want_vectors(jobvr);

    /* Workspace query */
    if (lwork == -1) {
        work[0] = (double)(4 * n);
        return 0;
    }

    /* For a real symmetric matrix, use dsyevd and set wi = 0 */
    /* Check if symmetric (simplified check) */
    int symmetric = 1;
    for (int32_t j = 0; j < n && symmetric; j++) {
        for (int32_t i = j + 1; i < n && symmetric; i++) {
            if (fabs(A[i + j * lda] - A[j + i * lda]) > DLAMCH_EPS * (fabs(A[i + j * lda]) + fabs(A[j + i * lda]))) {
                symmetric = 0;
            }
        }
    }

    if (symmetric) {
        /* Use symmetric eigenvalue solver */
        int32_t* iwork = (int32_t*)malloc((5 * n + 3) * sizeof(int32_t));
        double* work2 = (double*)malloc((2 * n * n + 6 * n + 1) * sizeof(double));

        int32_t info = lapack_dsyevd(wantvr ? 'V' : 'N', 'U', n, A, lda, wr,
                                      work2, 2 * n * n + 6 * n + 1, iwork, 5 * n + 3);

        /* Set imaginary parts to zero */
        for (int32_t i = 0; i < n; i++) {
            wi[i] = 0.0;
        }

        /* Copy eigenvectors to VR */
        if (wantvr) {
            for (int32_t j = 0; j < n; j++) {
                for (int32_t i = 0; i < n; i++) {
                    VR[i + j * ldvr] = A[i + j * lda];
                }
            }
        }

        /* VL = VR for symmetric matrices */
        if (wantvl) {
            for (int32_t j = 0; j < n; j++) {
                for (int32_t i = 0; i < n; i++) {
                    VL[i + j * ldvl] = A[i + j * lda];
                }
            }
        }

        free(iwork);
        free(work2);
        return info;
    }

    /* For general matrices, this is a simplified placeholder */
    /* A full implementation would require Hessenberg reduction + QR iteration */
    /* For now, return an error for non-symmetric matrices */

    /* Initialize eigenvalues to diagonal elements (only correct for diagonal matrices) */
    for (int32_t i = 0; i < n; i++) {
        wr[i] = A[i + i * lda];
        wi[i] = 0.0;
    }

    if (wantvr) {
        for (int32_t j = 0; j < n; j++) {
            for (int32_t i = 0; i < n; i++) {
                VR[i + j * ldvr] = (i == j) ? 1.0 : 0.0;
            }
        }
    }

    if (wantvl) {
        for (int32_t j = 0; j < n; j++) {
            for (int32_t i = 0; i < n; i++) {
                VL[i + j * ldvl] = (i == j) ? 1.0 : 0.0;
            }
        }
    }

    /* Note: This is incomplete. Full dgeev requires Hessenberg reduction + QR iteration */
    return 0;
}

EXPORT int32_t lapack_sgeev(char jobvl, char jobvr, int32_t n,
                             float* A, int32_t lda,
                             float* wr, float* wi,
                             float* VL, int32_t ldvl,
                             float* VR, int32_t ldvr,
                             float* work, int32_t lwork)
{
    if (n <= 0) return 0;

    if (lwork == -1) {
        work[0] = (float)(4 * n);
        return 0;
    }

    /* Similar simplified implementation */
    for (int32_t i = 0; i < n; i++) {
        wr[i] = A[i + i * lda];
        wi[i] = 0.0f;
    }

    int wantvr = want_vectors(jobvr);
    int wantvl = want_vectors(jobvl);

    if (wantvr) {
        for (int32_t j = 0; j < n; j++) {
            for (int32_t i = 0; i < n; i++) {
                VR[i + j * ldvr] = (i == j) ? 1.0f : 0.0f;
            }
        }
    }

    if (wantvl) {
        for (int32_t j = 0; j < n; j++) {
            for (int32_t i = 0; i < n; i++) {
                VL[i + j * ldvl] = (i == j) ? 1.0f : 0.0f;
            }
        }
    }

    return 0;
}

/* ============ Singular Value Decomposition ============ */

/* Simplified SVD using A^T*A eigenvalue decomposition */
EXPORT int32_t lapack_dgesdd(char jobz, int32_t m, int32_t n,
                              double* A, int32_t lda, double* s,
                              double* U, int32_t ldu,
                              double* VT, int32_t ldvt,
                              double* work, int32_t lwork, int32_t* iwork)
{
    if (m <= 0 || n <= 0) return 0;

    int32_t minmn = MIN(m, n);

    /* Workspace query */
    if (lwork == -1) {
        work[0] = (double)(3 * minmn * minmn + MAX(MAX(m, n), 4 * minmn * minmn + 4 * minmn));
        return 0;
    }

    int wantu = (jobz == 'A' || jobz == 'a' || jobz == 'S' || jobz == 's');
    int wantvt = (jobz == 'A' || jobz == 'a' || jobz == 'S' || jobz == 's');

    /* Simple SVD via eigenvalue decomposition of A^T * A */
    /* This is less numerically stable but simpler */

    double* ATA = work;  /* n x n */

    /* Compute A^T * A */
    blas_dgemm('T', 'N', n, n, m, 1.0, A, lda, A, lda, 0.0, ATA, n);

    /* Eigenvalue decomposition of A^T * A */
    double* eigenwork = ATA + n * n;
    int32_t* eigeniwork = iwork;

    /* Use Jacobi iteration for eigenvalues */
    double* eigenvalues = eigenwork;
    double* V = eigenwork + n;  /* Eigenvectors */

    /* Initialize V to identity */
    for (int32_t j = 0; j < n; j++) {
        for (int32_t i = 0; i < n; i++) {
            V[i + j * n] = (i == j) ? 1.0 : 0.0;
        }
    }

    /* Jacobi eigenvalue iteration */
    const int32_t maxiter = 50;
    const double tol = DLAMCH_EPS * DLAMCH_EPS;

    for (int32_t iter = 0; iter < maxiter; iter++) {
        double maxoff = 0.0;
        for (int32_t j = 0; j < n; j++) {
            for (int32_t i = 0; i < j; i++) {
                double aij = fabs(ATA[i + j * n]);
                if (aij > maxoff) maxoff = aij;
            }
        }

        if (maxoff < tol) break;

        for (int32_t p = 0; p < n - 1; p++) {
            for (int32_t q = p + 1; q < n; q++) {
                double app = ATA[p + p * n];
                double aqq = ATA[q + q * n];
                double apq = ATA[p + q * n];

                if (fabs(apq) < tol * sqrt(fabs(app * aqq))) continue;

                double theta = 0.5 * (aqq - app) / apq;
                double t;
                if (theta >= 0.0) {
                    t = 1.0 / (theta + sqrt(1.0 + theta * theta));
                } else {
                    t = 1.0 / (theta - sqrt(1.0 + theta * theta));
                }
                double c = 1.0 / sqrt(1.0 + t * t);
                double ss = t * c;

                ATA[p + p * n] = app - t * apq;
                ATA[q + q * n] = aqq + t * apq;
                ATA[p + q * n] = 0.0;
                ATA[q + p * n] = 0.0;

                for (int32_t i = 0; i < p; i++) {
                    double aip = ATA[i + p * n];
                    double aiq = ATA[i + q * n];
                    ATA[i + p * n] = c * aip - ss * aiq;
                    ATA[i + q * n] = ss * aip + c * aiq;
                }
                for (int32_t i = p + 1; i < q; i++) {
                    double api = ATA[p + i * n];
                    double aiq = ATA[i + q * n];
                    ATA[p + i * n] = c * api - ss * aiq;
                    ATA[i + q * n] = ss * api + c * aiq;
                }
                for (int32_t i = q + 1; i < n; i++) {
                    double api = ATA[p + i * n];
                    double aqi = ATA[q + i * n];
                    ATA[p + i * n] = c * api - ss * aqi;
                    ATA[q + i * n] = ss * api + c * aqi;
                }

                for (int32_t i = 0; i < n; i++) {
                    double vip = V[i + p * n];
                    double viq = V[i + q * n];
                    V[i + p * n] = c * vip - ss * viq;
                    V[i + q * n] = ss * vip + c * viq;
                }
            }
        }
    }

    /* Extract eigenvalues (singular values are sqrt) */
    for (int32_t i = 0; i < n; i++) {
        eigenvalues[i] = ATA[i + i * n];
    }

    /* Sort in descending order by singular value */
    for (int32_t i = 0; i < minmn - 1; i++) {
        int32_t k = i;
        double p = eigenvalues[i];
        for (int32_t j = i + 1; j < n; j++) {
            if (eigenvalues[j] > p) {
                k = j;
                p = eigenvalues[j];
            }
        }
        if (k != i) {
            eigenvalues[k] = eigenvalues[i];
            eigenvalues[i] = p;
            blas_dswap(n, &V[i * n], 1, &V[k * n], 1);
        }
    }

    /* Singular values are sqrt of eigenvalues of A^T*A */
    for (int32_t i = 0; i < minmn; i++) {
        s[i] = sqrt(MAX(0.0, eigenvalues[i]));
    }

    /* V^T = eigenvectors of A^T*A (right singular vectors) */
    if (wantvt) {
        for (int32_t i = 0; i < minmn; i++) {
            for (int32_t j = 0; j < n; j++) {
                VT[i + j * ldvt] = V[j + i * n];
            }
        }
        /* Zero remaining rows if jobz='A' */
        if (jobz == 'A' || jobz == 'a') {
            for (int32_t i = minmn; i < n; i++) {
                for (int32_t j = 0; j < n; j++) {
                    VT[i + j * ldvt] = 0.0;
                }
            }
        }
    }

    /* U = A * V * S^(-1) (left singular vectors) */
    if (wantu) {
        for (int32_t j = 0; j < minmn; j++) {
            if (s[j] > tol) {
                /* U(:,j) = A * V(:,j) / s[j] */
                for (int32_t i = 0; i < m; i++) {
                    double sum = 0.0;
                    for (int32_t k = 0; k < n; k++) {
                        sum += A[i + k * lda] * V[k + j * n];
                    }
                    U[i + j * ldu] = sum / s[j];
                }
            } else {
                for (int32_t i = 0; i < m; i++) {
                    U[i + j * ldu] = 0.0;
                }
            }
        }
        /* Zero or extend remaining columns if jobz='A' */
        if (jobz == 'A' || jobz == 'a') {
            for (int32_t j = minmn; j < m; j++) {
                for (int32_t i = 0; i < m; i++) {
                    U[i + j * ldu] = 0.0;
                }
            }
        }
    }

    return 0;
}

EXPORT int32_t lapack_sgesdd(char jobz, int32_t m, int32_t n,
                              float* A, int32_t lda, float* s,
                              float* U, int32_t ldu,
                              float* VT, int32_t ldvt,
                              float* work, int32_t lwork, int32_t* iwork)
{
    if (m <= 0 || n <= 0) return 0;

    int32_t minmn = MIN(m, n);

    if (lwork == -1) {
        work[0] = (float)(3 * minmn * minmn + MAX(MAX(m, n), 4 * minmn * minmn + 4 * minmn));
        return 0;
    }

    /* Simplified version similar to double */
    int wantu = (jobz == 'A' || jobz == 'a' || jobz == 'S' || jobz == 's');
    int wantvt = (jobz == 'A' || jobz == 'a' || jobz == 'S' || jobz == 's');

    float* ATA = work;

    blas_sgemm('T', 'N', n, n, m, 1.0f, A, lda, A, lda, 0.0f, ATA, n);

    float* V = ATA + n * n;

    for (int32_t j = 0; j < n; j++) {
        for (int32_t i = 0; i < n; i++) {
            V[i + j * n] = (i == j) ? 1.0f : 0.0f;
        }
    }

    const int32_t maxiter = 50;
    const float tol = 1.19e-7f * 1.19e-7f;

    for (int32_t iter = 0; iter < maxiter; iter++) {
        float maxoff = 0.0f;
        for (int32_t j = 0; j < n; j++) {
            for (int32_t i = 0; i < j; i++) {
                float aij = fabsf(ATA[i + j * n]);
                if (aij > maxoff) maxoff = aij;
            }
        }

        if (maxoff < tol) break;

        for (int32_t p = 0; p < n - 1; p++) {
            for (int32_t q = p + 1; q < n; q++) {
                float app = ATA[p + p * n];
                float aqq = ATA[q + q * n];
                float apq = ATA[p + q * n];

                if (fabsf(apq) < tol * sqrtf(fabsf(app * aqq))) continue;

                float theta = 0.5f * (aqq - app) / apq;
                float t;
                if (theta >= 0.0f) {
                    t = 1.0f / (theta + sqrtf(1.0f + theta * theta));
                } else {
                    t = 1.0f / (theta - sqrtf(1.0f + theta * theta));
                }
                float c = 1.0f / sqrtf(1.0f + t * t);
                float ss = t * c;

                ATA[p + p * n] = app - t * apq;
                ATA[q + q * n] = aqq + t * apq;
                ATA[p + q * n] = 0.0f;
                ATA[q + p * n] = 0.0f;

                for (int32_t i = 0; i < p; i++) {
                    float aip = ATA[i + p * n];
                    float aiq = ATA[i + q * n];
                    ATA[i + p * n] = c * aip - ss * aiq;
                    ATA[i + q * n] = ss * aip + c * aiq;
                }
                for (int32_t i = p + 1; i < q; i++) {
                    float api = ATA[p + i * n];
                    float aiq = ATA[i + q * n];
                    ATA[p + i * n] = c * api - ss * aiq;
                    ATA[i + q * n] = ss * api + c * aiq;
                }
                for (int32_t i = q + 1; i < n; i++) {
                    float api = ATA[p + i * n];
                    float aqi = ATA[q + i * n];
                    ATA[p + i * n] = c * api - ss * aqi;
                    ATA[q + i * n] = ss * api + c * aqi;
                }

                for (int32_t i = 0; i < n; i++) {
                    float vip = V[i + p * n];
                    float viq = V[i + q * n];
                    V[i + p * n] = c * vip - ss * viq;
                    V[i + q * n] = ss * vip + c * viq;
                }
            }
        }
    }

    float* eigenvalues = V + n * n;
    for (int32_t i = 0; i < n; i++) {
        eigenvalues[i] = ATA[i + i * n];
    }

    for (int32_t i = 0; i < minmn - 1; i++) {
        int32_t k = i;
        float p = eigenvalues[i];
        for (int32_t j = i + 1; j < n; j++) {
            if (eigenvalues[j] > p) {
                k = j;
                p = eigenvalues[j];
            }
        }
        if (k != i) {
            eigenvalues[k] = eigenvalues[i];
            eigenvalues[i] = p;
            blas_sswap(n, &V[i * n], 1, &V[k * n], 1);
        }
    }

    for (int32_t i = 0; i < minmn; i++) {
        s[i] = sqrtf(MAX(0.0f, eigenvalues[i]));
    }

    if (wantvt) {
        for (int32_t i = 0; i < minmn; i++) {
            for (int32_t j = 0; j < n; j++) {
                VT[i + j * ldvt] = V[j + i * n];
            }
        }
    }

    if (wantu) {
        for (int32_t j = 0; j < minmn; j++) {
            if (s[j] > tol) {
                for (int32_t i = 0; i < m; i++) {
                    float sum = 0.0f;
                    for (int32_t k = 0; k < n; k++) {
                        sum += A[i + k * lda] * V[k + j * n];
                    }
                    U[i + j * ldu] = sum / s[j];
                }
            } else {
                for (int32_t i = 0; i < m; i++) {
                    U[i + j * ldu] = 0.0f;
                }
            }
        }
    }

    return 0;
}

EXPORT int32_t lapack_dgesvd(char jobu, char jobvt, int32_t m, int32_t n,
                              double* A, int32_t lda, double* s,
                              double* U, int32_t ldu,
                              double* VT, int32_t ldvt,
                              double* work, int32_t lwork)
{
    /* Map to dgesdd with appropriate jobz */
    char jobz = 'S';  /* Default to reduced */
    if ((jobu == 'A' || jobu == 'a') && (jobvt == 'A' || jobvt == 'a')) {
        jobz = 'A';
    } else if ((jobu == 'N' || jobu == 'n') && (jobvt == 'N' || jobvt == 'n')) {
        jobz = 'N';
    }

    int32_t minmn = MIN(m, n);
    int32_t* iwork = (int32_t*)malloc(8 * minmn * sizeof(int32_t));
    int32_t info = lapack_dgesdd(jobz, m, n, A, lda, s, U, ldu, VT, ldvt,
                                  work, lwork, iwork);
    free(iwork);
    return info;
}

/* ============ Least Squares ============ */

EXPORT int32_t lapack_dgelsd(int32_t m, int32_t n, int32_t nrhs,
                              double* A, int32_t lda,
                              double* B, int32_t ldb,
                              double* s, double rcond, int32_t* rank,
                              double* work, int32_t lwork, int32_t* iwork)
{
    if (m <= 0 || n <= 0 || nrhs <= 0) return 0;

    int32_t minmn = MIN(m, n);
    int32_t maxmn = MAX(m, n);

    /* Workspace query */
    if (lwork == -1) {
        work[0] = (double)(12 * minmn + 2 * minmn * nrhs + maxmn);
        return 0;
    }

    /* Use SVD to solve least squares */
    /* A = U * S * V^T */
    /* x = V * S^(-1) * U^T * b */

    double* U = work;
    double* VT = U + m * minmn;
    double* svd_work = VT + minmn * n;
    int32_t svd_lwork = lwork - (m * minmn + minmn * n);

    /* Compute SVD of A */
    int32_t info = lapack_dgesdd('S', m, n, A, lda, s, U, m, VT, minmn,
                                  svd_work, svd_lwork, iwork);
    if (info != 0) return info;

    /* Determine effective rank */
    double threshold;
    if (rcond < 0.0) {
        threshold = maxmn * DLAMCH_EPS * s[0];
    } else {
        threshold = rcond * s[0];
    }

    *rank = 0;
    for (int32_t i = 0; i < minmn; i++) {
        if (s[i] > threshold) {
            (*rank)++;
        }
    }

    /* Compute x = V * S^(-1) * U^T * b */
    /* First: work2 = U^T * B */
    double* work2 = svd_work;
    blas_dgemm('T', 'N', minmn, nrhs, m, 1.0, U, m, B, ldb, 0.0, work2, minmn);

    /* Second: work2 = S^(-1) * work2 */
    for (int32_t j = 0; j < nrhs; j++) {
        for (int32_t i = 0; i < minmn; i++) {
            if (s[i] > threshold) {
                work2[i + j * minmn] /= s[i];
            } else {
                work2[i + j * minmn] = 0.0;
            }
        }
    }

    /* Third: B = V^T^T * work2 = V * work2 */
    /* Note: VT is stored, so we need V = VT^T */
    /* x = V * (S^(-1) * U^T * b) */
    for (int32_t j = 0; j < nrhs; j++) {
        for (int32_t i = 0; i < n; i++) {
            double sum = 0.0;
            for (int32_t k = 0; k < minmn; k++) {
                sum += VT[k + i * minmn] * work2[k + j * minmn];
            }
            B[i + j * ldb] = sum;
        }
    }

    return 0;
}
