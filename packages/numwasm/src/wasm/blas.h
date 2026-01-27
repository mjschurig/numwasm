#ifndef NUMJS_BLAS_H
#define NUMJS_BLAS_H

#include <stdint.h>
#include <stddef.h>

/*
 * NumJS BLAS - Basic Linear Algebra Subprograms
 *
 * A minimal implementation of BLAS routines for WebAssembly.
 * Follows BLAS naming conventions:
 *   - d prefix = double precision
 *   - s prefix = single precision
 *
 * Reference: https://www.netlib.org/blas/
 */

/* ============ Level 1 BLAS (Vector Operations) ============ */

/**
 * Dot product of two vectors.
 * result = sum(x[i] * y[i])
 *
 * @param n    Number of elements
 * @param x    First vector
 * @param incx Stride for x (usually 1)
 * @param y    Second vector
 * @param incy Stride for y (usually 1)
 * @return     Dot product
 */
double blas_ddot(int32_t n, const double* x, int32_t incx,
                 const double* y, int32_t incy);

float blas_sdot(int32_t n, const float* x, int32_t incx,
                const float* y, int32_t incy);

/**
 * Euclidean norm (L2 norm) of a vector.
 * result = sqrt(sum(x[i]^2))
 *
 * Uses scaled computation to avoid overflow/underflow.
 *
 * @param n    Number of elements
 * @param x    Vector
 * @param incx Stride for x
 * @return     L2 norm
 */
double blas_dnrm2(int32_t n, const double* x, int32_t incx);

float blas_snrm2(int32_t n, const float* x, int32_t incx);

/**
 * Sum of absolute values (L1 norm).
 * result = sum(|x[i]|)
 *
 * @param n    Number of elements
 * @param x    Vector
 * @param incx Stride for x
 * @return     L1 norm
 */
double blas_dasum(int32_t n, const double* x, int32_t incx);

float blas_sasum(int32_t n, const float* x, int32_t incx);

/**
 * Index of element with maximum absolute value.
 * Returns 0-based index (unlike FORTRAN BLAS which returns 1-based).
 *
 * @param n    Number of elements
 * @param x    Vector
 * @param incx Stride for x
 * @return     Index of max |x[i]|, or -1 if n <= 0
 */
int32_t blas_idamax(int32_t n, const double* x, int32_t incx);

int32_t blas_isamax(int32_t n, const float* x, int32_t incx);

/**
 * Scale a vector by a constant.
 * x = alpha * x
 *
 * @param n     Number of elements
 * @param alpha Scalar multiplier
 * @param x     Vector (modified in place)
 * @param incx  Stride for x
 */
void blas_dscal(int32_t n, double alpha, double* x, int32_t incx);

void blas_sscal(int32_t n, float alpha, float* x, int32_t incx);

/**
 * Vector addition with scaling.
 * y = alpha * x + y
 *
 * @param n     Number of elements
 * @param alpha Scalar multiplier for x
 * @param x     Source vector
 * @param incx  Stride for x
 * @param y     Destination vector (modified in place)
 * @param incy  Stride for y
 */
void blas_daxpy(int32_t n, double alpha, const double* x, int32_t incx,
                double* y, int32_t incy);

void blas_saxpy(int32_t n, float alpha, const float* x, int32_t incx,
                float* y, int32_t incy);

/**
 * Copy a vector.
 * y = x
 *
 * @param n    Number of elements
 * @param x    Source vector
 * @param incx Stride for x
 * @param y    Destination vector
 * @param incy Stride for y
 */
void blas_dcopy(int32_t n, const double* x, int32_t incx,
                double* y, int32_t incy);

void blas_scopy(int32_t n, const float* x, int32_t incx,
                float* y, int32_t incy);

/**
 * Swap two vectors.
 * x <-> y
 *
 * @param n    Number of elements
 * @param x    First vector (modified)
 * @param incx Stride for x
 * @param y    Second vector (modified)
 * @param incy Stride for y
 */
void blas_dswap(int32_t n, double* x, int32_t incx,
                double* y, int32_t incy);

void blas_sswap(int32_t n, float* x, int32_t incx,
                float* y, int32_t incy);

/* ============ Level 2 BLAS (Matrix-Vector Operations) ============ */

/**
 * General matrix-vector multiplication.
 * y = alpha * op(A) * x + beta * y
 *
 * @param trans 'N' = A, 'T' = A^T, 'C' = A^H (conjugate transpose)
 * @param m     Number of rows of A
 * @param n     Number of columns of A
 * @param alpha Scalar multiplier for A*x
 * @param A     Matrix A (column-major, m x n)
 * @param lda   Leading dimension of A (>= m)
 * @param x     Vector x
 * @param incx  Stride for x
 * @param beta  Scalar multiplier for y
 * @param y     Vector y (modified in place)
 * @param incy  Stride for y
 */
void blas_dgemv(char trans, int32_t m, int32_t n,
                double alpha, const double* A, int32_t lda,
                const double* x, int32_t incx,
                double beta, double* y, int32_t incy);

void blas_sgemv(char trans, int32_t m, int32_t n,
                float alpha, const float* A, int32_t lda,
                const float* x, int32_t incx,
                float beta, float* y, int32_t incy);

/**
 * Triangular matrix-vector multiplication.
 * x = op(A) * x
 *
 * @param uplo  'U' = upper triangular, 'L' = lower triangular
 * @param trans 'N' = A, 'T' = A^T, 'C' = A^H
 * @param diag  'N' = non-unit diagonal, 'U' = unit diagonal
 * @param n     Order of matrix A
 * @param A     Triangular matrix (column-major)
 * @param lda   Leading dimension of A
 * @param x     Vector (modified in place)
 * @param incx  Stride for x
 */
void blas_dtrmv(char uplo, char trans, char diag, int32_t n,
                const double* A, int32_t lda,
                double* x, int32_t incx);

/**
 * Triangular solve.
 * Solve op(A) * x = b for x
 *
 * @param uplo  'U' = upper triangular, 'L' = lower triangular
 * @param trans 'N' = A, 'T' = A^T, 'C' = A^H
 * @param diag  'N' = non-unit diagonal, 'U' = unit diagonal
 * @param n     Order of matrix A
 * @param A     Triangular matrix (column-major)
 * @param lda   Leading dimension of A
 * @param x     On entry: b. On exit: solution x (modified in place)
 * @param incx  Stride for x
 */
void blas_dtrsv(char uplo, char trans, char diag, int32_t n,
                const double* A, int32_t lda,
                double* x, int32_t incx);

/**
 * Rank-1 update.
 * A = alpha * x * y^T + A
 *
 * @param m     Number of rows of A
 * @param n     Number of columns of A
 * @param alpha Scalar multiplier
 * @param x     Vector x (length m)
 * @param incx  Stride for x
 * @param y     Vector y (length n)
 * @param incy  Stride for y
 * @param A     Matrix A (modified in place, column-major)
 * @param lda   Leading dimension of A
 */
void blas_dger(int32_t m, int32_t n, double alpha,
               const double* x, int32_t incx,
               const double* y, int32_t incy,
               double* A, int32_t lda);

/* ============ Level 3 BLAS (Matrix-Matrix Operations) ============ */

/**
 * General matrix-matrix multiplication.
 * C = alpha * op(A) * op(B) + beta * C
 *
 * This is the most critical BLAS routine for performance.
 *
 * @param transA 'N' = A, 'T' = A^T, 'C' = A^H
 * @param transB 'N' = B, 'T' = B^T, 'C' = B^H
 * @param m      Number of rows of op(A) and C
 * @param n      Number of columns of op(B) and C
 * @param k      Number of columns of op(A), rows of op(B)
 * @param alpha  Scalar multiplier for A*B
 * @param A      Matrix A (column-major)
 * @param lda    Leading dimension of A
 * @param B      Matrix B (column-major)
 * @param ldb    Leading dimension of B
 * @param beta   Scalar multiplier for C
 * @param C      Matrix C (modified in place, column-major)
 * @param ldc    Leading dimension of C
 */
void blas_dgemm(char transA, char transB,
                int32_t m, int32_t n, int32_t k,
                double alpha,
                const double* A, int32_t lda,
                const double* B, int32_t ldb,
                double beta,
                double* C, int32_t ldc);

void blas_sgemm(char transA, char transB,
                int32_t m, int32_t n, int32_t k,
                float alpha,
                const float* A, int32_t lda,
                const float* B, int32_t ldb,
                float beta,
                float* C, int32_t ldc);

/**
 * Triangular matrix-matrix multiplication.
 * B = alpha * op(A) * B   (side = 'L')
 * B = alpha * B * op(A)   (side = 'R')
 *
 * @param side   'L' = left multiply, 'R' = right multiply
 * @param uplo   'U' = upper triangular, 'L' = lower triangular
 * @param transA 'N' = A, 'T' = A^T, 'C' = A^H
 * @param diag   'N' = non-unit diagonal, 'U' = unit diagonal
 * @param m      Number of rows of B
 * @param n      Number of columns of B
 * @param alpha  Scalar multiplier
 * @param A      Triangular matrix (column-major)
 * @param lda    Leading dimension of A
 * @param B      Matrix B (modified in place, column-major)
 * @param ldb    Leading dimension of B
 */
void blas_dtrmm(char side, char uplo, char transA, char diag,
                int32_t m, int32_t n,
                double alpha,
                const double* A, int32_t lda,
                double* B, int32_t ldb);

/**
 * Triangular solve with multiple right-hand sides.
 * Solve op(A) * X = alpha * B   (side = 'L')
 * Solve X * op(A) = alpha * B   (side = 'R')
 * Solution X overwrites B.
 *
 * @param side   'L' = left solve, 'R' = right solve
 * @param uplo   'U' = upper triangular, 'L' = lower triangular
 * @param transA 'N' = A, 'T' = A^T, 'C' = A^H
 * @param diag   'N' = non-unit diagonal, 'U' = unit diagonal
 * @param m      Number of rows of B
 * @param n      Number of columns of B
 * @param alpha  Scalar multiplier
 * @param A      Triangular matrix (column-major)
 * @param lda    Leading dimension of A
 * @param B      On entry: RHS. On exit: solution X (column-major)
 * @param ldb    Leading dimension of B
 */
void blas_dtrsm(char side, char uplo, char transA, char diag,
                int32_t m, int32_t n,
                double alpha,
                const double* A, int32_t lda,
                double* B, int32_t ldb);

/**
 * Symmetric rank-k update.
 * C = alpha * A * A^T + beta * C   (trans = 'N')
 * C = alpha * A^T * A + beta * C   (trans = 'T')
 *
 * @param uplo  'U' = upper triangle of C, 'L' = lower triangle
 * @param trans 'N' = A * A^T, 'T' = A^T * A
 * @param n     Order of C
 * @param k     Columns of A (if trans='N') or rows (if trans='T')
 * @param alpha Scalar multiplier for A * A^T
 * @param A     Matrix A (column-major)
 * @param lda   Leading dimension of A
 * @param beta  Scalar multiplier for C
 * @param C     Symmetric matrix C (modified in place)
 * @param ldc   Leading dimension of C
 */
void blas_dsyrk(char uplo, char trans, int32_t n, int32_t k,
                double alpha, const double* A, int32_t lda,
                double beta, double* C, int32_t ldc);

#endif /* NUMJS_BLAS_H */
