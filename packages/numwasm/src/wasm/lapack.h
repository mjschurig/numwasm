#ifndef NUMJS_LAPACK_H
#define NUMJS_LAPACK_H

#include <stdint.h>
#include <stddef.h>

/*
 * NumJS LAPACK - Linear Algebra PACKage
 *
 * A minimal implementation of essential LAPACK routines for WebAssembly.
 * Follows LAPACK naming conventions:
 *   - d prefix = double precision real
 *   - s prefix = single precision real
 *   - z prefix = double precision complex
 *   - c prefix = single precision complex
 *
 * All matrices use column-major (Fortran) storage.
 * Return value: 0 = success, >0 = algorithm-specific error, <0 = invalid argument
 *
 * Reference: https://www.netlib.org/lapack/
 */

/* ============ LU Factorization ============ */

/**
 * Compute LU factorization with partial pivoting.
 * A = P * L * U
 *
 * On exit, L is stored below the diagonal (unit diagonal implied),
 * U is stored on and above the diagonal.
 *
 * @param m     Number of rows of A
 * @param n     Number of columns of A
 * @param A     On entry: matrix to factor (column-major, m x n)
 *              On exit: L and U factors
 * @param lda   Leading dimension of A (>= max(1, m))
 * @param ipiv  Output: pivot indices (size min(m,n))
 *              Row i was interchanged with row ipiv[i]
 * @return      0 on success
 *              >0 if U(i,i) = 0, factorization complete but U is singular
 */
int32_t lapack_dgetrf(int32_t m, int32_t n, double* A, int32_t lda,
                       int32_t* ipiv);

int32_t lapack_sgetrf(int32_t m, int32_t n, float* A, int32_t lda,
                       int32_t* ipiv);

/**
 * Solve a system of linear equations using LU factorization from dgetrf.
 * A * X = B  or  A^T * X = B
 *
 * @param trans 'N' = no transpose, 'T' = transpose, 'C' = conjugate transpose
 * @param n     Order of matrix A
 * @param nrhs  Number of right-hand sides (columns of B)
 * @param A     LU factorization from dgetrf (n x n)
 * @param lda   Leading dimension of A
 * @param ipiv  Pivot indices from dgetrf
 * @param B     On entry: right-hand side matrix B (n x nrhs)
 *              On exit: solution matrix X
 * @param ldb   Leading dimension of B
 * @return      0 on success
 */
int32_t lapack_dgetrs(char trans, int32_t n, int32_t nrhs,
                       const double* A, int32_t lda, const int32_t* ipiv,
                       double* B, int32_t ldb);

int32_t lapack_sgetrs(char trans, int32_t n, int32_t nrhs,
                       const float* A, int32_t lda, const int32_t* ipiv,
                       float* B, int32_t ldb);

/**
 * Compute the inverse of a matrix using LU factorization from dgetrf.
 *
 * @param n     Order of matrix A
 * @param A     On entry: LU factorization from dgetrf
 *              On exit: inverse of original matrix
 * @param lda   Leading dimension of A
 * @param ipiv  Pivot indices from dgetrf
 * @param work  Workspace array (size >= n)
 * @param lwork Size of work array. If lwork = -1, optimal size is returned in work[0]
 * @return      0 on success
 *              >0 if matrix is singular
 */
int32_t lapack_dgetri(int32_t n, double* A, int32_t lda, const int32_t* ipiv,
                       double* work, int32_t lwork);

int32_t lapack_sgetri(int32_t n, float* A, int32_t lda, const int32_t* ipiv,
                       float* work, int32_t lwork);

/* ============ Cholesky Factorization ============ */

/**
 * Compute Cholesky factorization of symmetric positive definite matrix.
 * A = L * L^T  (if uplo = 'L')
 * A = U^T * U  (if uplo = 'U')
 *
 * @param uplo  'U' = upper triangle stored, compute U
 *              'L' = lower triangle stored, compute L
 * @param n     Order of matrix A
 * @param A     On entry: symmetric positive definite matrix
 *              On exit: Cholesky factor L or U
 * @param lda   Leading dimension of A
 * @return      0 on success
 *              >0 if leading minor of order i is not positive definite
 */
int32_t lapack_dpotrf(char uplo, int32_t n, double* A, int32_t lda);

int32_t lapack_spotrf(char uplo, int32_t n, float* A, int32_t lda);

/**
 * Solve A * X = B using Cholesky factorization from dpotrf.
 *
 * @param uplo  'U' or 'L', must match dpotrf call
 * @param n     Order of matrix A
 * @param nrhs  Number of right-hand sides
 * @param A     Cholesky factorization from dpotrf
 * @param lda   Leading dimension of A
 * @param B     On entry: right-hand side B
 *              On exit: solution X
 * @param ldb   Leading dimension of B
 * @return      0 on success
 */
int32_t lapack_dpotrs(char uplo, int32_t n, int32_t nrhs,
                       const double* A, int32_t lda,
                       double* B, int32_t ldb);

/**
 * Compute inverse using Cholesky factorization from dpotrf.
 *
 * @param uplo  'U' or 'L', must match dpotrf call
 * @param n     Order of matrix A
 * @param A     On entry: Cholesky factor
 *              On exit: inverse (same triangle as input)
 * @param lda   Leading dimension of A
 * @return      0 on success
 */
int32_t lapack_dpotri(char uplo, int32_t n, double* A, int32_t lda);

/* ============ QR Factorization ============ */

/**
 * Compute QR factorization using Householder reflectors.
 * A = Q * R
 *
 * Q is represented as a product of elementary reflectors:
 * Q = H(1) * H(2) * ... * H(k), where k = min(m, n)
 * Each H(i) = I - tau[i] * v * v^T
 *
 * @param m     Number of rows of A
 * @param n     Number of columns of A
 * @param A     On entry: matrix to factor (m x n)
 *              On exit: R in upper triangle, Householder vectors below
 * @param lda   Leading dimension of A
 * @param tau   Output: scalar factors of reflectors (size min(m,n))
 * @param work  Workspace
 * @param lwork Size of workspace. If -1, optimal size returned in work[0]
 * @return      0 on success
 */
int32_t lapack_dgeqrf(int32_t m, int32_t n, double* A, int32_t lda,
                       double* tau, double* work, int32_t lwork);

int32_t lapack_sgeqrf(int32_t m, int32_t n, float* A, int32_t lda,
                       float* tau, float* work, int32_t lwork);

/**
 * Generate orthogonal matrix Q from Householder reflectors.
 *
 * @param m     Number of rows of Q
 * @param n     Number of columns of Q (n <= m)
 * @param k     Number of reflectors (k <= n)
 * @param A     On entry: Householder vectors from dgeqrf
 *              On exit: orthogonal matrix Q
 * @param lda   Leading dimension of A
 * @param tau   Scalar factors from dgeqrf
 * @param work  Workspace
 * @param lwork Size of workspace. If -1, optimal size returned in work[0]
 * @return      0 on success
 */
int32_t lapack_dorgqr(int32_t m, int32_t n, int32_t k,
                       double* A, int32_t lda, const double* tau,
                       double* work, int32_t lwork);

int32_t lapack_sorgqr(int32_t m, int32_t n, int32_t k,
                       float* A, int32_t lda, const float* tau,
                       float* work, int32_t lwork);

/**
 * Solve triangular system with single/multiple RHS.
 * op(A) * X = B  or  X * op(A) = B
 *
 * @param uplo  'U' = upper triangular, 'L' = lower triangular
 * @param trans 'N' = no transpose, 'T' = transpose
 * @param diag  'N' = non-unit diagonal, 'U' = unit diagonal
 * @param n     Order of A
 * @param nrhs  Number of right-hand sides
 * @param A     Triangular matrix (n x n)
 * @param lda   Leading dimension of A
 * @param B     On entry: RHS. On exit: solution (n x nrhs)
 * @param ldb   Leading dimension of B
 * @return      0 on success
 */
int32_t lapack_dtrtrs(char uplo, char trans, char diag,
                       int32_t n, int32_t nrhs,
                       const double* A, int32_t lda,
                       double* B, int32_t ldb);

/* ============ Eigenvalue Decomposition ============ */

/**
 * Compute eigenvalues and optionally eigenvectors of general matrix.
 * A * v = lambda * v
 *
 * For real matrices, eigenvalues may be complex. Complex conjugate pairs
 * are stored consecutively: wr[j]+i*wi[j] and wr[j+1]+i*wi[j+1] = wr[j]-i*wi[j]
 *
 * @param jobvl 'N' = no left eigenvectors, 'V' = compute left eigenvectors
 * @param jobvr 'N' = no right eigenvectors, 'V' = compute right eigenvectors
 * @param n     Order of matrix A
 * @param A     On entry: matrix A (destroyed on exit)
 * @param lda   Leading dimension of A
 * @param wr    Output: real parts of eigenvalues (size n)
 * @param wi    Output: imaginary parts of eigenvalues (size n)
 * @param VL    Output: left eigenvectors if jobvl='V' (n x n)
 * @param ldvl  Leading dimension of VL
 * @param VR    Output: right eigenvectors if jobvr='V' (n x n)
 * @param ldvr  Leading dimension of VR
 * @param work  Workspace
 * @param lwork Size of workspace. If -1, optimal size returned in work[0]
 * @return      0 on success
 *              >0 if QR algorithm failed to converge
 */
int32_t lapack_dgeev(char jobvl, char jobvr, int32_t n,
                      double* A, int32_t lda,
                      double* wr, double* wi,
                      double* VL, int32_t ldvl,
                      double* VR, int32_t ldvr,
                      double* work, int32_t lwork);

int32_t lapack_sgeev(char jobvl, char jobvr, int32_t n,
                      float* A, int32_t lda,
                      float* wr, float* wi,
                      float* VL, int32_t ldvl,
                      float* VR, int32_t ldvr,
                      float* work, int32_t lwork);

/**
 * Compute eigenvalues and eigenvectors of symmetric matrix.
 * A * v = lambda * v
 *
 * Uses divide-and-conquer algorithm (more efficient for large matrices).
 * Eigenvalues are always real and returned in ascending order.
 *
 * @param jobz  'N' = eigenvalues only, 'V' = eigenvalues and eigenvectors
 * @param uplo  'U' = upper triangle stored, 'L' = lower triangle stored
 * @param n     Order of matrix A
 * @param A     On entry: symmetric matrix (only uplo triangle used)
 *              On exit: eigenvectors (columns) if jobz='V'
 * @param lda   Leading dimension of A
 * @param w     Output: eigenvalues in ascending order (size n)
 * @param work  Workspace
 * @param lwork Size of workspace. If -1, optimal size returned in work[0]
 * @param iwork Integer workspace
 * @param liwork Size of integer workspace. If -1, optimal size returned in iwork[0]
 * @return      0 on success
 *              >0 if algorithm failed to converge
 */
int32_t lapack_dsyevd(char jobz, char uplo, int32_t n,
                       double* A, int32_t lda, double* w,
                       double* work, int32_t lwork,
                       int32_t* iwork, int32_t liwork);

int32_t lapack_ssyevd(char jobz, char uplo, int32_t n,
                       float* A, int32_t lda, float* w,
                       float* work, int32_t lwork,
                       int32_t* iwork, int32_t liwork);

/* ============ Singular Value Decomposition ============ */

/**
 * Compute SVD using divide-and-conquer algorithm.
 * A = U * S * V^T
 *
 * @param jobz  'A' = all m columns of U and all n rows of V^T
 *              'S' = first min(m,n) columns of U and rows of V^T
 *              'O' = overwrite A with U if m >= n, else with V^T
 *              'N' = no columns of U or rows of V^T computed
 * @param m     Number of rows of A
 * @param n     Number of columns of A
 * @param A     On entry: matrix to decompose (m x n)
 *              On exit: overwritten (depends on jobz)
 * @param lda   Leading dimension of A
 * @param s     Output: singular values in descending order (size min(m,n))
 * @param U     Output: left singular vectors (m x ucol, where ucol depends on jobz)
 * @param ldu   Leading dimension of U
 * @param VT    Output: right singular vectors transposed (vtrow x n)
 * @param ldvt  Leading dimension of VT
 * @param work  Workspace
 * @param lwork Size of workspace. If -1, optimal size returned in work[0]
 * @param iwork Integer workspace (size >= 8*min(m,n))
 * @return      0 on success
 *              >0 if algorithm did not converge
 */
int32_t lapack_dgesdd(char jobz, int32_t m, int32_t n,
                       double* A, int32_t lda, double* s,
                       double* U, int32_t ldu,
                       double* VT, int32_t ldvt,
                       double* work, int32_t lwork, int32_t* iwork);

int32_t lapack_sgesdd(char jobz, int32_t m, int32_t n,
                       float* A, int32_t lda, float* s,
                       float* U, int32_t ldu,
                       float* VT, int32_t ldvt,
                       float* work, int32_t lwork, int32_t* iwork);

/**
 * Compute SVD using standard algorithm (more robust but slower).
 * A = U * S * V^T
 *
 * @param jobu  'A' = all m columns of U
 *              'S' = first min(m,n) columns of U
 *              'O' = first min(m,n) columns overwrite A
 *              'N' = no columns of U
 * @param jobvt 'A' = all n rows of V^T
 *              'S' = first min(m,n) rows of V^T
 *              'O' = first min(m,n) rows overwrite A
 *              'N' = no rows of V^T
 */
int32_t lapack_dgesvd(char jobu, char jobvt, int32_t m, int32_t n,
                       double* A, int32_t lda, double* s,
                       double* U, int32_t ldu,
                       double* VT, int32_t ldvt,
                       double* work, int32_t lwork);

/* ============ Least Squares ============ */

/**
 * Solve overdetermined or underdetermined linear systems using SVD.
 * Minimize ||b - A*x||_2
 *
 * @param m     Number of rows of A
 * @param n     Number of columns of A
 * @param nrhs  Number of right-hand sides
 * @param A     Matrix A (m x n, destroyed)
 * @param lda   Leading dimension of A
 * @param B     On entry: right-hand side (max(m,n) x nrhs)
 *              On exit: solution x (first n rows if m >= n)
 * @param ldb   Leading dimension of B
 * @param s     Output: singular values of A in descending order
 * @param rcond Singular values below rcond*s[0] are treated as zero
 *              Use -1 for machine precision default
 * @param rank  Output: effective rank of A
 * @param work  Workspace
 * @param lwork Size of workspace. If -1, optimal size returned in work[0]
 * @param iwork Integer workspace
 * @return      0 on success
 *              >0 if SVD did not converge
 */
int32_t lapack_dgelsd(int32_t m, int32_t n, int32_t nrhs,
                       double* A, int32_t lda,
                       double* B, int32_t ldb,
                       double* s, double rcond, int32_t* rank,
                       double* work, int32_t lwork, int32_t* iwork);

/* ============ Schur Decomposition ============ */

/**
 * Compute Schur decomposition: A = Z * T * Z^T
 * where T is upper quasi-triangular and Z is orthogonal.
 *
 * @param jobvs 'N' = no Schur vectors, 'V' = compute Schur vectors
 * @param sort  'N' = no sorting of eigenvalues
 * @param n     Order of matrix A
 * @param A     On entry: matrix. On exit: Schur form T
 * @param lda   Leading dimension of A
 * @param sdim  Output: number of selected eigenvalues (0 if sort='N')
 * @param wr    Output: real parts of eigenvalues
 * @param wi    Output: imaginary parts of eigenvalues
 * @param VS    Output: Schur vectors if jobvs='V'
 * @param ldvs  Leading dimension of VS
 * @param work  Workspace
 * @param lwork Size of workspace. If -1, optimal size returned in work[0]
 * @return      0 on success
 */
int32_t lapack_dgees(char jobvs, char sort, int32_t n,
                      double* A, int32_t lda, int32_t* sdim,
                      double* wr, double* wi,
                      double* VS, int32_t ldvs,
                      double* work, int32_t lwork);

/* ============ Symmetric Indefinite Factorization ============ */

/**
 * Compute Bunch-Kaufman factorization of symmetric matrix.
 * A = P^T * L * D * L^T * P (lower) or A = P^T * U * D * U^T * P (upper)
 *
 * @param uplo  'U' = upper, 'L' = lower
 * @param n     Order of matrix A
 * @param A     On entry: symmetric matrix. On exit: L/U and D factors
 * @param lda   Leading dimension of A
 * @param ipiv  Output: pivot indices
 * @param work  Workspace
 * @param lwork Size of workspace. If -1, optimal size returned in work[0]
 * @return      0 on success
 */
int32_t lapack_dsytrf(char uplo, int32_t n, double* A, int32_t lda,
                       int32_t* ipiv, double* work, int32_t lwork);

/* ============ Triangular Inverse ============ */

/**
 * Compute inverse of triangular matrix.
 *
 * @param uplo  'U' = upper triangular, 'L' = lower triangular
 * @param diag  'N' = non-unit diagonal, 'U' = unit diagonal
 * @param n     Order of matrix A
 * @param A     On entry: triangular matrix
 *              On exit: inverse
 * @param lda   Leading dimension of A
 * @return      0 on success
 *              >0 if A(i,i) = 0 (singular)
 */
int32_t lapack_dtrtri(char uplo, char diag, int32_t n,
                       double* A, int32_t lda);

/* ============ Utility Functions ============ */

/**
 * Apply row permutation to matrix.
 * Used for applying pivot sequence from LU factorization.
 *
 * @param n     Number of columns of A
 * @param A     Matrix to permute (modified in place)
 * @param lda   Leading dimension of A
 * @param k1    First element of ipiv to use (0-based)
 * @param k2    Last element of ipiv to use (0-based)
 * @param ipiv  Pivot indices
 * @param incx  Increment for ipiv (1 = forward, -1 = backward)
 */
void lapack_dlaswp(int32_t n, double* A, int32_t lda,
                    int32_t k1, int32_t k2, const int32_t* ipiv, int32_t incx);

/**
 * Compute machine parameters.
 */
double lapack_dlamch(char cmach);

#endif /* NUMJS_LAPACK_H */
