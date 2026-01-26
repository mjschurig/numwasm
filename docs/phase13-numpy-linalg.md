# Phase 13: numpy.linalg Implementation Plan

Complete implementation roadmap for the NumJS-WASM linear algebra module, providing NumPy-compatible matrix operations with WebAssembly acceleration.

---

## ⚠️ Implementation Guideline

**The original NumPy source code is available in `/numpy`.**

Key reference files:
- `numpy/linalg/_linalg.py` - Main Python interface (114KB)
- `numpy/linalg/umath_linalg.cpp` - LAPACK ufunc wrappers (140KB)
- `numpy/linalg/lapack_lite/` - LAPACK routines (f2c translated)

Implementation should follow NumPy's algorithms and error handling for consistency.

---

## Current State (Pre-Phase 13)

```
src/wasm/
├── ndarray.h/c        # Core NDArray with views, slicing
├── dtype.h/c          # DType system
├── broadcast.h/c      # Broadcasting
├── indexing.h/c       # Index operations
├── pairwise_sum.h/c   # Accurate summation
└── logic.c            # Logical operations

src/ts/
├── NDArray.ts         # Core array class
├── types.ts           # Type definitions
├── dtype.ts           # Type utilities
├── broadcast.ts       # Broadcasting functions
├── indexing.ts        # Index operations
├── slice.ts           # Slicing utilities
├── iterators.ts       # Iterators
└── index.ts           # Public exports
```

**Existing Infrastructure Used:**
- NDArray with contiguity flags (C_CONTIGUOUS, F_CONTIGUOUS)
- DType system with Float32, Float64, Complex64, Complex128
- Broadcasting for shape compatibility
- View system for efficient memory sharing

---

## Phase 13 Dependency Tree

```
PHASE 13: NUMPY.LINALG
│
├── 13.1 BLAS Foundation (C/WASM)
│   ├── 13.1.1 Level 1 BLAS (vector operations)
│   │   ├── dot(x, y) → scalar
│   │   ├── nrm2(x) → scalar (L2 norm)
│   │   ├── asum(x) → scalar (L1 norm)
│   │   ├── iamax(x) → index (max absolute)
│   │   ├── scal(alpha, x) → void (scale)
│   │   ├── axpy(alpha, x, y) → void (y += alpha*x)
│   │   └── copy(x, y) → void
│   │
│   ├── 13.1.2 Level 2 BLAS (matrix-vector)
│   │   ├── gemv(A, x, y) → y = alpha*A*x + beta*y
│   │   ├── trmv(A, x) → x = A*x (triangular)
│   │   └── trsv(A, x) → solve A*x = b (triangular)
│   │
│   └── 13.1.3 Level 3 BLAS (matrix-matrix)
│       ├── gemm(A, B, C) → C = alpha*A*B + beta*C
│       ├── trmm(A, B) → B = A*B (triangular)
│       └── trsm(A, B) → solve A*X = B (triangular)
│
│   Dependencies: NDArray core, DType system
│
├── 13.2 LAPACK Decompositions (C/WASM)
│   ├── 13.2.1 LU Factorization
│   │   ├── getrf(A, ipiv) → PA = LU
│   │   ├── getrs(LU, ipiv, B) → solve using LU
│   │   └── getri(LU, ipiv) → inverse from LU
│   │
│   ├── 13.2.2 Cholesky Factorization
│   │   ├── potrf(A, uplo) → A = L*L^T or U^T*U
│   │   ├── potrs(L, B) → solve using Cholesky
│   │   └── potri(L) → inverse from Cholesky
│   │
│   ├── 13.2.3 QR Factorization
│   │   ├── geqrf(A, tau) → A = Q*R (Householder)
│   │   ├── orgqr(A, tau) → construct Q matrix
│   │   └── trtrs(R, B) → solve R*X = B
│   │
│   ├── 13.2.4 Eigenvalue Decomposition
│   │   ├── geev(A) → general eigenvalues/vectors
│   │   ├── syev/heev(A) → symmetric/Hermitian
│   │   └── syevd/heevd(A) → divide-and-conquer
│   │
│   └── 13.2.5 Singular Value Decomposition
│       ├── gesvd(A) → standard SVD
│       └── gesdd(A) → divide-and-conquer SVD
│
│   Dependencies: 13.1.* (BLAS)
│
├── 13.3 High-Level Linear Algebra (TypeScript + C)
│   ├── 13.3.1 Matrix Products
│   │   ├── dot(a, b) → various dot products
│   │   ├── vdot(a, b) → vector dot (flattened)
│   │   ├── inner(a, b) → inner product
│   │   ├── outer(a, b) → outer product
│   │   ├── matmul(a, b) → matrix multiplication
│   │   ├── tensordot(a, b, axes) → tensor contraction
│   │   ├── einsum(subscripts, *ops) → Einstein summation
│   │   ├── kron(a, b) → Kronecker product
│   │   └── multi_dot(arrays) → optimized chain multiply
│   │
│   ├── 13.3.2 Decompositions (user-facing)
│   │   ├── cholesky(a, upper) → L or U
│   │   ├── qr(a, mode) → Q, R
│   │   ├── svd(a, full_matrices, compute_uv) → U, S, Vh
│   │   └── svdvals(a) → singular values only
│   │
│   ├── 13.3.3 Eigenvalues
│   │   ├── eig(a) → eigenvalues, eigenvectors
│   │   ├── eigh(a, UPLO) → Hermitian eigenvalues
│   │   ├── eigvals(a) → eigenvalues only
│   │   └── eigvalsh(a, UPLO) → Hermitian eigenvalues only
│   │
│   ├── 13.3.4 Norms & Numbers
│   │   ├── norm(x, ord, axis) → vector/matrix norm
│   │   ├── matrix_norm(x, ord) → matrix norm
│   │   ├── vector_norm(x, ord, axis) → vector norm
│   │   ├── cond(x, p) → condition number
│   │   ├── det(a) → determinant
│   │   ├── slogdet(a) → sign and log determinant
│   │   ├── matrix_rank(A, tol) → numerical rank
│   │   └── trace(a, offset, axis1, axis2) → sum of diagonal
│   │
│   ├── 13.3.5 Solving & Inverting
│   │   ├── solve(a, b) → solve a*x = b
│   │   ├── tensorsolve(a, b, axes) → tensor equation
│   │   ├── lstsq(a, b, rcond) → least squares
│   │   ├── inv(a) → matrix inverse
│   │   ├── pinv(a, rcond) → pseudo-inverse
│   │   └── tensorinv(a, ind) → tensor inverse
│   │
│   └── 13.3.6 Matrix Operations
│       ├── matrix_power(a, n) → a^n
│       ├── diagonal(a, offset, axis1, axis2)  ← ALREADY EXISTS
│       ├── matrix_transpose(x) → swap last two axes
│       └── cross(a, b, axis) → cross product
│
│   Dependencies: 13.2.* (LAPACK)
│
└── 13.4 Error Handling & Utilities
    ├── 13.4.1 LinAlgError exception class
    ├── 13.4.2 Input validation helpers
    ├── 13.4.3 Type promotion for linalg
    └── 13.4.4 Result type containers (namedtuple-like)

    Dependencies: Core types
```

---

## Detailed Implementation Specifications

### 13.1 BLAS Foundation

#### 13.1.1 Level 1 BLAS (Vector Operations)

**File:** `src/wasm/blas.h` (new file)

```c
#ifndef NUMJS_BLAS_H
#define NUMJS_BLAS_H

#include "ndarray.h"

/* ============ Level 1 BLAS ============ */

/**
 * Dot product of two vectors.
 * x and y must have same length n.
 */
double blas_ddot(int32_t n, const double* x, int32_t incx,
                 const double* y, int32_t incy);

float blas_sdot(int32_t n, const float* x, int32_t incx,
                const float* y, int32_t incy);

/**
 * Complex dot product (conjugated): result = sum(conj(x[i]) * y[i])
 */
void blas_zdotc(int32_t n, const double* x, int32_t incx,
                const double* y, int32_t incy, double* result);

void blas_cdotc(int32_t n, const float* x, int32_t incx,
                const float* y, int32_t incy, float* result);

/**
 * L2 norm (Euclidean): sqrt(sum(|x[i]|^2))
 */
double blas_dnrm2(int32_t n, const double* x, int32_t incx);
float blas_snrm2(int32_t n, const float* x, int32_t incx);
double blas_dznrm2(int32_t n, const double* x, int32_t incx);  /* complex */

/**
 * L1 norm (sum of absolute values): sum(|x[i]|)
 */
double blas_dasum(int32_t n, const double* x, int32_t incx);
float blas_sasum(int32_t n, const float* x, int32_t incx);
double blas_dzasum(int32_t n, const double* x, int32_t incx);  /* complex */

/**
 * Index of maximum absolute value element (1-based for BLAS compat)
 */
int32_t blas_idamax(int32_t n, const double* x, int32_t incx);
int32_t blas_isamax(int32_t n, const float* x, int32_t incx);
int32_t blas_izamax(int32_t n, const double* x, int32_t incx);  /* complex */

/**
 * Scale vector: x = alpha * x
 */
void blas_dscal(int32_t n, double alpha, double* x, int32_t incx);
void blas_sscal(int32_t n, float alpha, float* x, int32_t incx);
void blas_zscal(int32_t n, const double* alpha, double* x, int32_t incx);

/**
 * Vector addition: y = alpha * x + y
 */
void blas_daxpy(int32_t n, double alpha, const double* x, int32_t incx,
                double* y, int32_t incy);
void blas_saxpy(int32_t n, float alpha, const float* x, int32_t incx,
                float* y, int32_t incy);
void blas_zaxpy(int32_t n, const double* alpha, const double* x, int32_t incx,
                double* y, int32_t incy);

/**
 * Vector copy: y = x
 */
void blas_dcopy(int32_t n, const double* x, int32_t incx,
                double* y, int32_t incy);
void blas_scopy(int32_t n, const float* x, int32_t incx,
                float* y, int32_t incy);
void blas_zcopy(int32_t n, const double* x, int32_t incx,
                double* y, int32_t incy);

/**
 * Swap vectors: x <-> y
 */
void blas_dswap(int32_t n, double* x, int32_t incx, double* y, int32_t incy);

#endif /* NUMJS_BLAS_H */
```

**File:** `src/wasm/blas.c` (new file)

```c
#include "blas.h"
#include <math.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/* ============ Level 1 BLAS Implementation ============ */

EXPORT double blas_ddot(int32_t n, const double* x, int32_t incx,
                         const double* y, int32_t incy)
{
    double result = 0.0;

    if (incx == 1 && incy == 1) {
        /* Optimized path: use pairwise summation for accuracy */
        /* Unroll loop for performance */
        int32_t m = n % 4;
        for (int32_t i = 0; i < m; i++) {
            result += x[i] * y[i];
        }
        for (int32_t i = m; i < n; i += 4) {
            result += x[i] * y[i] + x[i+1] * y[i+1] +
                      x[i+2] * y[i+2] + x[i+3] * y[i+3];
        }
    } else {
        /* Strided access */
        int32_t ix = 0, iy = 0;
        for (int32_t i = 0; i < n; i++) {
            result += x[ix] * y[iy];
            ix += incx;
            iy += incy;
        }
    }

    return result;
}

EXPORT double blas_dnrm2(int32_t n, const double* x, int32_t incx)
{
    /* Use scaled summation to avoid overflow/underflow */
    double scale = 0.0;
    double ssq = 1.0;

    int32_t ix = 0;
    for (int32_t i = 0; i < n; i++) {
        double absxi = fabs(x[ix]);
        if (absxi > 0.0) {
            if (scale < absxi) {
                double ratio = scale / absxi;
                ssq = 1.0 + ssq * ratio * ratio;
                scale = absxi;
            } else {
                double ratio = absxi / scale;
                ssq += ratio * ratio;
            }
        }
        ix += incx;
    }

    return scale * sqrt(ssq);
}

EXPORT void blas_daxpy(int32_t n, double alpha, const double* x, int32_t incx,
                        double* y, int32_t incy)
{
    if (alpha == 0.0) return;

    if (incx == 1 && incy == 1) {
        /* Unrolled loop */
        int32_t m = n % 4;
        for (int32_t i = 0; i < m; i++) {
            y[i] += alpha * x[i];
        }
        for (int32_t i = m; i < n; i += 4) {
            y[i]   += alpha * x[i];
            y[i+1] += alpha * x[i+1];
            y[i+2] += alpha * x[i+2];
            y[i+3] += alpha * x[i+3];
        }
    } else {
        int32_t ix = 0, iy = 0;
        for (int32_t i = 0; i < n; i++) {
            y[iy] += alpha * x[ix];
            ix += incx;
            iy += incy;
        }
    }
}

/* ... additional implementations ... */
```

---

#### 13.1.3 Level 3 BLAS (Matrix-Matrix Operations)

**File:** `src/wasm/blas.h` (additions)

```c
/* ============ Level 3 BLAS ============ */

/**
 * General matrix-matrix multiplication.
 * C = alpha * op(A) * op(B) + beta * C
 *
 * @param transA  'N' = A, 'T' = A^T, 'C' = A^H
 * @param transB  'N' = B, 'T' = B^T, 'C' = B^H
 * @param m       Rows of op(A) and C
 * @param n       Columns of op(B) and C
 * @param k       Columns of op(A), rows of op(B)
 * @param alpha   Scalar multiplier
 * @param A       Matrix A (lda x k for 'N', lda x m for 'T')
 * @param lda     Leading dimension of A
 * @param B       Matrix B (ldb x n for 'N', ldb x k for 'T')
 * @param ldb     Leading dimension of B
 * @param beta    Scalar multiplier for C
 * @param C       Matrix C (ldc x n)
 * @param ldc     Leading dimension of C
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

void blas_zgemm(char transA, char transB,
                int32_t m, int32_t n, int32_t k,
                const double* alpha,  /* complex */
                const double* A, int32_t lda,
                const double* B, int32_t ldb,
                const double* beta,   /* complex */
                double* C, int32_t ldc);

/**
 * Triangular matrix-matrix multiplication.
 * B = alpha * op(A) * B  or  B = alpha * B * op(A)
 *
 * @param side   'L' = left (A*B), 'R' = right (B*A)
 * @param uplo   'U' = upper triangular, 'L' = lower triangular
 * @param transA 'N' = A, 'T' = A^T, 'C' = A^H
 * @param diag   'N' = non-unit diagonal, 'U' = unit diagonal
 */
void blas_dtrmm(char side, char uplo, char transA, char diag,
                int32_t m, int32_t n,
                double alpha,
                const double* A, int32_t lda,
                double* B, int32_t ldb);

/**
 * Triangular solve with multiple right-hand sides.
 * Solve op(A) * X = alpha * B  or  X * op(A) = alpha * B
 * Solution overwrites B.
 */
void blas_dtrsm(char side, char uplo, char transA, char diag,
                int32_t m, int32_t n,
                double alpha,
                const double* A, int32_t lda,
                double* B, int32_t ldb);
```

**File:** `src/wasm/blas.c` (additions)

```c
EXPORT void blas_dgemm(char transA, char transB,
                        int32_t m, int32_t n, int32_t k,
                        double alpha,
                        const double* A, int32_t lda,
                        const double* B, int32_t ldb,
                        double beta,
                        double* C, int32_t ldc)
{
    /* Handle beta scaling of C */
    if (beta == 0.0) {
        for (int32_t j = 0; j < n; j++) {
            for (int32_t i = 0; i < m; i++) {
                C[i + j * ldc] = 0.0;
            }
        }
    } else if (beta != 1.0) {
        for (int32_t j = 0; j < n; j++) {
            for (int32_t i = 0; i < m; i++) {
                C[i + j * ldc] *= beta;
            }
        }
    }

    if (alpha == 0.0) return;

    /* Main multiplication loop */
    int noTransA = (transA == 'N' || transA == 'n');
    int noTransB = (transB == 'N' || transB == 'n');

    if (noTransA && noTransB) {
        /* C = alpha * A * B + beta * C */
        for (int32_t j = 0; j < n; j++) {
            for (int32_t l = 0; l < k; l++) {
                double temp = alpha * B[l + j * ldb];
                for (int32_t i = 0; i < m; i++) {
                    C[i + j * ldc] += temp * A[i + l * lda];
                }
            }
        }
    } else if (!noTransA && noTransB) {
        /* C = alpha * A^T * B + beta * C */
        for (int32_t j = 0; j < n; j++) {
            for (int32_t i = 0; i < m; i++) {
                double temp = 0.0;
                for (int32_t l = 0; l < k; l++) {
                    temp += A[l + i * lda] * B[l + j * ldb];
                }
                C[i + j * ldc] += alpha * temp;
            }
        }
    } else if (noTransA && !noTransB) {
        /* C = alpha * A * B^T + beta * C */
        for (int32_t j = 0; j < n; j++) {
            for (int32_t l = 0; l < k; l++) {
                double temp = alpha * B[j + l * ldb];
                for (int32_t i = 0; i < m; i++) {
                    C[i + j * ldc] += temp * A[i + l * lda];
                }
            }
        }
    } else {
        /* C = alpha * A^T * B^T + beta * C */
        for (int32_t j = 0; j < n; j++) {
            for (int32_t i = 0; i < m; i++) {
                double temp = 0.0;
                for (int32_t l = 0; l < k; l++) {
                    temp += A[l + i * lda] * B[j + l * ldb];
                }
                C[i + j * ldc] += alpha * temp;
            }
        }
    }
}
```

---

### 13.2 LAPACK Decompositions

#### 13.2.1 LU Factorization

**File:** `src/wasm/lapack.h` (new file)

```c
#ifndef NUMJS_LAPACK_H
#define NUMJS_LAPACK_H

#include "ndarray.h"

/* ============ LU Factorization ============ */

/**
 * Compute LU factorization with partial pivoting.
 * A = P * L * U
 *
 * On exit, A contains L (unit lower triangular) and U (upper triangular).
 * L is stored below diagonal, U is stored on and above diagonal.
 *
 * @param m     Number of rows
 * @param n     Number of columns
 * @param A     Input/output matrix (column-major, m x n)
 * @param lda   Leading dimension of A
 * @param ipiv  Output: pivot indices (size min(m,n))
 * @return      0 on success, >0 if U(i,i)=0 (singular), <0 for invalid arg
 */
int32_t lapack_dgetrf(int32_t m, int32_t n, double* A, int32_t lda,
                       int32_t* ipiv);

int32_t lapack_sgetrf(int32_t m, int32_t n, float* A, int32_t lda,
                       int32_t* ipiv);

int32_t lapack_zgetrf(int32_t m, int32_t n, double* A, int32_t lda,
                       int32_t* ipiv);  /* complex */

/**
 * Solve A*X = B using LU factorization from getrf.
 *
 * @param trans  'N' = A*X=B, 'T' = A^T*X=B, 'C' = A^H*X=B
 * @param n      Order of matrix A
 * @param nrhs   Number of right-hand sides (columns of B)
 * @param A      LU factorization from getrf
 * @param lda    Leading dimension of A
 * @param ipiv   Pivot indices from getrf
 * @param B      Input/output: RHS on entry, solution on exit
 * @param ldb    Leading dimension of B
 * @return       0 on success
 */
int32_t lapack_dgetrs(char trans, int32_t n, int32_t nrhs,
                       const double* A, int32_t lda,
                       const int32_t* ipiv,
                       double* B, int32_t ldb);

/**
 * Compute inverse of matrix using LU factorization from getrf.
 *
 * @param n      Order of matrix
 * @param A      Input: LU from getrf. Output: inverse
 * @param lda    Leading dimension
 * @param ipiv   Pivot indices from getrf
 * @param work   Workspace (size >= n)
 * @param lwork  Size of workspace
 * @return       0 on success, >0 if singular
 */
int32_t lapack_dgetri(int32_t n, double* A, int32_t lda,
                       const int32_t* ipiv, double* work, int32_t lwork);

/* ============ Cholesky Factorization ============ */

/**
 * Compute Cholesky factorization of symmetric positive-definite matrix.
 * A = L * L^T (lower) or A = U^T * U (upper)
 *
 * @param uplo   'U' = upper, 'L' = lower
 * @param n      Order of matrix
 * @param A      Input/output matrix
 * @param lda    Leading dimension
 * @return       0 on success, >0 if not positive-definite
 */
int32_t lapack_dpotrf(char uplo, int32_t n, double* A, int32_t lda);

int32_t lapack_spotrf(char uplo, int32_t n, float* A, int32_t lda);

int32_t lapack_zpotrf(char uplo, int32_t n, double* A, int32_t lda);  /* complex */

/**
 * Solve A*X = B using Cholesky factorization.
 */
int32_t lapack_dpotrs(char uplo, int32_t n, int32_t nrhs,
                       const double* A, int32_t lda,
                       double* B, int32_t ldb);

/* ============ QR Factorization ============ */

/**
 * Compute QR factorization using Householder reflectors.
 * A = Q * R
 *
 * @param m     Number of rows
 * @param n     Number of columns
 * @param A     Input/output: on exit contains R and Householder vectors
 * @param lda   Leading dimension
 * @param tau   Output: scalar factors of Householder reflectors (size min(m,n))
 * @param work  Workspace
 * @param lwork Workspace size (-1 for query)
 * @return      0 on success
 */
int32_t lapack_dgeqrf(int32_t m, int32_t n, double* A, int32_t lda,
                       double* tau, double* work, int32_t lwork);

/**
 * Generate orthogonal matrix Q from Householder reflectors.
 *
 * @param m     Number of rows of Q
 * @param n     Number of columns of Q
 * @param k     Number of reflectors
 * @param A     Input: Householder vectors from geqrf. Output: Q
 * @param lda   Leading dimension
 * @param tau   Scalar factors from geqrf
 * @param work  Workspace
 * @param lwork Workspace size
 */
int32_t lapack_dorgqr(int32_t m, int32_t n, int32_t k,
                       double* A, int32_t lda, const double* tau,
                       double* work, int32_t lwork);

int32_t lapack_zungqr(int32_t m, int32_t n, int32_t k,
                       double* A, int32_t lda, const double* tau,
                       double* work, int32_t lwork);  /* complex */

/* ============ Eigenvalue Decomposition ============ */

/**
 * Compute eigenvalues and optionally eigenvectors of general matrix.
 * A * v = lambda * v
 *
 * @param jobvl  'N' = no left eigenvectors, 'V' = compute left
 * @param jobvr  'N' = no right eigenvectors, 'V' = compute right
 * @param n      Order of matrix
 * @param A      Input/output matrix (destroyed)
 * @param lda    Leading dimension
 * @param wr     Output: real parts of eigenvalues (size n)
 * @param wi     Output: imaginary parts of eigenvalues (size n)
 * @param VL     Output: left eigenvectors (if jobvl='V')
 * @param ldvl   Leading dimension of VL
 * @param VR     Output: right eigenvectors (if jobvr='V')
 * @param ldvr   Leading dimension of VR
 * @param work   Workspace
 * @param lwork  Workspace size
 * @return       0 on success, >0 if QR algorithm failed
 */
int32_t lapack_dgeev(char jobvl, char jobvr, int32_t n,
                      double* A, int32_t lda,
                      double* wr, double* wi,
                      double* VL, int32_t ldvl,
                      double* VR, int32_t ldvr,
                      double* work, int32_t lwork);

int32_t lapack_zgeev(char jobvl, char jobvr, int32_t n,
                      double* A, int32_t lda,  /* complex */
                      double* w,               /* complex eigenvalues */
                      double* VL, int32_t ldvl,
                      double* VR, int32_t ldvr,
                      double* work, int32_t lwork,
                      double* rwork);

/**
 * Compute eigenvalues and eigenvectors of symmetric/Hermitian matrix.
 * Uses divide-and-conquer algorithm.
 *
 * @param jobz   'N' = eigenvalues only, 'V' = eigenvalues and eigenvectors
 * @param uplo   'U' = upper triangular, 'L' = lower triangular
 * @param n      Order of matrix
 * @param A      Input/output: on exit contains eigenvectors (if jobz='V')
 * @param lda    Leading dimension
 * @param w      Output: eigenvalues in ascending order
 * @param work   Workspace
 * @param lwork  Workspace size
 * @param iwork  Integer workspace
 * @param liwork Integer workspace size
 * @return       0 on success
 */
int32_t lapack_dsyevd(char jobz, char uplo, int32_t n,
                       double* A, int32_t lda, double* w,
                       double* work, int32_t lwork,
                       int32_t* iwork, int32_t liwork);

int32_t lapack_zheevd(char jobz, char uplo, int32_t n,
                       double* A, int32_t lda, double* w,  /* complex A, real w */
                       double* work, int32_t lwork,
                       double* rwork, int32_t lrwork,
                       int32_t* iwork, int32_t liwork);

/* ============ Singular Value Decomposition ============ */

/**
 * Compute SVD using divide-and-conquer algorithm.
 * A = U * S * V^H
 *
 * @param jobz   'A' = all columns of U and V^H
 *               'S' = first min(m,n) columns
 *               'O' = overwrite A with U or V^H
 *               'N' = no singular vectors
 * @param m      Number of rows
 * @param n      Number of columns
 * @param A      Input/output matrix
 * @param lda    Leading dimension of A
 * @param s      Output: singular values in descending order (size min(m,n))
 * @param U      Output: left singular vectors (if jobz != 'N')
 * @param ldu    Leading dimension of U
 * @param VT     Output: right singular vectors (transposed/conjugated)
 * @param ldvt   Leading dimension of VT
 * @param work   Workspace
 * @param lwork  Workspace size
 * @param iwork  Integer workspace (size >= 8*min(m,n))
 * @return       0 on success, >0 if algorithm did not converge
 */
int32_t lapack_dgesdd(char jobz, int32_t m, int32_t n,
                       double* A, int32_t lda, double* s,
                       double* U, int32_t ldu,
                       double* VT, int32_t ldvt,
                       double* work, int32_t lwork, int32_t* iwork);

int32_t lapack_zgesdd(char jobz, int32_t m, int32_t n,
                       double* A, int32_t lda, double* s,  /* complex A, real s */
                       double* U, int32_t ldu,
                       double* VT, int32_t ldvt,
                       double* work, int32_t lwork,
                       double* rwork, int32_t* iwork);

/* ============ Least Squares ============ */

/**
 * Solve overdetermined/underdetermined linear systems using SVD.
 * Minimize ||b - A*x||_2
 *
 * @param m      Number of rows
 * @param n      Number of columns
 * @param nrhs   Number of right-hand sides
 * @param A      Input matrix (destroyed)
 * @param lda    Leading dimension of A
 * @param B      Input/output: RHS/solution
 * @param ldb    Leading dimension of B
 * @param s      Output: singular values
 * @param rcond  Singular values below rcond*s[0] treated as zero
 * @param rank   Output: effective rank
 * @param work   Workspace
 * @param lwork  Workspace size
 * @param iwork  Integer workspace
 * @return       0 on success
 */
int32_t lapack_dgelsd(int32_t m, int32_t n, int32_t nrhs,
                       double* A, int32_t lda,
                       double* B, int32_t ldb,
                       double* s, double rcond, int32_t* rank,
                       double* work, int32_t lwork, int32_t* iwork);

#endif /* NUMJS_LAPACK_H */
```

---

### 13.3 High-Level Linear Algebra (TypeScript)

#### 13.3.1 Matrix Products

**File:** `src/ts/linalg.ts` (new file)

```typescript
import { NDArray } from './NDArray.js';
import { DType } from './types.js';
import { broadcastArrays, broadcastShapes } from './broadcast.js';

/**
 * NumJS linear algebra error.
 * Raised when matrix operations fail (singular matrix, non-convergence, etc.)
 */
export class LinAlgError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'LinAlgError';
  }
}

/* ============ Result Types ============ */

/**
 * Result of eigenvalue decomposition.
 */
export interface EigResult {
  eigenvalues: NDArray;   // Shape (..., M) - may be complex
  eigenvectors: NDArray;  // Shape (..., M, M) - right eigenvectors as columns
}

/**
 * Result of Hermitian eigenvalue decomposition.
 */
export interface EighResult {
  eigenvalues: NDArray;   // Shape (..., M) - always real, ascending
  eigenvectors: NDArray;  // Shape (..., M, M) - unitary/orthogonal
}

/**
 * Result of singular value decomposition.
 */
export interface SVDResult {
  U: NDArray;   // Shape (..., M, M) or (..., M, K) - left singular vectors
  S: NDArray;   // Shape (..., K) - singular values, descending
  Vh: NDArray;  // Shape (..., N, N) or (..., K, N) - right singular vectors
}

/**
 * Result of QR decomposition.
 */
export interface QRResult {
  Q: NDArray;  // Shape (..., M, K) or (..., M, M) - orthogonal matrix
  R: NDArray;  // Shape (..., K, N) or (..., M, N) - upper triangular
}

/**
 * Result of slogdet (sign and log-determinant).
 */
export interface SlogdetResult {
  sign: NDArray;      // Sign of determinant
  logabsdet: NDArray; // Log of absolute value of determinant
}

/**
 * Result of lstsq (least squares).
 */
export interface LstsqResult {
  x: NDArray;         // Least-squares solution
  residuals: NDArray; // Sum of squared residuals
  rank: number;       // Effective rank
  s: NDArray;         // Singular values
}

/* ============ Type Utilities ============ */

/**
 * Determine the computation dtype for linalg operations.
 * Integer types are promoted to float64.
 * Float types are preserved but must be float32 or float64.
 */
function _linalgType(arr: NDArray): DType {
  const dtype = arr.dtype;

  // Integer types -> float64
  if (dtype === DType.Bool || dtype === DType.Int8 || dtype === DType.Int16 ||
      dtype === DType.Int32 || dtype === DType.Int64 ||
      dtype === DType.UInt8 || dtype === DType.UInt16 ||
      dtype === DType.UInt32 || dtype === DType.UInt64) {
    return DType.Float64;
  }

  // Float16 -> Float32
  if (dtype === DType.Float16) {
    return DType.Float32;
  }

  // Float32, Float64, Complex64, Complex128 -> keep as is
  return dtype;
}

/**
 * Get the real dtype corresponding to a complex dtype.
 */
function _realType(dtype: DType): DType {
  if (dtype === DType.Complex64) return DType.Float32;
  if (dtype === DType.Complex128) return DType.Float64;
  return dtype;
}

/**
 * Check if dtype is complex.
 */
function _isComplex(dtype: DType): boolean {
  return dtype === DType.Complex64 || dtype === DType.Complex128;
}

/**
 * Assert array is at least 2D.
 */
function _assertStacked2d(arr: NDArray, name: string = 'array'): void {
  if (arr.ndim < 2) {
    throw new LinAlgError(
      `${name} must be at least 2-dimensional, got ${arr.ndim}-dimensional`
    );
  }
}

/**
 * Assert array is square in last two dimensions.
 */
function _assertSquare(arr: NDArray, name: string = 'array'): void {
  _assertStacked2d(arr, name);
  const m = arr.shape[arr.ndim - 2];
  const n = arr.shape[arr.ndim - 1];
  if (m !== n) {
    throw new LinAlgError(
      `${name} must be square, got shape (${arr.shape.join(', ')})`
    );
  }
}

/**
 * Assert array contains only finite values.
 */
function _assertFinite(arr: NDArray, name: string = 'array'): void {
  // Implementation would iterate and check for NaN/Inf
  // For now, a placeholder
}

/* ============ Matrix Products ============ */

/**
 * Dot product of two arrays.
 *
 * For 1-D arrays: inner product of vectors.
 * For 2-D arrays: matrix multiplication.
 * For N-D arrays: sum product over last axis of a and second-to-last of b.
 *
 * @param a - First array
 * @param b - Second array
 * @returns Dot product result
 *
 * @example
 * // Vector dot product
 * dot([1, 2, 3], [4, 5, 6])  // 32
 *
 * // Matrix multiplication
 * dot([[1, 2], [3, 4]], [[5, 6], [7, 8]])
 */
export function dot(a: NDArray, b: NDArray): NDArray {
  if (a.ndim === 0 || b.ndim === 0) {
    throw new LinAlgError('dot does not support 0-d arrays');
  }

  if (a.ndim === 1 && b.ndim === 1) {
    // Vector dot product
    if (a.shape[0] !== b.shape[0]) {
      throw new LinAlgError(
        `shapes (${a.shape[0]},) and (${b.shape[0]},) not aligned`
      );
    }
    return _vectorDot(a, b);
  }

  if (a.ndim === 2 && b.ndim === 2) {
    // Matrix multiplication
    return matmul(a, b);
  }

  // General case: sum product over last axis of a and second-to-last of b
  return _generalDot(a, b);
}

/**
 * Compute the dot product of two vectors (flattened inputs).
 */
export function vdot(a: NDArray, b: NDArray): number {
  const aFlat = a.ravel();
  const bFlat = b.ravel();

  if (aFlat.size !== bFlat.size) {
    throw new LinAlgError('vdot requires arrays of equal size');
  }

  // For complex arrays, conjugate first argument
  // For real arrays, just dot product
  // Returns scalar
  return _vectorDotScalar(aFlat, bFlat);
}

/**
 * Inner product of two arrays.
 * For 1-D arrays: dot product.
 * For N-D arrays: sum product over last axes.
 */
export function inner(a: NDArray, b: NDArray): NDArray {
  if (a.ndim === 0 || b.ndim === 0) {
    // Scalar case: just multiply
    return _scalarMultiply(a, b);
  }

  if (a.shape[a.ndim - 1] !== b.shape[b.ndim - 1]) {
    throw new LinAlgError(
      `shapes ${a.shape} and ${b.shape} not aligned: ` +
      `${a.shape[a.ndim - 1]} (dim ${a.ndim - 1}) != ${b.shape[b.ndim - 1]} (dim ${b.ndim - 1})`
    );
  }

  return _innerProduct(a, b);
}

/**
 * Compute the outer product of two vectors.
 */
export function outer(a: NDArray, b: NDArray): NDArray {
  const aFlat = a.ravel();
  const bFlat = b.ravel();

  // Result shape: (a.size, b.size)
  return _outerProduct(aFlat, bFlat);
}

/**
 * Matrix product of two arrays.
 *
 * Follows NumPy matmul semantics:
 * - 2D arrays: regular matrix multiplication
 * - N-D arrays: broadcast over batch dimensions, matmul on last two
 *
 * @param a - First array (..., M, K)
 * @param b - Second array (..., K, N)
 * @returns Product array (..., M, N)
 */
export function matmul(a: NDArray, b: NDArray): NDArray {
  if (a.ndim < 1 || b.ndim < 1) {
    throw new LinAlgError('matmul requires at least 1-dimensional arrays');
  }

  // Handle 1-D arrays by temporarily adding dimensions
  let a2 = a;
  let b2 = b;
  let removeFirst = false;
  let removeLast = false;

  if (a.ndim === 1) {
    a2 = a.reshape([1, a.shape[0]]);
    removeFirst = true;
  }
  if (b.ndim === 1) {
    b2 = b.reshape([b.shape[0], 1]);
    removeLast = true;
  }

  // Check inner dimensions match
  const k1 = a2.shape[a2.ndim - 1];
  const k2 = b2.shape[b2.ndim - 2];
  if (k1 !== k2) {
    throw new LinAlgError(
      `matmul: Input operand 1 has a mismatch in its core dimension 0, ` +
      `with gufunc signature (n?,k),(k,m?)->(n?,m?) (size ${k1} is different from ${k2})`
    );
  }

  // Call WASM implementation
  const result = _matmulImpl(a2, b2);

  // Remove temporarily added dimensions
  if (removeFirst && removeLast) {
    return result.reshape([]);
  } else if (removeFirst) {
    return result.reshape(result.shape.slice(0, -2).concat(result.shape[result.ndim - 1]));
  } else if (removeLast) {
    return result.reshape(result.shape.slice(0, -1));
  }

  return result;
}

/**
 * Compute tensor dot product along specified axes.
 *
 * @param a - First tensor
 * @param b - Second tensor
 * @param axes - Axes to sum over (default: 2 = last 2 of a, first 2 of b)
 */
export function tensordot(
  a: NDArray,
  b: NDArray,
  axes: number | [number[], number[]] = 2
): NDArray {
  let axesA: number[];
  let axesB: number[];

  if (typeof axes === 'number') {
    // axes=N means last N axes of a, first N axes of b
    axesA = [];
    axesB = [];
    for (let i = 0; i < axes; i++) {
      axesA.push(a.ndim - axes + i);
      axesB.push(i);
    }
  } else {
    [axesA, axesB] = axes;
  }

  return _tensordotImpl(a, b, axesA, axesB);
}

/**
 * Evaluates the Einstein summation convention on the operands.
 *
 * @param subscripts - Subscripts string (e.g., 'ij,jk->ik')
 * @param operands - Input arrays
 */
export function einsum(subscripts: string, ...operands: NDArray[]): NDArray {
  // Parse subscripts and determine contraction
  const { inputSubs, outputSubs, dimensions } = _parseEinsum(subscripts, operands);

  return _einsumImpl(inputSubs, outputSubs, dimensions, operands);
}

/**
 * Kronecker product of two arrays.
 */
export function kron(a: NDArray, b: NDArray): NDArray {
  // Result shape: element-wise multiply of shapes
  return _kronProduct(a, b);
}

/**
 * Compute the dot product of two or more arrays in a single function call,
 * while automatically selecting the fastest evaluation order.
 */
export function multi_dot(arrays: NDArray[]): NDArray {
  if (arrays.length < 2) {
    throw new LinAlgError('multi_dot requires at least 2 arrays');
  }
  if (arrays.length === 2) {
    return dot(arrays[0], arrays[1]);
  }

  // Use dynamic programming to find optimal parenthesization
  const order = _multiDotOptimalOrder(arrays);
  return _multiDotExecute(arrays, order);
}

/* ============ Decompositions ============ */

/**
 * Cholesky decomposition.
 *
 * Returns the lower or upper Cholesky factor of a Hermitian
 * positive-definite matrix.
 *
 * @param a - Hermitian positive-definite matrix (..., M, M)
 * @param upper - If true, return upper triangular factor U such that a = U^H @ U
 *                If false (default), return lower triangular L such that a = L @ L^H
 */
export function cholesky(a: NDArray, upper: boolean = false): NDArray {
  _assertSquare(a, 'a');

  const dtype = _linalgType(a);
  const workArr = a.astype(dtype).copy();

  const info = upper
    ? _choleskyUpper(workArr)
    : _choleskyLower(workArr);

  if (info > 0) {
    throw new LinAlgError(
      `Matrix is not positive definite - Cholesky decomposition cannot be computed`
    );
  }

  // Zero out the other triangle
  return upper ? _zeroLowerTriangle(workArr) : _zeroUpperTriangle(workArr);
}

/**
 * QR decomposition.
 *
 * Factor the matrix a as qr, where q is orthonormal and r is upper-triangular.
 *
 * @param a - Matrix to factor (..., M, N)
 * @param mode - 'reduced' (default), 'complete', 'r', or 'raw'
 */
export function qr(a: NDArray, mode: 'reduced' | 'complete' | 'r' | 'raw' = 'reduced'): QRResult | NDArray {
  _assertStacked2d(a, 'a');

  const m = a.shape[a.ndim - 2];
  const n = a.shape[a.ndim - 1];
  const k = Math.min(m, n);

  const dtype = _linalgType(a);
  const workArr = a.astype(dtype).copy();

  // Compute QR factorization
  const tau = _qrFactorize(workArr);

  if (mode === 'r') {
    // Return only R
    return _extractR(workArr, k);
  }

  if (mode === 'raw') {
    // Return (h, tau) - Householder reflectors
    return { Q: workArr, R: tau } as QRResult;
  }

  // Construct Q and R
  const rMatrix = _extractR(workArr, mode === 'complete' ? m : k);
  const qMatrix = _constructQ(workArr, tau, mode === 'complete' ? m : k);

  return { Q: qMatrix, R: rMatrix };
}

/**
 * Singular Value Decomposition.
 *
 * @param a - Matrix to decompose (..., M, N)
 * @param full_matrices - If true, U and Vh have shapes (..., M, M) and (..., N, N)
 *                        If false, shapes are (..., M, K) and (..., K, N) where K = min(M, N)
 * @param compute_uv - If true (default), compute U and Vh in addition to S
 * @param hermitian - If true, a is assumed Hermitian (eigenvalue-based algorithm)
 */
export function svd(
  a: NDArray,
  full_matrices: boolean = true,
  compute_uv: boolean = true,
  hermitian: boolean = false
): SVDResult | NDArray {
  _assertStacked2d(a, 'a');

  if (hermitian) {
    // Use eigenvalue decomposition for Hermitian matrices
    return _svdHermitian(a, full_matrices, compute_uv);
  }

  const dtype = _linalgType(a);
  const workArr = a.astype(dtype).copy();

  if (!compute_uv) {
    // Singular values only
    return _svdValuesOnly(workArr);
  }

  const { U, S, Vh } = full_matrices
    ? _svdFull(workArr)
    : _svdReduced(workArr);

  return { U, S, Vh };
}

/**
 * Return singular values only.
 */
export function svdvals(a: NDArray): NDArray {
  return svd(a, false, false) as NDArray;
}

/* ============ Eigenvalues ============ */

/**
 * Compute eigenvalues and right eigenvectors of a general matrix.
 *
 * @param a - Square matrix (..., M, M)
 * @returns EigResult with eigenvalues and eigenvectors
 */
export function eig(a: NDArray): EigResult {
  _assertSquare(a, 'a');

  const dtype = _linalgType(a);
  const workArr = a.astype(dtype).copy();

  const { eigenvalues, eigenvectors, info } = _eigGeneral(workArr);

  if (info > 0) {
    throw new LinAlgError(
      `Eigenvalue computation did not converge`
    );
  }

  return { eigenvalues, eigenvectors };
}

/**
 * Compute eigenvalues and eigenvectors of a Hermitian or real symmetric matrix.
 *
 * @param a - Hermitian/symmetric matrix (..., M, M)
 * @param UPLO - 'L' (default) uses lower triangle, 'U' uses upper triangle
 */
export function eigh(a: NDArray, UPLO: 'L' | 'U' = 'L'): EighResult {
  _assertSquare(a, 'a');

  const dtype = _linalgType(a);
  const workArr = a.astype(dtype).copy();

  const { eigenvalues, eigenvectors, info } = _eigHermitian(workArr, UPLO);

  if (info > 0) {
    throw new LinAlgError(
      `Eigenvalue computation did not converge`
    );
  }

  return { eigenvalues, eigenvectors };
}

/**
 * Compute eigenvalues of a general matrix.
 */
export function eigvals(a: NDArray): NDArray {
  const result = eig(a);
  return result.eigenvalues;
}

/**
 * Compute eigenvalues of a Hermitian or real symmetric matrix.
 */
export function eigvalsh(a: NDArray, UPLO: 'L' | 'U' = 'L'): NDArray {
  const result = eigh(a, UPLO);
  return result.eigenvalues;
}

/* ============ Norms & Numbers ============ */

/**
 * Matrix or vector norm.
 *
 * @param x - Input array
 * @param ord - Order of the norm (see below)
 * @param axis - Axis or axes along which to compute norm
 * @param keepdims - If true, axes are left with size 1
 *
 * Vector norms (ord):
 * - None or 'fro': Frobenius (L2)
 * - inf: max(|x|)
 * - -inf: min(|x|)
 * - 0: count of non-zero elements
 * - 1: sum(|x|)
 * - 2: sqrt(sum(|x|^2))
 * - other p: sum(|x|^p)^(1/p)
 *
 * Matrix norms (ord, axis must be 2-tuple):
 * - None or 'fro': Frobenius norm
 * - 'nuc': nuclear norm (sum of singular values)
 * - inf: max(sum(|row|))
 * - -inf: min(sum(|row|))
 * - 1: max(sum(|col|))
 * - -1: min(sum(|col|))
 * - 2: largest singular value
 * - -2: smallest singular value
 */
export function norm(
  x: NDArray,
  ord: number | 'fro' | 'nuc' | null = null,
  axis: number | [number, number] | null = null,
  keepdims: boolean = false
): NDArray | number {
  if (axis === null) {
    if (x.ndim === 1 || (ord === null && x.ndim === 2)) {
      // Vector norm or Frobenius
      return _vectorNorm(x.ravel(), ord ?? 2, keepdims);
    }
    // Treat as matrix norm over last two axes
    axis = [x.ndim - 2, x.ndim - 1] as [number, number];
  }

  if (typeof axis === 'number') {
    // Vector norm along single axis
    return _vectorNormAxis(x, ord ?? 2, axis, keepdims);
  }

  // Matrix norm over two axes
  return _matrixNorm(x, ord ?? 'fro', axis, keepdims);
}

/**
 * Compute matrix norm.
 */
export function matrix_norm(
  x: NDArray,
  ord: number | 'fro' | 'nuc' = 'fro',
  keepdims: boolean = false
): NDArray | number {
  _assertStacked2d(x, 'x');
  return norm(x, ord, [x.ndim - 2, x.ndim - 1], keepdims);
}

/**
 * Compute vector norm.
 */
export function vector_norm(
  x: NDArray,
  ord: number = 2,
  axis: number | null = null,
  keepdims: boolean = false
): NDArray | number {
  if (axis === null) {
    return _vectorNorm(x.ravel(), ord, keepdims);
  }
  return _vectorNormAxis(x, ord, axis, keepdims);
}

/**
 * Compute the condition number of a matrix.
 *
 * @param x - Matrix
 * @param p - Norm order (default: 2, using SVD)
 */
export function cond(x: NDArray, p: number | 'fro' | null = null): number {
  _assertStacked2d(x, 'x');

  if (p === null || p === 2 || p === -2) {
    // Use SVD
    const s = svdvals(x);
    const sArr = s.toArray();
    if (p === -2) {
      return sArr[sArr.length - 1] / sArr[0];
    }
    return sArr[0] / sArr[sArr.length - 1];
  }

  // Use norm
  const normX = norm(x, p) as number;
  const normXinv = norm(inv(x), p) as number;
  return normX * normXinv;
}

/**
 * Compute the determinant of an array.
 */
export function det(a: NDArray): NDArray | number {
  _assertSquare(a, 'a');

  const dtype = _linalgType(a);
  const workArr = a.astype(dtype).copy();

  return _determinant(workArr);
}

/**
 * Compute sign and (natural) logarithm of the determinant.
 */
export function slogdet(a: NDArray): SlogdetResult {
  _assertSquare(a, 'a');

  const dtype = _linalgType(a);
  const workArr = a.astype(dtype).copy();

  return _signLogDet(workArr);
}

/**
 * Return matrix rank using SVD method.
 *
 * @param A - Matrix
 * @param tol - Threshold for considering singular values as zero
 * @param hermitian - If true, use eigenvalue decomposition
 */
export function matrix_rank(
  A: NDArray,
  tol: number | null = null,
  hermitian: boolean = false
): number {
  _assertStacked2d(A, 'A');

  let s: NDArray;
  if (hermitian) {
    s = eigvalsh(A);
    // Take absolute value for counting
  } else {
    s = svdvals(A);
  }

  if (tol === null) {
    // Default tolerance: max(M, N) * eps * max(singular values)
    const m = A.shape[A.ndim - 2];
    const n = A.shape[A.ndim - 1];
    const eps = A.dtype === DType.Float32 ? 1.19e-7 : 2.22e-16;
    const sMax = Math.max(...s.toArray().map(Math.abs));
    tol = Math.max(m, n) * eps * sMax;
  }

  // Count singular values above tolerance
  return s.toArray().filter(v => Math.abs(v) > tol!).length;
}

/**
 * Return the sum along diagonals of the array.
 */
export function trace(
  a: NDArray,
  offset: number = 0,
  axis1: number = 0,
  axis2: number = 1
): NDArray | number {
  // Extract diagonal and sum
  const diag = a.diagonal(offset, axis1, axis2);
  return diag.sum();
}

/* ============ Solving & Inverting ============ */

/**
 * Solve a linear matrix equation: a @ x = b
 *
 * @param a - Coefficient matrix (..., M, M)
 * @param b - Ordinate values (..., M) or (..., M, K)
 */
export function solve(a: NDArray, b: NDArray): NDArray {
  _assertSquare(a, 'a');

  const m = a.shape[a.ndim - 1];

  // b can be 1D (single RHS) or 2D+ (multiple RHS)
  if (b.ndim >= 1 && b.shape[b.ndim - (b.ndim === a.ndim - 1 ? 1 : 2)] !== m) {
    throw new LinAlgError(
      `solve: last dimension of b (${b.shape[b.ndim - 1]}) must match ` +
      `matrix dimension (${m})`
    );
  }

  const dtype = _linalgType(a);
  const aWork = a.astype(dtype).copy();
  const bWork = b.astype(dtype).copy();

  const info = _solveLU(aWork, bWork);

  if (info > 0) {
    throw new LinAlgError('Singular matrix');
  }

  return bWork;
}

/**
 * Solve the tensor equation a x = b for x.
 */
export function tensorsolve(
  a: NDArray,
  b: NDArray,
  axes: number[] | null = null
): NDArray {
  // Reshape to matrix equation and solve
  return _tensorSolveImpl(a, b, axes);
}

/**
 * Return the least-squares solution to a linear matrix equation.
 *
 * @param a - Coefficient matrix (M, N)
 * @param b - Ordinate values (M,) or (M, K)
 * @param rcond - Cutoff for small singular values (default: N * machine_epsilon)
 */
export function lstsq(
  a: NDArray,
  b: NDArray,
  rcond: number | null = null
): LstsqResult {
  _assertStacked2d(a, 'a');

  const m = a.shape[a.ndim - 2];
  const n = a.shape[a.ndim - 1];

  if (rcond === null) {
    const eps = a.dtype === DType.Float32 ? 1.19e-7 : 2.22e-16;
    rcond = Math.max(m, n) * eps;
  }

  const dtype = _linalgType(a);
  const aWork = a.astype(dtype).copy();
  const bWork = b.astype(dtype).copy();

  const { x, residuals, rank, s, info } = _lstsqSVD(aWork, bWork, rcond);

  if (info > 0) {
    throw new LinAlgError('SVD computation did not converge in lstsq');
  }

  return { x, residuals, rank, s };
}

/**
 * Compute the (multiplicative) inverse of a matrix.
 */
export function inv(a: NDArray): NDArray {
  _assertSquare(a, 'a');

  const dtype = _linalgType(a);
  const workArr = a.astype(dtype).copy();

  const info = _invertMatrix(workArr);

  if (info > 0) {
    throw new LinAlgError('Singular matrix');
  }

  return workArr;
}

/**
 * Compute the (Moore-Penrose) pseudo-inverse of a matrix.
 *
 * @param a - Matrix to pseudo-invert
 * @param rcond - Cutoff for small singular values
 * @param hermitian - If true, assume a is Hermitian
 */
export function pinv(
  a: NDArray,
  rcond: number = 1e-15,
  hermitian: boolean = false
): NDArray {
  _assertStacked2d(a, 'a');

  if (hermitian) {
    return _pinvHermitian(a, rcond);
  }

  // SVD-based pseudo-inverse
  const { U, S, Vh } = svd(a, false, true) as SVDResult;

  // Compute reciprocal of singular values above threshold
  const sMax = Math.max(...S.toArray());
  const cutoff = rcond * sMax;

  // S_inv with zero for small singular values
  const sInv = S.toArray().map(s => s > cutoff ? 1 / s : 0);

  // pinv = Vh^H @ diag(s_inv) @ U^H
  return _computePinv(U, sInv, Vh);
}

/**
 * Compute the 'inverse' of an N-dimensional array.
 */
export function tensorinv(a: NDArray, ind: number = 2): NDArray {
  return _tensorInvImpl(a, ind);
}

/* ============ Matrix Operations ============ */

/**
 * Raise a square matrix to the (integer) power n.
 *
 * @param a - Square matrix
 * @param n - Exponent (can be negative for inverse)
 */
export function matrix_power(a: NDArray, n: number): NDArray {
  _assertSquare(a, 'a');

  if (!Number.isInteger(n)) {
    throw new LinAlgError('exponent must be an integer');
  }

  const m = a.shape[a.ndim - 1];

  if (n === 0) {
    // Return identity matrix
    return eye(m, m, 0, a.dtype);
  }

  if (n < 0) {
    a = inv(a);
    n = -n;
  }

  // Binary exponentiation
  let result = eye(m, m, 0, a.dtype);
  let base = a.copy();

  while (n > 0) {
    if (n & 1) {
      result = matmul(result, base);
    }
    base = matmul(base, base);
    n >>= 1;
  }

  return result;
}

/**
 * Permute the last two axes of an array.
 * Equivalent to x.swapaxes(-2, -1) for 2-D arrays.
 */
export function matrix_transpose(x: NDArray): NDArray {
  if (x.ndim < 2) {
    throw new LinAlgError('matrix_transpose requires at least 2-dimensional array');
  }
  return x.swapaxes(x.ndim - 2, x.ndim - 1);
}

/**
 * Return the cross product of two (arrays of) vectors.
 */
export function cross(
  a: NDArray,
  b: NDArray,
  axisa: number = -1,
  axisb: number = -1,
  axisc: number = -1,
  axis: number | null = null
): NDArray {
  if (axis !== null) {
    axisa = axisb = axisc = axis;
  }

  return _crossProduct(a, b, axisa, axisb, axisc);
}

/* ============ Internal Implementation Stubs ============ */

// These functions would call into WASM for actual computation

function _vectorDot(a: NDArray, b: NDArray): NDArray {
  // Call WASM blas_ddot
  throw new Error('Not implemented');
}

function _vectorDotScalar(a: NDArray, b: NDArray): number {
  // Call WASM blas_ddot for scalar result
  throw new Error('Not implemented');
}

function _scalarMultiply(a: NDArray, b: NDArray): NDArray {
  throw new Error('Not implemented');
}

function _innerProduct(a: NDArray, b: NDArray): NDArray {
  throw new Error('Not implemented');
}

function _outerProduct(a: NDArray, b: NDArray): NDArray {
  throw new Error('Not implemented');
}

function _generalDot(a: NDArray, b: NDArray): NDArray {
  throw new Error('Not implemented');
}

function _matmulImpl(a: NDArray, b: NDArray): NDArray {
  // Call WASM blas_dgemm
  throw new Error('Not implemented');
}

function _tensordotImpl(a: NDArray, b: NDArray, axesA: number[], axesB: number[]): NDArray {
  throw new Error('Not implemented');
}

function _parseEinsum(subscripts: string, operands: NDArray[]): any {
  throw new Error('Not implemented');
}

function _einsumImpl(inputSubs: any, outputSubs: any, dimensions: any, operands: NDArray[]): NDArray {
  throw new Error('Not implemented');
}

function _kronProduct(a: NDArray, b: NDArray): NDArray {
  throw new Error('Not implemented');
}

function _multiDotOptimalOrder(arrays: NDArray[]): number[][] {
  throw new Error('Not implemented');
}

function _multiDotExecute(arrays: NDArray[], order: number[][]): NDArray {
  throw new Error('Not implemented');
}

function _choleskyLower(a: NDArray): number {
  // Call WASM lapack_dpotrf with 'L'
  throw new Error('Not implemented');
}

function _choleskyUpper(a: NDArray): number {
  // Call WASM lapack_dpotrf with 'U'
  throw new Error('Not implemented');
}

function _zeroLowerTriangle(a: NDArray): NDArray {
  throw new Error('Not implemented');
}

function _zeroUpperTriangle(a: NDArray): NDArray {
  throw new Error('Not implemented');
}

function _qrFactorize(a: NDArray): NDArray {
  // Call WASM lapack_dgeqrf
  throw new Error('Not implemented');
}

function _extractR(a: NDArray, rows: number): NDArray {
  throw new Error('Not implemented');
}

function _constructQ(a: NDArray, tau: NDArray, cols: number): NDArray {
  // Call WASM lapack_dorgqr
  throw new Error('Not implemented');
}

function _svdHermitian(a: NDArray, full_matrices: boolean, compute_uv: boolean): SVDResult | NDArray {
  throw new Error('Not implemented');
}

function _svdValuesOnly(a: NDArray): NDArray {
  // Call WASM lapack_dgesdd with jobz='N'
  throw new Error('Not implemented');
}

function _svdFull(a: NDArray): SVDResult {
  // Call WASM lapack_dgesdd with jobz='A'
  throw new Error('Not implemented');
}

function _svdReduced(a: NDArray): SVDResult {
  // Call WASM lapack_dgesdd with jobz='S'
  throw new Error('Not implemented');
}

function _eigGeneral(a: NDArray): { eigenvalues: NDArray; eigenvectors: NDArray; info: number } {
  // Call WASM lapack_dgeev
  throw new Error('Not implemented');
}

function _eigHermitian(a: NDArray, UPLO: 'L' | 'U'): { eigenvalues: NDArray; eigenvectors: NDArray; info: number } {
  // Call WASM lapack_dsyevd or lapack_zheevd
  throw new Error('Not implemented');
}

function _vectorNorm(x: NDArray, ord: number, keepdims: boolean): NDArray | number {
  throw new Error('Not implemented');
}

function _vectorNormAxis(x: NDArray, ord: number, axis: number, keepdims: boolean): NDArray | number {
  throw new Error('Not implemented');
}

function _matrixNorm(x: NDArray, ord: number | 'fro' | 'nuc', axis: [number, number], keepdims: boolean): NDArray | number {
  throw new Error('Not implemented');
}

function _determinant(a: NDArray): NDArray | number {
  // Use LU factorization
  throw new Error('Not implemented');
}

function _signLogDet(a: NDArray): SlogdetResult {
  // Use LU factorization
  throw new Error('Not implemented');
}

function _solveLU(a: NDArray, b: NDArray): number {
  // Call WASM lapack_dgetrf + lapack_dgetrs
  throw new Error('Not implemented');
}

function _tensorSolveImpl(a: NDArray, b: NDArray, axes: number[] | null): NDArray {
  throw new Error('Not implemented');
}

function _lstsqSVD(a: NDArray, b: NDArray, rcond: number): { x: NDArray; residuals: NDArray; rank: number; s: NDArray; info: number } {
  // Call WASM lapack_dgelsd
  throw new Error('Not implemented');
}

function _invertMatrix(a: NDArray): number {
  // Call WASM lapack_dgetrf + lapack_dgetri
  throw new Error('Not implemented');
}

function _pinvHermitian(a: NDArray, rcond: number): NDArray {
  throw new Error('Not implemented');
}

function _computePinv(U: NDArray, sInv: number[], Vh: NDArray): NDArray {
  throw new Error('Not implemented');
}

function _tensorInvImpl(a: NDArray, ind: number): NDArray {
  throw new Error('Not implemented');
}

function _crossProduct(a: NDArray, b: NDArray, axisa: number, axisb: number, axisc: number): NDArray {
  throw new Error('Not implemented');
}

function eye(N: number, M: number, k: number, dtype: DType): NDArray {
  // Reference existing NDArray.eye implementation
  throw new Error('Not implemented');
}

/* ============ Exports ============ */

export const linalg = {
  // Error class
  LinAlgError,

  // Matrix products
  dot,
  vdot,
  inner,
  outer,
  matmul,
  tensordot,
  einsum,
  kron,
  multi_dot,

  // Decompositions
  cholesky,
  qr,
  svd,
  svdvals,

  // Eigenvalues
  eig,
  eigh,
  eigvals,
  eigvalsh,

  // Norms & Numbers
  norm,
  matrix_norm,
  vector_norm,
  cond,
  det,
  slogdet,
  matrix_rank,
  trace,

  // Solving & Inverting
  solve,
  tensorsolve,
  lstsq,
  inv,
  pinv,
  tensorinv,

  // Matrix Operations
  matrix_power,
  matrix_transpose,
  cross,
};
```

---

## File Changes Summary

### New Files to Create

```
src/wasm/
├── blas.h             # BLAS declarations (Level 1, 2, 3)
├── blas.c             # BLAS implementation
├── lapack.h           # LAPACK declarations
├── lapack.c           # LAPACK implementation (decompositions)
└── linalg.c           # High-level linalg WASM functions

src/ts/
└── linalg.ts          # Full linear algebra module
```

### Files to Modify

```
src/ts/types.ts
├── Add LinAlgError type
├── Add result interfaces (EigResult, SVDResult, etc.)
└── Add WASM function declarations for linalg

src/ts/index.ts
├── Export linalg module
├── Export individual functions (dot, matmul, solve, etc.)
└── Export LinAlgError

scripts/build-wasm.sh
├── Add blas.c to compilation
├── Add lapack.c to compilation
├── Add linalg.c to compilation
└── Add EXPORTED_FUNCTIONS for all linalg operations
```

---

## Build Script Updates

Add to `scripts/build-wasm.sh` EXPORTED_FUNCTIONS:

```bash
# BLAS Level 1
"_blas_ddot",
"_blas_sdot",
"_blas_zdotc",
"_blas_dnrm2",
"_blas_snrm2",
"_blas_dasum",
"_blas_idamax",
"_blas_dscal",
"_blas_daxpy",
"_blas_dcopy",
"_blas_dswap",

# BLAS Level 2
"_blas_dgemv",
"_blas_dtrmv",
"_blas_dtrsv",

# BLAS Level 3
"_blas_dgemm",
"_blas_sgemm",
"_blas_zgemm",
"_blas_dtrmm",
"_blas_dtrsm",

# LAPACK Decompositions
"_lapack_dgetrf",
"_lapack_dgetrs",
"_lapack_dgetri",
"_lapack_dpotrf",
"_lapack_dpotrs",
"_lapack_dgeqrf",
"_lapack_dorgqr",
"_lapack_dgeev",
"_lapack_dsyevd",
"_lapack_dgesdd",
"_lapack_dgelsd",

# High-level operations
"_linalg_matmul",
"_linalg_solve",
"_linalg_inv",
"_linalg_det",
"_linalg_slogdet",
"_linalg_norm",
"_linalg_eig",
"_linalg_svd",
"_linalg_qr",
"_linalg_cholesky"
```

Add new source files:

```bash
emcc \
    "$SRC_DIR/ndarray.c" \
    "$SRC_DIR/dtype.c" \
    "$SRC_DIR/pairwise_sum.c" \
    "$SRC_DIR/broadcast.c" \
    "$SRC_DIR/indexing.c" \
    "$SRC_DIR/blas.c" \
    "$SRC_DIR/lapack.c" \
    "$SRC_DIR/linalg.c" \
    -o "$OUT_DIR/numjs.cjs" \
    ...
```

---

## WasmModule Interface Updates

Add to `src/ts/types.ts`:

```typescript
// BLAS Level 1
_blas_ddot(n: number, xPtr: number, incx: number, yPtr: number, incy: number): number;
_blas_dnrm2(n: number, xPtr: number, incx: number): number;
_blas_daxpy(n: number, alpha: number, xPtr: number, incx: number, yPtr: number, incy: number): void;

// BLAS Level 3
_blas_dgemm(transA: number, transB: number, m: number, n: number, k: number,
            alpha: number, aPtr: number, lda: number, bPtr: number, ldb: number,
            beta: number, cPtr: number, ldc: number): void;

// LAPACK
_lapack_dgetrf(m: number, n: number, aPtr: number, lda: number, ipivPtr: number): number;
_lapack_dgetrs(trans: number, n: number, nrhs: number, aPtr: number, lda: number,
               ipivPtr: number, bPtr: number, ldb: number): number;
_lapack_dgetri(n: number, aPtr: number, lda: number, ipivPtr: number,
               workPtr: number, lwork: number): number;
_lapack_dpotrf(uplo: number, n: number, aPtr: number, lda: number): number;
_lapack_dgeqrf(m: number, n: number, aPtr: number, lda: number, tauPtr: number,
               workPtr: number, lwork: number): number;
_lapack_dorgqr(m: number, n: number, k: number, aPtr: number, lda: number,
               tauPtr: number, workPtr: number, lwork: number): number;
_lapack_dgeev(jobvl: number, jobvr: number, n: number, aPtr: number, lda: number,
              wrPtr: number, wiPtr: number, vlPtr: number, ldvl: number,
              vrPtr: number, ldvr: number, workPtr: number, lwork: number): number;
_lapack_dsyevd(jobz: number, uplo: number, n: number, aPtr: number, lda: number,
               wPtr: number, workPtr: number, lwork: number,
               iworkPtr: number, liwork: number): number;
_lapack_dgesdd(jobz: number, m: number, n: number, aPtr: number, lda: number,
               sPtr: number, uPtr: number, ldu: number, vtPtr: number, ldvt: number,
               workPtr: number, lwork: number, iworkPtr: number): number;
_lapack_dgelsd(m: number, n: number, nrhs: number, aPtr: number, lda: number,
               bPtr: number, ldb: number, sPtr: number, rcond: number,
               rankPtr: number, workPtr: number, lwork: number, iworkPtr: number): number;

// High-level
_linalg_matmul(aPtr: number, bPtr: number): number;
_linalg_solve(aPtr: number, bPtr: number): number;
_linalg_inv(aPtr: number): number;
```

---

## Implementation Order

```
Phase 13.1: BLAS Foundation (Weeks 1-2)
├── Week 1: Level 1 BLAS
│   ├── Day 1: ddot, sdot (vector dot product)
│   ├── Day 2: dnrm2, snrm2 (L2 norm with scaling)
│   ├── Day 3: dasum, idamax (L1 norm, max element)
│   ├── Day 4: dscal, daxpy (scale, axpy)
│   └── Day 5: dcopy, dswap + tests
│
└── Week 2: Level 2 & 3 BLAS
    ├── Day 1: dgemv (matrix-vector multiply)
    ├── Day 2: dtrmv, dtrsv (triangular operations)
    ├── Day 3: dgemm (matrix multiply - critical!)
    ├── Day 4: dtrmm, dtrsm (triangular matrix ops)
    └── Day 5: Complex variants (zgemm, etc.) + tests

Phase 13.2: LAPACK Decompositions (Weeks 3-5)
├── Week 3: LU & Cholesky
│   ├── Day 1: dgetrf (LU factorization)
│   ├── Day 2: dgetrs, dgetri (solve, inverse from LU)
│   ├── Day 3: dpotrf (Cholesky factorization)
│   ├── Day 4: dpotrs, dpotri (solve, inverse from Cholesky)
│   └── Day 5: Tests + error handling
│
├── Week 4: QR & Eigenvalues
│   ├── Day 1: dgeqrf (QR with Householder)
│   ├── Day 2: dorgqr (reconstruct Q)
│   ├── Day 3: dgeev (general eigenvalues)
│   ├── Day 4: dsyevd (symmetric eigenvalues, divide-conquer)
│   └── Day 5: Tests + complex variants
│
└── Week 5: SVD & Least Squares
    ├── Day 1: dgesdd (SVD divide-and-conquer, values only)
    ├── Day 2: dgesdd (full SVD with U, Vh)
    ├── Day 3: dgelsd (least squares via SVD)
    ├── Day 4: Complex variants (zgesdd, etc.)
    └── Day 5: Comprehensive tests

Phase 13.3: TypeScript High-Level API (Weeks 6-7)
├── Week 6: Core Functions
│   ├── Day 1: matmul, dot, vdot
│   ├── Day 2: inner, outer, tensordot
│   ├── Day 3: solve, inv, pinv
│   ├── Day 4: det, slogdet, matrix_rank
│   └── Day 5: norm variants
│
└── Week 7: Advanced Functions
    ├── Day 1: cholesky, qr
    ├── Day 2: svd, svdvals
    ├── Day 3: eig, eigh, eigvals, eigvalsh
    ├── Day 4: lstsq, matrix_power
    └── Day 5: einsum, kron, cross

Phase 13.4: Polish & Testing (Week 8)
├── Day 1: Error handling refinement
├── Day 2: Edge cases (empty arrays, 1x1, etc.)
├── Day 3: Numerical accuracy validation vs NumPy
├── Day 4: Performance optimization
└── Day 5: Documentation & examples
```

---

## Verification Plan

After Phase 13 completion, verify:

```bash
# Build
npm run build

# Run all tests
npm test

# Phase 13 specific tests:

# BLAS
✓ blas_ddot computes correct dot product
✓ blas_dnrm2 handles scale/overflow correctly
✓ blas_dgemm produces correct matrix product
✓ blas_dgemm handles transposed inputs

# LAPACK Decompositions
✓ lapack_dgetrf produces valid LU factors
✓ lapack_dpotrf fails on non-positive-definite
✓ lapack_dgeqrf produces orthogonal Q
✓ lapack_dgeev returns correct eigenvalues
✓ lapack_dsyevd eigenvalues are real and sorted
✓ lapack_dgesdd singular values are non-negative and sorted

# High-Level API
✓ matmul matches NumPy for various shapes
✓ solve(A, b) satisfies A @ x = b
✓ inv(A) @ A ≈ I
✓ det matches NumPy
✓ cholesky(A) @ cholesky(A).T ≈ A
✓ Q @ R ≈ A for QR decomposition
✓ U @ diag(S) @ Vh ≈ A for SVD
✓ eigenvectors satisfy A @ v = lambda * v
✓ lstsq minimizes ||b - A @ x||

# Edge Cases
✓ 1x1 matrices work correctly
✓ Empty arrays handled appropriately
✓ Singular matrices raise LinAlgError
✓ Non-square matrices work where applicable
✓ Batched operations broadcast correctly
```

Generate NumPy comparison vectors:

```python
# tests/python/generate_linalg_tests.py
import numpy as np
import json

np.random.seed(42)

tests = {
    "matmul": [
        {"a_shape": [3, 4], "b_shape": [4, 5], "result_shape": [3, 5]},
        {"a_shape": [2, 3, 4], "b_shape": [4, 5], "result_shape": [2, 3, 5]},
    ],
    "solve": [
        {"a": [[3, 1], [1, 2]], "b": [9, 8], "x": [2, 3]},
        {"a": [[1, 2, 3], [4, 5, 6], [7, 8, 10]], "b": [1, 2, 3]},
    ],
    "det": [
        {"a": [[1, 2], [3, 4]], "det": -2},
        {"a": [[1, 0, 0], [0, 2, 0], [0, 0, 3]], "det": 6},
    ],
    "eig": [
        {"a": [[1, 0], [0, 2]], "eigenvalues": [1, 2]},
        {"a": [[0, -1], [1, 0]], "eigenvalues_complex": True},
    ],
    "svd": [
        {"a": [[1, 0], [0, 2], [0, 0]], "singular_values": [2, 1]},
    ],
    "qr": [
        {"a_shape": [4, 3], "q_shape": [4, 3], "r_shape": [3, 3]},
    ],
    "cholesky": [
        {"a": [[4, 2], [2, 5]], "L": [[2, 0], [1, 2]]},
    ],
    "norm": [
        {"x": [1, 2, 3], "ord": 2, "result": 3.7416573867739413},
        {"x": [[1, 2], [3, 4]], "ord": "fro", "result": 5.477225575051661},
    ],
}

# Generate actual test data with NumPy
for test_type, cases in tests.items():
    for case in cases:
        if "a" in case and isinstance(case["a"], list):
            a = np.array(case["a"], dtype=np.float64)
            if test_type == "det":
                case["result"] = float(np.linalg.det(a))
            elif test_type == "solve":
                b = np.array(case["b"], dtype=np.float64)
                case["result"] = np.linalg.solve(a, b).tolist()
            elif test_type == "eig":
                w, v = np.linalg.eig(a)
                case["eigenvalues"] = w.tolist()
                case["eigenvectors"] = v.tolist()

with open("tests/fixtures/linalg_vectors.json", "w") as f:
    json.dump(tests, f, indent=2)
```

---

## Critical Dependencies for Later Phases

Phase 13 completion enables:

- **Phase 14 (numpy.fft)**: FFT can use matmul for DFT matrix approach
- **Phase 15 (numpy.random)**: Multivariate distributions need Cholesky decomposition
- **Statistics**: Covariance matrices, regression analysis
- **Machine Learning**: PCA (via SVD), linear regression (via lstsq)

Phase 13 should be implemented after completing:
- Phase 4 (Ufuncs) - needed for element-wise operations in linalg
- Phase 5 (Array Manipulation) - needed for reshape, transpose operations

---

## Performance Considerations

### Memory Layout

```
Column-Major (Fortran) vs Row-Major (C):
- LAPACK uses column-major (Fortran) layout
- NumJS uses row-major (C) layout by default
- Need transpose or layout conversion at boundaries
- Consider storing matrices column-major for linalg operations
```

### Workspace Management

```
LAPACK Workspace Pattern:
1. Query optimal workspace: call with lwork = -1
2. Allocate workspace: work = malloc(optimal_lwork * sizeof)
3. Execute: call with actual workspace
4. Free workspace

For WASM:
- Pre-allocate workspace pools
- Reuse buffers for repeated operations
- Consider lazy allocation for rarely-used functions
```

### Numerical Stability

```
Key Considerations:
- Use scaled norms (dnrm2) to avoid overflow
- Use slogdet instead of det for large matrices
- Use lstsq with appropriate rcond for rank-deficient systems
- Prefer divide-and-conquer algorithms (gesdd, syevd)
```

---

## API Compatibility Notes

### NumPy Differences

```typescript
// NumPy: returns tuple
// w, v = np.linalg.eig(a)

// NumJS: returns object (more TypeScript-friendly)
// const { eigenvalues, eigenvectors } = linalg.eig(a);
```

### Broadcasting Support

```typescript
// Batch operations on stacked matrices
const A = np.random.rand(10, 3, 3);  // 10 matrices of 3x3
const B = np.random.rand(10, 3, 2);  // 10 matrices of 3x2

const X = linalg.solve(A, B);  // Solves all 10 systems
// X.shape = [10, 3, 2]
```

### Complex Number Handling

```typescript
// Complex eigenvalues for non-symmetric matrices
const A = [[0, -1], [1, 0]];  // Rotation matrix
const { eigenvalues, eigenvectors } = linalg.eig(A);
// eigenvalues: [0+1j, 0-1j] (complex)
```
