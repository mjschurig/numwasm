/*
 * NumJS Linear Algebra - High-level WASM Interface
 *
 * Provides the interface between TypeScript and BLAS/LAPACK routines.
 * Handles NDArray memory layout conversion (C-order to Fortran-order).
 */

#include "ndarray.h"
#include "blas.h"
#include "lapack.h"
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/* Helper macros */
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

/* ============ Memory Layout Conversion ============ */

/**
 * Convert C-order (row-major) 2D array to Fortran-order (column-major).
 * Input: data[i * n + j] for row i, column j
 * Output: out[i + j * m] for row i, column j
 */
static void c_to_fortran_d(const double* data, double* out, int32_t m, int32_t n) {
    for (int32_t j = 0; j < n; j++) {
        for (int32_t i = 0; i < m; i++) {
            out[i + j * m] = data[i * n + j];
        }
    }
}

/**
 * Convert Fortran-order (column-major) to C-order (row-major).
 */
static void fortran_to_c_d(const double* data, double* out, int32_t m, int32_t n) {
    for (int32_t i = 0; i < m; i++) {
        for (int32_t j = 0; j < n; j++) {
            out[i * n + j] = data[i + j * m];
        }
    }
}

static void c_to_fortran_f(const float* data, float* out, int32_t m, int32_t n) {
    for (int32_t j = 0; j < n; j++) {
        for (int32_t i = 0; i < m; i++) {
            out[i + j * m] = data[i * n + j];
        }
    }
}

static void fortran_to_c_f(const float* data, float* out, int32_t m, int32_t n) {
    for (int32_t i = 0; i < m; i++) {
        for (int32_t j = 0; j < n; j++) {
            out[i * n + j] = data[i + j * m];
        }
    }
}

/* ============ Matrix Multiplication ============ */

/**
 * Matrix multiplication: C = A @ B
 *
 * @param a  Input NDArray A (..., M, K)
 * @param b  Input NDArray B (..., K, N)
 * @return   Result NDArray C (..., M, N), or NULL on error
 */
EXPORT NDArray* linalg_matmul(const NDArray* a, const NDArray* b) {
    if (!a || !b) return NULL;
    if (a->ndim < 2 || b->ndim < 2) return NULL;

    int32_t m = a->shape[a->ndim - 2];
    int32_t k1 = a->shape[a->ndim - 1];
    int32_t k2 = b->shape[b->ndim - 2];
    int32_t n = b->shape[b->ndim - 1];

    /* Check inner dimensions match */
    if (k1 != k2) return NULL;
    int32_t k = k1;

    /* For now, only handle 2D matrices */
    if (a->ndim != 2 || b->ndim != 2) {
        /* TODO: Handle batched matmul */
        return NULL;
    }

    /* Create result array */
    int32_t result_shape[2] = {m, n};
    NDArray* result = ndarray_empty(2, result_shape, a->dtype);
    if (!result) return NULL;

    /* Get data pointers */
    const double* A_data = (const double*)a->data;
    const double* B_data = (const double*)b->data;
    double* C_data = (double*)result->data;

    if (a->dtype == DTYPE_FLOAT64) {
        /* Allocate Fortran-order work arrays */
        double* A_f = (double*)malloc(m * k * sizeof(double));
        double* B_f = (double*)malloc(k * n * sizeof(double));
        double* C_f = (double*)malloc(m * n * sizeof(double));

        if (!A_f || !B_f || !C_f) {
            free(A_f); free(B_f); free(C_f);
            ndarray_free(result);
            return NULL;
        }

        /* Convert to Fortran order */
        c_to_fortran_d(A_data, A_f, m, k);
        c_to_fortran_d(B_data, B_f, k, n);

        /* Call BLAS dgemm: C = 1.0 * A * B + 0.0 * C */
        blas_dgemm('N', 'N', m, n, k, 1.0, A_f, m, B_f, k, 0.0, C_f, m);

        /* Convert back to C order */
        fortran_to_c_d(C_f, C_data, m, n);

        free(A_f);
        free(B_f);
        free(C_f);
    } else if (a->dtype == DTYPE_FLOAT32) {
        const float* A_data_f = (const float*)a->data;
        const float* B_data_f = (const float*)b->data;
        float* C_data_f = (float*)result->data;

        float* A_f = (float*)malloc(m * k * sizeof(float));
        float* B_f = (float*)malloc(k * n * sizeof(float));
        float* C_f = (float*)malloc(m * n * sizeof(float));

        if (!A_f || !B_f || !C_f) {
            free(A_f); free(B_f); free(C_f);
            ndarray_free(result);
            return NULL;
        }

        c_to_fortran_f(A_data_f, A_f, m, k);
        c_to_fortran_f(B_data_f, B_f, k, n);

        blas_sgemm('N', 'N', m, n, k, 1.0f, A_f, m, B_f, k, 0.0f, C_f, m);

        fortran_to_c_f(C_f, C_data_f, m, n);

        free(A_f);
        free(B_f);
        free(C_f);
    } else {
        ndarray_free(result);
        return NULL;
    }

    return result;
}

/* ============ Dot Product ============ */

/**
 * Dot product of two arrays.
 * For 1D arrays: inner product.
 * For 2D arrays: matrix multiplication.
 */
EXPORT NDArray* linalg_dot(const NDArray* a, const NDArray* b) {
    if (!a || !b) return NULL;

    /* 1D vectors: inner product */
    if (a->ndim == 1 && b->ndim == 1) {
        if (a->shape[0] != b->shape[0]) return NULL;

        int32_t n = a->shape[0];
        NDArray* result = ndarray_scalar(0.0, a->dtype);
        if (!result) return NULL;

        if (a->dtype == DTYPE_FLOAT64) {
            double dot = blas_ddot(n, (const double*)a->data, 1,
                                   (const double*)b->data, 1);
            *(double*)result->data = dot;
        } else if (a->dtype == DTYPE_FLOAT32) {
            float dot = blas_sdot(n, (const float*)a->data, 1,
                                  (const float*)b->data, 1);
            *(float*)result->data = dot;
        }

        return result;
    }

    /* 2D: matrix multiplication */
    if (a->ndim == 2 && b->ndim == 2) {
        return linalg_matmul(a, b);
    }

    /* Other cases: not implemented */
    return NULL;
}

/* ============ Linear Solve ============ */

/**
 * Solve linear system A * x = b
 *
 * @param a  Coefficient matrix A (N x N)
 * @param b  Right-hand side b (N) or (N x NRHS)
 * @return   Solution x, or NULL on error (singular matrix)
 */
EXPORT NDArray* linalg_solve(const NDArray* a, const NDArray* b) {
    if (!a || !b) return NULL;
    if (a->ndim != 2) return NULL;
    if (a->shape[0] != a->shape[1]) return NULL;  /* Must be square */

    int32_t n = a->shape[0];
    int32_t nrhs = 1;

    if (b->ndim == 1) {
        if (b->shape[0] != n) return NULL;
        nrhs = 1;
    } else if (b->ndim == 2) {
        if (b->shape[0] != n) return NULL;
        nrhs = b->shape[1];
    } else {
        return NULL;
    }

    /* Create result array (copy of b) */
    NDArray* x = ndarray_copy(b);
    if (!x) return NULL;

    /* Allocate work arrays */
    int32_t* ipiv = (int32_t*)malloc(n * sizeof(int32_t));
    if (!ipiv) {
        ndarray_free(x);
        return NULL;
    }

    if (a->dtype == DTYPE_FLOAT64) {
        /* Copy A and convert to Fortran order */
        double* A_f = (double*)malloc(n * n * sizeof(double));
        double* B_f = (double*)malloc(n * nrhs * sizeof(double));

        if (!A_f || !B_f) {
            free(ipiv); free(A_f); free(B_f);
            ndarray_free(x);
            return NULL;
        }

        c_to_fortran_d((const double*)a->data, A_f, n, n);
        if (b->ndim == 1) {
            memcpy(B_f, b->data, n * sizeof(double));
        } else {
            c_to_fortran_d((const double*)b->data, B_f, n, nrhs);
        }

        /* LU factorization */
        int32_t info = lapack_dgetrf(n, n, A_f, n, ipiv);
        if (info > 0) {
            /* Singular matrix */
            free(ipiv); free(A_f); free(B_f);
            ndarray_free(x);
            return NULL;
        }

        /* Solve */
        info = lapack_dgetrs('N', n, nrhs, A_f, n, ipiv, B_f, n);

        /* Convert result back */
        if (b->ndim == 1) {
            memcpy(x->data, B_f, n * sizeof(double));
        } else {
            fortran_to_c_d(B_f, (double*)x->data, n, nrhs);
        }

        free(A_f);
        free(B_f);
    } else if (a->dtype == DTYPE_FLOAT32) {
        float* A_f = (float*)malloc(n * n * sizeof(float));
        float* B_f = (float*)malloc(n * nrhs * sizeof(float));

        if (!A_f || !B_f) {
            free(ipiv); free(A_f); free(B_f);
            ndarray_free(x);
            return NULL;
        }

        c_to_fortran_f((const float*)a->data, A_f, n, n);
        if (b->ndim == 1) {
            memcpy(B_f, b->data, n * sizeof(float));
        } else {
            c_to_fortran_f((const float*)b->data, B_f, n, nrhs);
        }

        int32_t info = lapack_sgetrf(n, n, A_f, n, ipiv);
        if (info > 0) {
            free(ipiv); free(A_f); free(B_f);
            ndarray_free(x);
            return NULL;
        }

        info = lapack_sgetrs('N', n, nrhs, A_f, n, ipiv, B_f, n);

        if (b->ndim == 1) {
            memcpy(x->data, B_f, n * sizeof(float));
        } else {
            fortran_to_c_f(B_f, (float*)x->data, n, nrhs);
        }

        free(A_f);
        free(B_f);
    } else {
        free(ipiv);
        ndarray_free(x);
        return NULL;
    }

    free(ipiv);
    return x;
}

/* ============ Matrix Inverse ============ */

/**
 * Compute matrix inverse.
 *
 * @param a  Square matrix A (N x N)
 * @return   Inverse A^(-1), or NULL on error (singular matrix)
 */
EXPORT NDArray* linalg_inv(const NDArray* a) {
    if (!a) return NULL;
    if (a->ndim != 2) return NULL;
    if (a->shape[0] != a->shape[1]) return NULL;

    int32_t n = a->shape[0];

    /* Create result array */
    NDArray* result = ndarray_copy(a);
    if (!result) return NULL;

    int32_t* ipiv = (int32_t*)malloc(n * sizeof(int32_t));
    if (!ipiv) {
        ndarray_free(result);
        return NULL;
    }

    if (a->dtype == DTYPE_FLOAT64) {
        double* A_f = (double*)malloc(n * n * sizeof(double));
        double* work = (double*)malloc(n * sizeof(double));

        if (!A_f || !work) {
            free(ipiv); free(A_f); free(work);
            ndarray_free(result);
            return NULL;
        }

        c_to_fortran_d((const double*)a->data, A_f, n, n);

        /* LU factorization */
        int32_t info = lapack_dgetrf(n, n, A_f, n, ipiv);
        if (info > 0) {
            free(ipiv); free(A_f); free(work);
            ndarray_free(result);
            return NULL;
        }

        /* Compute inverse */
        info = lapack_dgetri(n, A_f, n, ipiv, work, n);
        if (info > 0) {
            free(ipiv); free(A_f); free(work);
            ndarray_free(result);
            return NULL;
        }

        fortran_to_c_d(A_f, (double*)result->data, n, n);

        free(A_f);
        free(work);
    } else if (a->dtype == DTYPE_FLOAT32) {
        float* A_f = (float*)malloc(n * n * sizeof(float));
        float* work = (float*)malloc(n * sizeof(float));

        if (!A_f || !work) {
            free(ipiv); free(A_f); free(work);
            ndarray_free(result);
            return NULL;
        }

        c_to_fortran_f((const float*)a->data, A_f, n, n);

        int32_t info = lapack_sgetrf(n, n, A_f, n, ipiv);
        if (info > 0) {
            free(ipiv); free(A_f); free(work);
            ndarray_free(result);
            return NULL;
        }

        info = lapack_sgetri(n, A_f, n, ipiv, work, n);
        if (info > 0) {
            free(ipiv); free(A_f); free(work);
            ndarray_free(result);
            return NULL;
        }

        fortran_to_c_f(A_f, (float*)result->data, n, n);

        free(A_f);
        free(work);
    } else {
        free(ipiv);
        ndarray_free(result);
        return NULL;
    }

    free(ipiv);
    return result;
}

/* ============ Determinant ============ */

/**
 * Compute matrix determinant using LU factorization.
 *
 * @param a  Square matrix A (N x N)
 * @return   Scalar array containing determinant, or NULL on error
 */
EXPORT NDArray* linalg_det(const NDArray* a) {
    if (!a) return NULL;
    if (a->ndim != 2) return NULL;
    if (a->shape[0] != a->shape[1]) return NULL;

    int32_t n = a->shape[0];

    NDArray* result = ndarray_scalar(0.0, a->dtype);
    if (!result) return NULL;

    int32_t* ipiv = (int32_t*)malloc(n * sizeof(int32_t));
    if (!ipiv) {
        ndarray_free(result);
        return NULL;
    }

    if (a->dtype == DTYPE_FLOAT64) {
        double* A_f = (double*)malloc(n * n * sizeof(double));
        if (!A_f) {
            free(ipiv);
            ndarray_free(result);
            return NULL;
        }

        c_to_fortran_d((const double*)a->data, A_f, n, n);

        /* LU factorization */
        int32_t info = lapack_dgetrf(n, n, A_f, n, ipiv);

        /* Determinant = product of diagonal of U * sign from pivots */
        double det = 1.0;
        int32_t num_swaps = 0;

        for (int32_t i = 0; i < n; i++) {
            det *= A_f[i + i * n];
            if (ipiv[i] != i) {
                num_swaps++;
            }
        }

        if (num_swaps % 2 == 1) {
            det = -det;
        }

        *(double*)result->data = det;

        free(A_f);
    } else if (a->dtype == DTYPE_FLOAT32) {
        float* A_f = (float*)malloc(n * n * sizeof(float));
        if (!A_f) {
            free(ipiv);
            ndarray_free(result);
            return NULL;
        }

        c_to_fortran_f((const float*)a->data, A_f, n, n);

        int32_t info = lapack_sgetrf(n, n, A_f, n, ipiv);

        float det = 1.0f;
        int32_t num_swaps = 0;

        for (int32_t i = 0; i < n; i++) {
            det *= A_f[i + i * n];
            if (ipiv[i] != i) {
                num_swaps++;
            }
        }

        if (num_swaps % 2 == 1) {
            det = -det;
        }

        *(float*)result->data = det;

        free(A_f);
    } else {
        free(ipiv);
        ndarray_free(result);
        return NULL;
    }

    free(ipiv);
    return result;
}

/* ============ Cholesky Decomposition ============ */

/**
 * Cholesky decomposition for positive definite matrices.
 *
 * @param a     Symmetric positive definite matrix A (N x N)
 * @param upper If non-zero, return upper triangular U (A = U^T U)
 *              Otherwise return lower triangular L (A = L L^T)
 * @return      Cholesky factor, or NULL on error
 */
EXPORT NDArray* linalg_cholesky(const NDArray* a, int32_t upper) {
    if (!a) return NULL;
    if (a->ndim != 2) return NULL;
    if (a->shape[0] != a->shape[1]) return NULL;

    int32_t n = a->shape[0];

    NDArray* result = ndarray_copy(a);
    if (!result) return NULL;

    if (a->dtype == DTYPE_FLOAT64) {
        double* A_f = (double*)malloc(n * n * sizeof(double));
        if (!A_f) {
            ndarray_free(result);
            return NULL;
        }

        c_to_fortran_d((const double*)a->data, A_f, n, n);

        int32_t info = lapack_dpotrf(upper ? 'U' : 'L', n, A_f, n);
        if (info > 0) {
            /* Not positive definite */
            free(A_f);
            ndarray_free(result);
            return NULL;
        }

        /* Zero out the other triangle */
        if (upper) {
            for (int32_t j = 0; j < n; j++) {
                for (int32_t i = j + 1; i < n; i++) {
                    A_f[i + j * n] = 0.0;
                }
            }
        } else {
            for (int32_t j = 0; j < n; j++) {
                for (int32_t i = 0; i < j; i++) {
                    A_f[i + j * n] = 0.0;
                }
            }
        }

        fortran_to_c_d(A_f, (double*)result->data, n, n);

        free(A_f);
    } else if (a->dtype == DTYPE_FLOAT32) {
        float* A_f = (float*)malloc(n * n * sizeof(float));
        if (!A_f) {
            ndarray_free(result);
            return NULL;
        }

        c_to_fortran_f((const float*)a->data, A_f, n, n);

        int32_t info = lapack_spotrf(upper ? 'U' : 'L', n, A_f, n);
        if (info > 0) {
            free(A_f);
            ndarray_free(result);
            return NULL;
        }

        if (upper) {
            for (int32_t j = 0; j < n; j++) {
                for (int32_t i = j + 1; i < n; i++) {
                    A_f[i + j * n] = 0.0f;
                }
            }
        } else {
            for (int32_t j = 0; j < n; j++) {
                for (int32_t i = 0; i < j; i++) {
                    A_f[i + j * n] = 0.0f;
                }
            }
        }

        fortran_to_c_f(A_f, (float*)result->data, n, n);

        free(A_f);
    } else {
        ndarray_free(result);
        return NULL;
    }

    return result;
}

/* ============ QR Decomposition ============ */

/**
 * QR decomposition result structure.
 */
typedef struct {
    NDArray* Q;
    NDArray* R;
} QRResult;

/**
 * QR decomposition.
 *
 * @param a  Matrix A (M x N)
 * @return   QRResult with Q (M x K) and R (K x N) where K = min(M, N),
 *           or NULL on error
 */
EXPORT QRResult* linalg_qr(const NDArray* a) {
    if (!a) return NULL;
    if (a->ndim != 2) return NULL;

    int32_t m = a->shape[0];
    int32_t n = a->shape[1];
    int32_t k = MIN(m, n);

    QRResult* result = (QRResult*)malloc(sizeof(QRResult));
    if (!result) return NULL;

    /* Create Q and R arrays */
    int32_t q_shape[2] = {m, k};
    int32_t r_shape[2] = {k, n};

    result->Q = ndarray_empty(2, q_shape, a->dtype);
    result->R = ndarray_empty(2, r_shape, a->dtype);

    if (!result->Q || !result->R) {
        if (result->Q) ndarray_free(result->Q);
        if (result->R) ndarray_free(result->R);
        free(result);
        return NULL;
    }

    if (a->dtype == DTYPE_FLOAT64) {
        double* A_f = (double*)malloc(m * n * sizeof(double));
        double* tau = (double*)malloc(k * sizeof(double));
        double* work = (double*)malloc(n * sizeof(double));

        if (!A_f || !tau || !work) {
            free(A_f); free(tau); free(work);
            ndarray_free(result->Q);
            ndarray_free(result->R);
            free(result);
            return NULL;
        }

        c_to_fortran_d((const double*)a->data, A_f, m, n);

        /* QR factorization */
        int32_t info = lapack_dgeqrf(m, n, A_f, m, tau, work, n);

        /* Extract R (upper triangular part) */
        double* R_f = (double*)calloc(k * n, sizeof(double));
        for (int32_t j = 0; j < n; j++) {
            for (int32_t i = 0; i <= MIN(j, k - 1); i++) {
                R_f[i + j * k] = A_f[i + j * m];
            }
        }
        fortran_to_c_d(R_f, (double*)result->R->data, k, n);
        free(R_f);

        /* Generate Q */
        info = lapack_dorgqr(m, k, k, A_f, m, tau, work, n);

        /* Copy Q (first k columns) */
        double* Q_f = (double*)malloc(m * k * sizeof(double));
        for (int32_t j = 0; j < k; j++) {
            for (int32_t i = 0; i < m; i++) {
                Q_f[i + j * m] = A_f[i + j * m];
            }
        }
        fortran_to_c_d(Q_f, (double*)result->Q->data, m, k);
        free(Q_f);

        free(A_f);
        free(tau);
        free(work);
    } else if (a->dtype == DTYPE_FLOAT32) {
        float* A_f = (float*)malloc(m * n * sizeof(float));
        float* tau = (float*)malloc(k * sizeof(float));
        float* work = (float*)malloc(n * sizeof(float));

        if (!A_f || !tau || !work) {
            free(A_f); free(tau); free(work);
            ndarray_free(result->Q);
            ndarray_free(result->R);
            free(result);
            return NULL;
        }

        c_to_fortran_f((const float*)a->data, A_f, m, n);

        int32_t info = lapack_sgeqrf(m, n, A_f, m, tau, work, n);

        float* R_f = (float*)calloc(k * n, sizeof(float));
        for (int32_t j = 0; j < n; j++) {
            for (int32_t i = 0; i <= MIN(j, k - 1); i++) {
                R_f[i + j * k] = A_f[i + j * m];
            }
        }
        fortran_to_c_f(R_f, (float*)result->R->data, k, n);
        free(R_f);

        info = lapack_sorgqr(m, k, k, A_f, m, tau, work, n);

        float* Q_f = (float*)malloc(m * k * sizeof(float));
        for (int32_t j = 0; j < k; j++) {
            for (int32_t i = 0; i < m; i++) {
                Q_f[i + j * m] = A_f[i + j * m];
            }
        }
        fortran_to_c_f(Q_f, (float*)result->Q->data, m, k);
        free(Q_f);

        free(A_f);
        free(tau);
        free(work);
    } else {
        ndarray_free(result->Q);
        ndarray_free(result->R);
        free(result);
        return NULL;
    }

    return result;
}

EXPORT void linalg_qr_free(QRResult* result) {
    if (result) {
        if (result->Q) ndarray_free(result->Q);
        if (result->R) ndarray_free(result->R);
        free(result);
    }
}

EXPORT NDArray* linalg_qr_get_q(QRResult* result) {
    return result ? result->Q : NULL;
}

EXPORT NDArray* linalg_qr_get_r(QRResult* result) {
    return result ? result->R : NULL;
}

/* ============ Eigenvalue Decomposition ============ */

typedef struct {
    NDArray* eigenvalues;   /* Real parts (or all if complex) */
    NDArray* eigenvalues_i; /* Imaginary parts (may be NULL) */
    NDArray* eigenvectors;
} EigResult;

/**
 * Compute eigenvalues and eigenvectors of a square matrix.
 */
EXPORT EigResult* linalg_eig(const NDArray* a) {
    if (!a) return NULL;
    if (a->ndim != 2) return NULL;
    if (a->shape[0] != a->shape[1]) return NULL;

    int32_t n = a->shape[0];

    EigResult* result = (EigResult*)malloc(sizeof(EigResult));
    if (!result) return NULL;

    int32_t ev_shape[1] = {n};
    int32_t vec_shape[2] = {n, n};

    result->eigenvalues = ndarray_empty(1, ev_shape, a->dtype);
    result->eigenvalues_i = ndarray_empty(1, ev_shape, a->dtype);
    result->eigenvectors = ndarray_empty(2, vec_shape, a->dtype);

    if (!result->eigenvalues || !result->eigenvalues_i || !result->eigenvectors) {
        if (result->eigenvalues) ndarray_free(result->eigenvalues);
        if (result->eigenvalues_i) ndarray_free(result->eigenvalues_i);
        if (result->eigenvectors) ndarray_free(result->eigenvectors);
        free(result);
        return NULL;
    }

    if (a->dtype == DTYPE_FLOAT64) {
        double* A_f = (double*)malloc(n * n * sizeof(double));
        double* wr = (double*)malloc(n * sizeof(double));
        double* wi = (double*)malloc(n * sizeof(double));
        double* VR = (double*)malloc(n * n * sizeof(double));
        double* work = (double*)malloc(4 * n * sizeof(double));

        if (!A_f || !wr || !wi || !VR || !work) {
            free(A_f); free(wr); free(wi); free(VR); free(work);
            ndarray_free(result->eigenvalues);
            ndarray_free(result->eigenvalues_i);
            ndarray_free(result->eigenvectors);
            free(result);
            return NULL;
        }

        c_to_fortran_d((const double*)a->data, A_f, n, n);

        int32_t info = lapack_dgeev('N', 'V', n, A_f, n, wr, wi,
                                     NULL, 1, VR, n, work, 4 * n);

        /* Copy eigenvalues */
        memcpy(result->eigenvalues->data, wr, n * sizeof(double));
        memcpy(result->eigenvalues_i->data, wi, n * sizeof(double));

        /* Copy eigenvectors */
        fortran_to_c_d(VR, (double*)result->eigenvectors->data, n, n);

        free(A_f); free(wr); free(wi); free(VR); free(work);
    } else if (a->dtype == DTYPE_FLOAT32) {
        float* A_f = (float*)malloc(n * n * sizeof(float));
        float* wr = (float*)malloc(n * sizeof(float));
        float* wi = (float*)malloc(n * sizeof(float));
        float* VR = (float*)malloc(n * n * sizeof(float));
        float* work = (float*)malloc(4 * n * sizeof(float));

        if (!A_f || !wr || !wi || !VR || !work) {
            free(A_f); free(wr); free(wi); free(VR); free(work);
            ndarray_free(result->eigenvalues);
            ndarray_free(result->eigenvalues_i);
            ndarray_free(result->eigenvectors);
            free(result);
            return NULL;
        }

        c_to_fortran_f((const float*)a->data, A_f, n, n);

        int32_t info = lapack_sgeev('N', 'V', n, A_f, n, wr, wi,
                                     NULL, 1, VR, n, work, 4 * n);

        memcpy(result->eigenvalues->data, wr, n * sizeof(float));
        memcpy(result->eigenvalues_i->data, wi, n * sizeof(float));

        fortran_to_c_f(VR, (float*)result->eigenvectors->data, n, n);

        free(A_f); free(wr); free(wi); free(VR); free(work);
    } else {
        ndarray_free(result->eigenvalues);
        ndarray_free(result->eigenvalues_i);
        ndarray_free(result->eigenvectors);
        free(result);
        return NULL;
    }

    return result;
}

EXPORT void linalg_eig_free(EigResult* result) {
    if (result) {
        if (result->eigenvalues) ndarray_free(result->eigenvalues);
        if (result->eigenvalues_i) ndarray_free(result->eigenvalues_i);
        if (result->eigenvectors) ndarray_free(result->eigenvectors);
        free(result);
    }
}

EXPORT NDArray* linalg_eig_get_values(EigResult* result) {
    return result ? result->eigenvalues : NULL;
}

EXPORT NDArray* linalg_eig_get_values_imag(EigResult* result) {
    return result ? result->eigenvalues_i : NULL;
}

EXPORT NDArray* linalg_eig_get_vectors(EigResult* result) {
    return result ? result->eigenvectors : NULL;
}

/* ============ SVD ============ */

typedef struct {
    NDArray* U;
    NDArray* S;
    NDArray* Vh;
} SVDResult;

/**
 * Singular Value Decomposition: A = U @ diag(S) @ Vh
 */
EXPORT SVDResult* linalg_svd(const NDArray* a, int32_t full_matrices) {
    if (!a) return NULL;
    if (a->ndim != 2) return NULL;

    int32_t m = a->shape[0];
    int32_t n = a->shape[1];
    int32_t k = MIN(m, n);

    SVDResult* result = (SVDResult*)malloc(sizeof(SVDResult));
    if (!result) return NULL;

    int32_t u_cols = full_matrices ? m : k;
    int32_t vh_rows = full_matrices ? n : k;

    int32_t u_shape[2] = {m, u_cols};
    int32_t s_shape[1] = {k};
    int32_t vh_shape[2] = {vh_rows, n};

    result->U = ndarray_empty(2, u_shape, a->dtype);
    result->S = ndarray_empty(1, s_shape, a->dtype);
    result->Vh = ndarray_empty(2, vh_shape, a->dtype);

    if (!result->U || !result->S || !result->Vh) {
        if (result->U) ndarray_free(result->U);
        if (result->S) ndarray_free(result->S);
        if (result->Vh) ndarray_free(result->Vh);
        free(result);
        return NULL;
    }

    char jobz = full_matrices ? 'A' : 'S';

    if (a->dtype == DTYPE_FLOAT64) {
        double* A_f = (double*)malloc(m * n * sizeof(double));
        double* U_f = (double*)malloc(m * u_cols * sizeof(double));
        double* VT_f = (double*)malloc(vh_rows * n * sizeof(double));
        double* s = (double*)malloc(k * sizeof(double));
        int32_t lwork = 4 * k * k + 7 * k + MAX(m, n);
        double* work = (double*)malloc(lwork * sizeof(double));
        int32_t* iwork = (int32_t*)malloc(8 * k * sizeof(int32_t));

        if (!A_f || !U_f || !VT_f || !s || !work || !iwork) {
            free(A_f); free(U_f); free(VT_f); free(s); free(work); free(iwork);
            ndarray_free(result->U);
            ndarray_free(result->S);
            ndarray_free(result->Vh);
            free(result);
            return NULL;
        }

        c_to_fortran_d((const double*)a->data, A_f, m, n);

        int32_t info = lapack_dgesdd(jobz, m, n, A_f, m, s, U_f, m, VT_f, vh_rows,
                                      work, lwork, iwork);

        /* Copy results */
        memcpy(result->S->data, s, k * sizeof(double));
        fortran_to_c_d(U_f, (double*)result->U->data, m, u_cols);
        fortran_to_c_d(VT_f, (double*)result->Vh->data, vh_rows, n);

        free(A_f); free(U_f); free(VT_f); free(s); free(work); free(iwork);
    } else if (a->dtype == DTYPE_FLOAT32) {
        float* A_f = (float*)malloc(m * n * sizeof(float));
        float* U_f = (float*)malloc(m * u_cols * sizeof(float));
        float* VT_f = (float*)malloc(vh_rows * n * sizeof(float));
        float* s = (float*)malloc(k * sizeof(float));
        int32_t lwork = 4 * k * k + 7 * k + MAX(m, n);
        float* work = (float*)malloc(lwork * sizeof(float));
        int32_t* iwork = (int32_t*)malloc(8 * k * sizeof(int32_t));

        if (!A_f || !U_f || !VT_f || !s || !work || !iwork) {
            free(A_f); free(U_f); free(VT_f); free(s); free(work); free(iwork);
            ndarray_free(result->U);
            ndarray_free(result->S);
            ndarray_free(result->Vh);
            free(result);
            return NULL;
        }

        c_to_fortran_f((const float*)a->data, A_f, m, n);

        int32_t info = lapack_sgesdd(jobz, m, n, A_f, m, s, U_f, m, VT_f, vh_rows,
                                      work, lwork, iwork);

        memcpy(result->S->data, s, k * sizeof(float));
        fortran_to_c_f(U_f, (float*)result->U->data, m, u_cols);
        fortran_to_c_f(VT_f, (float*)result->Vh->data, vh_rows, n);

        free(A_f); free(U_f); free(VT_f); free(s); free(work); free(iwork);
    } else {
        ndarray_free(result->U);
        ndarray_free(result->S);
        ndarray_free(result->Vh);
        free(result);
        return NULL;
    }

    return result;
}

EXPORT void linalg_svd_free(SVDResult* result) {
    if (result) {
        if (result->U) ndarray_free(result->U);
        if (result->S) ndarray_free(result->S);
        if (result->Vh) ndarray_free(result->Vh);
        free(result);
    }
}

EXPORT NDArray* linalg_svd_get_u(SVDResult* result) {
    return result ? result->U : NULL;
}

EXPORT NDArray* linalg_svd_get_s(SVDResult* result) {
    return result ? result->S : NULL;
}

EXPORT NDArray* linalg_svd_get_vh(SVDResult* result) {
    return result ? result->Vh : NULL;
}

/* ============ Vector/Matrix Norm ============ */

/**
 * Compute norm of an array.
 *
 * @param a     Input array
 * @param ord   Norm order: 0=Frobenius, 1=L1, 2=L2, -1=inf, -2=-inf
 * @return      Scalar array with norm value
 */
EXPORT NDArray* linalg_norm(const NDArray* a, int32_t ord) {
    if (!a) return NULL;

    NDArray* result = ndarray_scalar(0.0, a->dtype);
    if (!result) return NULL;

    int32_t size = a->size;

    if (a->dtype == DTYPE_FLOAT64) {
        const double* data = (const double*)a->data;
        double norm_val = 0.0;

        switch (ord) {
            case 2:  /* L2 / Frobenius (default) */
            case 0:
                norm_val = blas_dnrm2(size, data, 1);
                break;
            case 1:  /* L1 norm */
                norm_val = blas_dasum(size, data, 1);
                break;
            case -1: /* Inf norm (max abs) */
                {
                    int32_t idx = blas_idamax(size, data, 1);
                    norm_val = fabs(data[idx]);
                }
                break;
            case -2: /* -Inf norm (min abs) */
                norm_val = fabs(data[0]);
                for (int32_t i = 1; i < size; i++) {
                    double absval = fabs(data[i]);
                    if (absval < norm_val) norm_val = absval;
                }
                break;
            default: /* p-norm */
                {
                    double sum = 0.0;
                    double p = (double)ord;
                    for (int32_t i = 0; i < size; i++) {
                        sum += pow(fabs(data[i]), p);
                    }
                    norm_val = pow(sum, 1.0 / p);
                }
        }

        *(double*)result->data = norm_val;
    } else if (a->dtype == DTYPE_FLOAT32) {
        const float* data = (const float*)a->data;
        float norm_val = 0.0f;

        switch (ord) {
            case 2:
            case 0:
                norm_val = blas_snrm2(size, data, 1);
                break;
            case 1:
                norm_val = blas_sasum(size, data, 1);
                break;
            case -1:
                {
                    int32_t idx = blas_isamax(size, data, 1);
                    norm_val = fabsf(data[idx]);
                }
                break;
            case -2:
                norm_val = fabsf(data[0]);
                for (int32_t i = 1; i < size; i++) {
                    float absval = fabsf(data[i]);
                    if (absval < norm_val) norm_val = absval;
                }
                break;
            default:
                {
                    float sum = 0.0f;
                    float p = (float)ord;
                    for (int32_t i = 0; i < size; i++) {
                        sum += powf(fabsf(data[i]), p);
                    }
                    norm_val = powf(sum, 1.0f / p);
                }
        }

        *(float*)result->data = norm_val;
    } else {
        ndarray_free(result);
        return NULL;
    }

    return result;
}

/* ============ Cholesky Solve ============ */

/**
 * Solve linear system using Cholesky factorization.
 * Solves A * x = b where c is the Cholesky factor of A.
 *
 * @param c      Cholesky factor from linalg_cholesky (N x N)
 * @param b      Right-hand side (N,) or (N x NRHS)
 * @param lower  Non-zero if c is lower triangular (A = L L^T)
 * @return       Solution x, or NULL on error
 */
EXPORT NDArray* linalg_cholesky_solve(const NDArray* c, const NDArray* b, int32_t lower) {
    if (!c || !b) return NULL;
    if (c->ndim != 2 || c->shape[0] != c->shape[1]) return NULL;

    int32_t n = c->shape[0];
    int32_t nrhs;
    int32_t is_vector = (b->ndim == 1);

    if (is_vector) {
        if (b->shape[0] != n) return NULL;
        nrhs = 1;
    } else if (b->ndim == 2) {
        if (b->shape[0] != n) return NULL;
        nrhs = b->shape[1];
    } else {
        return NULL;
    }

    /* Create result array (same shape as b) */
    NDArray* x = ndarray_copy(b);
    if (!x) return NULL;

    if (c->dtype == DTYPE_FLOAT64) {
        double* C_f = (double*)malloc(n * n * sizeof(double));
        double* B_f = (double*)malloc(n * nrhs * sizeof(double));

        if (!C_f || !B_f) {
            free(C_f); free(B_f);
            ndarray_free(x);
            return NULL;
        }

        c_to_fortran_d((const double*)c->data, C_f, n, n);

        if (is_vector) {
            memcpy(B_f, b->data, n * sizeof(double));
        } else {
            c_to_fortran_d((const double*)b->data, B_f, n, nrhs);
        }

        int32_t info = lapack_dpotrs(lower ? 'L' : 'U', n, nrhs, C_f, n, B_f, n);

        if (info != 0) {
            free(C_f); free(B_f);
            ndarray_free(x);
            return NULL;
        }

        if (is_vector) {
            memcpy(x->data, B_f, n * sizeof(double));
        } else {
            fortran_to_c_d(B_f, (double*)x->data, n, nrhs);
        }

        free(C_f);
        free(B_f);
    } else if (c->dtype == DTYPE_FLOAT32) {
        float* C_f = (float*)malloc(n * n * sizeof(float));
        float* B_f = (float*)malloc(n * nrhs * sizeof(float));

        if (!C_f || !B_f) {
            free(C_f); free(B_f);
            ndarray_free(x);
            return NULL;
        }

        c_to_fortran_f((const float*)c->data, C_f, n, n);

        if (is_vector) {
            memcpy(B_f, b->data, n * sizeof(float));
        } else {
            c_to_fortran_f((const float*)b->data, B_f, n, nrhs);
        }

        /* Note: would need lapack_spotrs for float, using double conversion for now */
        /* TODO: add spotrs to lapack.c */
        free(C_f); free(B_f);
        ndarray_free(x);
        return NULL;  /* Float not yet supported */
    } else {
        ndarray_free(x);
        return NULL;
    }

    return x;
}

/* ============ LU Decomposition ============ */

/**
 * LU decomposition result structure.
 */
typedef struct {
    NDArray* P;  /* Permutation matrix (M x M) */
    NDArray* L;  /* Lower triangular (M x K) where K = min(M, N) */
    NDArray* U;  /* Upper triangular (K x N) */
} LUResult;

/**
 * LU decomposition with partial pivoting.
 * A = P @ L @ U where P is permutation, L is lower triangular with unit diagonal,
 * and U is upper triangular.
 *
 * @param a  Matrix A (M x N)
 * @return   LUResult with P, L, U, or NULL on error
 */
EXPORT LUResult* linalg_lu(const NDArray* a) {
    if (!a) return NULL;
    if (a->ndim != 2) return NULL;

    int32_t m = a->shape[0];
    int32_t n = a->shape[1];
    int32_t k = MIN(m, n);

    LUResult* result = (LUResult*)malloc(sizeof(LUResult));
    if (!result) return NULL;

    /* Create output arrays */
    int32_t p_shape[2] = {m, m};
    int32_t l_shape[2] = {m, k};
    int32_t u_shape[2] = {k, n};

    result->P = ndarray_empty(2, p_shape, DTYPE_FLOAT64);
    result->L = ndarray_empty(2, l_shape, a->dtype);
    result->U = ndarray_empty(2, u_shape, a->dtype);

    if (!result->P || !result->L || !result->U) {
        if (result->P) ndarray_free(result->P);
        if (result->L) ndarray_free(result->L);
        if (result->U) ndarray_free(result->U);
        free(result);
        return NULL;
    }

    if (a->dtype == DTYPE_FLOAT64) {
        double* A_f = (double*)malloc(m * n * sizeof(double));
        int32_t* ipiv = (int32_t*)malloc(k * sizeof(int32_t));

        if (!A_f || !ipiv) {
            free(A_f); free(ipiv);
            ndarray_free(result->P);
            ndarray_free(result->L);
            ndarray_free(result->U);
            free(result);
            return NULL;
        }

        c_to_fortran_d((const double*)a->data, A_f, m, n);

        int32_t info = lapack_dgetrf(m, n, A_f, m, ipiv);
        /* info > 0 means singular, but we still return the decomposition */

        /* Extract L (lower with unit diagonal) */
        double* L_f = (double*)calloc(m * k, sizeof(double));
        for (int32_t j = 0; j < k; j++) {
            L_f[j + j * m] = 1.0;  /* Unit diagonal */
            for (int32_t i = j + 1; i < m; i++) {
                L_f[i + j * m] = A_f[i + j * m];
            }
        }
        fortran_to_c_d(L_f, (double*)result->L->data, m, k);
        free(L_f);

        /* Extract U (upper triangular) */
        double* U_f = (double*)calloc(k * n, sizeof(double));
        for (int32_t j = 0; j < n; j++) {
            for (int32_t i = 0; i <= MIN(j, k - 1); i++) {
                U_f[i + j * k] = A_f[i + j * m];
            }
        }
        fortran_to_c_d(U_f, (double*)result->U->data, k, n);
        free(U_f);

        /* Build permutation matrix from pivot indices */
        /* ipiv is 1-based from LAPACK; ipiv[i] means row i was swapped with row ipiv[i]-1 */
        double* P_data = (double*)result->P->data;
        memset(P_data, 0, m * m * sizeof(double));

        /* Start with identity */
        int32_t* perm = (int32_t*)malloc(m * sizeof(int32_t));
        for (int32_t i = 0; i < m; i++) perm[i] = i;

        /* Apply pivot swaps */
        for (int32_t i = 0; i < k; i++) {
            int32_t swap_row = ipiv[i] - 1;  /* Convert to 0-based */
            if (swap_row != i) {
                int32_t tmp = perm[i];
                perm[i] = perm[swap_row];
                perm[swap_row] = tmp;
            }
        }

        /* Build P matrix: P[perm[i], i] = 1 means row perm[i] of A becomes row i of P@A */
        /* Actually P[i, perm[i]] = 1 for P @ A = L @ U */
        for (int32_t i = 0; i < m; i++) {
            P_data[i * m + perm[i]] = 1.0;
        }

        free(perm);
        free(A_f);
        free(ipiv);
    } else if (a->dtype == DTYPE_FLOAT32) {
        float* A_f = (float*)malloc(m * n * sizeof(float));
        int32_t* ipiv = (int32_t*)malloc(k * sizeof(int32_t));

        if (!A_f || !ipiv) {
            free(A_f); free(ipiv);
            ndarray_free(result->P);
            ndarray_free(result->L);
            ndarray_free(result->U);
            free(result);
            return NULL;
        }

        c_to_fortran_f((const float*)a->data, A_f, m, n);

        int32_t info = lapack_sgetrf(m, n, A_f, m, ipiv);

        /* Extract L */
        float* L_f = (float*)calloc(m * k, sizeof(float));
        for (int32_t j = 0; j < k; j++) {
            L_f[j + j * m] = 1.0f;
            for (int32_t i = j + 1; i < m; i++) {
                L_f[i + j * m] = A_f[i + j * m];
            }
        }
        fortran_to_c_f(L_f, (float*)result->L->data, m, k);
        free(L_f);

        /* Extract U */
        float* U_f = (float*)calloc(k * n, sizeof(float));
        for (int32_t j = 0; j < n; j++) {
            for (int32_t i = 0; i <= MIN(j, k - 1); i++) {
                U_f[i + j * k] = A_f[i + j * m];
            }
        }
        fortran_to_c_f(U_f, (float*)result->U->data, k, n);
        free(U_f);

        /* Build permutation matrix (always double for consistency) */
        double* P_data = (double*)result->P->data;
        memset(P_data, 0, m * m * sizeof(double));

        int32_t* perm = (int32_t*)malloc(m * sizeof(int32_t));
        for (int32_t i = 0; i < m; i++) perm[i] = i;

        for (int32_t i = 0; i < k; i++) {
            int32_t swap_row = ipiv[i] - 1;
            if (swap_row != i) {
                int32_t tmp = perm[i];
                perm[i] = perm[swap_row];
                perm[swap_row] = tmp;
            }
        }

        for (int32_t i = 0; i < m; i++) {
            P_data[i * m + perm[i]] = 1.0;
        }

        free(perm);
        free(A_f);
        free(ipiv);
    } else {
        ndarray_free(result->P);
        ndarray_free(result->L);
        ndarray_free(result->U);
        free(result);
        return NULL;
    }

    return result;
}

EXPORT void linalg_lu_free(LUResult* result) {
    if (result) {
        if (result->P) ndarray_free(result->P);
        if (result->L) ndarray_free(result->L);
        if (result->U) ndarray_free(result->U);
        free(result);
    }
}

EXPORT NDArray* linalg_lu_get_p(LUResult* result) {
    return result ? result->P : NULL;
}

EXPORT NDArray* linalg_lu_get_l(LUResult* result) {
    return result ? result->L : NULL;
}

EXPORT NDArray* linalg_lu_get_u(LUResult* result) {
    return result ? result->U : NULL;
}

/* ============ Matrix Rank ============ */

/**
 * Compute numerical rank of a matrix using SVD.
 *
 * @param a    Input matrix (M x N)
 * @param tol  Tolerance for singular values. Values <= tol * max(s) are considered zero.
 *             Use negative value for default: max(M,N) * eps * max(s)
 * @return     Numerical rank, or -1 on error
 */
EXPORT int32_t linalg_matrix_rank(const NDArray* a, double tol) {
    if (!a) return -1;
    if (a->ndim != 2) return -1;

    int32_t m = a->shape[0];
    int32_t n = a->shape[1];
    int32_t k = MIN(m, n);

    if (k == 0) return 0;

    /* Compute SVD to get singular values only (using jobz='N' would be more efficient,
       but we reuse the existing SVD implementation) */
    if (a->dtype == DTYPE_FLOAT64) {
        double* A_f = (double*)malloc(m * n * sizeof(double));
        double* s = (double*)malloc(k * sizeof(double));
        int32_t lwork = 4 * k * k + 7 * k + MAX(m, n);
        double* work = (double*)malloc(lwork * sizeof(double));
        int32_t* iwork = (int32_t*)malloc(8 * k * sizeof(int32_t));

        /* Minimal U and VT for 'N' mode - not needed but LAPACK may require valid pointers */
        double* U_f = (double*)malloc(1 * sizeof(double));
        double* VT_f = (double*)malloc(1 * sizeof(double));

        if (!A_f || !s || !work || !iwork || !U_f || !VT_f) {
            free(A_f); free(s); free(work); free(iwork); free(U_f); free(VT_f);
            return -1;
        }

        c_to_fortran_d((const double*)a->data, A_f, m, n);

        /* Use jobz='N' - no U or VT computed, just singular values */
        int32_t info = lapack_dgesdd('N', m, n, A_f, m, s, U_f, 1, VT_f, 1,
                                      work, lwork, iwork);

        if (info != 0) {
            free(A_f); free(s); free(work); free(iwork); free(U_f); free(VT_f);
            return -1;
        }

        /* Determine tolerance */
        double max_s = s[0];  /* Singular values are in descending order */
        double eps = lapack_dlamch('E');  /* Machine epsilon */
        double actual_tol = (tol < 0) ? MAX(m, n) * eps * max_s : tol * max_s;

        /* Count singular values above tolerance */
        int32_t rank = 0;
        for (int32_t i = 0; i < k; i++) {
            if (s[i] > actual_tol) {
                rank++;
            }
        }

        free(A_f); free(s); free(work); free(iwork); free(U_f); free(VT_f);
        return rank;
    } else if (a->dtype == DTYPE_FLOAT32) {
        float* A_f = (float*)malloc(m * n * sizeof(float));
        float* s = (float*)malloc(k * sizeof(float));
        int32_t lwork = 4 * k * k + 7 * k + MAX(m, n);
        float* work = (float*)malloc(lwork * sizeof(float));
        int32_t* iwork = (int32_t*)malloc(8 * k * sizeof(int32_t));
        float* U_f = (float*)malloc(1 * sizeof(float));
        float* VT_f = (float*)malloc(1 * sizeof(float));

        if (!A_f || !s || !work || !iwork || !U_f || !VT_f) {
            free(A_f); free(s); free(work); free(iwork); free(U_f); free(VT_f);
            return -1;
        }

        c_to_fortran_f((const float*)a->data, A_f, m, n);

        int32_t info = lapack_sgesdd('N', m, n, A_f, m, s, U_f, 1, VT_f, 1,
                                      work, lwork, iwork);

        if (info != 0) {
            free(A_f); free(s); free(work); free(iwork); free(U_f); free(VT_f);
            return -1;
        }

        float max_s = s[0];
        float eps = (float)lapack_dlamch('E');
        float actual_tol = (tol < 0) ? MAX(m, n) * eps * max_s : (float)tol * max_s;

        int32_t rank = 0;
        for (int32_t i = 0; i < k; i++) {
            if (s[i] > actual_tol) {
                rank++;
            }
        }

        free(A_f); free(s); free(work); free(iwork); free(U_f); free(VT_f);
        return rank;
    }

    return -1;
}

/* ============ Pseudo-Inverse (pinv) ============ */

/**
 * Compute Moore-Penrose pseudo-inverse using SVD.
 * pinv(A) = V @ diag(1/s) @ U^T, with small singular values zeroed.
 *
 * @param a      Input matrix (M x N)
 * @param rcond  Relative condition number. Singular values s[i] < rcond * max(s)
 *               are treated as zero. Use negative value for default.
 * @return       Pseudo-inverse (N x M), or NULL on error
 */
EXPORT NDArray* linalg_pinv(const NDArray* a, double rcond) {
    if (!a) return NULL;
    if (a->ndim != 2) return NULL;

    int32_t m = a->shape[0];
    int32_t n = a->shape[1];
    int32_t k = MIN(m, n);

    if (m == 0 || n == 0) {
        /* Return empty (n, m) array */
        int32_t result_shape[2] = {n, m};
        return ndarray_empty(2, result_shape, a->dtype);
    }

    /* Create result array (N x M) */
    int32_t result_shape[2] = {n, m};
    NDArray* result = ndarray_empty(2, result_shape, a->dtype);
    if (!result) return NULL;

    if (a->dtype == DTYPE_FLOAT64) {
        /* Compute full SVD: A = U @ diag(s) @ Vh */
        double* A_f = (double*)malloc(m * n * sizeof(double));
        double* U_f = (double*)malloc(m * k * sizeof(double));
        double* VT_f = (double*)malloc(k * n * sizeof(double));
        double* s = (double*)malloc(k * sizeof(double));
        int32_t lwork = 4 * k * k + 7 * k + MAX(m, n);
        double* work = (double*)malloc(lwork * sizeof(double));
        int32_t* iwork = (int32_t*)malloc(8 * k * sizeof(int32_t));

        if (!A_f || !U_f || !VT_f || !s || !work || !iwork) {
            free(A_f); free(U_f); free(VT_f); free(s); free(work); free(iwork);
            ndarray_free(result);
            return NULL;
        }

        c_to_fortran_d((const double*)a->data, A_f, m, n);

        /* Use jobz='S' for reduced SVD */
        int32_t info = lapack_dgesdd('S', m, n, A_f, m, s, U_f, m, VT_f, k,
                                      work, lwork, iwork);

        if (info != 0) {
            free(A_f); free(U_f); free(VT_f); free(s); free(work); free(iwork);
            ndarray_free(result);
            return NULL;
        }

        /* Determine cutoff for singular values */
        double eps = lapack_dlamch('E');
        double cutoff = (rcond < 0) ? MAX(m, n) * eps * s[0] : rcond * s[0];

        /* Compute pinv = V @ diag(1/s) @ U^T = VT^T @ diag(1/s) @ U^T
         * First: scale columns of U by 1/s[i] (or 0 if s[i] <= cutoff) */
        for (int32_t i = 0; i < k; i++) {
            double scale = (s[i] > cutoff) ? 1.0 / s[i] : 0.0;
            for (int32_t j = 0; j < m; j++) {
                U_f[j + i * m] *= scale;  /* Scale column i of U */
            }
        }

        /* Now compute pinv = VT^T @ U_scaled^T = V @ U_scaled^T
         * Result is (n x m): pinv[i,j] = sum_l V[i,l] * U_scaled[j,l]
         *                              = sum_l VT[l,i] * U_scaled[j,l]
         * In Fortran order: pinv_f[i + j*n] = sum_l VT_f[l + i*k] * U_f[j + l*m]
         */
        double* pinv_f = (double*)calloc(n * m, sizeof(double));
        if (!pinv_f) {
            free(A_f); free(U_f); free(VT_f); free(s); free(work); free(iwork);
            ndarray_free(result);
            return NULL;
        }

        /* pinv = V @ U_scaled^T where V = VT^T
         * Using BLAS: pinv = alpha * VT^T @ U^T + beta * pinv
         * VT is k x n, VT^T is n x k
         * U_scaled is m x k, U_scaled^T is k x m
         * Result: n x m
         */
        blas_dgemm('T', 'T', n, m, k, 1.0, VT_f, k, U_f, m, 0.0, pinv_f, n);

        fortran_to_c_d(pinv_f, (double*)result->data, n, m);

        free(pinv_f);
        free(A_f); free(U_f); free(VT_f); free(s); free(work); free(iwork);
    } else if (a->dtype == DTYPE_FLOAT32) {
        float* A_f = (float*)malloc(m * n * sizeof(float));
        float* U_f = (float*)malloc(m * k * sizeof(float));
        float* VT_f = (float*)malloc(k * n * sizeof(float));
        float* s = (float*)malloc(k * sizeof(float));
        int32_t lwork = 4 * k * k + 7 * k + MAX(m, n);
        float* work = (float*)malloc(lwork * sizeof(float));
        int32_t* iwork = (int32_t*)malloc(8 * k * sizeof(int32_t));

        if (!A_f || !U_f || !VT_f || !s || !work || !iwork) {
            free(A_f); free(U_f); free(VT_f); free(s); free(work); free(iwork);
            ndarray_free(result);
            return NULL;
        }

        c_to_fortran_f((const float*)a->data, A_f, m, n);

        int32_t info = lapack_sgesdd('S', m, n, A_f, m, s, U_f, m, VT_f, k,
                                      work, lwork, iwork);

        if (info != 0) {
            free(A_f); free(U_f); free(VT_f); free(s); free(work); free(iwork);
            ndarray_free(result);
            return NULL;
        }

        float eps = (float)lapack_dlamch('E');
        float cutoff = (rcond < 0) ? MAX(m, n) * eps * s[0] : (float)rcond * s[0];

        for (int32_t i = 0; i < k; i++) {
            float scale = (s[i] > cutoff) ? 1.0f / s[i] : 0.0f;
            for (int32_t j = 0; j < m; j++) {
                U_f[j + i * m] *= scale;
            }
        }

        float* pinv_f = (float*)calloc(n * m, sizeof(float));
        if (!pinv_f) {
            free(A_f); free(U_f); free(VT_f); free(s); free(work); free(iwork);
            ndarray_free(result);
            return NULL;
        }

        blas_sgemm('T', 'T', n, m, k, 1.0f, VT_f, k, U_f, m, 0.0f, pinv_f, n);

        fortran_to_c_f(pinv_f, (float*)result->data, n, m);

        free(pinv_f);
        free(A_f); free(U_f); free(VT_f); free(s); free(work); free(iwork);
    } else {
        ndarray_free(result);
        return NULL;
    }

    return result;
}

/* ============ Condition Number ============ */

/**
 * Compute condition number of a matrix.
 *
 * @param a    Input matrix (M x N)
 * @param p    Norm type:
 *             2: 2-norm (default, ratio of largest to smallest singular value)
 *            -2: ratio of smallest to largest singular value
 *             Other values not yet supported (would need dgecon)
 * @return     Condition number, or -1.0 on error, or INFINITY if singular
 */
EXPORT double linalg_cond(const NDArray* a, int32_t p) {
    if (!a) return -1.0;
    if (a->ndim != 2) return -1.0;

    int32_t m = a->shape[0];
    int32_t n = a->shape[1];
    int32_t k = MIN(m, n);

    if (k == 0) return 0.0;

    /* For 2-norm and -2-norm, use SVD */
    if (p == 2 || p == -2) {
        if (a->dtype == DTYPE_FLOAT64) {
            double* A_f = (double*)malloc(m * n * sizeof(double));
            double* s = (double*)malloc(k * sizeof(double));
            int32_t lwork = 4 * k * k + 7 * k + MAX(m, n);
            double* work = (double*)malloc(lwork * sizeof(double));
            int32_t* iwork = (int32_t*)malloc(8 * k * sizeof(int32_t));
            double* U_f = (double*)malloc(1 * sizeof(double));
            double* VT_f = (double*)malloc(1 * sizeof(double));

            if (!A_f || !s || !work || !iwork || !U_f || !VT_f) {
                free(A_f); free(s); free(work); free(iwork); free(U_f); free(VT_f);
                return -1.0;
            }

            c_to_fortran_d((const double*)a->data, A_f, m, n);

            int32_t info = lapack_dgesdd('N', m, n, A_f, m, s, U_f, 1, VT_f, 1,
                                          work, lwork, iwork);

            if (info != 0) {
                free(A_f); free(s); free(work); free(iwork); free(U_f); free(VT_f);
                return -1.0;
            }

            double cond;
            double max_s = s[0];
            double min_s = s[k - 1];

            if (min_s == 0.0) {
                cond = INFINITY;
            } else if (p == 2) {
                cond = max_s / min_s;
            } else {  /* p == -2 */
                cond = min_s / max_s;
            }

            free(A_f); free(s); free(work); free(iwork); free(U_f); free(VT_f);
            return cond;
        } else if (a->dtype == DTYPE_FLOAT32) {
            float* A_f = (float*)malloc(m * n * sizeof(float));
            float* s = (float*)malloc(k * sizeof(float));
            int32_t lwork = 4 * k * k + 7 * k + MAX(m, n);
            float* work = (float*)malloc(lwork * sizeof(float));
            int32_t* iwork = (int32_t*)malloc(8 * k * sizeof(int32_t));
            float* U_f = (float*)malloc(1 * sizeof(float));
            float* VT_f = (float*)malloc(1 * sizeof(float));

            if (!A_f || !s || !work || !iwork || !U_f || !VT_f) {
                free(A_f); free(s); free(work); free(iwork); free(U_f); free(VT_f);
                return -1.0;
            }

            c_to_fortran_f((const float*)a->data, A_f, m, n);

            int32_t info = lapack_sgesdd('N', m, n, A_f, m, s, U_f, 1, VT_f, 1,
                                          work, lwork, iwork);

            if (info != 0) {
                free(A_f); free(s); free(work); free(iwork); free(U_f); free(VT_f);
                return -1.0;
            }

            double cond;
            float max_s = s[0];
            float min_s = s[k - 1];

            if (min_s == 0.0f) {
                cond = INFINITY;
            } else if (p == 2) {
                cond = (double)max_s / (double)min_s;
            } else {
                cond = (double)min_s / (double)max_s;
            }

            free(A_f); free(s); free(work); free(iwork); free(U_f); free(VT_f);
            return cond;
        }
    }

    /* Other norms not yet supported */
    return -1.0;
}

/* ============ Least Squares ============ */

/**
 * Least squares result structure.
 */
typedef struct {
    NDArray* x;         /* Solution (N,) or (N x NRHS) */
    NDArray* residuals; /* Residuals (NRHS,) for overdetermined, empty otherwise */
    int32_t rank;       /* Effective rank of matrix A */
    NDArray* s;         /* Singular values (min(M,N),) */
} LstsqResult;

/**
 * Solve least squares problem: minimize ||b - A @ x||_2.
 *
 * @param a      Matrix A (M x N)
 * @param b      Right-hand side (M,) or (M x NRHS)
 * @param rcond  Cutoff for small singular values. Values s[i] < rcond * max(s)
 *               are treated as zero. Use negative for machine precision default.
 * @return       LstsqResult with x, residuals, rank, s, or NULL on error
 */
EXPORT LstsqResult* linalg_lstsq(const NDArray* a, const NDArray* b, double rcond) {
    if (!a || !b) return NULL;
    if (a->ndim != 2) return NULL;

    int32_t m = a->shape[0];
    int32_t n = a->shape[1];
    int32_t k = MIN(m, n);
    int32_t nrhs;
    int32_t is_vector = (b->ndim == 1);

    if (is_vector) {
        if (b->shape[0] != m) return NULL;
        nrhs = 1;
    } else if (b->ndim == 2) {
        if (b->shape[0] != m) return NULL;
        nrhs = b->shape[1];
    } else {
        return NULL;
    }

    LstsqResult* result = (LstsqResult*)malloc(sizeof(LstsqResult));
    if (!result) return NULL;

    /* Create output arrays */
    int32_t x_shape[2] = {n, nrhs};
    int32_t s_shape[1] = {k};

    if (is_vector) {
        int32_t x_shape_1d[1] = {n};
        result->x = ndarray_empty(1, x_shape_1d, a->dtype);
    } else {
        result->x = ndarray_empty(2, x_shape, a->dtype);
    }
    result->s = ndarray_empty(1, s_shape, a->dtype);
    result->rank = 0;

    /* Residuals only for overdetermined systems (m > n) with full rank */
    if (m > n) {
        int32_t res_shape[1] = {nrhs};
        result->residuals = ndarray_empty(1, res_shape, a->dtype);
    } else {
        int32_t res_shape[1] = {0};
        result->residuals = ndarray_empty(1, res_shape, a->dtype);
    }

    if (!result->x || !result->s || !result->residuals) {
        if (result->x) ndarray_free(result->x);
        if (result->s) ndarray_free(result->s);
        if (result->residuals) ndarray_free(result->residuals);
        free(result);
        return NULL;
    }

    if (a->dtype == DTYPE_FLOAT64) {
        double* A_f = (double*)malloc(m * n * sizeof(double));
        double* B_f = (double*)malloc(MAX(m, n) * nrhs * sizeof(double));
        double* s = (double*)malloc(k * sizeof(double));

        /* Query workspace size */
        double work_query;
        int32_t iwork_query;
        int32_t dummy_rank;
        lapack_dgelsd(m, n, nrhs, NULL, m, NULL, MAX(m, n), NULL, rcond, &dummy_rank,
                      &work_query, -1, &iwork_query);
        int32_t lwork = (int32_t)work_query + 100;  /* Add some padding */
        double* work = (double*)malloc(lwork * sizeof(double));
        int32_t* iwork = (int32_t*)malloc(8 * k * sizeof(int32_t));

        if (!A_f || !B_f || !s || !work || !iwork) {
            free(A_f); free(B_f); free(s); free(work); free(iwork);
            ndarray_free(result->x);
            ndarray_free(result->s);
            ndarray_free(result->residuals);
            free(result);
            return NULL;
        }

        /* Convert A to Fortran order */
        c_to_fortran_d((const double*)a->data, A_f, m, n);

        /* Copy b to B (extended to max(m,n) rows for solution) */
        memset(B_f, 0, MAX(m, n) * nrhs * sizeof(double));
        if (is_vector) {
            memcpy(B_f, b->data, m * sizeof(double));
        } else {
            c_to_fortran_d((const double*)b->data, B_f, m, nrhs);
        }

        /* Store original b for residual computation */
        double* b_orig = NULL;
        if (m > n) {
            b_orig = (double*)malloc(m * nrhs * sizeof(double));
            if (b_orig) {
                memcpy(b_orig, B_f, m * nrhs * sizeof(double));
            }
        }

        /* Solve least squares */
        int32_t info = lapack_dgelsd(m, n, nrhs, A_f, m, B_f, MAX(m, n), s, rcond,
                                      &result->rank, work, lwork, iwork);

        if (info != 0) {
            free(A_f); free(B_f); free(s); free(work); free(iwork);
            if (b_orig) free(b_orig);
            ndarray_free(result->x);
            ndarray_free(result->s);
            ndarray_free(result->residuals);
            free(result);
            return NULL;
        }

        /* Copy solution (first n rows of B) */
        if (is_vector) {
            memcpy(result->x->data, B_f, n * sizeof(double));
        } else {
            /* B_f has solution in first n rows (Fortran order) */
            double* x_temp = (double*)malloc(n * nrhs * sizeof(double));
            for (int32_t j = 0; j < nrhs; j++) {
                for (int32_t i = 0; i < n; i++) {
                    x_temp[i + j * n] = B_f[i + j * MAX(m, n)];
                }
            }
            fortran_to_c_d(x_temp, (double*)result->x->data, n, nrhs);
            free(x_temp);
        }

        /* Copy singular values */
        memcpy(result->s->data, s, k * sizeof(double));

        /* Compute residuals for overdetermined systems with full rank */
        if (m > n && result->rank == n && b_orig) {
            /* residuals[j] = ||b[:,j] - A @ x[:,j]||^2 */
            /* Using the fact that dgelsd leaves residuals in rows n+1 to m of B */
            double* res_data = (double*)result->residuals->data;
            for (int32_t j = 0; j < nrhs; j++) {
                double sum_sq = 0.0;
                for (int32_t i = n; i < m; i++) {
                    double r = B_f[i + j * MAX(m, n)];
                    sum_sq += r * r;
                }
                res_data[j] = sum_sq;
            }
        }

        free(A_f); free(B_f); free(s); free(work); free(iwork);
        if (b_orig) free(b_orig);
    } else if (a->dtype == DTYPE_FLOAT32) {
        /* Float version - convert to double, solve, convert back */
        /* TODO: implement sgelsd for native float support */
        ndarray_free(result->x);
        ndarray_free(result->s);
        ndarray_free(result->residuals);
        free(result);
        return NULL;  /* Float not yet supported */
    } else {
        ndarray_free(result->x);
        ndarray_free(result->s);
        ndarray_free(result->residuals);
        free(result);
        return NULL;
    }

    return result;
}

EXPORT void linalg_lstsq_free(LstsqResult* result) {
    if (result) {
        if (result->x) ndarray_free(result->x);
        if (result->residuals) ndarray_free(result->residuals);
        if (result->s) ndarray_free(result->s);
        free(result);
    }
}

EXPORT NDArray* linalg_lstsq_get_x(LstsqResult* result) {
    return result ? result->x : NULL;
}

EXPORT NDArray* linalg_lstsq_get_residuals(LstsqResult* result) {
    return result ? result->residuals : NULL;
}

EXPORT int32_t linalg_lstsq_get_rank(LstsqResult* result) {
    return result ? result->rank : -1;
}

EXPORT NDArray* linalg_lstsq_get_s(LstsqResult* result) {
    return result ? result->s : NULL;
}

/* ============ Schur Decomposition ============ */

/**
 * Schur decomposition result structure.
 */
typedef struct {
    NDArray* T;  /* Schur form (upper quasi-triangular, N x N) */
    NDArray* Z;  /* Orthogonal matrix (N x N) */
} SchurResult;

/**
 * Compute Schur decomposition: A = Z @ T @ Z^T
 * where T is upper quasi-triangular (real Schur form) and Z is orthogonal.
 *
 * Real eigenvalues appear as 1x1 blocks on the diagonal of T.
 * Complex conjugate pairs appear as 2x2 blocks.
 *
 * @param a  Square matrix A (N x N)
 * @return   SchurResult with T and Z, or NULL on error
 */
EXPORT SchurResult* linalg_schur(const NDArray* a) {
    if (!a) return NULL;
    if (a->ndim != 2) return NULL;
    if (a->shape[0] != a->shape[1]) return NULL;

    int32_t n = a->shape[0];

    SchurResult* result = (SchurResult*)malloc(sizeof(SchurResult));
    if (!result) return NULL;

    int32_t shape[2] = {n, n};
    result->T = ndarray_empty(2, shape, a->dtype);
    result->Z = ndarray_empty(2, shape, a->dtype);

    if (!result->T || !result->Z) {
        if (result->T) ndarray_free(result->T);
        if (result->Z) ndarray_free(result->Z);
        free(result);
        return NULL;
    }

    if (a->dtype == DTYPE_FLOAT64) {
        double* A_f = (double*)malloc(n * n * sizeof(double));
        double* VS_f = (double*)malloc(n * n * sizeof(double));
        double* wr = (double*)malloc(n * sizeof(double));
        double* wi = (double*)malloc(n * sizeof(double));

        /* Query workspace */
        double work_query;
        lapack_dgees('V', 'N', n, NULL, n, NULL, NULL, NULL, NULL, n, &work_query, -1);
        int32_t lwork = (int32_t)work_query + 100;
        double* work = (double*)malloc(lwork * sizeof(double));

        if (!A_f || !VS_f || !wr || !wi || !work) {
            free(A_f); free(VS_f); free(wr); free(wi); free(work);
            ndarray_free(result->T);
            ndarray_free(result->Z);
            free(result);
            return NULL;
        }

        c_to_fortran_d((const double*)a->data, A_f, n, n);

        int32_t sdim;
        int32_t info = lapack_dgees('V', 'N', n, A_f, n, &sdim, wr, wi, VS_f, n,
                                     work, lwork);

        if (info != 0) {
            free(A_f); free(VS_f); free(wr); free(wi); free(work);
            ndarray_free(result->T);
            ndarray_free(result->Z);
            free(result);
            return NULL;
        }

        /* A_f now contains Schur form T, VS_f contains Z */
        fortran_to_c_d(A_f, (double*)result->T->data, n, n);
        fortran_to_c_d(VS_f, (double*)result->Z->data, n, n);

        free(A_f); free(VS_f); free(wr); free(wi); free(work);
    } else if (a->dtype == DTYPE_FLOAT32) {
        /* Float version not yet implemented */
        ndarray_free(result->T);
        ndarray_free(result->Z);
        free(result);
        return NULL;
    } else {
        ndarray_free(result->T);
        ndarray_free(result->Z);
        free(result);
        return NULL;
    }

    return result;
}

EXPORT void linalg_schur_free(SchurResult* result) {
    if (result) {
        if (result->T) ndarray_free(result->T);
        if (result->Z) ndarray_free(result->Z);
        free(result);
    }
}

EXPORT NDArray* linalg_schur_get_t(SchurResult* result) {
    return result ? result->T : NULL;
}

EXPORT NDArray* linalg_schur_get_z(SchurResult* result) {
    return result ? result->Z : NULL;
}

/* ============ LDL Decomposition ============ */

/**
 * LDL decomposition result structure.
 */
typedef struct {
    NDArray* L;     /* Lower triangular with unit diagonal (N x N) */
    NDArray* D;     /* Block diagonal as 1D array (N,) - diagonal elements */
    NDArray* perm;  /* Permutation indices (N,) */
} LDLResult;

/**
 * Compute LDL decomposition for symmetric matrix.
 * A = P^T @ L @ D @ L^T @ P
 * where L is lower triangular with unit diagonal, D is diagonal,
 * and P is a permutation matrix.
 *
 * @param a      Symmetric matrix A (N x N)
 * @param lower  Non-zero to use lower triangle (default), 0 for upper
 * @return       LDLResult with L, D, perm, or NULL on error
 */
EXPORT LDLResult* linalg_ldl(const NDArray* a, int32_t lower) {
    if (!a) return NULL;
    if (a->ndim != 2) return NULL;
    if (a->shape[0] != a->shape[1]) return NULL;

    int32_t n = a->shape[0];

    LDLResult* result = (LDLResult*)malloc(sizeof(LDLResult));
    if (!result) return NULL;

    int32_t l_shape[2] = {n, n};
    int32_t d_shape[1] = {n};
    int32_t perm_shape[1] = {n};

    result->L = ndarray_empty(2, l_shape, a->dtype);
    result->D = ndarray_empty(1, d_shape, a->dtype);
    result->perm = ndarray_empty(1, perm_shape, DTYPE_INT32);

    if (!result->L || !result->D || !result->perm) {
        if (result->L) ndarray_free(result->L);
        if (result->D) ndarray_free(result->D);
        if (result->perm) ndarray_free(result->perm);
        free(result);
        return NULL;
    }

    if (a->dtype == DTYPE_FLOAT64) {
        double* A_f = (double*)malloc(n * n * sizeof(double));
        int32_t* ipiv = (int32_t*)malloc(n * sizeof(int32_t));

        /* Query workspace */
        double work_query;
        lapack_dsytrf(lower ? 'L' : 'U', n, NULL, n, NULL, &work_query, -1);
        int32_t lwork = (int32_t)work_query + 100;
        double* work = (double*)malloc(lwork * sizeof(double));

        if (!A_f || !ipiv || !work) {
            free(A_f); free(ipiv); free(work);
            ndarray_free(result->L);
            ndarray_free(result->D);
            ndarray_free(result->perm);
            free(result);
            return NULL;
        }

        c_to_fortran_d((const double*)a->data, A_f, n, n);

        int32_t info = lapack_dsytrf(lower ? 'L' : 'U', n, A_f, n, ipiv, work, lwork);

        if (info < 0) {
            free(A_f); free(ipiv); free(work);
            ndarray_free(result->L);
            ndarray_free(result->D);
            ndarray_free(result->perm);
            free(result);
            return NULL;
        }

        /* Extract L, D, and permutation from the factored matrix */
        double* L_data = (double*)result->L->data;
        double* D_data = (double*)result->D->data;
        int32_t* perm_data = (int32_t*)result->perm->data;

        /* Initialize L as identity */
        memset(L_data, 0, n * n * sizeof(double));
        for (int32_t i = 0; i < n; i++) {
            L_data[i * n + i] = 1.0;
        }

        /* Extract L and D from factored A_f */
        if (lower) {
            /* L is stored below diagonal, D on diagonal */
            for (int32_t j = 0; j < n; j++) {
                /* Diagonal element goes to D */
                D_data[j] = A_f[j + j * n];

                /* Below diagonal goes to L (in C-order: L[i,j] = L_data[i*n + j]) */
                for (int32_t i = j + 1; i < n; i++) {
                    L_data[i * n + j] = A_f[i + j * n];
                }
            }
        } else {
            /* Upper triangle stored - transpose to get L */
            for (int32_t j = 0; j < n; j++) {
                D_data[j] = A_f[j + j * n];
                for (int32_t i = 0; i < j; i++) {
                    L_data[j * n + i] = A_f[i + j * n];
                }
            }
        }

        /* Convert pivot indices to 0-based permutation */
        for (int32_t i = 0; i < n; i++) {
            perm_data[i] = (ipiv[i] > 0) ? ipiv[i] - 1 : -ipiv[i] - 1;
        }

        free(A_f); free(ipiv); free(work);
    } else if (a->dtype == DTYPE_FLOAT32) {
        /* Float version not yet implemented */
        ndarray_free(result->L);
        ndarray_free(result->D);
        ndarray_free(result->perm);
        free(result);
        return NULL;
    } else {
        ndarray_free(result->L);
        ndarray_free(result->D);
        ndarray_free(result->perm);
        free(result);
        return NULL;
    }

    return result;
}

EXPORT void linalg_ldl_free(LDLResult* result) {
    if (result) {
        if (result->L) ndarray_free(result->L);
        if (result->D) ndarray_free(result->D);
        if (result->perm) ndarray_free(result->perm);
        free(result);
    }
}

EXPORT NDArray* linalg_ldl_get_l(LDLResult* result) {
    return result ? result->L : NULL;
}

EXPORT NDArray* linalg_ldl_get_d(LDLResult* result) {
    return result ? result->D : NULL;
}

EXPORT NDArray* linalg_ldl_get_perm(LDLResult* result) {
    return result ? result->perm : NULL;
}
