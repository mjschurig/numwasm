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
        int32_t lwork = 3 * k * k + MAX(MAX(m, n), 4 * k * k + 4 * k);
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
        int32_t lwork = 3 * k * k + MAX(MAX(m, n), 4 * k * k + 4 * k);
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
