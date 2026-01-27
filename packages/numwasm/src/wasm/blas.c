/*
 * NumJS BLAS - Basic Linear Algebra Subprograms Implementation
 *
 * A minimal but correct implementation of essential BLAS routines
 * optimized for WebAssembly.
 *
 * Memory layout: Column-major (Fortran order) for LAPACK compatibility.
 */

#include "blas.h"
#include <math.h>
#include <string.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/* Helper macros */
#define MIN(a, b) ((a) < (b) ? (a) : (b))
#define MAX(a, b) ((a) > (b) ? (a) : (b))

/* Character comparison helpers */
static inline int is_trans(char c) {
    return c == 'T' || c == 't' || c == 'C' || c == 'c';
}

static inline int is_upper(char c) {
    return c == 'U' || c == 'u';
}

static inline int is_unit(char c) {
    return c == 'U' || c == 'u';
}

static inline int is_left(char c) {
    return c == 'L' || c == 'l';
}

/* ============ Level 1 BLAS Implementation ============ */

EXPORT double blas_ddot(int32_t n, const double* x, int32_t incx,
                         const double* y, int32_t incy)
{
    if (n <= 0) return 0.0;

    double result = 0.0;

    if (incx == 1 && incy == 1) {
        /* Optimized path: contiguous access with loop unrolling */
        int32_t m = n % 5;
        for (int32_t i = 0; i < m; i++) {
            result += x[i] * y[i];
        }
        for (int32_t i = m; i < n; i += 5) {
            result += x[i] * y[i]
                    + x[i+1] * y[i+1]
                    + x[i+2] * y[i+2]
                    + x[i+3] * y[i+3]
                    + x[i+4] * y[i+4];
        }
    } else {
        /* Strided access */
        int32_t ix = (incx > 0) ? 0 : (1 - n) * incx;
        int32_t iy = (incy > 0) ? 0 : (1 - n) * incy;
        for (int32_t i = 0; i < n; i++) {
            result += x[ix] * y[iy];
            ix += incx;
            iy += incy;
        }
    }

    return result;
}

EXPORT float blas_sdot(int32_t n, const float* x, int32_t incx,
                        const float* y, int32_t incy)
{
    if (n <= 0) return 0.0f;

    float result = 0.0f;

    if (incx == 1 && incy == 1) {
        int32_t m = n % 5;
        for (int32_t i = 0; i < m; i++) {
            result += x[i] * y[i];
        }
        for (int32_t i = m; i < n; i += 5) {
            result += x[i] * y[i]
                    + x[i+1] * y[i+1]
                    + x[i+2] * y[i+2]
                    + x[i+3] * y[i+3]
                    + x[i+4] * y[i+4];
        }
    } else {
        int32_t ix = (incx > 0) ? 0 : (1 - n) * incx;
        int32_t iy = (incy > 0) ? 0 : (1 - n) * incy;
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
    if (n <= 0 || incx <= 0) return 0.0;

    /*
     * Use scaled computation to avoid overflow/underflow.
     * Algorithm: compute scale and sum of squares separately.
     * result = scale * sqrt(ssq) where ssq = sum((x[i]/scale)^2)
     */
    double scale = 0.0;
    double ssq = 1.0;

    int32_t ix = 0;
    for (int32_t i = 0; i < n; i++) {
        double absxi = fabs(x[ix]);
        if (absxi > 0.0) {
            if (scale < absxi) {
                /* Re-scale: ssq = 1 + ssq * (scale/absxi)^2 */
                double ratio = scale / absxi;
                ssq = 1.0 + ssq * ratio * ratio;
                scale = absxi;
            } else {
                /* Accumulate: ssq += (absxi/scale)^2 */
                double ratio = absxi / scale;
                ssq += ratio * ratio;
            }
        }
        ix += incx;
    }

    return scale * sqrt(ssq);
}

EXPORT float blas_snrm2(int32_t n, const float* x, int32_t incx)
{
    if (n <= 0 || incx <= 0) return 0.0f;

    float scale = 0.0f;
    float ssq = 1.0f;

    int32_t ix = 0;
    for (int32_t i = 0; i < n; i++) {
        float absxi = fabsf(x[ix]);
        if (absxi > 0.0f) {
            if (scale < absxi) {
                float ratio = scale / absxi;
                ssq = 1.0f + ssq * ratio * ratio;
                scale = absxi;
            } else {
                float ratio = absxi / scale;
                ssq += ratio * ratio;
            }
        }
        ix += incx;
    }

    return scale * sqrtf(ssq);
}

EXPORT double blas_dasum(int32_t n, const double* x, int32_t incx)
{
    if (n <= 0 || incx <= 0) return 0.0;

    double result = 0.0;

    if (incx == 1) {
        /* Unrolled loop for contiguous access */
        int32_t m = n % 6;
        for (int32_t i = 0; i < m; i++) {
            result += fabs(x[i]);
        }
        for (int32_t i = m; i < n; i += 6) {
            result += fabs(x[i]) + fabs(x[i+1]) + fabs(x[i+2])
                    + fabs(x[i+3]) + fabs(x[i+4]) + fabs(x[i+5]);
        }
    } else {
        int32_t ix = 0;
        for (int32_t i = 0; i < n; i++) {
            result += fabs(x[ix]);
            ix += incx;
        }
    }

    return result;
}

EXPORT float blas_sasum(int32_t n, const float* x, int32_t incx)
{
    if (n <= 0 || incx <= 0) return 0.0f;

    float result = 0.0f;

    if (incx == 1) {
        int32_t m = n % 6;
        for (int32_t i = 0; i < m; i++) {
            result += fabsf(x[i]);
        }
        for (int32_t i = m; i < n; i += 6) {
            result += fabsf(x[i]) + fabsf(x[i+1]) + fabsf(x[i+2])
                    + fabsf(x[i+3]) + fabsf(x[i+4]) + fabsf(x[i+5]);
        }
    } else {
        int32_t ix = 0;
        for (int32_t i = 0; i < n; i++) {
            result += fabsf(x[ix]);
            ix += incx;
        }
    }

    return result;
}

EXPORT int32_t blas_idamax(int32_t n, const double* x, int32_t incx)
{
    if (n <= 0 || incx <= 0) return -1;
    if (n == 1) return 0;

    int32_t result = 0;
    double dmax = fabs(x[0]);

    if (incx == 1) {
        for (int32_t i = 1; i < n; i++) {
            double absxi = fabs(x[i]);
            if (absxi > dmax) {
                result = i;
                dmax = absxi;
            }
        }
    } else {
        int32_t ix = incx;
        for (int32_t i = 1; i < n; i++) {
            double absxi = fabs(x[ix]);
            if (absxi > dmax) {
                result = i;
                dmax = absxi;
            }
            ix += incx;
        }
    }

    return result;
}

EXPORT int32_t blas_isamax(int32_t n, const float* x, int32_t incx)
{
    if (n <= 0 || incx <= 0) return -1;
    if (n == 1) return 0;

    int32_t result = 0;
    float smax = fabsf(x[0]);

    if (incx == 1) {
        for (int32_t i = 1; i < n; i++) {
            float absxi = fabsf(x[i]);
            if (absxi > smax) {
                result = i;
                smax = absxi;
            }
        }
    } else {
        int32_t ix = incx;
        for (int32_t i = 1; i < n; i++) {
            float absxi = fabsf(x[ix]);
            if (absxi > smax) {
                result = i;
                smax = absxi;
            }
            ix += incx;
        }
    }

    return result;
}

EXPORT void blas_dscal(int32_t n, double alpha, double* x, int32_t incx)
{
    if (n <= 0 || incx <= 0) return;

    if (incx == 1) {
        /* Unrolled loop */
        int32_t m = n % 5;
        for (int32_t i = 0; i < m; i++) {
            x[i] *= alpha;
        }
        for (int32_t i = m; i < n; i += 5) {
            x[i] *= alpha;
            x[i+1] *= alpha;
            x[i+2] *= alpha;
            x[i+3] *= alpha;
            x[i+4] *= alpha;
        }
    } else {
        int32_t ix = 0;
        for (int32_t i = 0; i < n; i++) {
            x[ix] *= alpha;
            ix += incx;
        }
    }
}

EXPORT void blas_sscal(int32_t n, float alpha, float* x, int32_t incx)
{
    if (n <= 0 || incx <= 0) return;

    if (incx == 1) {
        int32_t m = n % 5;
        for (int32_t i = 0; i < m; i++) {
            x[i] *= alpha;
        }
        for (int32_t i = m; i < n; i += 5) {
            x[i] *= alpha;
            x[i+1] *= alpha;
            x[i+2] *= alpha;
            x[i+3] *= alpha;
            x[i+4] *= alpha;
        }
    } else {
        int32_t ix = 0;
        for (int32_t i = 0; i < n; i++) {
            x[ix] *= alpha;
            ix += incx;
        }
    }
}

EXPORT void blas_daxpy(int32_t n, double alpha, const double* x, int32_t incx,
                        double* y, int32_t incy)
{
    if (n <= 0 || alpha == 0.0) return;

    if (incx == 1 && incy == 1) {
        /* Unrolled loop for contiguous access */
        int32_t m = n % 4;
        for (int32_t i = 0; i < m; i++) {
            y[i] += alpha * x[i];
        }
        for (int32_t i = m; i < n; i += 4) {
            y[i] += alpha * x[i];
            y[i+1] += alpha * x[i+1];
            y[i+2] += alpha * x[i+2];
            y[i+3] += alpha * x[i+3];
        }
    } else {
        int32_t ix = (incx > 0) ? 0 : (1 - n) * incx;
        int32_t iy = (incy > 0) ? 0 : (1 - n) * incy;
        for (int32_t i = 0; i < n; i++) {
            y[iy] += alpha * x[ix];
            ix += incx;
            iy += incy;
        }
    }
}

EXPORT void blas_saxpy(int32_t n, float alpha, const float* x, int32_t incx,
                        float* y, int32_t incy)
{
    if (n <= 0 || alpha == 0.0f) return;

    if (incx == 1 && incy == 1) {
        int32_t m = n % 4;
        for (int32_t i = 0; i < m; i++) {
            y[i] += alpha * x[i];
        }
        for (int32_t i = m; i < n; i += 4) {
            y[i] += alpha * x[i];
            y[i+1] += alpha * x[i+1];
            y[i+2] += alpha * x[i+2];
            y[i+3] += alpha * x[i+3];
        }
    } else {
        int32_t ix = (incx > 0) ? 0 : (1 - n) * incx;
        int32_t iy = (incy > 0) ? 0 : (1 - n) * incy;
        for (int32_t i = 0; i < n; i++) {
            y[iy] += alpha * x[ix];
            ix += incx;
            iy += incy;
        }
    }
}

EXPORT void blas_dcopy(int32_t n, const double* x, int32_t incx,
                        double* y, int32_t incy)
{
    if (n <= 0) return;

    if (incx == 1 && incy == 1) {
        memcpy(y, x, n * sizeof(double));
    } else {
        int32_t ix = (incx > 0) ? 0 : (1 - n) * incx;
        int32_t iy = (incy > 0) ? 0 : (1 - n) * incy;
        for (int32_t i = 0; i < n; i++) {
            y[iy] = x[ix];
            ix += incx;
            iy += incy;
        }
    }
}

EXPORT void blas_scopy(int32_t n, const float* x, int32_t incx,
                        float* y, int32_t incy)
{
    if (n <= 0) return;

    if (incx == 1 && incy == 1) {
        memcpy(y, x, n * sizeof(float));
    } else {
        int32_t ix = (incx > 0) ? 0 : (1 - n) * incx;
        int32_t iy = (incy > 0) ? 0 : (1 - n) * incy;
        for (int32_t i = 0; i < n; i++) {
            y[iy] = x[ix];
            ix += incx;
            iy += incy;
        }
    }
}

EXPORT void blas_dswap(int32_t n, double* x, int32_t incx,
                        double* y, int32_t incy)
{
    if (n <= 0) return;

    if (incx == 1 && incy == 1) {
        for (int32_t i = 0; i < n; i++) {
            double temp = x[i];
            x[i] = y[i];
            y[i] = temp;
        }
    } else {
        int32_t ix = (incx > 0) ? 0 : (1 - n) * incx;
        int32_t iy = (incy > 0) ? 0 : (1 - n) * incy;
        for (int32_t i = 0; i < n; i++) {
            double temp = x[ix];
            x[ix] = y[iy];
            y[iy] = temp;
            ix += incx;
            iy += incy;
        }
    }
}

EXPORT void blas_sswap(int32_t n, float* x, int32_t incx,
                        float* y, int32_t incy)
{
    if (n <= 0) return;

    if (incx == 1 && incy == 1) {
        for (int32_t i = 0; i < n; i++) {
            float temp = x[i];
            x[i] = y[i];
            y[i] = temp;
        }
    } else {
        int32_t ix = (incx > 0) ? 0 : (1 - n) * incx;
        int32_t iy = (incy > 0) ? 0 : (1 - n) * incy;
        for (int32_t i = 0; i < n; i++) {
            float temp = x[ix];
            x[ix] = y[iy];
            y[iy] = temp;
            ix += incx;
            iy += incy;
        }
    }
}

/* ============ Level 2 BLAS Implementation ============ */

EXPORT void blas_dgemv(char trans, int32_t m, int32_t n,
                        double alpha, const double* A, int32_t lda,
                        const double* x, int32_t incx,
                        double beta, double* y, int32_t incy)
{
    if (m <= 0 || n <= 0) return;

    int32_t lenx, leny;
    if (!is_trans(trans)) {
        /* y = alpha * A * x + beta * y */
        lenx = n;
        leny = m;
    } else {
        /* y = alpha * A^T * x + beta * y */
        lenx = m;
        leny = n;
    }

    /* Scale y by beta */
    if (beta != 1.0) {
        if (incy == 1) {
            if (beta == 0.0) {
                memset(y, 0, leny * sizeof(double));
            } else {
                for (int32_t i = 0; i < leny; i++) {
                    y[i] *= beta;
                }
            }
        } else {
            int32_t iy = (incy > 0) ? 0 : (1 - leny) * incy;
            if (beta == 0.0) {
                for (int32_t i = 0; i < leny; i++) {
                    y[iy] = 0.0;
                    iy += incy;
                }
            } else {
                for (int32_t i = 0; i < leny; i++) {
                    y[iy] *= beta;
                    iy += incy;
                }
            }
        }
    }

    if (alpha == 0.0) return;

    if (!is_trans(trans)) {
        /* y = alpha * A * x + y */
        /* Column-major: A[i,j] = A[i + j*lda] */
        int32_t jx = (incx > 0) ? 0 : (1 - n) * incx;
        for (int32_t j = 0; j < n; j++) {
            double temp = alpha * x[jx];
            int32_t iy = (incy > 0) ? 0 : (1 - m) * incy;
            for (int32_t i = 0; i < m; i++) {
                y[iy] += temp * A[i + j * lda];
                iy += incy;
            }
            jx += incx;
        }
    } else {
        /* y = alpha * A^T * x + y */
        int32_t jy = (incy > 0) ? 0 : (1 - n) * incy;
        for (int32_t j = 0; j < n; j++) {
            double temp = 0.0;
            int32_t ix = (incx > 0) ? 0 : (1 - m) * incx;
            for (int32_t i = 0; i < m; i++) {
                temp += A[i + j * lda] * x[ix];
                ix += incx;
            }
            y[jy] += alpha * temp;
            jy += incy;
        }
    }
}

EXPORT void blas_sgemv(char trans, int32_t m, int32_t n,
                        float alpha, const float* A, int32_t lda,
                        const float* x, int32_t incx,
                        float beta, float* y, int32_t incy)
{
    if (m <= 0 || n <= 0) return;

    int32_t lenx, leny;
    if (!is_trans(trans)) {
        lenx = n;
        leny = m;
    } else {
        lenx = m;
        leny = n;
    }

    if (beta != 1.0f) {
        if (incy == 1) {
            if (beta == 0.0f) {
                memset(y, 0, leny * sizeof(float));
            } else {
                for (int32_t i = 0; i < leny; i++) {
                    y[i] *= beta;
                }
            }
        } else {
            int32_t iy = (incy > 0) ? 0 : (1 - leny) * incy;
            if (beta == 0.0f) {
                for (int32_t i = 0; i < leny; i++) {
                    y[iy] = 0.0f;
                    iy += incy;
                }
            } else {
                for (int32_t i = 0; i < leny; i++) {
                    y[iy] *= beta;
                    iy += incy;
                }
            }
        }
    }

    if (alpha == 0.0f) return;

    if (!is_trans(trans)) {
        int32_t jx = (incx > 0) ? 0 : (1 - n) * incx;
        for (int32_t j = 0; j < n; j++) {
            float temp = alpha * x[jx];
            int32_t iy = (incy > 0) ? 0 : (1 - m) * incy;
            for (int32_t i = 0; i < m; i++) {
                y[iy] += temp * A[i + j * lda];
                iy += incy;
            }
            jx += incx;
        }
    } else {
        int32_t jy = (incy > 0) ? 0 : (1 - n) * incy;
        for (int32_t j = 0; j < n; j++) {
            float temp = 0.0f;
            int32_t ix = (incx > 0) ? 0 : (1 - m) * incx;
            for (int32_t i = 0; i < m; i++) {
                temp += A[i + j * lda] * x[ix];
                ix += incx;
            }
            y[jy] += alpha * temp;
            jy += incy;
        }
    }
}

EXPORT void blas_dtrmv(char uplo, char trans, char diag, int32_t n,
                        const double* A, int32_t lda,
                        double* x, int32_t incx)
{
    if (n <= 0) return;

    int nounit = !is_unit(diag);

    if (!is_trans(trans)) {
        /* x = A * x */
        if (is_upper(uplo)) {
            /* Upper triangular */
            int32_t jx = (incx > 0) ? 0 : (1 - n) * incx;
            for (int32_t j = 0; j < n; j++) {
                double temp = x[jx];
                int32_t ix = (incx > 0) ? 0 : (1 - n) * incx;
                for (int32_t i = 0; i < j; i++) {
                    x[ix] += temp * A[i + j * lda];
                    ix += incx;
                }
                if (nounit) {
                    x[jx] *= A[j + j * lda];
                }
                jx += incx;
            }
        } else {
            /* Lower triangular */
            int32_t jx = (incx > 0) ? (n - 1) * incx : 0;
            for (int32_t j = n - 1; j >= 0; j--) {
                double temp = x[jx];
                int32_t ix = (incx > 0) ? (n - 1) * incx : 0;
                for (int32_t i = n - 1; i > j; i--) {
                    x[ix] += temp * A[i + j * lda];
                    ix -= incx;
                }
                if (nounit) {
                    x[jx] *= A[j + j * lda];
                }
                jx -= incx;
            }
        }
    } else {
        /* x = A^T * x */
        if (is_upper(uplo)) {
            int32_t jx = (incx > 0) ? (n - 1) * incx : 0;
            for (int32_t j = n - 1; j >= 0; j--) {
                double temp = x[jx];
                if (nounit) {
                    temp *= A[j + j * lda];
                }
                int32_t ix = jx - incx;
                for (int32_t i = j - 1; i >= 0; i--) {
                    temp += A[i + j * lda] * x[ix];
                    ix -= incx;
                }
                x[jx] = temp;
                jx -= incx;
            }
        } else {
            int32_t jx = (incx > 0) ? 0 : (1 - n) * incx;
            for (int32_t j = 0; j < n; j++) {
                double temp = x[jx];
                if (nounit) {
                    temp *= A[j + j * lda];
                }
                int32_t ix = jx + incx;
                for (int32_t i = j + 1; i < n; i++) {
                    temp += A[i + j * lda] * x[ix];
                    ix += incx;
                }
                x[jx] = temp;
                jx += incx;
            }
        }
    }
}

EXPORT void blas_dtrsv(char uplo, char trans, char diag, int32_t n,
                        const double* A, int32_t lda,
                        double* x, int32_t incx)
{
    if (n <= 0) return;

    int nounit = !is_unit(diag);

    if (!is_trans(trans)) {
        /* Solve A * x = b */
        if (is_upper(uplo)) {
            /* Upper triangular: back substitution */
            int32_t jx = (incx > 0) ? (n - 1) * incx : 0;
            for (int32_t j = n - 1; j >= 0; j--) {
                if (nounit) {
                    x[jx] /= A[j + j * lda];
                }
                double temp = x[jx];
                int32_t ix = jx - incx;
                for (int32_t i = j - 1; i >= 0; i--) {
                    x[ix] -= temp * A[i + j * lda];
                    ix -= incx;
                }
                jx -= incx;
            }
        } else {
            /* Lower triangular: forward substitution */
            int32_t jx = (incx > 0) ? 0 : (1 - n) * incx;
            for (int32_t j = 0; j < n; j++) {
                if (nounit) {
                    x[jx] /= A[j + j * lda];
                }
                double temp = x[jx];
                int32_t ix = jx + incx;
                for (int32_t i = j + 1; i < n; i++) {
                    x[ix] -= temp * A[i + j * lda];
                    ix += incx;
                }
                jx += incx;
            }
        }
    } else {
        /* Solve A^T * x = b */
        if (is_upper(uplo)) {
            int32_t jx = (incx > 0) ? 0 : (1 - n) * incx;
            for (int32_t j = 0; j < n; j++) {
                double temp = x[jx];
                int32_t ix = (incx > 0) ? 0 : (1 - n) * incx;
                for (int32_t i = 0; i < j; i++) {
                    temp -= A[i + j * lda] * x[ix];
                    ix += incx;
                }
                if (nounit) {
                    temp /= A[j + j * lda];
                }
                x[jx] = temp;
                jx += incx;
            }
        } else {
            int32_t jx = (incx > 0) ? (n - 1) * incx : 0;
            for (int32_t j = n - 1; j >= 0; j--) {
                double temp = x[jx];
                int32_t ix = (incx > 0) ? (n - 1) * incx : 0;
                for (int32_t i = n - 1; i > j; i--) {
                    temp -= A[i + j * lda] * x[ix];
                    ix -= incx;
                }
                if (nounit) {
                    temp /= A[j + j * lda];
                }
                x[jx] = temp;
                jx -= incx;
            }
        }
    }
}

EXPORT void blas_dger(int32_t m, int32_t n, double alpha,
                       const double* x, int32_t incx,
                       const double* y, int32_t incy,
                       double* A, int32_t lda)
{
    if (m <= 0 || n <= 0 || alpha == 0.0) return;

    int32_t jy = (incy > 0) ? 0 : (1 - n) * incy;
    for (int32_t j = 0; j < n; j++) {
        double temp = alpha * y[jy];
        int32_t ix = (incx > 0) ? 0 : (1 - m) * incx;
        for (int32_t i = 0; i < m; i++) {
            A[i + j * lda] += x[ix] * temp;
            ix += incx;
        }
        jy += incy;
    }
}

/* ============ Level 3 BLAS Implementation ============ */

EXPORT void blas_dgemm(char transA, char transB,
                        int32_t m, int32_t n, int32_t k,
                        double alpha,
                        const double* A, int32_t lda,
                        const double* B, int32_t ldb,
                        double beta,
                        double* C, int32_t ldc)
{
    if (m <= 0 || n <= 0) return;

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

    if (alpha == 0.0 || k <= 0) return;

    int notransA = !is_trans(transA);
    int notransB = !is_trans(transB);

    if (notransA && notransB) {
        /* C = alpha * A * B + C */
        /* A is m x k, B is k x n, C is m x n */
        for (int32_t j = 0; j < n; j++) {
            for (int32_t l = 0; l < k; l++) {
                double temp = alpha * B[l + j * ldb];
                for (int32_t i = 0; i < m; i++) {
                    C[i + j * ldc] += temp * A[i + l * lda];
                }
            }
        }
    } else if (!notransA && notransB) {
        /* C = alpha * A^T * B + C */
        /* A^T is m x k (A is k x m), B is k x n */
        for (int32_t j = 0; j < n; j++) {
            for (int32_t i = 0; i < m; i++) {
                double temp = 0.0;
                for (int32_t l = 0; l < k; l++) {
                    temp += A[l + i * lda] * B[l + j * ldb];
                }
                C[i + j * ldc] += alpha * temp;
            }
        }
    } else if (notransA && !notransB) {
        /* C = alpha * A * B^T + C */
        /* A is m x k, B^T is k x n (B is n x k) */
        for (int32_t j = 0; j < n; j++) {
            for (int32_t l = 0; l < k; l++) {
                double temp = alpha * B[j + l * ldb];
                for (int32_t i = 0; i < m; i++) {
                    C[i + j * ldc] += temp * A[i + l * lda];
                }
            }
        }
    } else {
        /* C = alpha * A^T * B^T + C */
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

EXPORT void blas_sgemm(char transA, char transB,
                        int32_t m, int32_t n, int32_t k,
                        float alpha,
                        const float* A, int32_t lda,
                        const float* B, int32_t ldb,
                        float beta,
                        float* C, int32_t ldc)
{
    if (m <= 0 || n <= 0) return;

    if (beta == 0.0f) {
        for (int32_t j = 0; j < n; j++) {
            for (int32_t i = 0; i < m; i++) {
                C[i + j * ldc] = 0.0f;
            }
        }
    } else if (beta != 1.0f) {
        for (int32_t j = 0; j < n; j++) {
            for (int32_t i = 0; i < m; i++) {
                C[i + j * ldc] *= beta;
            }
        }
    }

    if (alpha == 0.0f || k <= 0) return;

    int notransA = !is_trans(transA);
    int notransB = !is_trans(transB);

    if (notransA && notransB) {
        for (int32_t j = 0; j < n; j++) {
            for (int32_t l = 0; l < k; l++) {
                float temp = alpha * B[l + j * ldb];
                for (int32_t i = 0; i < m; i++) {
                    C[i + j * ldc] += temp * A[i + l * lda];
                }
            }
        }
    } else if (!notransA && notransB) {
        for (int32_t j = 0; j < n; j++) {
            for (int32_t i = 0; i < m; i++) {
                float temp = 0.0f;
                for (int32_t l = 0; l < k; l++) {
                    temp += A[l + i * lda] * B[l + j * ldb];
                }
                C[i + j * ldc] += alpha * temp;
            }
        }
    } else if (notransA && !notransB) {
        for (int32_t j = 0; j < n; j++) {
            for (int32_t l = 0; l < k; l++) {
                float temp = alpha * B[j + l * ldb];
                for (int32_t i = 0; i < m; i++) {
                    C[i + j * ldc] += temp * A[i + l * lda];
                }
            }
        }
    } else {
        for (int32_t j = 0; j < n; j++) {
            for (int32_t i = 0; i < m; i++) {
                float temp = 0.0f;
                for (int32_t l = 0; l < k; l++) {
                    temp += A[l + i * lda] * B[j + l * ldb];
                }
                C[i + j * ldc] += alpha * temp;
            }
        }
    }
}

EXPORT void blas_dtrmm(char side, char uplo, char transA, char diag,
                        int32_t m, int32_t n,
                        double alpha,
                        const double* A, int32_t lda,
                        double* B, int32_t ldb)
{
    if (m <= 0 || n <= 0) return;

    if (alpha == 0.0) {
        for (int32_t j = 0; j < n; j++) {
            for (int32_t i = 0; i < m; i++) {
                B[i + j * ldb] = 0.0;
            }
        }
        return;
    }

    int lside = is_left(side);
    int upper = is_upper(uplo);
    int nounit = !is_unit(diag);
    int notrans = !is_trans(transA);

    if (lside) {
        /* B = alpha * op(A) * B */
        if (notrans) {
            if (upper) {
                for (int32_t j = 0; j < n; j++) {
                    for (int32_t k = 0; k < m; k++) {
                        double temp = alpha * B[k + j * ldb];
                        for (int32_t i = 0; i < k; i++) {
                            B[i + j * ldb] += temp * A[i + k * lda];
                        }
                        if (nounit) {
                            temp *= A[k + k * lda];
                        }
                        B[k + j * ldb] = temp;
                    }
                }
            } else {
                for (int32_t j = 0; j < n; j++) {
                    for (int32_t k = m - 1; k >= 0; k--) {
                        double temp = alpha * B[k + j * ldb];
                        if (nounit) {
                            temp *= A[k + k * lda];
                        }
                        B[k + j * ldb] = temp;
                        for (int32_t i = k + 1; i < m; i++) {
                            B[i + j * ldb] += temp * A[i + k * lda];
                        }
                    }
                }
            }
        } else {
            if (upper) {
                for (int32_t j = 0; j < n; j++) {
                    for (int32_t i = m - 1; i >= 0; i--) {
                        double temp = B[i + j * ldb];
                        if (nounit) {
                            temp *= A[i + i * lda];
                        }
                        for (int32_t k = 0; k < i; k++) {
                            temp += A[k + i * lda] * B[k + j * ldb];
                        }
                        B[i + j * ldb] = alpha * temp;
                    }
                }
            } else {
                for (int32_t j = 0; j < n; j++) {
                    for (int32_t i = 0; i < m; i++) {
                        double temp = B[i + j * ldb];
                        if (nounit) {
                            temp *= A[i + i * lda];
                        }
                        for (int32_t k = i + 1; k < m; k++) {
                            temp += A[k + i * lda] * B[k + j * ldb];
                        }
                        B[i + j * ldb] = alpha * temp;
                    }
                }
            }
        }
    } else {
        /* B = alpha * B * op(A) */
        if (notrans) {
            if (upper) {
                for (int32_t j = n - 1; j >= 0; j--) {
                    double temp = alpha;
                    if (nounit) {
                        temp *= A[j + j * lda];
                    }
                    for (int32_t i = 0; i < m; i++) {
                        B[i + j * ldb] *= temp;
                    }
                    for (int32_t k = 0; k < j; k++) {
                        double akj = alpha * A[k + j * lda];
                        for (int32_t i = 0; i < m; i++) {
                            B[i + j * ldb] += akj * B[i + k * ldb];
                        }
                    }
                }
            } else {
                for (int32_t j = 0; j < n; j++) {
                    double temp = alpha;
                    if (nounit) {
                        temp *= A[j + j * lda];
                    }
                    for (int32_t i = 0; i < m; i++) {
                        B[i + j * ldb] *= temp;
                    }
                    for (int32_t k = j + 1; k < n; k++) {
                        double akj = alpha * A[k + j * lda];
                        for (int32_t i = 0; i < m; i++) {
                            B[i + j * ldb] += akj * B[i + k * ldb];
                        }
                    }
                }
            }
        } else {
            if (upper) {
                for (int32_t k = 0; k < n; k++) {
                    for (int32_t j = 0; j < k; j++) {
                        double ajk = alpha * A[j + k * lda];
                        for (int32_t i = 0; i < m; i++) {
                            B[i + j * ldb] += ajk * B[i + k * ldb];
                        }
                    }
                    double temp = alpha;
                    if (nounit) {
                        temp *= A[k + k * lda];
                    }
                    for (int32_t i = 0; i < m; i++) {
                        B[i + k * ldb] *= temp;
                    }
                }
            } else {
                for (int32_t k = n - 1; k >= 0; k--) {
                    for (int32_t j = k + 1; j < n; j++) {
                        double ajk = alpha * A[j + k * lda];
                        for (int32_t i = 0; i < m; i++) {
                            B[i + j * ldb] += ajk * B[i + k * ldb];
                        }
                    }
                    double temp = alpha;
                    if (nounit) {
                        temp *= A[k + k * lda];
                    }
                    for (int32_t i = 0; i < m; i++) {
                        B[i + k * ldb] *= temp;
                    }
                }
            }
        }
    }
}

EXPORT void blas_dtrsm(char side, char uplo, char transA, char diag,
                        int32_t m, int32_t n,
                        double alpha,
                        const double* A, int32_t lda,
                        double* B, int32_t ldb)
{
    if (m <= 0 || n <= 0) return;

    if (alpha == 0.0) {
        for (int32_t j = 0; j < n; j++) {
            for (int32_t i = 0; i < m; i++) {
                B[i + j * ldb] = 0.0;
            }
        }
        return;
    }

    int lside = is_left(side);
    int upper = is_upper(uplo);
    int nounit = !is_unit(diag);
    int notrans = !is_trans(transA);

    if (lside) {
        /* Solve op(A) * X = alpha * B */
        if (notrans) {
            if (upper) {
                /* Upper triangular, no transpose: back substitution */
                for (int32_t j = 0; j < n; j++) {
                    for (int32_t i = 0; i < m; i++) {
                        B[i + j * ldb] *= alpha;
                    }
                    for (int32_t k = m - 1; k >= 0; k--) {
                        if (nounit) {
                            B[k + j * ldb] /= A[k + k * lda];
                        }
                        double temp = B[k + j * ldb];
                        for (int32_t i = 0; i < k; i++) {
                            B[i + j * ldb] -= temp * A[i + k * lda];
                        }
                    }
                }
            } else {
                /* Lower triangular, no transpose: forward substitution */
                for (int32_t j = 0; j < n; j++) {
                    for (int32_t i = 0; i < m; i++) {
                        B[i + j * ldb] *= alpha;
                    }
                    for (int32_t k = 0; k < m; k++) {
                        if (nounit) {
                            B[k + j * ldb] /= A[k + k * lda];
                        }
                        double temp = B[k + j * ldb];
                        for (int32_t i = k + 1; i < m; i++) {
                            B[i + j * ldb] -= temp * A[i + k * lda];
                        }
                    }
                }
            }
        } else {
            /* Transposed cases */
            if (upper) {
                for (int32_t j = 0; j < n; j++) {
                    for (int32_t i = 0; i < m; i++) {
                        double temp = alpha * B[i + j * ldb];
                        for (int32_t k = 0; k < i; k++) {
                            temp -= A[k + i * lda] * B[k + j * ldb];
                        }
                        if (nounit) {
                            temp /= A[i + i * lda];
                        }
                        B[i + j * ldb] = temp;
                    }
                }
            } else {
                for (int32_t j = 0; j < n; j++) {
                    for (int32_t i = m - 1; i >= 0; i--) {
                        double temp = alpha * B[i + j * ldb];
                        for (int32_t k = i + 1; k < m; k++) {
                            temp -= A[k + i * lda] * B[k + j * ldb];
                        }
                        if (nounit) {
                            temp /= A[i + i * lda];
                        }
                        B[i + j * ldb] = temp;
                    }
                }
            }
        }
    } else {
        /* Solve X * op(A) = alpha * B */
        if (notrans) {
            if (upper) {
                for (int32_t j = 0; j < n; j++) {
                    for (int32_t i = 0; i < m; i++) {
                        B[i + j * ldb] *= alpha;
                    }
                    for (int32_t k = 0; k < j; k++) {
                        double akj = A[k + j * lda];
                        for (int32_t i = 0; i < m; i++) {
                            B[i + j * ldb] -= akj * B[i + k * ldb];
                        }
                    }
                    if (nounit) {
                        double ajj = 1.0 / A[j + j * lda];
                        for (int32_t i = 0; i < m; i++) {
                            B[i + j * ldb] *= ajj;
                        }
                    }
                }
            } else {
                for (int32_t j = n - 1; j >= 0; j--) {
                    for (int32_t i = 0; i < m; i++) {
                        B[i + j * ldb] *= alpha;
                    }
                    for (int32_t k = j + 1; k < n; k++) {
                        double akj = A[k + j * lda];
                        for (int32_t i = 0; i < m; i++) {
                            B[i + j * ldb] -= akj * B[i + k * ldb];
                        }
                    }
                    if (nounit) {
                        double ajj = 1.0 / A[j + j * lda];
                        for (int32_t i = 0; i < m; i++) {
                            B[i + j * ldb] *= ajj;
                        }
                    }
                }
            }
        } else {
            if (upper) {
                for (int32_t k = n - 1; k >= 0; k--) {
                    if (nounit) {
                        double akk = 1.0 / A[k + k * lda];
                        for (int32_t i = 0; i < m; i++) {
                            B[i + k * ldb] *= akk;
                        }
                    }
                    for (int32_t j = 0; j < k; j++) {
                        double ajk = A[j + k * lda];
                        for (int32_t i = 0; i < m; i++) {
                            B[i + j * ldb] -= ajk * B[i + k * ldb];
                        }
                    }
                    for (int32_t i = 0; i < m; i++) {
                        B[i + k * ldb] *= alpha;
                    }
                }
            } else {
                for (int32_t k = 0; k < n; k++) {
                    if (nounit) {
                        double akk = 1.0 / A[k + k * lda];
                        for (int32_t i = 0; i < m; i++) {
                            B[i + k * ldb] *= akk;
                        }
                    }
                    for (int32_t j = k + 1; j < n; j++) {
                        double ajk = A[j + k * lda];
                        for (int32_t i = 0; i < m; i++) {
                            B[i + j * ldb] -= ajk * B[i + k * ldb];
                        }
                    }
                    for (int32_t i = 0; i < m; i++) {
                        B[i + k * ldb] *= alpha;
                    }
                }
            }
        }
    }
}

EXPORT void blas_dsyrk(char uplo, char trans, int32_t n, int32_t k,
                        double alpha, const double* A, int32_t lda,
                        double beta, double* C, int32_t ldc)
{
    if (n <= 0) return;

    int upper = is_upper(uplo);
    int notrans = !is_trans(trans);

    /* Scale C by beta */
    if (beta == 0.0) {
        if (upper) {
            for (int32_t j = 0; j < n; j++) {
                for (int32_t i = 0; i <= j; i++) {
                    C[i + j * ldc] = 0.0;
                }
            }
        } else {
            for (int32_t j = 0; j < n; j++) {
                for (int32_t i = j; i < n; i++) {
                    C[i + j * ldc] = 0.0;
                }
            }
        }
    } else if (beta != 1.0) {
        if (upper) {
            for (int32_t j = 0; j < n; j++) {
                for (int32_t i = 0; i <= j; i++) {
                    C[i + j * ldc] *= beta;
                }
            }
        } else {
            for (int32_t j = 0; j < n; j++) {
                for (int32_t i = j; i < n; i++) {
                    C[i + j * ldc] *= beta;
                }
            }
        }
    }

    if (alpha == 0.0 || k <= 0) return;

    if (notrans) {
        /* C = alpha * A * A^T + beta * C */
        if (upper) {
            for (int32_t j = 0; j < n; j++) {
                for (int32_t l = 0; l < k; l++) {
                    double temp = alpha * A[j + l * lda];
                    for (int32_t i = 0; i <= j; i++) {
                        C[i + j * ldc] += temp * A[i + l * lda];
                    }
                }
            }
        } else {
            for (int32_t j = 0; j < n; j++) {
                for (int32_t l = 0; l < k; l++) {
                    double temp = alpha * A[j + l * lda];
                    for (int32_t i = j; i < n; i++) {
                        C[i + j * ldc] += temp * A[i + l * lda];
                    }
                }
            }
        }
    } else {
        /* C = alpha * A^T * A + beta * C */
        if (upper) {
            for (int32_t j = 0; j < n; j++) {
                for (int32_t i = 0; i <= j; i++) {
                    double temp = 0.0;
                    for (int32_t l = 0; l < k; l++) {
                        temp += A[l + i * lda] * A[l + j * lda];
                    }
                    C[i + j * ldc] += alpha * temp;
                }
            }
        } else {
            for (int32_t j = 0; j < n; j++) {
                for (int32_t i = j; i < n; i++) {
                    double temp = 0.0;
                    for (int32_t l = 0; l < k; l++) {
                        temp += A[l + i * lda] * A[l + j * lda];
                    }
                    C[i + j * ldc] += alpha * temp;
                }
            }
        }
    }
}
