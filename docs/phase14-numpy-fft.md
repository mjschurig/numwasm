# Phase 14: numpy.fft Implementation Plan

Complete implementation roadmap for the NumJS-WASM Fast Fourier Transform module, providing NumPy-compatible FFT operations with WebAssembly acceleration.

---

## Implementation Guideline

**The original NumPy source code is available in `/numpy`.**

Key reference files:
- `numpy/fft/_pocketfft.py` - Main Python FFT interface (~1,693 lines)
- `numpy/fft/_pocketfft_umath.cpp` - C++/Cython pocketfft wrapper
- `numpy/fft/_helper.py` - Helper functions (fftfreq, fftshift)
- `numpy/fft/tests/test_pocketfft.py` - Comprehensive test suite

Implementation should follow NumPy's algorithms and normalization conventions for consistency.

---

## Current State (Pre-Phase 14)

```
src/wasm/
├── ndarray.h/c        # Core NDArray with views, slicing
├── dtype.h/c          # DType system (Float32, Float64, Complex64, Complex128)
├── broadcast.h/c      # Broadcasting
├── indexing.h/c       # Index operations
├── ufunc.h/c          # Ufunc infrastructure
├── ufunc_unary.c      # Unary operations (exp, sin, cos, etc.)
├── ufunc_binary.c     # Binary operations
├── manipulation.c     # Array manipulation (concat, roll, etc.)
├── statistics.c       # Statistical operations
└── setops.c           # Set operations

src/ts/
├── NDArray.ts         # Core array class
├── types.ts           # Type definitions
├── dtype.ts           # Type utilities
├── broadcast.ts       # Broadcasting functions
├── indexing.ts        # Index operations
├── manipulation.ts    # Array manipulation
├── ufunc.ts           # Ufunc TypeScript layer
└── index.ts           # Public exports
```

**Required Infrastructure:**
- Complex number support (Complex64, Complex128) - already exists
- Array manipulation (`roll` for fftshift) - already exists
- Element-wise operations (multiply, divide, exp) - already exists

---

## Phase 14 Dependency Tree

```
PHASE 14: NUMPY.FFT
│
├── 14.1 Core FFT Algorithm (C/WASM)
│   ├── 14.1.1 Twiddle Factor Generation
│   │   ├── fft_twiddle_factors(n) → complex array
│   │   └── Caching for repeated sizes
│   │
│   ├── 14.1.2 Cooley-Tukey Radix-2 FFT
│   │   ├── fft_radix2_dit(x, n) → Decimation-in-time
│   │   ├── fft_radix2_dif(x, n) → Decimation-in-frequency
│   │   └── Bit-reversal permutation
│   │
│   ├── 14.1.3 Mixed-Radix FFT
│   │   ├── fft_radix_4(x, n) → Radix-4 butterfly
│   │   ├── fft_radix_3(x, n) → Radix-3 butterfly
│   │   ├── fft_radix_5(x, n) → Radix-5 butterfly
│   │   └── fft_composite(x, n) → Factor-based decomposition
│   │
│   ├── 14.1.4 Bluestein's Algorithm (Chirp-Z)
│   │   ├── fft_bluestein(x, n) → Arbitrary size via convolution
│   │   └── Used for prime-sized transforms
│   │
│   └── 14.1.5 Real FFT Optimizations
│       ├── rfft_real_to_complex(x, n) → Half-spectrum
│       ├── irfft_complex_to_real(x, n) → Hermitian reconstruction
│       └── Hermitian symmetry exploitation
│
│   Dependencies: NDArray core, Complex DType
│
├── 14.2 1D FFT Functions (TypeScript + C)
│   ├── 14.2.1 Complex-to-Complex
│   │   ├── fft(a, n, axis, norm, out)
│   │   └── ifft(a, n, axis, norm, out)
│   │
│   ├── 14.2.2 Real-to-Complex
│   │   ├── rfft(a, n, axis, norm, out) → n//2+1 complex
│   │   └── irfft(a, n, axis, norm, out) → n real
│   │
│   └── 14.2.3 Hermitian Transforms
│       ├── hfft(a, n, axis, norm, out) → real output
│       └── ihfft(a, n, axis, norm, out) → complex output
│
│   Dependencies: 14.1.* (Core FFT)
│
├── 14.3 Multi-Dimensional FFT (TypeScript + C)
│   ├── 14.3.1 2D Transforms
│   │   ├── fft2(a, s, axes, norm, out)
│   │   ├── ifft2(a, s, axes, norm, out)
│   │   ├── rfft2(a, s, axes, norm, out)
│   │   └── irfft2(a, s, axes, norm, out)
│   │
│   └── 14.3.2 N-Dimensional Transforms
│       ├── fftn(a, s, axes, norm, out)
│       ├── ifftn(a, s, axes, norm, out)
│       ├── rfftn(a, s, axes, norm, out)
│       └── irfftn(a, s, axes, norm, out)
│
│   Dependencies: 14.2.* (1D FFT)
│
└── 14.4 Helper Functions (TypeScript)
    ├── 14.4.1 Frequency Arrays
    │   ├── fftfreq(n, d, device) → frequency bins
    │   └── rfftfreq(n, d, device) → positive frequencies
    │
    └── 14.4.2 Shifting Functions
        ├── fftshift(x, axes) → zero-frequency to center
        └── ifftshift(x, axes) → inverse shift

    Dependencies: Core NDArray, manipulation (roll)
```

---

## Mathematical Background

### Discrete Fourier Transform (DFT)

**Forward DFT:**
```
A[k] = Σ(m=0 to n-1) a[m] * W_n^(mk)    for k = 0, 1, ..., n-1

where W_n = exp(-2πi/n) is the primitive nth root of unity
```

**Inverse DFT:**
```
a[m] = (1/n) * Σ(k=0 to n-1) A[k] * W_n^(-mk)    for m = 0, 1, ..., n-1
```

### Normalization Modes

| Mode | Forward | Inverse |
|------|---------|---------|
| `"backward"` (default) | 1 | 1/n |
| `"ortho"` | 1/√n | 1/√n |
| `"forward"` | 1/n | 1 |

### Real FFT (Hermitian Symmetry)

For real input `a`, the DFT output satisfies:
```
A[n-k] = conj(A[k])
```

Therefore, only `n//2 + 1` unique complex values need to be stored:
- DC component: `A[0]` (always real)
- Positive frequencies: `A[1], A[2], ..., A[n//2-1]`
- Nyquist frequency: `A[n//2]` (real for even n)

---

## Detailed Implementation Specifications

### 14.1 Core FFT Algorithm

#### 14.1.1 Twiddle Factors

**File:** `src/wasm/fft.h`

```c
#ifndef NUMJS_FFT_H
#define NUMJS_FFT_H

#include "ndarray.h"

/* ============ Twiddle Factor Generation ============ */

/**
 * Compute twiddle factors W_n^k = exp(-2πik/n) for k = 0, 1, ..., n/2-1
 *
 * @param n     FFT size
 * @param out   Output array (complex, size n/2)
 */
void fft_twiddle_factors(int32_t n, double* out);

/**
 * Compute inverse twiddle factors W_n^(-k) = exp(2πik/n)
 */
void fft_twiddle_factors_inv(int32_t n, double* out);

/* ============ Bit Reversal ============ */

/**
 * Compute bit-reversed index for radix-2 FFT.
 *
 * @param i     Input index
 * @param log2n log2(n) where n is FFT size
 * @return      Bit-reversed index
 */
int32_t fft_bit_reverse(int32_t i, int32_t log2n);

/**
 * Perform in-place bit-reversal permutation on complex array.
 *
 * @param data  Complex array (interleaved real/imag)
 * @param n     Number of complex elements
 */
void fft_bit_reverse_permute(double* data, int32_t n);

#endif /* NUMJS_FFT_H */
```

**File:** `src/wasm/fft.c`

```c
#include "fft.h"
#include <math.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/* ============ Twiddle Factor Implementation ============ */

EXPORT void fft_twiddle_factors(int32_t n, double* out)
{
    double theta = -2.0 * M_PI / (double)n;
    for (int32_t k = 0; k < n / 2; k++) {
        double angle = theta * k;
        out[2 * k]     = cos(angle);  /* Real part */
        out[2 * k + 1] = sin(angle);  /* Imaginary part */
    }
}

EXPORT void fft_twiddle_factors_inv(int32_t n, double* out)
{
    double theta = 2.0 * M_PI / (double)n;
    for (int32_t k = 0; k < n / 2; k++) {
        double angle = theta * k;
        out[2 * k]     = cos(angle);
        out[2 * k + 1] = sin(angle);
    }
}

/* ============ Bit Reversal Implementation ============ */

EXPORT int32_t fft_bit_reverse(int32_t i, int32_t log2n)
{
    int32_t rev = 0;
    for (int32_t j = 0; j < log2n; j++) {
        rev = (rev << 1) | (i & 1);
        i >>= 1;
    }
    return rev;
}

EXPORT void fft_bit_reverse_permute(double* data, int32_t n)
{
    int32_t log2n = 0;
    int32_t temp = n;
    while (temp > 1) {
        temp >>= 1;
        log2n++;
    }

    for (int32_t i = 0; i < n; i++) {
        int32_t j = fft_bit_reverse(i, log2n);
        if (j > i) {
            /* Swap complex elements */
            double temp_re = data[2 * i];
            double temp_im = data[2 * i + 1];
            data[2 * i]     = data[2 * j];
            data[2 * i + 1] = data[2 * j + 1];
            data[2 * j]     = temp_re;
            data[2 * j + 1] = temp_im;
        }
    }
}
```

---

#### 14.1.2 Cooley-Tukey Radix-2 FFT

**File:** `src/wasm/fft.h` (additions)

```c
/* ============ Radix-2 FFT ============ */

/**
 * In-place radix-2 decimation-in-time FFT.
 * Input must be power-of-2 length.
 *
 * @param data      Complex array (interleaved real/imag, size 2*n)
 * @param n         Number of complex elements (must be power of 2)
 * @param inverse   0 for forward FFT, 1 for inverse FFT
 * @return          0 on success, -1 if n is not power of 2
 */
int32_t fft_radix2(double* data, int32_t n, int32_t inverse);

/**
 * Check if n is a power of 2.
 */
int32_t fft_is_power_of_2(int32_t n);
```

**File:** `src/wasm/fft.c` (additions)

```c
EXPORT int32_t fft_is_power_of_2(int32_t n)
{
    return (n > 0) && ((n & (n - 1)) == 0);
}

EXPORT int32_t fft_radix2(double* data, int32_t n, int32_t inverse)
{
    if (!fft_is_power_of_2(n)) {
        return -1;  /* Error: n must be power of 2 */
    }

    if (n <= 1) return 0;

    /* Bit-reversal permutation */
    fft_bit_reverse_permute(data, n);

    /* Direction: -1 for forward, +1 for inverse */
    double sign = inverse ? 1.0 : -1.0;

    /* Cooley-Tukey iterative algorithm */
    for (int32_t len = 2; len <= n; len <<= 1) {
        double theta = sign * 2.0 * M_PI / len;
        double wpr = cos(theta);
        double wpi = sin(theta);

        for (int32_t i = 0; i < n; i += len) {
            double wr = 1.0;
            double wi = 0.0;

            for (int32_t j = 0; j < len / 2; j++) {
                int32_t idx1 = i + j;
                int32_t idx2 = i + j + len / 2;

                /* Butterfly operation */
                double t_re = wr * data[2 * idx2]     - wi * data[2 * idx2 + 1];
                double t_im = wr * data[2 * idx2 + 1] + wi * data[2 * idx2];

                data[2 * idx2]     = data[2 * idx1]     - t_re;
                data[2 * idx2 + 1] = data[2 * idx1 + 1] - t_im;
                data[2 * idx1]     = data[2 * idx1]     + t_re;
                data[2 * idx1 + 1] = data[2 * idx1 + 1] + t_im;

                /* Update twiddle factor */
                double temp = wr;
                wr = wr * wpr - wi * wpi;
                wi = temp * wpi + wi * wpr;
            }
        }
    }

    return 0;
}
```

---

#### 14.1.3 Mixed-Radix FFT

**File:** `src/wasm/fft.h` (additions)

```c
/* ============ Mixed-Radix FFT ============ */

/**
 * Factorize n into prime factors suitable for FFT.
 * Preferred factors: 2, 3, 4, 5
 *
 * @param n         Number to factorize
 * @param factors   Output array for factors (max 32)
 * @return          Number of factors
 */
int32_t fft_factorize(int32_t n, int32_t* factors);

/**
 * Mixed-radix FFT for composite sizes.
 * Handles sizes with factors 2, 3, 4, 5.
 *
 * @param data      Complex array
 * @param n         Number of complex elements
 * @param inverse   0 for forward, 1 for inverse
 * @return          0 on success
 */
int32_t fft_mixed_radix(double* data, int32_t n, int32_t inverse);

/**
 * Radix-4 butterfly operation.
 */
void fft_radix4_butterfly(double* data, int32_t n, int32_t inverse);

/**
 * Radix-3 butterfly operation.
 */
void fft_radix3_butterfly(double* data, int32_t n, int32_t inverse);
```

**File:** `src/wasm/fft.c` (additions)

```c
EXPORT int32_t fft_factorize(int32_t n, int32_t* factors)
{
    int32_t count = 0;

    /* Extract factors of 4 (2^2) for efficiency */
    while (n % 4 == 0 && count < 32) {
        factors[count++] = 4;
        n /= 4;
    }

    /* Extract remaining factor of 2 */
    if (n % 2 == 0 && count < 32) {
        factors[count++] = 2;
        n /= 2;
    }

    /* Extract factors of 3 */
    while (n % 3 == 0 && count < 32) {
        factors[count++] = 3;
        n /= 3;
    }

    /* Extract factors of 5 */
    while (n % 5 == 0 && count < 32) {
        factors[count++] = 5;
        n /= 5;
    }

    /* Remaining odd factors */
    for (int32_t f = 7; f * f <= n; f += 2) {
        while (n % f == 0 && count < 32) {
            factors[count++] = f;
            n /= f;
        }
    }

    if (n > 1 && count < 32) {
        factors[count++] = n;  /* Remaining prime */
    }

    return count;
}

/* Constants for radix-3 */
static const double FFT_SQRT3_2 = 0.86602540378443864676;  /* sqrt(3)/2 */

EXPORT void fft_radix3_butterfly(double* data, int32_t n, int32_t inverse)
{
    double sign = inverse ? 1.0 : -1.0;
    int32_t m = n / 3;

    double theta = sign * 2.0 * M_PI / 3.0;
    double w1_re = cos(theta);
    double w1_im = sin(theta);
    double w2_re = cos(2.0 * theta);
    double w2_im = sin(2.0 * theta);

    for (int32_t k = 0; k < m; k++) {
        /* Load three inputs */
        double a0_re = data[2 * k];
        double a0_im = data[2 * k + 1];
        double a1_re = data[2 * (k + m)];
        double a1_im = data[2 * (k + m) + 1];
        double a2_re = data[2 * (k + 2 * m)];
        double a2_im = data[2 * (k + 2 * m) + 1];

        /* Apply twiddle factors and combine */
        double t1_re = a1_re * w1_re - a1_im * w1_im;
        double t1_im = a1_re * w1_im + a1_im * w1_re;
        double t2_re = a2_re * w2_re - a2_im * w2_im;
        double t2_im = a2_re * w2_im + a2_im * w2_re;

        /* Output */
        data[2 * k]             = a0_re + t1_re + t2_re;
        data[2 * k + 1]         = a0_im + t1_im + t2_im;
        data[2 * (k + m)]       = a0_re - 0.5 * (t1_re + t2_re) + sign * FFT_SQRT3_2 * (t1_im - t2_im);
        data[2 * (k + m) + 1]   = a0_im - 0.5 * (t1_im + t2_im) - sign * FFT_SQRT3_2 * (t1_re - t2_re);
        data[2 * (k + 2*m)]     = a0_re - 0.5 * (t1_re + t2_re) - sign * FFT_SQRT3_2 * (t1_im - t2_im);
        data[2 * (k + 2*m) + 1] = a0_im - 0.5 * (t1_im + t2_im) + sign * FFT_SQRT3_2 * (t1_re - t2_re);
    }
}
```

---

#### 14.1.4 Bluestein's Algorithm

**File:** `src/wasm/fft.h` (additions)

```c
/* ============ Bluestein's Algorithm (Chirp-Z Transform) ============ */

/**
 * FFT for arbitrary size using Bluestein's algorithm.
 * Converts arbitrary-size FFT to power-of-2 convolution.
 *
 * @param data      Complex input/output array
 * @param n         Number of complex elements (any positive integer)
 * @param inverse   0 for forward, 1 for inverse
 * @param work      Workspace (size >= 2 * next_power_of_2(2*n-1) complex)
 * @return          0 on success
 */
int32_t fft_bluestein(double* data, int32_t n, int32_t inverse, double* work);

/**
 * Find the smallest power of 2 >= n.
 */
int32_t fft_next_power_of_2(int32_t n);
```

**File:** `src/wasm/fft.c` (additions)

```c
EXPORT int32_t fft_next_power_of_2(int32_t n)
{
    int32_t p = 1;
    while (p < n) p <<= 1;
    return p;
}

EXPORT int32_t fft_bluestein(double* data, int32_t n, int32_t inverse, double* work)
{
    if (n <= 0) return -1;
    if (n == 1) return 0;

    /* Convolution size: smallest power of 2 >= 2*n - 1 */
    int32_t m = fft_next_power_of_2(2 * n - 1);

    double sign = inverse ? 1.0 : -1.0;
    double theta = sign * M_PI / n;

    /* Pointers into workspace */
    double* chirp = work;              /* Chirp signal: m complex */
    double* conv_a = work + 2 * m;     /* Convolution array a: m complex */
    double* conv_b = work + 4 * m;     /* Convolution array b: m complex */

    /* Initialize chirp signal: w[k] = exp(sign * i * pi * k^2 / n) */
    for (int32_t k = 0; k < n; k++) {
        double angle = theta * k * k;
        chirp[2 * k]     = cos(angle);
        chirp[2 * k + 1] = sin(angle);
    }

    /* Build convolution array a: a[k] = x[k] * conj(chirp[k]) */
    for (int32_t k = 0; k < n; k++) {
        double x_re = data[2 * k];
        double x_im = data[2 * k + 1];
        double c_re = chirp[2 * k];
        double c_im = -chirp[2 * k + 1];  /* Conjugate */
        conv_a[2 * k]     = x_re * c_re - x_im * c_im;
        conv_a[2 * k + 1] = x_re * c_im + x_im * c_re;
    }
    /* Zero-pad conv_a */
    for (int32_t k = n; k < m; k++) {
        conv_a[2 * k]     = 0.0;
        conv_a[2 * k + 1] = 0.0;
    }

    /* Build convolution array b: b[k] = chirp[k] for k < n, chirp[m-k] for k > m-n */
    for (int32_t k = 0; k < m; k++) {
        conv_b[2 * k]     = 0.0;
        conv_b[2 * k + 1] = 0.0;
    }
    for (int32_t k = 0; k < n; k++) {
        conv_b[2 * k]     = chirp[2 * k];
        conv_b[2 * k + 1] = chirp[2 * k + 1];
    }
    for (int32_t k = 1; k < n; k++) {
        conv_b[2 * (m - k)]     = chirp[2 * k];
        conv_b[2 * (m - k) + 1] = chirp[2 * k + 1];
    }

    /* Circular convolution via FFT */
    fft_radix2(conv_a, m, 0);
    fft_radix2(conv_b, m, 0);

    /* Point-wise multiply */
    for (int32_t k = 0; k < m; k++) {
        double a_re = conv_a[2 * k];
        double a_im = conv_a[2 * k + 1];
        double b_re = conv_b[2 * k];
        double b_im = conv_b[2 * k + 1];
        conv_a[2 * k]     = a_re * b_re - a_im * b_im;
        conv_a[2 * k + 1] = a_re * b_im + a_im * b_re;
    }

    /* Inverse FFT */
    fft_radix2(conv_a, m, 1);

    /* Scale and multiply by chirp */
    double scale = 1.0 / m;
    for (int32_t k = 0; k < n; k++) {
        double c_re = conv_a[2 * k] * scale;
        double c_im = conv_a[2 * k + 1] * scale;
        double w_re = chirp[2 * k];
        double w_im = -chirp[2 * k + 1];  /* Conjugate */
        data[2 * k]     = c_re * w_re - c_im * w_im;
        data[2 * k + 1] = c_re * w_im + c_im * w_re;
    }

    return 0;
}
```

---

#### 14.1.5 Real FFT Optimizations

**File:** `src/wasm/fft.h` (additions)

```c
/* ============ Real FFT ============ */

/**
 * Real-to-complex FFT exploiting Hermitian symmetry.
 * Output is n/2 + 1 complex values.
 *
 * @param real_in   Real input array (size n)
 * @param out       Complex output (size n/2 + 1, interleaved)
 * @param n         Number of real input samples
 * @param work      Workspace (size >= n complex)
 * @return          0 on success
 */
int32_t fft_rfft(const double* real_in, double* out, int32_t n, double* work);

/**
 * Complex-to-real inverse FFT.
 * Input is n/2 + 1 complex values with Hermitian symmetry.
 *
 * @param complex_in  Complex input (size n/2 + 1)
 * @param out         Real output (size n)
 * @param n           Number of real output samples
 * @param work        Workspace (size >= n complex)
 * @return            0 on success
 */
int32_t fft_irfft(const double* complex_in, double* out, int32_t n, double* work);

/**
 * Hermitian FFT: complex Hermitian-symmetric input to real output.
 */
int32_t fft_hfft(const double* hermitian_in, double* out, int32_t n, double* work);

/**
 * Inverse Hermitian FFT: real input to complex Hermitian-symmetric output.
 */
int32_t fft_ihfft(const double* real_in, double* out, int32_t n, double* work);
```

**File:** `src/wasm/fft.c` (additions)

```c
EXPORT int32_t fft_rfft(const double* real_in, double* out, int32_t n, double* work)
{
    /* Pack real input as complex with zero imaginary */
    for (int32_t i = 0; i < n; i++) {
        work[2 * i]     = real_in[i];
        work[2 * i + 1] = 0.0;
    }

    /* Perform complex FFT */
    int32_t result;
    if (fft_is_power_of_2(n)) {
        result = fft_radix2(work, n, 0);
    } else {
        /* Use Bluestein for non-power-of-2 */
        int32_t m = fft_next_power_of_2(2 * n - 1);
        double* bluestein_work = work + 2 * n;
        result = fft_bluestein(work, n, 0, bluestein_work);
    }

    if (result != 0) return result;

    /* Copy first n/2 + 1 values (Hermitian symmetry) */
    int32_t out_size = n / 2 + 1;
    for (int32_t i = 0; i < out_size; i++) {
        out[2 * i]     = work[2 * i];
        out[2 * i + 1] = work[2 * i + 1];
    }

    return 0;
}

EXPORT int32_t fft_irfft(const double* complex_in, double* out, int32_t n, double* work)
{
    int32_t in_size = n / 2 + 1;

    /* Reconstruct full complex spectrum using Hermitian symmetry */
    /* X[k] = conj(X[n-k]) for k > n/2 */
    for (int32_t i = 0; i < in_size; i++) {
        work[2 * i]     = complex_in[2 * i];
        work[2 * i + 1] = complex_in[2 * i + 1];
    }

    for (int32_t i = in_size; i < n; i++) {
        int32_t j = n - i;
        work[2 * i]     =  complex_in[2 * j];      /* Real part */
        work[2 * i + 1] = -complex_in[2 * j + 1];  /* Conjugate: negate imag */
    }

    /* Perform inverse FFT */
    int32_t result;
    if (fft_is_power_of_2(n)) {
        result = fft_radix2(work, n, 1);
    } else {
        int32_t m = fft_next_power_of_2(2 * n - 1);
        double* bluestein_work = work + 2 * n;
        result = fft_bluestein(work, n, 1, bluestein_work);
    }

    if (result != 0) return result;

    /* Extract real parts and scale */
    double scale = 1.0 / n;
    for (int32_t i = 0; i < n; i++) {
        out[i] = work[2 * i] * scale;
    }

    return 0;
}
```

---

### 14.2 1D FFT Functions (TypeScript)

**File:** `src/ts/fft.ts`

```typescript
import { NDArray } from './NDArray.js';
import { DType } from './types.js';
import { getWasmModule } from './wasm-loader.js';

/* ============ Types ============ */

/**
 * FFT normalization modes.
 * - "backward": No scaling on forward, 1/n on inverse (default)
 * - "ortho": 1/sqrt(n) on both forward and inverse
 * - "forward": 1/n on forward, no scaling on inverse
 */
export type FFTNorm = 'backward' | 'ortho' | 'forward' | null;

/* ============ Internal Utilities ============ */

/**
 * Get normalization factor for FFT.
 */
function _getNormFactor(n: number, norm: FFTNorm, inverse: boolean): number {
  const effectiveNorm = norm ?? 'backward';

  switch (effectiveNorm) {
    case 'backward':
      return inverse ? 1 / n : 1;
    case 'ortho':
      return 1 / Math.sqrt(n);
    case 'forward':
      return inverse ? 1 : 1 / n;
    default:
      throw new Error(`Unknown normalization mode: ${norm}`);
  }
}

/**
 * Ensure array is complex type for FFT operations.
 */
function _asComplex(arr: NDArray): NDArray {
  if (arr.dtype === DType.Complex64 || arr.dtype === DType.Complex128) {
    return arr;
  }
  // Promote to complex128 for maximum precision
  return arr.astype(DType.Complex128);
}

/**
 * Apply 1D FFT along specified axis using WASM.
 */
function _fft1d(
  arr: NDArray,
  n: number | null,
  axis: number,
  inverse: boolean,
  norm: FFTNorm
): NDArray {
  // Normalize axis
  const ndim = arr.ndim;
  if (axis < 0) axis += ndim;
  if (axis < 0 || axis >= ndim) {
    throw new Error(`axis ${axis} is out of bounds for array of dimension ${ndim}`);
  }

  // Determine FFT size
  const inputSize = arr.shape[axis];
  const fftSize = n ?? inputSize;

  // Prepare input: pad or truncate along axis
  let input = _asComplex(arr);
  if (fftSize !== inputSize) {
    input = _resizeAxis(input, axis, fftSize);
  }

  // Get normalization factor
  const normFactor = _getNormFactor(fftSize, norm, inverse);

  // Apply FFT using WASM
  const result = _applyFFTAxis(input, axis, inverse);

  // Apply normalization
  if (normFactor !== 1) {
    return result.multiply(normFactor);
  }

  return result;
}

/**
 * Resize array along specified axis (pad with zeros or truncate).
 */
function _resizeAxis(arr: NDArray, axis: number, newSize: number): NDArray {
  const oldSize = arr.shape[axis];

  if (newSize === oldSize) {
    return arr;
  }

  if (newSize < oldSize) {
    // Truncate
    const slices = arr.shape.map((s, i) =>
      i === axis ? [0, newSize] : [0, s]
    );
    return arr.slice(...slices);
  }

  // Pad with zeros
  const newShape = [...arr.shape];
  newShape[axis] = newSize;
  const result = NDArray.zeros(newShape, arr.dtype);

  // Copy existing data
  // ... (implementation details)

  return result;
}

/**
 * Apply FFT along specified axis using WASM module.
 */
function _applyFFTAxis(arr: NDArray, axis: number, inverse: boolean): NDArray {
  const wasm = getWasmModule();

  // Move axis to last position for efficient processing
  const axes = [...Array(arr.ndim).keys()];
  const targetAxis = arr.ndim - 1;
  if (axis !== targetAxis) {
    // Swap axis to end
    [axes[axis], axes[targetAxis]] = [axes[targetAxis], axes[axis]];
    arr = arr.transpose(axes);
  }

  const fftSize = arr.shape[arr.ndim - 1];
  const batchSize = arr.size / fftSize;

  // Allocate WASM memory
  const dataPtr = wasm._malloc(arr.size * 16);  // 16 bytes per complex128
  const workSize = fftSize * 2 * 8 * 6;  // Workspace for Bluestein
  const workPtr = wasm._malloc(workSize);

  try {
    // Copy data to WASM
    const data = arr.toTypedArray();
    wasm.HEAPF64.set(data, dataPtr / 8);

    // Process each batch
    for (let b = 0; b < batchSize; b++) {
      const offset = dataPtr + b * fftSize * 16;
      wasm._fft_complex(offset, fftSize, inverse ? 1 : 0, workPtr);
    }

    // Copy result back
    const resultData = new Float64Array(arr.size * 2);
    resultData.set(wasm.HEAPF64.subarray(dataPtr / 8, dataPtr / 8 + arr.size * 2));

    let result = NDArray.fromTypedArray(resultData, arr.shape, DType.Complex128);

    // Restore original axis order
    if (axis !== targetAxis) {
      const inverseAxes = [...axes];
      [inverseAxes[axis], inverseAxes[targetAxis]] = [inverseAxes[targetAxis], inverseAxes[axis]];
      result = result.transpose(inverseAxes);
    }

    return result;
  } finally {
    wasm._free(dataPtr);
    wasm._free(workPtr);
  }
}

/* ============ Public API ============ */

/**
 * Compute the one-dimensional discrete Fourier Transform.
 *
 * @param a - Input array (can be complex or real)
 * @param n - Length of the transformed axis. If n < a.shape[axis], input is truncated.
 *            If n > a.shape[axis], input is zero-padded. Default: a.shape[axis]
 * @param axis - Axis over which to compute the FFT. Default: -1 (last axis)
 * @param norm - Normalization mode: "backward", "ortho", or "forward"
 * @param out - Optional output array (not yet supported)
 * @returns Complex array containing the FFT result
 *
 * @example
 * // Simple FFT of a signal
 * const signal = np.array([1, 2, 3, 4]);
 * const spectrum = np.fft.fft(signal);
 *
 * @example
 * // Zero-padded FFT for higher frequency resolution
 * const spectrum = np.fft.fft(signal, 8);
 */
export function fft(
  a: NDArray,
  n: number | null = null,
  axis: number = -1,
  norm: FFTNorm = null,
  out: NDArray | null = null
): NDArray {
  if (out !== null) {
    throw new Error('out parameter not yet supported');
  }
  return _fft1d(a, n, axis, false, norm);
}

/**
 * Compute the one-dimensional inverse discrete Fourier Transform.
 *
 * @param a - Input array
 * @param n - Length of the transformed axis
 * @param axis - Axis over which to compute the inverse FFT. Default: -1
 * @param norm - Normalization mode
 * @param out - Optional output array (not yet supported)
 * @returns Complex array containing the inverse FFT result
 *
 * @example
 * // Round-trip: ifft(fft(x)) ≈ x
 * const x = np.array([1, 2, 3, 4]);
 * const y = np.fft.ifft(np.fft.fft(x));
 * // y ≈ x (within numerical precision)
 */
export function ifft(
  a: NDArray,
  n: number | null = null,
  axis: number = -1,
  norm: FFTNorm = null,
  out: NDArray | null = null
): NDArray {
  if (out !== null) {
    throw new Error('out parameter not yet supported');
  }
  return _fft1d(a, n, axis, true, norm);
}

/**
 * Compute the one-dimensional discrete Fourier Transform for real input.
 *
 * This function computes the one-dimensional n-point DFT of a real-valued
 * array. The output is Hermitian-symmetric, so only the positive frequencies
 * are returned: output has shape (..., n//2 + 1).
 *
 * @param a - Input array (real-valued)
 * @param n - Number of points in the FFT. If n < a.shape[axis], input is truncated.
 *            If n > a.shape[axis], input is zero-padded.
 * @param axis - Axis over which to compute the FFT. Default: -1
 * @param norm - Normalization mode
 * @param out - Optional output array (not yet supported)
 * @returns Complex array with shape (..., n//2 + 1)
 *
 * @example
 * const signal = np.array([1, 2, 3, 4]);
 * const spectrum = np.fft.rfft(signal);
 * // spectrum.shape = [3]  (4//2 + 1 = 3)
 */
export function rfft(
  a: NDArray,
  n: number | null = null,
  axis: number = -1,
  norm: FFTNorm = null,
  out: NDArray | null = null
): NDArray {
  if (out !== null) {
    throw new Error('out parameter not yet supported');
  }

  // Normalize axis
  const ndim = a.ndim;
  let normAxis = axis < 0 ? axis + ndim : axis;
  if (normAxis < 0 || normAxis >= ndim) {
    throw new Error(`axis ${axis} is out of bounds for array of dimension ${ndim}`);
  }

  const fftSize = n ?? a.shape[normAxis];
  const normFactor = _getNormFactor(fftSize, norm, false);

  // Use WASM rfft implementation
  const result = _applyRFFTAxis(a, normAxis, fftSize, false);

  if (normFactor !== 1) {
    return result.multiply(normFactor);
  }

  return result;
}

/**
 * Compute the inverse of the one-dimensional discrete Fourier Transform for real input.
 *
 * This function computes the inverse of the one-dimensional n-point DFT of
 * Hermitian-symmetric input (from rfft). The input should have shape (..., n//2 + 1).
 *
 * @param a - Input array (complex, Hermitian-symmetric)
 * @param n - Length of the output. If n is not given, it is determined from the
 *            input shape: n = 2 * (a.shape[axis] - 1)
 * @param axis - Axis over which to compute the inverse FFT. Default: -1
 * @param norm - Normalization mode
 * @param out - Optional output array (not yet supported)
 * @returns Real array with shape (..., n)
 *
 * @example
 * // Round-trip: irfft(rfft(x)) ≈ x
 * const x = np.array([1, 2, 3, 4]);
 * const y = np.fft.irfft(np.fft.rfft(x));
 * // y ≈ x (within numerical precision)
 */
export function irfft(
  a: NDArray,
  n: number | null = null,
  axis: number = -1,
  norm: FFTNorm = null,
  out: NDArray | null = null
): NDArray {
  if (out !== null) {
    throw new Error('out parameter not yet supported');
  }

  // Normalize axis
  const ndim = a.ndim;
  let normAxis = axis < 0 ? axis + ndim : axis;
  if (normAxis < 0 || normAxis >= ndim) {
    throw new Error(`axis ${axis} is out of bounds for array of dimension ${ndim}`);
  }

  // Determine output size from input
  const inputSize = a.shape[normAxis];
  const outputSize = n ?? 2 * (inputSize - 1);

  const normFactor = _getNormFactor(outputSize, norm, true);

  // Use WASM irfft implementation
  const result = _applyRFFTAxis(a, normAxis, outputSize, true);

  if (normFactor !== 1) {
    return result.multiply(normFactor);
  }

  return result;
}

/**
 * Compute the FFT of a signal that has Hermitian symmetry (real spectrum).
 *
 * @param a - Input array with Hermitian symmetry
 * @param n - Length of the transformed axis
 * @param axis - Axis over which to compute the FFT. Default: -1
 * @param norm - Normalization mode
 * @param out - Optional output array (not yet supported)
 * @returns Real array
 */
export function hfft(
  a: NDArray,
  n: number | null = null,
  axis: number = -1,
  norm: FFTNorm = null,
  out: NDArray | null = null
): NDArray {
  if (out !== null) {
    throw new Error('out parameter not yet supported');
  }

  // hfft(a) = irfft(conj(a))
  const conj_a = a.conj();
  return irfft(conj_a, n, axis, norm);
}

/**
 * Compute the inverse FFT of a signal that has Hermitian symmetry.
 *
 * @param a - Input array (real-valued)
 * @param n - Length of the transformed axis
 * @param axis - Axis over which to compute the inverse FFT. Default: -1
 * @param norm - Normalization mode
 * @param out - Optional output array (not yet supported)
 * @returns Complex array with Hermitian symmetry
 */
export function ihfft(
  a: NDArray,
  n: number | null = null,
  axis: number = -1,
  norm: FFTNorm = null,
  out: NDArray | null = null
): NDArray {
  if (out !== null) {
    throw new Error('out parameter not yet supported');
  }

  // ihfft(a) = conj(rfft(a))
  const rfft_a = rfft(a, n, axis, norm);
  return rfft_a.conj();
}

/* Internal helper for real FFT */
function _applyRFFTAxis(
  arr: NDArray,
  axis: number,
  fftSize: number,
  inverse: boolean
): NDArray {
  // Implementation using WASM fft_rfft / fft_irfft
  throw new Error('Not implemented');
}
```

---

### 14.3 Multi-Dimensional FFT

**File:** `src/ts/fft.ts` (additions)

```typescript
/* ============ N-Dimensional FFT ============ */

/**
 * Internal function to apply FFT along multiple axes.
 */
function _fftn_impl(
  a: NDArray,
  s: number[] | null,
  axes: number[] | null,
  inverse: boolean,
  norm: FFTNorm
): NDArray {
  const ndim = a.ndim;

  // Default axes: all axes
  let effectiveAxes = axes ?? [...Array(ndim).keys()];

  // Normalize negative axes
  effectiveAxes = effectiveAxes.map(ax => ax < 0 ? ax + ndim : ax);

  // Default sizes: original shape along axes
  const effectiveS = s ?? effectiveAxes.map(ax => a.shape[ax]);

  if (effectiveAxes.length !== effectiveS.length) {
    throw new Error('Shape and axes must have the same length');
  }

  // Apply 1D FFT along each axis in reverse order (for efficiency)
  let result = a;
  for (let i = effectiveAxes.length - 1; i >= 0; i--) {
    const ax = effectiveAxes[i];
    const n = effectiveS[i];
    result = _fft1d(result, n, ax, inverse, norm);
  }

  return result;
}

/**
 * Compute the 2-dimensional discrete Fourier Transform.
 *
 * @param a - Input array
 * @param s - Shape (length of each axis) of the output. Default: shape of a
 * @param axes - Axes over which to compute the FFT. Default: (-2, -1)
 * @param norm - Normalization mode
 * @param out - Optional output array (not yet supported)
 * @returns Complex array containing the 2D FFT result
 *
 * @example
 * const image = np.random.rand(64, 64);
 * const spectrum = np.fft.fft2(image);
 */
export function fft2(
  a: NDArray,
  s: [number, number] | null = null,
  axes: [number, number] = [-2, -1],
  norm: FFTNorm = null,
  out: NDArray | null = null
): NDArray {
  if (out !== null) {
    throw new Error('out parameter not yet supported');
  }
  return _fftn_impl(a, s, axes, false, norm);
}

/**
 * Compute the 2-dimensional inverse discrete Fourier Transform.
 */
export function ifft2(
  a: NDArray,
  s: [number, number] | null = null,
  axes: [number, number] = [-2, -1],
  norm: FFTNorm = null,
  out: NDArray | null = null
): NDArray {
  if (out !== null) {
    throw new Error('out parameter not yet supported');
  }
  return _fftn_impl(a, s, axes, true, norm);
}

/**
 * Compute the 2-dimensional FFT of a real array.
 *
 * @param a - Input array (real-valued)
 * @param s - Shape of the output
 * @param axes - Axes over which to compute the FFT. Default: (-2, -1)
 * @param norm - Normalization mode
 * @param out - Optional output array (not yet supported)
 * @returns Complex array with Hermitian symmetry
 */
export function rfft2(
  a: NDArray,
  s: [number, number] | null = null,
  axes: [number, number] = [-2, -1],
  norm: FFTNorm = null,
  out: NDArray | null = null
): NDArray {
  if (out !== null) {
    throw new Error('out parameter not yet supported');
  }
  return rfftn(a, s, axes, norm);
}

/**
 * Compute the inverse of the 2-dimensional FFT of real input.
 */
export function irfft2(
  a: NDArray,
  s: [number, number] | null = null,
  axes: [number, number] = [-2, -1],
  norm: FFTNorm = null,
  out: NDArray | null = null
): NDArray {
  if (out !== null) {
    throw new Error('out parameter not yet supported');
  }
  return irfftn(a, s, axes, norm);
}

/**
 * Compute the N-dimensional discrete Fourier Transform.
 *
 * @param a - Input array
 * @param s - Shape of the output along transformed axes
 * @param axes - Axes over which to compute the FFT. Default: all axes
 * @param norm - Normalization mode
 * @param out - Optional output array (not yet supported)
 * @returns Complex array containing the N-dimensional FFT result
 *
 * @example
 * const volume = np.random.rand(32, 32, 32);
 * const spectrum = np.fft.fftn(volume);
 */
export function fftn(
  a: NDArray,
  s: number[] | null = null,
  axes: number[] | null = null,
  norm: FFTNorm = null,
  out: NDArray | null = null
): NDArray {
  if (out !== null) {
    throw new Error('out parameter not yet supported');
  }
  return _fftn_impl(a, s, axes, false, norm);
}

/**
 * Compute the N-dimensional inverse discrete Fourier Transform.
 */
export function ifftn(
  a: NDArray,
  s: number[] | null = null,
  axes: number[] | null = null,
  norm: FFTNorm = null,
  out: NDArray | null = null
): NDArray {
  if (out !== null) {
    throw new Error('out parameter not yet supported');
  }
  return _fftn_impl(a, s, axes, true, norm);
}

/**
 * Compute the N-dimensional discrete Fourier Transform for real input.
 *
 * The output has Hermitian symmetry along the last transformed axis.
 * Output shape along the last axis: s[-1]//2 + 1
 *
 * @param a - Input array (real-valued)
 * @param s - Shape of the output along transformed axes
 * @param axes - Axes over which to compute the FFT. Default: all axes
 * @param norm - Normalization mode
 * @param out - Optional output array (not yet supported)
 * @returns Complex array with shape (..., s[-1]//2 + 1) along last transformed axis
 */
export function rfftn(
  a: NDArray,
  s: number[] | null = null,
  axes: number[] | null = null,
  norm: FFTNorm = null,
  out: NDArray | null = null
): NDArray {
  if (out !== null) {
    throw new Error('out parameter not yet supported');
  }

  const ndim = a.ndim;

  // Default axes: all axes
  let effectiveAxes = axes ?? [...Array(ndim).keys()];
  effectiveAxes = effectiveAxes.map(ax => ax < 0 ? ax + ndim : ax);

  // Default sizes
  const effectiveS = s ?? effectiveAxes.map(ax => a.shape[ax]);

  if (effectiveAxes.length === 0) {
    return a.astype(DType.Complex128);
  }

  // Apply real FFT on last axis
  const lastAxis = effectiveAxes[effectiveAxes.length - 1];
  const lastN = effectiveS[effectiveS.length - 1];
  let result = rfft(a, lastN, lastAxis, norm);

  // Apply complex FFT on remaining axes
  for (let i = effectiveAxes.length - 2; i >= 0; i--) {
    const ax = effectiveAxes[i];
    const n = effectiveS[i];
    result = fft(result, n, ax, norm);
  }

  return result;
}

/**
 * Compute the inverse of the N-dimensional FFT of real input.
 */
export function irfftn(
  a: NDArray,
  s: number[] | null = null,
  axes: number[] | null = null,
  norm: FFTNorm = null,
  out: NDArray | null = null
): NDArray {
  if (out !== null) {
    throw new Error('out parameter not yet supported');
  }

  const ndim = a.ndim;

  // Default axes: all axes
  let effectiveAxes = axes ?? [...Array(ndim).keys()];
  effectiveAxes = effectiveAxes.map(ax => ax < 0 ? ax + ndim : ax);

  // For irfftn, we need to determine s from input if not provided
  // The last axis of input has size s[-1]//2 + 1
  let effectiveS: number[];
  if (s !== null) {
    effectiveS = s;
  } else {
    effectiveS = effectiveAxes.map((ax, i) => {
      if (i === effectiveAxes.length - 1) {
        // Last axis: n//2 + 1 -> n = 2 * (size - 1)
        return 2 * (a.shape[ax] - 1);
      }
      return a.shape[ax];
    });
  }

  if (effectiveAxes.length === 0) {
    // Return real part
    return a.real();
  }

  // Apply complex inverse FFT on all but last axis
  let result = a;
  for (let i = 0; i < effectiveAxes.length - 1; i++) {
    const ax = effectiveAxes[i];
    const n = effectiveS[i];
    result = ifft(result, n, ax, norm);
  }

  // Apply real inverse FFT on last axis
  const lastAxis = effectiveAxes[effectiveAxes.length - 1];
  const lastN = effectiveS[effectiveS.length - 1];
  result = irfft(result, lastN, lastAxis, norm);

  return result;
}
```

---

### 14.4 Helper Functions

**File:** `src/ts/fft.ts` (additions)

```typescript
/* ============ Helper Functions ============ */

/**
 * Return the Discrete Fourier Transform sample frequencies.
 *
 * Returns an array of length n with the sample frequencies:
 * f = [0, 1, ..., n/2-1, -n/2, ..., -1] / (d*n)   if n is even
 * f = [0, 1, ..., (n-1)/2, -(n-1)/2, ..., -1] / (d*n)   if n is odd
 *
 * @param n - Window length
 * @param d - Sample spacing (default: 1.0)
 * @param device - Not used, for NumPy API compatibility
 * @returns Array of sample frequencies
 *
 * @example
 * const freqs = np.fft.fftfreq(8, 0.1);
 * // freqs = [0, 1.25, 2.5, 3.75, -5, -3.75, -2.5, -1.25]
 */
export function fftfreq(n: number, d: number = 1.0, device: null = null): NDArray {
  if (n <= 0) {
    throw new Error('n must be positive');
  }

  const freq = new Float64Array(n);
  const val = 1.0 / (n * d);

  const N = Math.floor((n - 1) / 2) + 1;  // Number of positive frequencies

  // Positive frequencies: 0, 1, ..., N-1
  for (let i = 0; i < N; i++) {
    freq[i] = i * val;
  }

  // Negative frequencies: -N, -(N-1), ..., -1 (if n is even, starts at -n/2)
  for (let i = N; i < n; i++) {
    freq[i] = (i - n) * val;
  }

  return NDArray.fromTypedArray(freq, [n], DType.Float64);
}

/**
 * Return the Discrete Fourier Transform sample frequencies for rfft.
 *
 * Returns an array of length n//2 + 1 with the positive sample frequencies:
 * f = [0, 1, ..., n//2] / (d*n)
 *
 * @param n - Window length
 * @param d - Sample spacing (default: 1.0)
 * @param device - Not used, for NumPy API compatibility
 * @returns Array of positive sample frequencies
 *
 * @example
 * const freqs = np.fft.rfftfreq(8, 0.1);
 * // freqs = [0, 1.25, 2.5, 3.75, 5]
 */
export function rfftfreq(n: number, d: number = 1.0, device: null = null): NDArray {
  if (n <= 0) {
    throw new Error('n must be positive');
  }

  const outSize = Math.floor(n / 2) + 1;
  const freq = new Float64Array(outSize);
  const val = 1.0 / (n * d);

  for (let i = 0; i < outSize; i++) {
    freq[i] = i * val;
  }

  return NDArray.fromTypedArray(freq, [outSize], DType.Float64);
}

/**
 * Shift the zero-frequency component to the center of the spectrum.
 *
 * This function swaps half-spaces for all axes listed (or all axes by default).
 *
 * @param x - Input array
 * @param axes - Axes over which to shift. Default: all axes
 * @returns Shifted array
 *
 * @example
 * const freqs = np.fft.fftfreq(4);
 * // freqs = [0, 0.25, -0.5, -0.25]
 * const shifted = np.fft.fftshift(freqs);
 * // shifted = [-0.5, -0.25, 0, 0.25]
 */
export function fftshift(x: NDArray, axes: number | number[] | null = null): NDArray {
  const ndim = x.ndim;

  // Default: all axes
  let effectiveAxes: number[];
  if (axes === null) {
    effectiveAxes = [...Array(ndim).keys()];
  } else if (typeof axes === 'number') {
    effectiveAxes = [axes];
  } else {
    effectiveAxes = axes;
  }

  // Normalize negative axes
  effectiveAxes = effectiveAxes.map(ax => ax < 0 ? ax + ndim : ax);

  // Compute shift amounts: n // 2 for each axis
  const shift = effectiveAxes.map(ax => Math.floor(x.shape[ax] / 2));

  // Use roll to perform the shift
  return x.roll(shift, effectiveAxes);
}

/**
 * The inverse of fftshift.
 *
 * Undoes the effect of fftshift.
 *
 * @param x - Input array
 * @param axes - Axes over which to shift. Default: all axes
 * @returns Shifted array
 *
 * @example
 * const x = np.fft.fftfreq(4);
 * const y = np.fft.fftshift(x);
 * const z = np.fft.ifftshift(y);
 * // z ≈ x
 */
export function ifftshift(x: NDArray, axes: number | number[] | null = null): NDArray {
  const ndim = x.ndim;

  // Default: all axes
  let effectiveAxes: number[];
  if (axes === null) {
    effectiveAxes = [...Array(ndim).keys()];
  } else if (typeof axes === 'number') {
    effectiveAxes = [axes];
  } else {
    effectiveAxes = axes;
  }

  // Normalize negative axes
  effectiveAxes = effectiveAxes.map(ax => ax < 0 ? ax + ndim : ax);

  // Compute shift amounts: -(n // 2) which is equivalent to (n - n//2) = (n+1)//2
  const shift = effectiveAxes.map(ax => -Math.floor(x.shape[ax] / 2));

  // Use roll to perform the shift
  return x.roll(shift, effectiveAxes);
}

/* ============ Module Export ============ */

export const fftModule = {
  // 1D transforms
  fft,
  ifft,
  rfft,
  irfft,
  hfft,
  ihfft,

  // 2D transforms
  fft2,
  ifft2,
  rfft2,
  irfft2,

  // N-D transforms
  fftn,
  ifftn,
  rfftn,
  irfftn,

  // Helper functions
  fftfreq,
  rfftfreq,
  fftshift,
  ifftshift,
};
```

---

## File Changes Summary

### New Files to Create

```
src/wasm/
├── fft.h              # FFT declarations (core algorithms)
└── fft.c              # FFT implementation (Cooley-Tukey, Bluestein, etc.)

src/ts/
└── fft.ts             # Full FFT module (TypeScript layer)
```

### Files to Modify

```
src/ts/types.ts
├── Add FFTNorm type
└── Add WASM function declarations for FFT

src/ts/index.ts
├── Export fft module
├── Export individual functions (fft, ifft, rfft, etc.)
└── Export helper functions (fftfreq, fftshift, etc.)

scripts/build-wasm.sh
├── Add fft.c to compilation
└── Add EXPORTED_FUNCTIONS for all FFT operations
```

---

## Build Script Updates

Add to `scripts/build-wasm.sh` source files:

```bash
emcc \
    "$SRC_DIR/ndarray.c" \
    "$SRC_DIR/dtype.c" \
    "$SRC_DIR/pairwise_sum.c" \
    "$SRC_DIR/broadcast.c" \
    "$SRC_DIR/indexing.c" \
    "$SRC_DIR/logic.c" \
    "$SRC_DIR/manipulation.c" \
    "$SRC_DIR/ufunc.c" \
    "$SRC_DIR/ufunc_unary.c" \
    "$SRC_DIR/ufunc_binary.c" \
    "$SRC_DIR/sorting.c" \
    "$SRC_DIR/searching.c" \
    "$SRC_DIR/statistics.c" \
    "$SRC_DIR/setops.c" \
    "$SRC_DIR/fft.c" \
    -o "$OUT_DIR/numjs.cjs" \
    ...
```

Add to EXPORTED_FUNCTIONS:

```bash
# FFT Core
"_fft_twiddle_factors",
"_fft_twiddle_factors_inv",
"_fft_bit_reverse",
"_fft_bit_reverse_permute",
"_fft_is_power_of_2",
"_fft_next_power_of_2",
"_fft_radix2",
"_fft_factorize",
"_fft_mixed_radix",
"_fft_bluestein",

# Real FFT
"_fft_rfft",
"_fft_irfft",
"_fft_hfft",
"_fft_ihfft",

# High-level (1D along axis)
"_fft_complex_axis",
"_fft_real_axis"
```

---

## WasmModule Interface Updates

Add to `src/ts/types.ts`:

```typescript
// FFT Core Operations
_fft_twiddle_factors(n: number, outPtr: number): void;
_fft_radix2(dataPtr: number, n: number, inverse: number): number;
_fft_bluestein(dataPtr: number, n: number, inverse: number, workPtr: number): number;
_fft_is_power_of_2(n: number): number;
_fft_next_power_of_2(n: number): number;

// Real FFT
_fft_rfft(realInPtr: number, outPtr: number, n: number, workPtr: number): number;
_fft_irfft(complexInPtr: number, outPtr: number, n: number, workPtr: number): number;

// High-level operations
_fft_complex_axis(arrPtr: number, axis: number, n: number, inverse: number, workPtr: number): number;
```

---

## Implementation Order

```
Phase 14.1: Core FFT Algorithm (Week 1-2)
├── Day 1-2: Twiddle factors and bit-reversal
│   ├── fft_twiddle_factors() - complex exponentials
│   ├── fft_bit_reverse() - index reversal
│   └── fft_bit_reverse_permute() - in-place permutation
│
├── Day 3-4: Radix-2 Cooley-Tukey FFT
│   ├── fft_radix2() - iterative implementation
│   ├── Forward and inverse transforms
│   └── Unit tests for power-of-2 sizes
│
├── Day 5-6: Mixed-Radix FFT
│   ├── fft_factorize() - prime factorization
│   ├── fft_radix3_butterfly() - factor 3
│   ├── fft_radix4_butterfly() - factor 4
│   └── fft_mixed_radix() - composite sizes
│
└── Day 7: Bluestein's Algorithm
    ├── fft_bluestein() - arbitrary size FFT
    ├── Convolution-based implementation
    └── Tests for prime sizes

Phase 14.2: 1D FFT Functions (Week 2-3)
├── Day 1-2: Complex FFT
│   ├── fft() - forward transform
│   ├── ifft() - inverse transform
│   ├── Axis handling
│   └── Normalization modes
│
├── Day 3-4: Real FFT
│   ├── fft_rfft() - WASM implementation
│   ├── fft_irfft() - WASM implementation
│   ├── rfft() - TypeScript wrapper
│   ├── irfft() - TypeScript wrapper
│   └── Hermitian symmetry tests
│
└── Day 5: Hermitian FFT
    ├── hfft() - Hermitian to real
    ├── ihfft() - real to Hermitian
    └── Round-trip tests

Phase 14.3: Multi-Dimensional FFT (Week 3-4)
├── Day 1-2: 2D FFT
│   ├── fft2() - 2D forward
│   ├── ifft2() - 2D inverse
│   ├── rfft2() - 2D real forward
│   └── irfft2() - 2D real inverse
│
└── Day 3-4: N-Dimensional FFT
    ├── fftn() - ND forward
    ├── ifftn() - ND inverse
    ├── rfftn() - ND real forward
    ├── irfftn() - ND real inverse
    └── Batch operation tests

Phase 14.4: Helper Functions & Polish (Week 4)
├── Day 1: Frequency Functions
│   ├── fftfreq() - sample frequencies
│   └── rfftfreq() - positive frequencies
│
├── Day 2: Shift Functions
│   ├── fftshift() - center zero-frequency
│   └── ifftshift() - undo shift
│
├── Day 3-4: Testing & Validation
│   ├── NumPy comparison tests
│   ├── Edge cases (empty, 1-element, large)
│   └── Performance benchmarks
│
└── Day 5: Documentation & Examples
    ├── API documentation
    ├── Usage examples
    └── Integration tests
```

---

## Verification Plan

After Phase 14 completion, verify:

```bash
# Build
npm run build

# Run all tests
npm test

# Phase 14 specific tests:

# Core FFT
✓ fft_radix2 produces correct DFT for power-of-2 sizes
✓ fft_bluestein handles arbitrary sizes (prime, composite)
✓ Forward and inverse are conjugates: ifft(fft(x)) ≈ x
✓ Twiddle factors satisfy W_n^n = 1

# 1D Transforms
✓ fft matches NumPy for various sizes and dtypes
✓ ifft matches NumPy for various normalization modes
✓ rfft produces n//2 + 1 complex outputs
✓ irfft reconstructs original real signal
✓ Round-trip: irfft(rfft(x)) ≈ x for real x
✓ hfft/ihfft handle Hermitian symmetry correctly

# Multi-Dimensional
✓ fft2 matches NumPy for 2D arrays
✓ fftn handles arbitrary dimensions
✓ rfftn/irfftn handle real input correctly
✓ Axis specification works correctly

# Helper Functions
✓ fftfreq produces correct frequency bins
✓ rfftfreq produces correct positive frequencies
✓ fftshift centers zero-frequency component
✓ ifftshift undoes fftshift: ifftshift(fftshift(x)) = x

# Edge Cases
✓ n=1 FFT returns input
✓ Empty arrays handled appropriately
✓ Large transforms (n > 10000) complete successfully
✓ All normalization modes ("backward", "ortho", "forward")
```

Generate NumPy comparison vectors:

```python
# tests/python/generate_fft_tests.py
import numpy as np
import json

np.random.seed(42)

tests = {
    "fft": [
        {"input": [1, 2, 3, 4], "n": None, "norm": None},
        {"input": [1, 2, 3, 4, 5], "n": None, "norm": None},  # Non-power-of-2
        {"input": [1, 2, 3, 4], "n": 8, "norm": None},  # Zero-padded
        {"input": [1, 2, 3, 4], "n": None, "norm": "ortho"},  # Orthonormal
    ],
    "rfft": [
        {"input": [1, 2, 3, 4], "n": None, "norm": None},
        {"input": [1, 2, 3, 4, 5], "n": None, "norm": None},
    ],
    "fft2": [
        {"input_shape": [4, 4], "s": None, "norm": None},
        {"input_shape": [3, 5], "s": None, "norm": None},
    ],
    "fftfreq": [
        {"n": 8, "d": 1.0},
        {"n": 7, "d": 0.1},
    ],
    "fftshift": [
        {"input": [0, 1, 2, 3, -4, -3, -2, -1]},
        {"input_shape": [4, 4]},
    ],
}

# Generate actual test data with NumPy
for test_type, cases in tests.items():
    for case in cases:
        if "input" in case:
            a = np.array(case["input"], dtype=np.float64)
        elif "input_shape" in case:
            a = np.random.rand(*case["input_shape"])
            case["input"] = a.tolist()

        if test_type == "fft":
            result = np.fft.fft(a, n=case.get("n"), norm=case.get("norm"))
            case["result_real"] = result.real.tolist()
            case["result_imag"] = result.imag.tolist()
        elif test_type == "rfft":
            result = np.fft.rfft(a, n=case.get("n"), norm=case.get("norm"))
            case["result_real"] = result.real.tolist()
            case["result_imag"] = result.imag.tolist()
        elif test_type == "fft2":
            result = np.fft.fft2(a, s=case.get("s"), norm=case.get("norm"))
            case["result_real"] = result.real.tolist()
            case["result_imag"] = result.imag.tolist()
        elif test_type == "fftfreq":
            result = np.fft.fftfreq(case["n"], case["d"])
            case["result"] = result.tolist()
        elif test_type == "fftshift":
            a = np.array(case["input"]) if "input" in case else np.random.rand(*case["input_shape"])
            result = np.fft.fftshift(a)
            case["result"] = result.tolist()

with open("tests/fixtures/fft_vectors.json", "w") as f:
    json.dump(tests, f, indent=2)
```

---

## Critical Dependencies

Phase 14 completion enables:

- **Signal Processing**: Filtering, spectral analysis, convolution via FFT
- **Image Processing**: Frequency domain filtering, compression
- **Phase 15 (numpy.random)**: Some distributions use FFT-based algorithms
- **Scientific Computing**: Solving PDEs, spectral methods

Phase 14 should be implemented after completing:
- Phase 4 (Ufuncs) - needed for element-wise complex operations
- Phase 5 (Array Manipulation) - needed for roll, reshape operations
- Complex number support in NDArray

---

## Performance Considerations

### Algorithm Selection Strategy

```
Size Classification:
├── Power of 2 (n = 2^k)
│   └── Use radix-2 Cooley-Tukey (fastest)
│
├── Smooth numbers (n = 2^a * 3^b * 5^c)
│   └── Use mixed-radix FFT
│
├── Products of small primes
│   └── Use mixed-radix with higher radices
│
└── Prime or large prime factors
    └── Use Bluestein's algorithm
```

### Memory Management

```
Workspace Requirements:
├── Radix-2: In-place (no extra memory)
├── Mixed-radix: O(n) workspace
├── Bluestein: O(m) where m = next_power_of_2(2n-1) ≈ 4n
└── Real FFT: O(n) for unpacking
```

### WASM-Specific Optimizations

```
Optimizations:
├── Pre-compute twiddle factors (cache for repeated sizes)
├── Use WASM SIMD when available
├── Batch small FFTs together
├── Minimize WASM ↔ JS data transfers
└── Consider WebWorkers for large transforms
```

---

## API Compatibility Notes

### NumPy Differences

```typescript
// NumPy: Returns ndarray
// fft.fft(a)

// NumJS: Returns NDArray with identical semantics
// fft(a)
```

### Normalization Convention

```
NumPy default ("backward"):
- Forward: no scaling
- Inverse: divide by n

This matches the mathematical DFT definition.
```

### Data Type Handling

```typescript
// Input promotion:
// - Real (float32/float64) → stays real for rfft
// - Integer types → promoted to float64
// - Complex stays complex

// Output types:
// - fft/ifft: complex128
// - rfft: complex128
// - irfft: float64
```

---

## Test Categories

### Correctness Tests

1. **Identity tests**: `ifft(fft(x)) ≈ x`
2. **Parseval's theorem**: `sum(|x|^2) = sum(|FFT(x)|^2) / n`
3. **Linearity**: `FFT(ax + by) = a*FFT(x) + b*FFT(y)`
4. **Shift theorem**: Time shift ↔ phase shift
5. **Convolution theorem**: `FFT(x * y) = FFT(x) · FFT(y)`

### Edge Cases

1. **n = 1**: Should return input unchanged
2. **n = 0 or empty**: Should handle gracefully
3. **Large n**: Performance and memory handling
4. **Non-contiguous arrays**: Proper stride handling
5. **Various dtypes**: float32, float64, complex64, complex128

### Numerical Accuracy

1. **Compare to NumPy**: Maximum relative error < 1e-10
2. **Test known transforms**: DC signal, single frequency, delta function
3. **Accumulation errors**: Long sequences of operations
