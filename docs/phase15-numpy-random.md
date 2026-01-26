# Phase 15: numpy.random Implementation Plan

Complete implementation roadmap for the NumJS-WASM random number generation module, providing NumPy-compatible random sampling with WebAssembly acceleration.

---

## ⚠️ Implementation Guideline

**The original NumPy source code is available in `/numpy`.**

Key reference files:
- `numpy/random/_generator.pyx` - Generator class with 40+ distribution methods (~5,086 lines)
- `numpy/random/bit_generator.pyx` - BitGenerator base class, SeedSequence (~1,500 lines)
- `numpy/random/_pcg64.pyx` - PCG64 BitGenerator wrapper
- `numpy/random/src/pcg64/pcg64.c` - PCG64 C implementation
- `numpy/random/src/distributions/distributions.c` - Distribution algorithms (Ziggurat, etc.)
- `numpy/random/src/distributions/ziggurat_constants.h` - Pre-computed Ziggurat tables

Implementation should follow NumPy's algorithms and error handling for consistency.

---

## Current State (Pre-Phase 15)

```
src/wasm/
├── ndarray.h/c        # Core NDArray with views, slicing
├── dtype.h/c          # DType system
├── broadcast.h/c      # Broadcasting
├── indexing.h/c       # Index operations
├── ufunc.h/c          # Ufunc infrastructure
├── statistics.h/c     # Statistics (mean, var, std)
├── sorting.h/c        # Sorting algorithms
└── blas.h/c           # BLAS operations (if Phase 13 complete)

src/ts/
├── NDArray.ts         # Core array class
├── types.ts           # Type definitions
├── ufunc.ts           # Ufunc TypeScript layer
├── statistics.ts      # Statistics functions
├── sorting.ts         # Sorting functions
└── index.ts           # Public exports
```

**Existing Infrastructure Used:**
- NDArray with all dtypes (Float32, Float64, Int32, Int64, etc.)
- Broadcasting for parameter arrays
- Statistics module for multivariate distributions (Cholesky from linalg)

---

## Phase 15 Dependency Tree

```
PHASE 15: NUMPY.RANDOM
│
├── 15.1 BitGenerator Infrastructure (C/WASM)
│   ├── 15.1.1 bitgen_t Interface
│   │   ├── next_uint64() → uint64
│   │   ├── next_uint32() → uint32
│   │   ├── next_double() → [0, 1) double
│   │   └── next_raw() → raw bits
│   │
│   ├── 15.1.2 PCG64 Implementation (Primary)
│   │   ├── pcg64_state structure
│   │   ├── pcg64_seed(state, seed, inc)
│   │   ├── pcg64_next64(state) → uint64
│   │   ├── pcg64_next32(state) → uint32
│   │   ├── pcg64_advance(state, delta)
│   │   └── pcg64_jump(state)
│   │
│   ├── 15.1.3 Xoshiro256** Implementation (Alternative)
│   │   ├── xoshiro256_state structure
│   │   ├── xoshiro256_seed(state, seed)
│   │   ├── xoshiro256_next64(state) → uint64
│   │   └── xoshiro256_jump(state)
│   │
│   └── 15.1.4 SeedSequence
│       ├── seed_seq_init(entropy, spawn_key)
│       ├── seed_seq_generate(n_words) → uint32[]
│       └── seed_seq_spawn(n_children)
│
│   Dependencies: None (standalone)
│
├── 15.2 Distribution Algorithms (C/WASM)
│   ├── 15.2.1 Uniform Distributions
│   │   ├── random_uniform() → [0, 1)
│   │   ├── random_uniform_range(low, high)
│   │   └── random_integers(low, high)
│   │
│   ├── 15.2.2 Normal/Gaussian (Ziggurat Algorithm)
│   │   ├── random_standard_normal() → N(0, 1)
│   │   ├── random_normal(loc, scale)
│   │   └── Ziggurat tables (256 entries)
│   │
│   ├── 15.2.3 Exponential (Ziggurat Algorithm)
│   │   ├── random_standard_exponential() → Exp(1)
│   │   ├── random_exponential(scale)
│   │   └── Ziggurat tables (256 entries)
│   │
│   ├── 15.2.4 Gamma Distribution
│   │   ├── random_standard_gamma(shape)
│   │   ├── random_gamma(shape, scale)
│   │   └── Marsaglia-Tsang method
│   │
│   ├── 15.2.5 Other Continuous Distributions
│   │   ├── random_beta(a, b)
│   │   ├── random_chisquare(df)
│   │   ├── random_f(dfnum, dfden)
│   │   ├── random_t(df)
│   │   ├── random_cauchy()
│   │   ├── random_pareto(a)
│   │   ├── random_weibull(a)
│   │   ├── random_power(a)
│   │   ├── random_laplace(loc, scale)
│   │   ├── random_gumbel(loc, scale)
│   │   ├── random_logistic(loc, scale)
│   │   ├── random_lognormal(mean, sigma)
│   │   ├── random_rayleigh(scale)
│   │   ├── random_wald(mean, scale)
│   │   ├── random_triangular(left, mode, right)
│   │   └── random_vonmises(mu, kappa)
│   │
│   ├── 15.2.6 Discrete Distributions
│   │   ├── random_binomial(n, p)
│   │   ├── random_negative_binomial(n, p)
│   │   ├── random_poisson(lam)
│   │   ├── random_geometric(p)
│   │   ├── random_hypergeometric(ngood, nbad, nsample)
│   │   ├── random_logseries(p)
│   │   └── random_zipf(a)
│   │
│   └── 15.2.7 Bounded Integers
│       ├── random_bounded_uint64(off, rng, mask)
│       ├── random_bounded_uint32(off, rng, mask)
│       └── Lemire's algorithm for unbiased bounded integers
│
│   Dependencies: 15.1.* (BitGenerator)
│
├── 15.3 Generator Class (TypeScript)
│   ├── 15.3.1 Core Generator
│   │   ├── constructor(bitGenerator)
│   │   ├── spawn(n_children) → Generator[]
│   │   ├── random(size?, dtype?) → uniform [0, 1)
│   │   └── bytes(length) → Uint8Array
│   │
│   ├── 15.3.2 Integer Methods
│   │   ├── integers(low, high?, size?, dtype?, endpoint?)
│   │   └── choice(a, size?, replace?, p?, axis?, shuffle?)
│   │
│   ├── 15.3.3 Continuous Distributions
│   │   ├── uniform(low, high, size?)
│   │   ├── standard_normal(size?, dtype?)
│   │   ├── normal(loc, scale, size?)
│   │   ├── standard_exponential(size?, dtype?, method?)
│   │   ├── exponential(scale, size?)
│   │   ├── standard_gamma(shape, size?, dtype?)
│   │   ├── gamma(shape, scale, size?)
│   │   ├── beta(a, b, size?)
│   │   ├── chisquare(df, size?)
│   │   ├── f(dfnum, dfden, size?)
│   │   ├── standard_t(df, size?)
│   │   ├── standard_cauchy(size?)
│   │   ├── pareto(a, size?)
│   │   ├── weibull(a, size?)
│   │   ├── power(a, size?)
│   │   ├── laplace(loc, scale, size?)
│   │   ├── gumbel(loc, scale, size?)
│   │   ├── logistic(loc, scale, size?)
│   │   ├── lognormal(mean, sigma, size?)
│   │   ├── rayleigh(scale, size?)
│   │   ├── wald(mean, scale, size?)
│   │   ├── triangular(left, mode, right, size?)
│   │   └── vonmises(mu, kappa, size?)
│   │
│   ├── 15.3.4 Discrete Distributions
│   │   ├── binomial(n, p, size?)
│   │   ├── negative_binomial(n, p, size?)
│   │   ├── poisson(lam, size?)
│   │   ├── geometric(p, size?)
│   │   ├── hypergeometric(ngood, nbad, nsample, size?)
│   │   ├── logseries(p, size?)
│   │   └── zipf(a, size?)
│   │
│   ├── 15.3.5 Multivariate Distributions
│   │   ├── multivariate_normal(mean, cov, size?, method?)
│   │   ├── multinomial(n, pvals, size?)
│   │   ├── dirichlet(alpha, size?)
│   │   └── multivariate_hypergeometric(colors, nsample, size?, method?)
│   │
│   └── 15.3.6 Permutation Methods
│       ├── shuffle(x, axis?)
│       ├── permutation(x, axis?)
│       └── permuted(x, axis?, out?)
│
│   Dependencies: 15.2.* (Distributions)
│
├── 15.4 BitGenerator Classes (TypeScript)
│   ├── 15.4.1 BitGenerator Base Class
│   │   ├── abstract next_uint64() → bigint
│   │   ├── abstract next_uint32() → number
│   │   ├── abstract next_double() → number
│   │   ├── getState() → object
│   │   ├── setState(state)
│   │   └── spawn(n_children) → BitGenerator[]
│   │
│   ├── 15.4.2 PCG64 Class
│   │   ├── constructor(seed?)
│   │   ├── jumped(jumps?) → PCG64
│   │   └── advance(delta)
│   │
│   ├── 15.4.3 Xoshiro256 Class
│   │   ├── constructor(seed?)
│   │   └── jumped(jumps?) → Xoshiro256
│   │
│   └── 15.4.4 SeedSequence Class
│       ├── constructor(entropy?, spawn_key?, pool_size?)
│       ├── generate_state(n_words, dtype?) → Uint32Array
│       └── spawn(n_children) → SeedSequence[]
│
│   Dependencies: 15.1.* (C implementations)
│
└── 15.5 Module API & Utilities
    ├── 15.5.1 default_rng(seed?) → Generator
    ├── 15.5.2 Legacy RandomState (optional compatibility)
    └── 15.5.3 Convenience functions (random, randn, randint)

    Dependencies: 15.3.*, 15.4.*
```

---

## Detailed Implementation Specifications

### 15.1 BitGenerator Infrastructure

#### 15.1.1 bitgen_t Interface

**File:** `src/wasm/random/bitgen.h` (new file)

```c
#ifndef NUMJS_BITGEN_H
#define NUMJS_BITGEN_H

#include <stdint.h>

/**
 * BitGenerator interface structure.
 * All bit generators must implement these function pointers.
 */
typedef struct {
    void *state;                                    /* Generator-specific state */
    uint64_t (*next_uint64)(void *state);          /* Generate 64-bit unsigned */
    uint32_t (*next_uint32)(void *state);          /* Generate 32-bit unsigned */
    double (*next_double)(void *state);            /* Generate [0, 1) double */
    uint64_t (*next_raw)(void *state);             /* Raw output (same as uint64) */
} bitgen_t;

/**
 * Convert uint64 to double in [0, 1).
 * Uses upper 53 bits for IEEE 754 double precision.
 */
static inline double uint64_to_double(uint64_t x) {
    return (x >> 11) * (1.0 / 9007199254740992.0);  /* 2^53 */
}

/**
 * Convert uint32 to float in [0, 1).
 * Uses upper 24 bits for IEEE 754 single precision.
 */
static inline float uint32_to_float(uint32_t x) {
    return (x >> 8) * (1.0f / 16777216.0f);  /* 2^24 */
}

#endif /* NUMJS_BITGEN_H */
```

#### 15.1.2 PCG64 Implementation

**File:** `src/wasm/random/pcg64.h` (new file)

```c
#ifndef NUMJS_PCG64_H
#define NUMJS_PCG64_H

#include "bitgen.h"

/**
 * 128-bit unsigned integer for PCG state.
 * Uses struct for platforms without native __uint128_t.
 */
typedef struct {
    uint64_t high;
    uint64_t low;
} pcg128_t;

/**
 * PCG64 state structure.
 * Uses Linear Congruential Generator with permutation output.
 */
typedef struct {
    pcg128_t state;    /* Current state */
    pcg128_t inc;      /* Increment (stream selector) */
} pcg_state_t;

/**
 * PCG64 wrapper with buffered 32-bit output.
 */
typedef struct {
    pcg_state_t *pcg_state;
    int has_uint32;     /* Flag: buffered 32-bit value available */
    uint32_t uinteger;  /* Buffered 32-bit value */
} pcg64_state;

/* PCG64 constants */
#define PCG_DEFAULT_MULTIPLIER_HIGH 2549297995355413924ULL
#define PCG_DEFAULT_MULTIPLIER_LOW  4865540595714422341ULL

/* Function declarations */
void pcg64_seed(pcg64_state *state, uint64_t seed_high, uint64_t seed_low,
                uint64_t inc_high, uint64_t inc_low);
uint64_t pcg64_next64(pcg64_state *state);
uint32_t pcg64_next32(pcg64_state *state);
double pcg64_next_double(pcg64_state *state);
void pcg64_advance(pcg64_state *state, uint64_t delta_high, uint64_t delta_low);
void pcg64_get_state(pcg64_state *state, uint64_t *out);
void pcg64_set_state(pcg64_state *state, uint64_t *in);

/* Fill bitgen_t interface */
void pcg64_init_bitgen(bitgen_t *bitgen, pcg64_state *state);

#endif /* NUMJS_PCG64_H */
```

**File:** `src/wasm/random/pcg64.c` (new file)

```c
#include "pcg64.h"
#include <stdlib.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/* 128-bit multiplication */
static inline pcg128_t pcg128_mult(pcg128_t a, pcg128_t b) {
    pcg128_t result;
    uint64_t lo_lo = (a.low & 0xFFFFFFFF) * (b.low & 0xFFFFFFFF);
    uint64_t hi_lo = (a.low >> 32) * (b.low & 0xFFFFFFFF);
    uint64_t lo_hi = (a.low & 0xFFFFFFFF) * (b.low >> 32);
    uint64_t hi_hi = (a.low >> 32) * (b.low >> 32);

    uint64_t cross = (lo_lo >> 32) + (hi_lo & 0xFFFFFFFF) + lo_hi;
    result.high = (hi_lo >> 32) + (cross >> 32) + hi_hi +
                  a.high * b.low + a.low * b.high;
    result.low = (cross << 32) | (lo_lo & 0xFFFFFFFF);
    return result;
}

/* 128-bit addition */
static inline pcg128_t pcg128_add(pcg128_t a, pcg128_t b) {
    pcg128_t result;
    result.low = a.low + b.low;
    result.high = a.high + b.high + (result.low < a.low);
    return result;
}

/* PCG output function: XSL-RR (xorshift low, random rotate) */
static inline uint64_t pcg_output_xsl_rr_128_64(pcg128_t state) {
    uint64_t xorshifted = state.high ^ state.low;
    int rot = state.high >> 58;
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 63));
}

/* Single step of the LCG */
static inline void pcg_step(pcg_state_t *state) {
    pcg128_t mult = {PCG_DEFAULT_MULTIPLIER_HIGH, PCG_DEFAULT_MULTIPLIER_LOW};
    state->state = pcg128_add(pcg128_mult(state->state, mult), state->inc);
}

EXPORT void pcg64_seed(pcg64_state *state, uint64_t seed_high, uint64_t seed_low,
                        uint64_t inc_high, uint64_t inc_low) {
    if (!state->pcg_state) {
        state->pcg_state = (pcg_state_t *)malloc(sizeof(pcg_state_t));
    }

    /* Ensure increment is odd */
    state->pcg_state->inc.high = inc_high;
    state->pcg_state->inc.low = (inc_low << 1) | 1;

    /* Initialize state */
    state->pcg_state->state.high = 0;
    state->pcg_state->state.low = 0;
    pcg_step(state->pcg_state);

    state->pcg_state->state.high += seed_high;
    state->pcg_state->state.low += seed_low;
    pcg_step(state->pcg_state);

    state->has_uint32 = 0;
    state->uinteger = 0;
}

EXPORT uint64_t pcg64_next64(pcg64_state *state) {
    uint64_t result = pcg_output_xsl_rr_128_64(state->pcg_state->state);
    pcg_step(state->pcg_state);
    return result;
}

EXPORT uint32_t pcg64_next32(pcg64_state *state) {
    if (state->has_uint32) {
        state->has_uint32 = 0;
        return state->uinteger;
    }
    uint64_t next = pcg64_next64(state);
    state->has_uint32 = 1;
    state->uinteger = (uint32_t)(next >> 32);
    return (uint32_t)next;
}

EXPORT double pcg64_next_double(pcg64_state *state) {
    return uint64_to_double(pcg64_next64(state));
}

/* Advance by delta steps using O(log n) algorithm */
EXPORT void pcg64_advance(pcg64_state *state, uint64_t delta_high, uint64_t delta_low) {
    pcg128_t delta = {delta_high, delta_low};
    pcg128_t cur_mult = {PCG_DEFAULT_MULTIPLIER_HIGH, PCG_DEFAULT_MULTIPLIER_LOW};
    pcg128_t cur_plus = state->pcg_state->inc;
    pcg128_t acc_mult = {0, 1};
    pcg128_t acc_plus = {0, 0};

    while (delta.high > 0 || delta.low > 0) {
        if (delta.low & 1) {
            acc_mult = pcg128_mult(acc_mult, cur_mult);
            acc_plus = pcg128_add(pcg128_mult(acc_plus, cur_mult), cur_plus);
        }
        pcg128_t one = {0, 1};
        cur_plus = pcg128_mult(pcg128_add(cur_mult, one), cur_plus);
        cur_mult = pcg128_mult(cur_mult, cur_mult);

        /* delta >>= 1 */
        delta.low = (delta.low >> 1) | (delta.high << 63);
        delta.high >>= 1;
    }

    state->pcg_state->state = pcg128_add(
        pcg128_mult(acc_mult, state->pcg_state->state),
        acc_plus
    );
}

EXPORT void pcg64_get_state(pcg64_state *state, uint64_t *out) {
    out[0] = state->pcg_state->state.high;
    out[1] = state->pcg_state->state.low;
    out[2] = state->pcg_state->inc.high;
    out[3] = state->pcg_state->inc.low;
    out[4] = state->has_uint32;
    out[5] = state->uinteger;
}

EXPORT void pcg64_set_state(pcg64_state *state, uint64_t *in) {
    if (!state->pcg_state) {
        state->pcg_state = (pcg_state_t *)malloc(sizeof(pcg_state_t));
    }
    state->pcg_state->state.high = in[0];
    state->pcg_state->state.low = in[1];
    state->pcg_state->inc.high = in[2];
    state->pcg_state->inc.low = in[3];
    state->has_uint32 = (int)in[4];
    state->uinteger = (uint32_t)in[5];
}

/* Initialize bitgen_t interface */
static uint64_t _pcg64_next64_wrapper(void *state) {
    return pcg64_next64((pcg64_state *)state);
}

static uint32_t _pcg64_next32_wrapper(void *state) {
    return pcg64_next32((pcg64_state *)state);
}

static double _pcg64_next_double_wrapper(void *state) {
    return pcg64_next_double((pcg64_state *)state);
}

EXPORT void pcg64_init_bitgen(bitgen_t *bitgen, pcg64_state *state) {
    bitgen->state = state;
    bitgen->next_uint64 = _pcg64_next64_wrapper;
    bitgen->next_uint32 = _pcg64_next32_wrapper;
    bitgen->next_double = _pcg64_next_double_wrapper;
    bitgen->next_raw = _pcg64_next64_wrapper;
}
```

#### 15.1.3 SeedSequence Implementation

**File:** `src/wasm/random/seedseq.h` (new file)

```c
#ifndef NUMJS_SEEDSEQ_H
#define NUMJS_SEEDSEQ_H

#include <stdint.h>

/**
 * SeedSequence for reproducible entropy management.
 * Based on Melissa E. O'Neill's seed_seq design.
 */
typedef struct {
    uint32_t *pool;       /* Entropy pool */
    int32_t pool_size;    /* Pool size in uint32 words */
    int32_t n_children;   /* Number of children spawned */
} seed_sequence_t;

/* Mixing constants */
#define SEED_SEQ_INIT_A  0x43b0d7e5u
#define SEED_SEQ_MULT_A  0x931e8875u
#define SEED_SEQ_INIT_B  0x8b51f9ddu
#define SEED_SEQ_MULT_B  0x58f38dadu

/* Function declarations */
void seed_seq_init(seed_sequence_t *seq, uint32_t *entropy, int32_t entropy_len,
                   uint32_t *spawn_key, int32_t spawn_key_len, int32_t pool_size);
void seed_seq_generate(seed_sequence_t *seq, uint32_t *out, int32_t n_words);
void seed_seq_free(seed_sequence_t *seq);

#endif /* NUMJS_SEEDSEQ_H */
```

**File:** `src/wasm/random/seedseq.c` (new file)

```c
#include "seedseq.h"
#include <stdlib.h>
#include <string.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/* Hash function for mixing */
static inline uint32_t hashmix(uint32_t value, uint32_t *hash) {
    value ^= *hash;
    *hash *= SEED_SEQ_MULT_A;
    value *= SEED_SEQ_MULT_A;
    value ^= value >> 16;
    return value;
}

EXPORT void seed_seq_init(seed_sequence_t *seq, uint32_t *entropy, int32_t entropy_len,
                           uint32_t *spawn_key, int32_t spawn_key_len, int32_t pool_size) {
    seq->pool_size = pool_size > 0 ? pool_size : 4;
    seq->pool = (uint32_t *)calloc(seq->pool_size, sizeof(uint32_t));
    seq->n_children = 0;

    /* Mix entropy into pool */
    uint32_t hash = SEED_SEQ_INIT_A;
    for (int32_t i = 0; i < entropy_len; i++) {
        int32_t idx = i % seq->pool_size;
        seq->pool[idx] ^= hashmix(entropy[i], &hash);
    }

    /* Mix spawn_key into pool */
    hash = SEED_SEQ_INIT_B;
    for (int32_t i = 0; i < spawn_key_len; i++) {
        int32_t idx = i % seq->pool_size;
        seq->pool[idx] ^= hashmix(spawn_key[i], &hash);
    }

    /* Final mixing pass to ensure all bits affect all others */
    for (int32_t i = 0; i < seq->pool_size; i++) {
        hash = SEED_SEQ_INIT_A;
        for (int32_t j = 0; j < seq->pool_size; j++) {
            seq->pool[j] = hashmix(seq->pool[j], &hash);
        }
    }
}

EXPORT void seed_seq_generate(seed_sequence_t *seq, uint32_t *out, int32_t n_words) {
    /* Generate output words using pool as seed */
    uint32_t hash = SEED_SEQ_INIT_A;

    for (int32_t i = 0; i < n_words; i++) {
        uint32_t value = seq->pool[i % seq->pool_size];
        value ^= (uint32_t)i;
        out[i] = hashmix(value, &hash);
    }
}

EXPORT void seed_seq_free(seed_sequence_t *seq) {
    if (seq->pool) {
        free(seq->pool);
        seq->pool = NULL;
    }
}
```

---

### 15.2 Distribution Algorithms

#### 15.2.1 Ziggurat Constants

**File:** `src/wasm/random/ziggurat_constants.h` (new file)

```c
#ifndef NUMJS_ZIGGURAT_CONSTANTS_H
#define NUMJS_ZIGGURAT_CONSTANTS_H

/*
 * Ziggurat constants for normal and exponential distributions.
 * Pre-computed tables with 256 entries for efficiency.
 * ~98.9% of samples return on first try.
 */

/* Normal distribution constants */
#define ZIGGURAT_NOR_R      3.6541528853610088
#define ZIGGURAT_NOR_INV_R  0.27366123732975828

/* Exponential distribution constants */
#define ZIGGURAT_EXP_R      7.6971174701310497

/* Normal distribution tables (256 entries each) */
extern const double ki_double[256];  /* Threshold table */
extern const double wi_double[256];  /* Width table */
extern const double fi_double[256];  /* f(x) table */

extern const float ki_float[256];
extern const float wi_float[256];
extern const float fi_float[256];

/* Exponential distribution tables */
extern const double ke_double[256];
extern const double we_double[256];
extern const double fe_double[256];

extern const float ke_float[256];
extern const float we_float[256];
extern const float fe_float[256];

#endif /* NUMJS_ZIGGURAT_CONSTANTS_H */
```

#### 15.2.2 Distribution Implementation

**File:** `src/wasm/random/distributions.h` (new file)

```c
#ifndef NUMJS_DISTRIBUTIONS_H
#define NUMJS_DISTRIBUTIONS_H

#include "bitgen.h"
#include <stdint.h>

/* ============ Uniform Distributions ============ */

double random_standard_uniform(bitgen_t *state);
float random_standard_uniform_f(bitgen_t *state);
void random_uniform_fill(bitgen_t *state, int64_t cnt, double *out);
void random_uniform_fill_f(bitgen_t *state, int64_t cnt, float *out);

/* ============ Normal Distribution (Ziggurat) ============ */

double random_standard_normal(bitgen_t *state);
float random_standard_normal_f(bitgen_t *state);
void random_standard_normal_fill(bitgen_t *state, int64_t cnt, double *out);
void random_standard_normal_fill_f(bitgen_t *state, int64_t cnt, float *out);

double random_normal(bitgen_t *state, double loc, double scale);

/* ============ Exponential Distribution (Ziggurat) ============ */

double random_standard_exponential(bitgen_t *state);
float random_standard_exponential_f(bitgen_t *state);
void random_standard_exponential_fill(bitgen_t *state, int64_t cnt, double *out);
void random_standard_exponential_fill_f(bitgen_t *state, int64_t cnt, float *out);

/* Inverse method (for alternative) */
void random_standard_exponential_inv_fill(bitgen_t *state, int64_t cnt, double *out);

double random_exponential(bitgen_t *state, double scale);

/* ============ Gamma Distribution ============ */

double random_standard_gamma(bitgen_t *state, double shape);
float random_standard_gamma_f(bitgen_t *state, float shape);
void random_standard_gamma_fill(bitgen_t *state, int64_t cnt, double shape, double *out);

double random_gamma(bitgen_t *state, double shape, double scale);

/* ============ Beta Distribution ============ */

double random_beta(bitgen_t *state, double a, double b);
void random_beta_fill(bitgen_t *state, int64_t cnt, double a, double b, double *out);

/* ============ Chi-Square Distribution ============ */

double random_chisquare(bitgen_t *state, double df);
double random_noncentral_chisquare(bitgen_t *state, double df, double nonc);

/* ============ F Distribution ============ */

double random_f(bitgen_t *state, double dfnum, double dfden);
double random_noncentral_f(bitgen_t *state, double dfnum, double dfden, double nonc);

/* ============ Student's t Distribution ============ */

double random_standard_t(bitgen_t *state, double df);

/* ============ Other Continuous Distributions ============ */

double random_standard_cauchy(bitgen_t *state);
double random_pareto(bitgen_t *state, double a);
double random_weibull(bitgen_t *state, double a);
double random_power(bitgen_t *state, double a);
double random_laplace(bitgen_t *state, double loc, double scale);
double random_gumbel(bitgen_t *state, double loc, double scale);
double random_logistic(bitgen_t *state, double loc, double scale);
double random_lognormal(bitgen_t *state, double mean, double sigma);
double random_rayleigh(bitgen_t *state, double scale);
double random_wald(bitgen_t *state, double mean, double scale);
double random_triangular(bitgen_t *state, double left, double mode, double right);
double random_vonmises(bitgen_t *state, double mu, double kappa);

/* ============ Discrete Distributions ============ */

int64_t random_binomial(bitgen_t *state, double p, int64_t n);
int64_t random_negative_binomial(bitgen_t *state, double n, double p);
int64_t random_poisson(bitgen_t *state, double lam);
int64_t random_geometric(bitgen_t *state, double p);
int64_t random_hypergeometric(bitgen_t *state, int64_t ngood, int64_t nbad, int64_t nsample);
int64_t random_logseries(bitgen_t *state, double p);
int64_t random_zipf(bitgen_t *state, double a);

/* ============ Bounded Integers ============ */

uint64_t random_bounded_uint64(bitgen_t *state, uint64_t off, uint64_t rng, uint64_t mask);
uint32_t random_bounded_uint32(bitgen_t *state, uint32_t off, uint32_t rng, uint32_t mask);
int64_t random_integers(bitgen_t *state, int64_t low, int64_t high);

/* ============ Binomial State (for caching) ============ */

typedef struct {
    int initialized;
    int64_t n;
    double p;
    double q;
    double r;
    double p_ratio;
    double psave;
    int64_t nsave;
    /* ... additional cached values ... */
} binomial_t;

int64_t random_binomial_btpe(bitgen_t *state, binomial_t *binomial, int64_t n, double p);
int64_t random_binomial_inversion(bitgen_t *state, int64_t n, double p);

#endif /* NUMJS_DISTRIBUTIONS_H */
```

**File:** `src/wasm/random/distributions.c` (new file - partial)

```c
#include "distributions.h"
#include "ziggurat_constants.h"
#include <math.h>
#include <float.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/* ============ Uniform Distributions ============ */

EXPORT double random_standard_uniform(bitgen_t *state) {
    return state->next_double(state->state);
}

EXPORT float random_standard_uniform_f(bitgen_t *state) {
    return (state->next_uint32(state->state) >> 8) * (1.0f / 16777216.0f);
}

EXPORT void random_uniform_fill(bitgen_t *state, int64_t cnt, double *out) {
    for (int64_t i = 0; i < cnt; i++) {
        out[i] = state->next_double(state->state);
    }
}

/* ============ Normal Distribution (Ziggurat) ============ */

EXPORT double random_standard_normal(bitgen_t *state) {
    uint64_t r;
    int sign;
    uint64_t rabs;
    int idx;
    double x, xx, yy;

    for (;;) {
        r = state->next_uint64(state->state);
        idx = r & 0xff;
        r >>= 8;
        sign = r & 0x1;
        rabs = (r >> 1) & 0x000fffffffffffffULL;
        x = rabs * wi_double[idx];

        if (sign & 0x1) x = -x;

        if (rabs < ki_double[idx]) {
            return x;  /* 99.3% of the time return here */
        }

        if (idx == 0) {
            /* Tail */
            for (;;) {
                xx = -ZIGGURAT_NOR_INV_R * log(1.0 - state->next_double(state->state));
                yy = -log(1.0 - state->next_double(state->state));
                if (yy + yy > xx * xx) {
                    return ((rabs >> 8) & 0x1) ? -(ZIGGURAT_NOR_R + xx)
                                               : ZIGGURAT_NOR_R + xx;
                }
            }
        } else {
            if (((fi_double[idx - 1] - fi_double[idx]) * state->next_double(state->state) +
                 fi_double[idx]) < exp(-0.5 * x * x)) {
                return x;
            }
        }
    }
}

EXPORT void random_standard_normal_fill(bitgen_t *state, int64_t cnt, double *out) {
    for (int64_t i = 0; i < cnt; i++) {
        out[i] = random_standard_normal(state);
    }
}

EXPORT double random_normal(bitgen_t *state, double loc, double scale) {
    return loc + scale * random_standard_normal(state);
}

/* ============ Exponential Distribution (Ziggurat) ============ */

static double standard_exponential_unlikely(bitgen_t *state, uint8_t idx, double x) {
    if (idx == 0) {
        return ZIGGURAT_EXP_R - log(1.0 - state->next_double(state->state));
    } else if ((fe_double[idx - 1] - fe_double[idx]) * state->next_double(state->state) +
               fe_double[idx] < exp(-x)) {
        return x;
    } else {
        return random_standard_exponential(state);
    }
}

EXPORT double random_standard_exponential(bitgen_t *state) {
    uint64_t ri;
    uint8_t idx;
    double x;

    ri = state->next_uint64(state->state);
    ri >>= 3;
    idx = ri & 0xFF;
    ri >>= 8;
    x = ri * we_double[idx];

    if (ri < ke_double[idx]) {
        return x;  /* 98.9% of the time return here */
    }
    return standard_exponential_unlikely(state, idx, x);
}

EXPORT void random_standard_exponential_fill(bitgen_t *state, int64_t cnt, double *out) {
    for (int64_t i = 0; i < cnt; i++) {
        out[i] = random_standard_exponential(state);
    }
}

EXPORT double random_exponential(bitgen_t *state, double scale) {
    return scale * random_standard_exponential(state);
}

/* ============ Gamma Distribution (Marsaglia-Tsang) ============ */

EXPORT double random_standard_gamma(bitgen_t *state, double shape) {
    double b, c, U, V, X, Y;

    if (shape == 1.0) {
        return random_standard_exponential(state);
    } else if (shape == 0.0) {
        return 0.0;
    } else if (shape < 1.0) {
        for (;;) {
            U = state->next_double(state->state);
            V = random_standard_exponential(state);
            if (U <= 1.0 - shape) {
                X = pow(U, 1.0 / shape);
                if (X <= V) return X;
            } else {
                Y = -log((1.0 - U) / shape);
                X = pow(1.0 - shape + shape * Y, 1.0 / shape);
                if (X <= V + Y) return X;
            }
        }
    } else {
        /* Marsaglia-Tsang method for shape > 1 */
        b = shape - 1.0 / 3.0;
        c = 1.0 / sqrt(9.0 * b);
        for (;;) {
            do {
                X = random_standard_normal(state);
                V = 1.0 + c * X;
            } while (V <= 0.0);

            V = V * V * V;
            U = state->next_double(state->state);
            if (U < 1.0 - 0.0331 * (X * X) * (X * X)) {
                return b * V;
            }
            if (log(U) < 0.5 * X * X + b * (1.0 - V + log(V))) {
                return b * V;
            }
        }
    }
}

EXPORT double random_gamma(bitgen_t *state, double shape, double scale) {
    return scale * random_standard_gamma(state, shape);
}

/* ============ Beta Distribution ============ */

EXPORT double random_beta(bitgen_t *state, double a, double b) {
    double Ga, Gb;

    if (a <= 1.0 && b <= 1.0) {
        /* Use Johnk's algorithm for small parameters */
        double U, V, X, Y;
        for (;;) {
            U = state->next_double(state->state);
            V = state->next_double(state->state);
            X = pow(U, 1.0 / a);
            Y = pow(V, 1.0 / b);
            if (X + Y <= 1.0) {
                if (X + Y > 0) {
                    return X / (X + Y);
                } else {
                    double logX = log(U) / a;
                    double logY = log(V) / b;
                    double logM = logX > logY ? logX : logY;
                    logX -= logM;
                    logY -= logM;
                    return exp(logX - log(exp(logX) + exp(logY)));
                }
            }
        }
    } else {
        /* Use gamma variates */
        Ga = random_standard_gamma(state, a);
        Gb = random_standard_gamma(state, b);
        return Ga / (Ga + Gb);
    }
}

/* ============ Chi-Square Distribution ============ */

EXPORT double random_chisquare(bitgen_t *state, double df) {
    return 2.0 * random_standard_gamma(state, df / 2.0);
}

/* ============ Student's t Distribution ============ */

EXPORT double random_standard_t(bitgen_t *state, double df) {
    double N = random_standard_normal(state);
    double G = random_standard_gamma(state, df / 2.0);
    return sqrt(df / 2.0) * N / sqrt(G);
}

/* ============ Other Continuous Distributions ============ */

EXPORT double random_standard_cauchy(bitgen_t *state) {
    return random_standard_normal(state) / random_standard_normal(state);
}

EXPORT double random_pareto(bitgen_t *state, double a) {
    return exp(random_standard_exponential(state) / a) - 1.0;
}

EXPORT double random_weibull(bitgen_t *state, double a) {
    if (a == 0.0) return 0.0;
    return pow(random_standard_exponential(state), 1.0 / a);
}

EXPORT double random_power(bitgen_t *state, double a) {
    return pow(1.0 - exp(-random_standard_exponential(state)), 1.0 / a);
}

EXPORT double random_laplace(bitgen_t *state, double loc, double scale) {
    double U = state->next_double(state->state);
    if (U < 0.5) {
        U = loc + scale * log(2.0 * U);
    } else {
        U = loc - scale * log(2.0 * (1.0 - U));
    }
    return U;
}

EXPORT double random_gumbel(bitgen_t *state, double loc, double scale) {
    double U = 1.0 - state->next_double(state->state);
    return loc - scale * log(-log(U));
}

EXPORT double random_logistic(bitgen_t *state, double loc, double scale) {
    double U = state->next_double(state->state);
    return loc + scale * log(U / (1.0 - U));
}

EXPORT double random_lognormal(bitgen_t *state, double mean, double sigma) {
    return exp(random_normal(state, mean, sigma));
}

EXPORT double random_rayleigh(bitgen_t *state, double scale) {
    return scale * sqrt(2.0 * random_standard_exponential(state));
}

EXPORT double random_wald(bitgen_t *state, double mean, double scale) {
    double mu_2l = mean / (2.0 * scale);
    double Y = random_standard_normal(state);
    Y = mean * Y * Y;
    double X = mean + mu_2l * (Y - sqrt(4.0 * scale * Y + Y * Y));
    double U = state->next_double(state->state);
    if (U <= mean / (mean + X)) {
        return X;
    } else {
        return mean * mean / X;
    }
}

EXPORT double random_triangular(bitgen_t *state, double left, double mode, double right) {
    double base = right - left;
    double leftbase = mode - left;
    double ratio = leftbase / base;
    double U = state->next_double(state->state);
    if (U <= ratio) {
        return left + sqrt(U * base * leftbase);
    } else {
        return right - sqrt((1.0 - U) * base * (right - mode));
    }
}

/* ============ Discrete Distributions ============ */

EXPORT int64_t random_poisson(bitgen_t *state, double lam) {
    if (lam >= 10.0) {
        /* Use PTRS algorithm for large lambda */
        double slam = sqrt(lam);
        double loglam = log(lam);
        double b = 0.931 + 2.53 * slam;
        double a = -0.059 + 0.02483 * b;
        double invalpha = 1.1239 + 1.1328 / (b - 3.4);
        double vr = 0.9277 - 3.6224 / (b - 2.0);

        for (;;) {
            double U, V, us, k;
            U = state->next_double(state->state) - 0.5;
            V = state->next_double(state->state);
            us = 0.5 - fabs(U);
            k = floor((2.0 * a / us + b) * U + lam + 0.43);

            if (us >= 0.07 && V <= vr) {
                return (int64_t)k;
            }
            if (k < 0 || (us < 0.013 && V > us)) {
                continue;
            }
            if (log(V) + log(invalpha) - log(a / (us * us) + b) <=
                -lam + k * loglam - lgamma(k + 1)) {
                return (int64_t)k;
            }
        }
    } else {
        /* Direct method for small lambda */
        double enlam = exp(-lam);
        int64_t k = 0;
        double prod = 1.0;

        for (;;) {
            prod *= state->next_double(state->state);
            if (prod > enlam) {
                k++;
            } else {
                return k;
            }
        }
    }
}

EXPORT int64_t random_geometric(bitgen_t *state, double p) {
    return (int64_t)ceil(log(1.0 - state->next_double(state->state)) / log(1.0 - p));
}

/* ============ Bounded Integers (Lemire's Algorithm) ============ */

EXPORT uint64_t random_bounded_uint64(bitgen_t *state, uint64_t off, uint64_t rng, uint64_t mask) {
    uint64_t val;
    if (rng == 0) {
        return off;
    } else if ((rng & (rng + 1)) == 0) {
        /* Power of 2 - simple masking */
        return off + (state->next_uint64(state->state) & rng);
    } else {
        /* Lemire's algorithm for unbiased bounded integers */
        uint64_t threshold = (-rng) % rng;
        for (;;) {
            val = state->next_uint64(state->state);
            if (val >= threshold) {
                return off + (val % (rng + 1));
            }
        }
    }
}

EXPORT int64_t random_integers(bitgen_t *state, int64_t low, int64_t high) {
    uint64_t rng = (uint64_t)(high - low);
    return low + (int64_t)random_bounded_uint64(state, 0, rng, 0);
}
```

---

### 15.3 Generator Class (TypeScript)

**File:** `src/ts/random.ts` (new file)

```typescript
import { NDArray } from './NDArray.js';
import { DType } from './types.js';

/**
 * Result type for multivariate normal.
 */
export type SizeType = number | number[] | null;

/**
 * BitGenerator abstract base class.
 * All bit generators must implement these methods.
 */
export abstract class BitGenerator {
  protected _wasmState: number = 0;  // Pointer to WASM state

  abstract next_uint64(): bigint;
  abstract next_uint32(): number;
  abstract next_double(): number;

  /**
   * Get the current state as a serializable object.
   */
  abstract getState(): object;

  /**
   * Restore state from a serializable object.
   */
  abstract setState(state: object): void;

  /**
   * Create independent child BitGenerators.
   */
  abstract spawn(n_children: number): BitGenerator[];

  /**
   * Get the WASM state pointer for C calls.
   */
  get wasmStatePtr(): number {
    return this._wasmState;
  }
}

/**
 * SeedSequence for reproducible entropy management.
 */
export class SeedSequence {
  private _entropy: Uint32Array;
  private _spawnKey: number[];
  private _poolSize: number;
  private _nChildrenSpawned: number = 0;

  constructor(
    entropy?: number | number[] | Uint32Array | null,
    spawnKey: number[] = [],
    poolSize: number = 4
  ) {
    this._spawnKey = spawnKey;
    this._poolSize = poolSize;

    if (entropy === null || entropy === undefined) {
      // Use crypto.getRandomValues for OS entropy
      this._entropy = new Uint32Array(poolSize);
      if (typeof crypto !== 'undefined' && crypto.getRandomValues) {
        crypto.getRandomValues(this._entropy);
      } else {
        // Fallback: use Date.now() + Math.random()
        for (let i = 0; i < poolSize; i++) {
          this._entropy[i] = (Date.now() * Math.random()) >>> 0;
        }
      }
    } else if (typeof entropy === 'number') {
      // Single integer seed
      this._entropy = new Uint32Array([entropy >>> 0, (entropy / 0x100000000) >>> 0]);
    } else if (entropy instanceof Uint32Array) {
      this._entropy = entropy;
    } else {
      // Array of numbers
      this._entropy = new Uint32Array(entropy);
    }
  }

  /**
   * Generate state words for seeding a BitGenerator.
   */
  generateState(nWords: number, dtype: 'uint32' | 'uint64' = 'uint32'): Uint32Array | BigUint64Array {
    const wasm = getWasmModule();
    const entropyPtr = wasm._malloc(this._entropy.length * 4);
    const spawnKeyPtr = wasm._malloc(this._spawnKey.length * 4);
    const outPtr = wasm._malloc(nWords * 4);

    try {
      // Copy entropy to WASM
      wasm.HEAPU32.set(this._entropy, entropyPtr >> 2);
      for (let i = 0; i < this._spawnKey.length; i++) {
        wasm.HEAPU32[(spawnKeyPtr >> 2) + i] = this._spawnKey[i];
      }

      // Call WASM seed_seq_generate
      wasm._seed_seq_generate_words(
        entropyPtr, this._entropy.length,
        spawnKeyPtr, this._spawnKey.length,
        this._poolSize,
        outPtr, nWords
      );

      // Copy result
      const result = new Uint32Array(nWords);
      for (let i = 0; i < nWords; i++) {
        result[i] = wasm.HEAPU32[(outPtr >> 2) + i];
      }

      if (dtype === 'uint64') {
        const result64 = new BigUint64Array(nWords / 2);
        for (let i = 0; i < result64.length; i++) {
          result64[i] = BigInt(result[i * 2]) | (BigInt(result[i * 2 + 1]) << 32n);
        }
        return result64;
      }

      return result;
    } finally {
      wasm._free(entropyPtr);
      wasm._free(spawnKeyPtr);
      wasm._free(outPtr);
    }
  }

  /**
   * Spawn n_children independent SeedSequences.
   */
  spawn(nChildren: number): SeedSequence[] {
    const children: SeedSequence[] = [];
    for (let i = 0; i < nChildren; i++) {
      children.push(new SeedSequence(
        this._entropy,
        [...this._spawnKey, this._nChildrenSpawned + i],
        this._poolSize
      ));
    }
    this._nChildrenSpawned += nChildren;
    return children;
  }

  get entropy(): Uint32Array {
    return this._entropy;
  }

  get spawnKey(): number[] {
    return this._spawnKey;
  }
}

/**
 * PCG64 BitGenerator - the default for NumJS.
 * 128-bit state with 64-bit output using XSL-RR permutation.
 */
export class PCG64 extends BitGenerator {
  private _seedSequence: SeedSequence | null = null;

  constructor(seed?: number | number[] | SeedSequence | null) {
    super();

    const wasm = getWasmModule();
    this._wasmState = wasm._pcg64_create();

    if (seed instanceof SeedSequence) {
      this._seedSequence = seed;
      const state = seed.generateState(4, 'uint32');
      this._initFromState(state);
    } else {
      this._seedSequence = new SeedSequence(seed);
      const state = this._seedSequence.generateState(4, 'uint32');
      this._initFromState(state);
    }
  }

  private _initFromState(state: Uint32Array): void {
    const wasm = getWasmModule();
    wasm._pcg64_seed(
      this._wasmState,
      state[0], state[1],  // seed_high, seed_low
      state[2], state[3]   // inc_high, inc_low
    );
  }

  next_uint64(): bigint {
    const wasm = getWasmModule();
    const low = wasm._pcg64_next64_low(this._wasmState);
    const high = wasm._pcg64_next64_high(this._wasmState);
    return BigInt(low) | (BigInt(high) << 32n);
  }

  next_uint32(): number {
    const wasm = getWasmModule();
    return wasm._pcg64_next32(this._wasmState);
  }

  next_double(): number {
    const wasm = getWasmModule();
    return wasm._pcg64_next_double(this._wasmState);
  }

  /**
   * Advance the state by delta steps.
   */
  advance(delta: bigint): this {
    const wasm = getWasmModule();
    const low = Number(delta & 0xFFFFFFFFn);
    const high = Number((delta >> 32n) & 0xFFFFFFFFn);
    wasm._pcg64_advance(this._wasmState, high, low);
    return this;
  }

  /**
   * Return a jumped copy of the generator.
   * Equivalent to 2^128 calls to next_uint64().
   */
  jumped(jumps: number = 1): PCG64 {
    const newGen = new PCG64(this._seedSequence);
    newGen.setState(this.getState());
    // Jump by 2^128 * jumps
    for (let i = 0; i < jumps; i++) {
      newGen.advance(1n << 128n);
    }
    return newGen;
  }

  getState(): object {
    const wasm = getWasmModule();
    const statePtr = wasm._malloc(6 * 8);
    try {
      wasm._pcg64_get_state(this._wasmState, statePtr);
      return {
        bit_generator: 'PCG64',
        state: {
          state: BigInt(wasm.HEAPU32[statePtr >> 2]) |
                 (BigInt(wasm.HEAPU32[(statePtr >> 2) + 1]) << 32n),
          inc: BigInt(wasm.HEAPU32[(statePtr >> 2) + 2]) |
               (BigInt(wasm.HEAPU32[(statePtr >> 2) + 3]) << 32n),
        },
        has_uint32: wasm.HEAPU32[(statePtr >> 2) + 4],
        uinteger: wasm.HEAPU32[(statePtr >> 2) + 5],
      };
    } finally {
      wasm._free(statePtr);
    }
  }

  setState(state: object): void {
    const wasm = getWasmModule();
    const s = state as any;
    const statePtr = wasm._malloc(6 * 8);
    try {
      const st = BigInt(s.state.state);
      const inc = BigInt(s.state.inc);
      wasm.HEAPU32[statePtr >> 2] = Number(st & 0xFFFFFFFFn);
      wasm.HEAPU32[(statePtr >> 2) + 1] = Number((st >> 32n) & 0xFFFFFFFFn);
      wasm.HEAPU32[(statePtr >> 2) + 2] = Number(inc & 0xFFFFFFFFn);
      wasm.HEAPU32[(statePtr >> 2) + 3] = Number((inc >> 32n) & 0xFFFFFFFFn);
      wasm.HEAPU32[(statePtr >> 2) + 4] = s.has_uint32;
      wasm.HEAPU32[(statePtr >> 2) + 5] = s.uinteger;
      wasm._pcg64_set_state(this._wasmState, statePtr);
    } finally {
      wasm._free(statePtr);
    }
  }

  spawn(nChildren: number): PCG64[] {
    if (!this._seedSequence) {
      throw new Error('Cannot spawn from a BitGenerator without a SeedSequence');
    }
    const childSeqs = this._seedSequence.spawn(nChildren);
    return childSeqs.map(seq => new PCG64(seq));
  }
}

/**
 * Generator - main interface for random number generation.
 * Follows NumPy's Generator API.
 */
export class Generator {
  private _bitGenerator: BitGenerator;
  private _wasmBitgen: number = 0;  // bitgen_t pointer

  constructor(bitGenerator?: BitGenerator) {
    this._bitGenerator = bitGenerator ?? new PCG64();
    this._initWasmBitgen();
  }

  private _initWasmBitgen(): void {
    const wasm = getWasmModule();
    this._wasmBitgen = wasm._bitgen_create(this._bitGenerator.wasmStatePtr);
  }

  /**
   * Access the underlying BitGenerator.
   */
  get bitGenerator(): BitGenerator {
    return this._bitGenerator;
  }

  /**
   * Spawn n_children independent Generators.
   */
  spawn(nChildren: number): Generator[] {
    const childBitGens = this._bitGenerator.spawn(nChildren);
    return childBitGens.map(bg => new Generator(bg));
  }

  /* ============ Uniform Distributions ============ */

  /**
   * Return random floats in the half-open interval [0.0, 1.0).
   */
  random(size?: SizeType, dtype: DType = DType.Float64): NDArray | number {
    if (size === null || size === undefined) {
      return this._bitGenerator.next_double();
    }

    const shape = typeof size === 'number' ? [size] : size;
    const result = NDArray.empty(shape, dtype);
    const wasm = getWasmModule();

    if (dtype === DType.Float64) {
      wasm._random_uniform_fill_f64(this._wasmBitgen, result.size, result.dataPtr);
    } else if (dtype === DType.Float32) {
      wasm._random_uniform_fill_f32(this._wasmBitgen, result.size, result.dataPtr);
    } else {
      throw new Error(`Unsupported dtype for random(): ${dtype}`);
    }

    return result;
  }

  /**
   * Return random bytes.
   */
  bytes(length: number): Uint8Array {
    const result = new Uint8Array(length);
    for (let i = 0; i < length; i += 4) {
      const val = this._bitGenerator.next_uint32();
      for (let j = 0; j < 4 && i + j < length; j++) {
        result[i + j] = (val >> (j * 8)) & 0xFF;
      }
    }
    return result;
  }

  /* ============ Integer Methods ============ */

  /**
   * Return random integers from low (inclusive) to high (exclusive).
   */
  integers(
    low: number,
    high?: number | null,
    size?: SizeType,
    dtype: DType = DType.Int64,
    endpoint: boolean = false
  ): NDArray | number {
    // Handle single argument case: integers(high) -> [0, high)
    if (high === null || high === undefined) {
      high = low;
      low = 0;
    }

    if (endpoint) {
      high = high + 1;
    }

    if (size === null || size === undefined) {
      const wasm = getWasmModule();
      return wasm._random_integers(this._wasmBitgen, low, high - 1);
    }

    const shape = typeof size === 'number' ? [size] : size;
    const result = NDArray.empty(shape, dtype);
    const wasm = getWasmModule();

    wasm._random_integers_fill(this._wasmBitgen, low, high - 1, result.size, result.dataPtr, dtype);

    return result;
  }

  /**
   * Generates a random sample from a given array.
   */
  choice(
    a: number | NDArray,
    size?: SizeType,
    replace: boolean = true,
    p?: NDArray | number[] | null,
    axis: number = 0,
    shuffle: boolean = true
  ): NDArray | number {
    const population = typeof a === 'number' ? a : a.size;

    if (!replace && size !== null && size !== undefined) {
      const sizeNum = typeof size === 'number' ? size : size.reduce((a, b) => a * b, 1);
      if (sizeNum > population) {
        throw new Error('Cannot take a larger sample than population when replace=false');
      }
    }

    if (size === null || size === undefined) {
      // Return single element
      const idx = this.integers(0, population);
      if (typeof a === 'number') {
        return idx as number;
      }
      return a.flat.get(idx as number);
    }

    const shape = typeof size === 'number' ? [size] : size;
    const totalSize = shape.reduce((a, b) => a * b, 1);

    if (replace) {
      // With replacement: simple random indices
      const indices = this.integers(0, population, totalSize);
      if (typeof a === 'number') {
        return indices as NDArray;
      }
      // Use take to select elements
      return a.take(indices as NDArray);
    } else {
      // Without replacement: Fisher-Yates shuffle subset
      const indices = new Int32Array(population);
      for (let i = 0; i < population; i++) indices[i] = i;

      for (let i = 0; i < totalSize; i++) {
        const j = i + (this.integers(0, population - i) as number);
        [indices[i], indices[j]] = [indices[j], indices[i]];
      }

      const selectedIndices = NDArray.fromArray(Array.from(indices.slice(0, totalSize)));
      if (typeof a === 'number') {
        return selectedIndices.reshape(shape);
      }
      return a.take(selectedIndices).reshape(shape);
    }
  }

  /* ============ Continuous Distributions ============ */

  /**
   * Draw samples from a uniform distribution.
   */
  uniform(low: number = 0.0, high: number = 1.0, size?: SizeType): NDArray | number {
    if (size === null || size === undefined) {
      return low + (high - low) * this._bitGenerator.next_double();
    }

    const result = this.random(size, DType.Float64) as NDArray;
    // Scale: result = low + (high - low) * result
    const wasm = getWasmModule();
    wasm._array_scale_shift(result.dataPtr, result.size, high - low, low);
    return result;
  }

  /**
   * Draw samples from a standard Normal distribution (mean=0, stdev=1).
   */
  standard_normal(size?: SizeType, dtype: DType = DType.Float64): NDArray | number {
    if (size === null || size === undefined) {
      const wasm = getWasmModule();
      return wasm._random_standard_normal(this._wasmBitgen);
    }

    const shape = typeof size === 'number' ? [size] : size;
    const result = NDArray.empty(shape, dtype);
    const wasm = getWasmModule();

    if (dtype === DType.Float64) {
      wasm._random_standard_normal_fill_f64(this._wasmBitgen, result.size, result.dataPtr);
    } else {
      wasm._random_standard_normal_fill_f32(this._wasmBitgen, result.size, result.dataPtr);
    }

    return result;
  }

  /**
   * Draw random samples from a normal (Gaussian) distribution.
   */
  normal(loc: number = 0.0, scale: number = 1.0, size?: SizeType): NDArray | number {
    if (scale < 0) {
      throw new Error('scale must be non-negative');
    }

    if (size === null || size === undefined) {
      return loc + scale * (this.standard_normal() as number);
    }

    const result = this.standard_normal(size) as NDArray;
    const wasm = getWasmModule();
    wasm._array_scale_shift(result.dataPtr, result.size, scale, loc);
    return result;
  }

  /**
   * Draw samples from a standard exponential distribution.
   */
  standard_exponential(
    size?: SizeType,
    dtype: DType = DType.Float64,
    method: 'zig' | 'inv' = 'zig'
  ): NDArray | number {
    if (size === null || size === undefined) {
      const wasm = getWasmModule();
      return wasm._random_standard_exponential(this._wasmBitgen);
    }

    const shape = typeof size === 'number' ? [size] : size;
    const result = NDArray.empty(shape, dtype);
    const wasm = getWasmModule();

    if (method === 'zig') {
      wasm._random_standard_exponential_fill_f64(this._wasmBitgen, result.size, result.dataPtr);
    } else {
      wasm._random_standard_exponential_inv_fill_f64(this._wasmBitgen, result.size, result.dataPtr);
    }

    return result;
  }

  /**
   * Draw samples from an exponential distribution.
   */
  exponential(scale: number = 1.0, size?: SizeType): NDArray | number {
    if (scale < 0) {
      throw new Error('scale must be non-negative');
    }

    if (size === null || size === undefined) {
      return scale * (this.standard_exponential() as number);
    }

    const result = this.standard_exponential(size) as NDArray;
    const wasm = getWasmModule();
    wasm._array_scale(result.dataPtr, result.size, scale);
    return result;
  }

  /**
   * Draw samples from a Gamma distribution.
   */
  gamma(shape: number, scale: number = 1.0, size?: SizeType): NDArray | number {
    if (shape < 0) {
      throw new Error('shape must be non-negative');
    }
    if (scale < 0) {
      throw new Error('scale must be non-negative');
    }

    if (size === null || size === undefined) {
      const wasm = getWasmModule();
      return wasm._random_gamma(this._wasmBitgen, shape, scale);
    }

    const sizeArr = typeof size === 'number' ? [size] : size;
    const result = NDArray.empty(sizeArr, DType.Float64);
    const wasm = getWasmModule();
    wasm._random_gamma_fill(this._wasmBitgen, shape, scale, result.size, result.dataPtr);
    return result;
  }

  /**
   * Draw samples from a Beta distribution.
   */
  beta(a: number, b: number, size?: SizeType): NDArray | number {
    if (a <= 0) throw new Error('a must be positive');
    if (b <= 0) throw new Error('b must be positive');

    if (size === null || size === undefined) {
      const wasm = getWasmModule();
      return wasm._random_beta(this._wasmBitgen, a, b);
    }

    const sizeArr = typeof size === 'number' ? [size] : size;
    const result = NDArray.empty(sizeArr, DType.Float64);
    const wasm = getWasmModule();
    wasm._random_beta_fill(this._wasmBitgen, a, b, result.size, result.dataPtr);
    return result;
  }

  /**
   * Draw samples from a chi-square distribution.
   */
  chisquare(df: number, size?: SizeType): NDArray | number {
    if (df <= 0) throw new Error('df must be positive');
    return this.gamma(df / 2.0, 2.0, size);
  }

  /**
   * Draw samples from a Student's t distribution.
   */
  standard_t(df: number, size?: SizeType): NDArray | number {
    if (df <= 0) throw new Error('df must be positive');

    if (size === null || size === undefined) {
      const wasm = getWasmModule();
      return wasm._random_standard_t(this._wasmBitgen, df);
    }

    const sizeArr = typeof size === 'number' ? [size] : size;
    const result = NDArray.empty(sizeArr, DType.Float64);
    const wasm = getWasmModule();
    wasm._random_standard_t_fill(this._wasmBitgen, df, result.size, result.dataPtr);
    return result;
  }

  /**
   * Draw samples from a standard Cauchy distribution.
   */
  standard_cauchy(size?: SizeType): NDArray | number {
    if (size === null || size === undefined) {
      const wasm = getWasmModule();
      return wasm._random_standard_cauchy(this._wasmBitgen);
    }

    const sizeArr = typeof size === 'number' ? [size] : size;
    const result = NDArray.empty(sizeArr, DType.Float64);
    const wasm = getWasmModule();
    wasm._random_standard_cauchy_fill(this._wasmBitgen, result.size, result.dataPtr);
    return result;
  }

  /* ============ Discrete Distributions ============ */

  /**
   * Draw samples from a binomial distribution.
   */
  binomial(n: number, p: number, size?: SizeType): NDArray | number {
    if (n < 0) throw new Error('n must be non-negative');
    if (p < 0 || p > 1) throw new Error('p must be in [0, 1]');

    if (size === null || size === undefined) {
      const wasm = getWasmModule();
      return wasm._random_binomial(this._wasmBitgen, n, p);
    }

    const sizeArr = typeof size === 'number' ? [size] : size;
    const result = NDArray.empty(sizeArr, DType.Int64);
    const wasm = getWasmModule();
    wasm._random_binomial_fill(this._wasmBitgen, n, p, result.size, result.dataPtr);
    return result;
  }

  /**
   * Draw samples from a Poisson distribution.
   */
  poisson(lam: number = 1.0, size?: SizeType): NDArray | number {
    if (lam < 0) throw new Error('lam must be non-negative');

    if (size === null || size === undefined) {
      const wasm = getWasmModule();
      return wasm._random_poisson(this._wasmBitgen, lam);
    }

    const sizeArr = typeof size === 'number' ? [size] : size;
    const result = NDArray.empty(sizeArr, DType.Int64);
    const wasm = getWasmModule();
    wasm._random_poisson_fill(this._wasmBitgen, lam, result.size, result.dataPtr);
    return result;
  }

  /**
   * Draw samples from a geometric distribution.
   */
  geometric(p: number, size?: SizeType): NDArray | number {
    if (p <= 0 || p > 1) throw new Error('p must be in (0, 1]');

    if (size === null || size === undefined) {
      const wasm = getWasmModule();
      return wasm._random_geometric(this._wasmBitgen, p);
    }

    const sizeArr = typeof size === 'number' ? [size] : size;
    const result = NDArray.empty(sizeArr, DType.Int64);
    const wasm = getWasmModule();
    wasm._random_geometric_fill(this._wasmBitgen, p, result.size, result.dataPtr);
    return result;
  }

  /* ============ Multivariate Distributions ============ */

  /**
   * Draw random samples from a multivariate normal distribution.
   */
  multivariate_normal(
    mean: NDArray | number[],
    cov: NDArray | number[][],
    size?: SizeType,
    method: 'svd' | 'eigh' | 'cholesky' = 'svd'
  ): NDArray {
    const meanArr = mean instanceof NDArray ? mean : NDArray.fromArray(mean);
    const covArr = cov instanceof NDArray ? cov : NDArray.fromArray(cov);

    const d = meanArr.size;
    if (covArr.shape[0] !== d || covArr.shape[1] !== d) {
      throw new Error('cov must be a square matrix with size matching mean');
    }

    // Compute decomposition of covariance matrix
    let A: NDArray;
    if (method === 'cholesky') {
      // A = cholesky(cov), then X = mean + A @ Z
      const { linalg } = await import('./linalg.js');
      A = linalg.cholesky(covArr);
    } else if (method === 'eigh') {
      const { linalg } = await import('./linalg.js');
      const { eigenvalues, eigenvectors } = linalg.eigh(covArr);
      // A = V @ sqrt(D)
      const sqrtD = eigenvalues.map(Math.sqrt);
      A = linalg.matmul(eigenvectors, NDArray.diag(sqrtD));
    } else {
      // SVD method (most robust)
      const { linalg } = await import('./linalg.js');
      const { U, S, Vh } = linalg.svd(covArr) as any;
      const sqrtS = S.map(Math.sqrt);
      A = linalg.matmul(U, NDArray.diag(sqrtS));
    }

    // Determine output shape
    let outShape: number[];
    if (size === null || size === undefined) {
      outShape = [d];
    } else if (typeof size === 'number') {
      outShape = [size, d];
    } else {
      outShape = [...size, d];
    }

    const totalSamples = outShape.slice(0, -1).reduce((a, b) => a * b, 1);

    // Generate standard normal samples
    const Z = this.standard_normal([totalSamples, d]) as NDArray;

    // Transform: X = mean + Z @ A.T
    const { linalg } = await import('./linalg.js');
    const result = linalg.matmul(Z, A.T);

    // Add mean
    // result += mean (broadcasted)
    for (let i = 0; i < totalSamples; i++) {
      for (let j = 0; j < d; j++) {
        result.setFlat(i * d + j, result.getFlat(i * d + j) + meanArr.getFlat(j));
      }
    }

    return result.reshape(outShape);
  }

  /**
   * Draw samples from a multinomial distribution.
   */
  multinomial(n: number, pvals: NDArray | number[], size?: SizeType): NDArray {
    const p = pvals instanceof NDArray ? pvals.toArray() : pvals;
    const k = p.length;

    // Normalize probabilities
    const pSum = p.reduce((a, b) => a + b, 0);
    const pNorm = p.map(v => v / pSum);

    // Determine output shape
    let outShape: number[];
    if (size === null || size === undefined) {
      outShape = [k];
    } else if (typeof size === 'number') {
      outShape = [size, k];
    } else {
      outShape = [...size, k];
    }

    const totalSamples = outShape.slice(0, -1).reduce((a, b) => a * b, 1);
    const result = NDArray.zeros(outShape, DType.Int64);

    // For each sample, use conditional binomials
    for (let s = 0; s < totalSamples; s++) {
      let remaining = n;
      let pRemaining = 1.0;

      for (let i = 0; i < k - 1 && remaining > 0; i++) {
        const pCond = pNorm[i] / pRemaining;
        const count = this.binomial(remaining, pCond) as number;
        result.setFlat(s * k + i, count);
        remaining -= count;
        pRemaining -= pNorm[i];
      }
      result.setFlat(s * k + k - 1, remaining);
    }

    return result;
  }

  /**
   * Draw samples from a Dirichlet distribution.
   */
  dirichlet(alpha: NDArray | number[], size?: SizeType): NDArray {
    const a = alpha instanceof NDArray ? alpha.toArray() : alpha;
    const k = a.length;

    // Determine output shape
    let outShape: number[];
    if (size === null || size === undefined) {
      outShape = [k];
    } else if (typeof size === 'number') {
      outShape = [size, k];
    } else {
      outShape = [...size, k];
    }

    const totalSamples = outShape.slice(0, -1).reduce((a, b) => a * b, 1);
    const result = NDArray.empty(outShape, DType.Float64);

    // Dirichlet = normalized gamma variates
    for (let s = 0; s < totalSamples; s++) {
      let sum = 0;
      for (let i = 0; i < k; i++) {
        const g = this.gamma(a[i], 1.0) as number;
        result.setFlat(s * k + i, g);
        sum += g;
      }
      // Normalize
      for (let i = 0; i < k; i++) {
        result.setFlat(s * k + i, result.getFlat(s * k + i) / sum);
      }
    }

    return result;
  }

  /* ============ Permutation Methods ============ */

  /**
   * Modify an array in-place by shuffling its contents.
   */
  shuffle(x: NDArray, axis: number = 0): void {
    const n = x.shape[axis];

    if (axis === 0 && x.ndim === 1) {
      // Simple 1D shuffle
      for (let i = n - 1; i > 0; i--) {
        const j = this.integers(0, i + 1) as number;
        const temp = x.getFlat(i);
        x.setFlat(i, x.getFlat(j));
        x.setFlat(j, temp);
      }
    } else {
      // Shuffle along axis using index swapping
      for (let i = n - 1; i > 0; i--) {
        const j = this.integers(0, i + 1) as number;
        // Swap slices at i and j along axis
        const sliceI = x.slice([...Array(axis).fill(':'), i]);
        const sliceJ = x.slice([...Array(axis).fill(':'), j]);
        const temp = sliceI.copy();
        sliceI.set(sliceJ);
        sliceJ.set(temp);
      }
    }
  }

  /**
   * Randomly permute a sequence, or return a permuted range.
   */
  permutation(x: number | NDArray, axis: number = 0): NDArray {
    if (typeof x === 'number') {
      const result = NDArray.arange(0, x, 1, DType.Int64);
      this.shuffle(result);
      return result;
    } else {
      const result = x.copy();
      this.shuffle(result, axis);
      return result;
    }
  }

  /**
   * Randomly permute x along axis, returning a new array.
   */
  permuted(x: NDArray, axis: number | null = null, out?: NDArray): NDArray {
    if (axis === null) {
      // Permute flattened array
      const flat = x.ravel();
      const result = out ?? flat.copy();
      this.shuffle(result);
      return result;
    }

    const result = out ?? x.copy();
    this.shuffle(result, axis);
    return result;
  }
}

/* ============ Module-level Functions ============ */

let _defaultRng: Generator | null = null;

/**
 * Construct a new Generator with the specified BitGenerator.
 * If seed is None, fresh entropy will be pulled from the OS.
 */
export function default_rng(seed?: number | number[] | SeedSequence | BitGenerator | null): Generator {
  if (seed instanceof BitGenerator) {
    return new Generator(seed);
  }
  if (seed instanceof Generator) {
    return seed;
  }
  return new Generator(new PCG64(seed));
}

// Placeholder for WASM module reference
function getWasmModule(): any {
  // This would be imported from wasm-loader.ts
  throw new Error('WASM module not initialized');
}

/* ============ Exports ============ */

export const random = {
  // Classes
  Generator,
  BitGenerator,
  PCG64,
  SeedSequence,

  // Module function
  default_rng,
};
```

---

## File Changes Summary

### New Files to Create

```
src/wasm/random/
├── bitgen.h           # BitGenerator interface
├── pcg64.h            # PCG64 declarations
├── pcg64.c            # PCG64 implementation
├── xoshiro256.h       # Xoshiro256** declarations (optional)
├── xoshiro256.c       # Xoshiro256** implementation (optional)
├── seedseq.h          # SeedSequence declarations
├── seedseq.c          # SeedSequence implementation
├── distributions.h    # Distribution function declarations
├── distributions.c    # Distribution implementations
├── ziggurat_constants.h  # Pre-computed Ziggurat tables
└── ziggurat_constants.c  # Ziggurat table data

src/ts/
└── random.ts          # Full random module (Generator, BitGenerator, etc.)
```

### Files to Modify

```
src/ts/types.ts
├── Add random-related type exports
└── Add WASM function declarations for random

src/ts/index.ts
├── Export random module
├── Export Generator, PCG64, SeedSequence classes
└── Export default_rng function

scripts/build-wasm.sh
├── Add random/*.c to compilation
└── Add EXPORTED_FUNCTIONS for all random operations
```

---

## Build Script Updates

Add to `scripts/build-wasm.sh` EXPORTED_FUNCTIONS:

```bash
# BitGenerator Infrastructure
"_pcg64_create",
"_pcg64_seed",
"_pcg64_next64_low",
"_pcg64_next64_high",
"_pcg64_next32",
"_pcg64_next_double",
"_pcg64_advance",
"_pcg64_get_state",
"_pcg64_set_state",
"_bitgen_create",

# SeedSequence
"_seed_seq_generate_words",

# Uniform Distributions
"_random_uniform_fill_f64",
"_random_uniform_fill_f32",

# Normal Distribution
"_random_standard_normal",
"_random_standard_normal_fill_f64",
"_random_standard_normal_fill_f32",

# Exponential Distribution
"_random_standard_exponential",
"_random_standard_exponential_fill_f64",
"_random_standard_exponential_inv_fill_f64",

# Gamma Distribution
"_random_gamma",
"_random_gamma_fill",

# Other Continuous
"_random_beta",
"_random_beta_fill",
"_random_standard_t",
"_random_standard_t_fill",
"_random_standard_cauchy",
"_random_standard_cauchy_fill",

# Discrete Distributions
"_random_binomial",
"_random_binomial_fill",
"_random_poisson",
"_random_poisson_fill",
"_random_geometric",
"_random_geometric_fill",

# Bounded Integers
"_random_integers",
"_random_integers_fill",

# Utility
"_array_scale",
"_array_scale_shift"
```

Add new source files:

```bash
emcc \
    "$SRC_DIR/ndarray.c" \
    "$SRC_DIR/dtype.c" \
    "$SRC_DIR/broadcast.c" \
    "$SRC_DIR/indexing.c" \
    "$SRC_DIR/ufunc.c" \
    "$SRC_DIR/statistics.c" \
    "$SRC_DIR/random/pcg64.c" \
    "$SRC_DIR/random/seedseq.c" \
    "$SRC_DIR/random/distributions.c" \
    "$SRC_DIR/random/ziggurat_constants.c" \
    -o "$OUT_DIR/numjs.cjs" \
    ...
```

---

## WasmModule Interface Updates

Add to `src/ts/types.ts`:

```typescript
// BitGenerator
_pcg64_create(): number;
_pcg64_seed(state: number, seedHigh: number, seedLow: number, incHigh: number, incLow: number): void;
_pcg64_next64_low(state: number): number;
_pcg64_next64_high(state: number): number;
_pcg64_next32(state: number): number;
_pcg64_next_double(state: number): number;
_pcg64_advance(state: number, deltaHigh: number, deltaLow: number): void;
_pcg64_get_state(state: number, outPtr: number): void;
_pcg64_set_state(state: number, inPtr: number): void;
_bitgen_create(bitgenState: number): number;

// SeedSequence
_seed_seq_generate_words(entropyPtr: number, entropyLen: number,
                         spawnKeyPtr: number, spawnKeyLen: number,
                         poolSize: number, outPtr: number, nWords: number): void;

// Distributions
_random_standard_normal(bitgen: number): number;
_random_standard_normal_fill_f64(bitgen: number, n: number, outPtr: number): void;
_random_standard_normal_fill_f32(bitgen: number, n: number, outPtr: number): void;
_random_standard_exponential(bitgen: number): number;
_random_standard_exponential_fill_f64(bitgen: number, n: number, outPtr: number): void;
_random_gamma(bitgen: number, shape: number, scale: number): number;
_random_gamma_fill(bitgen: number, shape: number, scale: number, n: number, outPtr: number): void;
_random_beta(bitgen: number, a: number, b: number): number;
_random_beta_fill(bitgen: number, a: number, b: number, n: number, outPtr: number): void;
_random_binomial(bitgen: number, n: number, p: number): number;
_random_binomial_fill(bitgen: number, n: number, p: number, cnt: number, outPtr: number): void;
_random_poisson(bitgen: number, lam: number): number;
_random_poisson_fill(bitgen: number, lam: number, cnt: number, outPtr: number): void;
_random_integers(bitgen: number, low: number, high: number): number;
_random_integers_fill(bitgen: number, low: number, high: number, cnt: number, outPtr: number, dtype: number): void;

// Array utilities
_array_scale(ptr: number, n: number, scale: number): void;
_array_scale_shift(ptr: number, n: number, scale: number, shift: number): void;
```

---

## Implementation Order

```
Phase 15.1: BitGenerator Infrastructure (Week 1)
├── Day 1: bitgen_t interface, uint64_to_double conversion
├── Day 2: PCG64 state structure, seed function
├── Day 3: PCG64 next64, next32, next_double
├── Day 4: PCG64 advance (jump-ahead algorithm)
├── Day 5: PCG64 state get/set, tests

Phase 15.2: SeedSequence & Core Distributions (Week 2)
├── Day 1: SeedSequence mixing algorithm
├── Day 2: SeedSequence generate_state, spawn
├── Day 3: Ziggurat tables generation/validation
├── Day 4: standard_normal (Ziggurat)
├── Day 5: standard_exponential (Ziggurat) + tests

Phase 15.3: Continuous Distributions (Week 3)
├── Day 1: Gamma distribution (Marsaglia-Tsang)
├── Day 2: Beta distribution (Johnk's + gamma)
├── Day 3: Chi-square, F, Student's t
├── Day 4: Cauchy, Pareto, Weibull, Power
├── Day 5: Laplace, Gumbel, Logistic, Lognormal, Rayleigh, Wald, Triangular

Phase 15.4: Discrete Distributions (Week 4)
├── Day 1: Bounded integers (Lemire's algorithm)
├── Day 2: Binomial (BTPE algorithm)
├── Day 3: Poisson (PTRS algorithm)
├── Day 4: Geometric, Negative Binomial
├── Day 5: Hypergeometric, Logseries, Zipf

Phase 15.5: TypeScript Generator Class (Week 5)
├── Day 1: BitGenerator base class, PCG64 TypeScript wrapper
├── Day 2: SeedSequence TypeScript implementation
├── Day 3: Generator class - random, integers, uniform
├── Day 4: Generator class - normal, exponential, gamma, beta
├── Day 5: Generator class - discrete distributions

Phase 15.6: Advanced Features (Week 6)
├── Day 1: choice() with/without replacement
├── Day 2: shuffle(), permutation(), permuted()
├── Day 3: multivariate_normal()
├── Day 4: multinomial(), dirichlet()
├── Day 5: Integration tests, edge cases

Phase 15.7: Polish & Testing (Week 7)
├── Day 1: Performance optimization
├── Day 2: Statistical validation (chi-square tests)
├── Day 3: NumPy comparison tests
├── Day 4: Documentation
└── Day 5: Final integration
```

---

## Verification Plan

After Phase 15 completion, verify:

```bash
# Build
npm run build

# Run all tests
npm test

# Phase 15 specific tests:

# BitGenerator
✓ PCG64 produces consistent sequence from same seed
✓ PCG64 state can be saved and restored
✓ PCG64 advance produces correct offset
✓ SeedSequence produces different states from different seeds
✓ SeedSequence spawn creates independent streams

# Distributions - Statistical Tests
✓ uniform() samples are in [0, 1)
✓ normal() mean ≈ 0, std ≈ 1 for large samples
✓ exponential() mean ≈ scale
✓ gamma() mean ≈ shape * scale
✓ beta() mean ≈ a / (a + b)
✓ binomial() mean ≈ n * p
✓ poisson() mean ≈ lambda

# Generator Methods
✓ integers() returns values in [low, high)
✓ choice() without replacement never repeats
✓ shuffle() modifies array in-place
✓ permutation() returns new shuffled array
✓ multivariate_normal() covariance matches input

# Edge Cases
✓ Empty size returns scalar
✓ Zero-length arrays handled
✓ Extreme parameters (very small/large) work
✓ dtype parameter respected
```

Generate NumPy comparison vectors:

```python
# tests/python/generate_random_tests.py
import numpy as np
import json

# Use fixed seed for reproducibility
rng = np.random.default_rng(12345)

tests = {
    "pcg64_sequence": {
        "seed": 12345,
        "first_10_uint64": [rng.integers(0, 2**64, dtype=np.uint64) for _ in range(10)],
    },
    "uniform": [
        {"size": 1000, "mean": 0.5, "std": 0.2887},  # sqrt(1/12)
    ],
    "normal": [
        {"loc": 0, "scale": 1, "size": 10000},
        {"loc": 5, "scale": 2, "size": 10000},
    ],
    "exponential": [
        {"scale": 1, "size": 10000, "expected_mean": 1},
        {"scale": 2.5, "size": 10000, "expected_mean": 2.5},
    ],
    "gamma": [
        {"shape": 2, "scale": 1, "size": 10000, "expected_mean": 2},
        {"shape": 0.5, "scale": 2, "size": 10000, "expected_mean": 1},
    ],
    "beta": [
        {"a": 2, "b": 5, "size": 10000, "expected_mean": 2/7},
    ],
    "binomial": [
        {"n": 10, "p": 0.5, "size": 10000, "expected_mean": 5},
        {"n": 100, "p": 0.1, "size": 10000, "expected_mean": 10},
    ],
    "poisson": [
        {"lam": 5, "size": 10000, "expected_mean": 5},
        {"lam": 100, "size": 10000, "expected_mean": 100},
    ],
    "multivariate_normal": [
        {
            "mean": [0, 0],
            "cov": [[1, 0.5], [0.5, 1]],
            "size": 10000,
        },
    ],
}

with open("tests/fixtures/random_vectors.json", "w") as f:
    json.dump(tests, f, indent=2, default=lambda x: x.tolist() if hasattr(x, 'tolist') else x)
```

---

## Critical Dependencies

Phase 15 completion enables:

- **Machine Learning**: Random initialization, stochastic gradient descent, dropout
- **Monte Carlo Simulations**: Financial modeling, physics simulations
- **Statistical Testing**: Bootstrap, permutation tests
- **Data Augmentation**: Random transforms, noise injection

Phase 15 should be implemented after completing:
- Phase 4 (Ufuncs) - needed for element-wise operations
- Phase 6 (Statistics) - provides mean, var for validation
- Phase 13 (numpy.linalg) - needed for multivariate_normal (Cholesky)

---

## Performance Considerations

### Ziggurat Algorithm

```
Performance Characteristics:
- 98.9% of samples return on first iteration
- Pre-computed tables (256 entries): ~4KB memory
- Avoids expensive exp(), log() calls in hot path
- Float32 and Float64 variants
```

### PCG64 Optimization

```
Key Optimizations:
- 128-bit state allows 2^128 period
- XSL-RR output function is fast (shift, rotate, xor)
- Buffered 32-bit output from 64-bit generation
- Jump-ahead in O(log n) steps
```

### Memory Layout

```
Array Generation:
- Direct fill into NDArray data buffer
- Avoid intermediate allocations
- WASM linear memory for all operations
- Batch generation reduces call overhead
```

### Threading Considerations

```
Thread Safety:
- Each Generator has independent BitGenerator
- SeedSequence.spawn() creates non-overlapping streams
- No shared state between spawned generators
- Consider Web Workers for parallel generation
```

---

## API Compatibility Notes

### NumPy Differences

```typescript
// NumPy: Generator is the primary interface
// rng = np.random.default_rng(seed)
// samples = rng.random(size)

// NumJS: Same pattern
// const rng = random.default_rng(seed);
// const samples = rng.random(size);
```

### Reproducibility

```typescript
// Same seed produces same sequence
const rng1 = random.default_rng(12345);
const rng2 = random.default_rng(12345);

rng1.random(10);  // [0.22, 0.87, ...]
rng2.random(10);  // [0.22, 0.87, ...] - identical
```

### Parallel Streams

```typescript
// Spawn independent generators for parallel work
const rng = random.default_rng(12345);
const [rng1, rng2, rng3] = rng.spawn(3);

// Each produces independent, non-overlapping sequence
// Safe for Web Workers
```

---

## Algorithm References

| Distribution | Algorithm | Reference |
|--------------|-----------|-----------|
| Normal | Ziggurat | Marsaglia & Tsang (2000) |
| Exponential | Ziggurat | Marsaglia & Tsang (2000) |
| Gamma | Marsaglia-Tsang | Marsaglia & Tsang (2000) |
| Beta | Johnk's + Gamma | Johnk (1964) |
| Binomial | BTPE | Kachitvichyanukul & Schmeiser (1988) |
| Poisson | PTRS | Hörmann (1993) |
| Bounded Int | Lemire's | Lemire (2019) |
| PCG64 | PCG Family | O'Neill (2014) |

---

## Ziggurat Table Generation

The Ziggurat tables can be generated using the following algorithm:

```python
# Generate Ziggurat tables for normal distribution
import numpy as np

def generate_ziggurat_normal(n=256):
    """Generate Ziggurat tables for N(0,1)."""
    # Tail cutoff
    r = 3.6541528853610088

    # Area of each rectangle
    v = 9.91256303526217e-3

    # Build tables
    x = np.zeros(n + 1)
    x[n] = r
    x[n-1] = r

    for i in range(n-2, 0, -1):
        x[i] = np.sqrt(-2.0 * np.log(v / x[i+1] + np.exp(-0.5 * x[i+1]**2)))

    x[0] = 0

    # k table (threshold)
    k = np.zeros(n, dtype=np.uint64)
    for i in range(n):
        k[i] = int((x[i+1] / x[i]) * 2**53)

    # w table (width)
    w = np.zeros(n)
    for i in range(n):
        w[i] = x[i] / 2**53

    # f table (f(x))
    f = np.zeros(n)
    for i in range(n):
        f[i] = np.exp(-0.5 * x[i]**2)

    return k, w, f

k, w, f = generate_ziggurat_normal()
print(f"const uint64_t ki_double[256] = {{{', '.join(map(str, k))}}};")
```
