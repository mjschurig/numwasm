/**
 * NumJS Random - PCG64 BitGenerator Implementation
 *
 * PCG64 (Permuted Congruential Generator) with 128-bit state.
 * Uses XSL-RR (xorshift low, random rotate) output function.
 */

#include "pcg64.h"
#include <stdlib.h>
#include <string.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/* ============ 128-bit Arithmetic ============ */

/**
 * Multiply two 128-bit integers using schoolbook multiplication.
 * Result = a * b (mod 2^128)
 */
pcg128_t pcg128_mult(pcg128_t a, pcg128_t b) {
    pcg128_t result;

    /* Split into 32-bit pieces for multiplication */
    uint64_t a_lo_lo = a.low & 0xFFFFFFFFULL;
    uint64_t a_lo_hi = a.low >> 32;
    uint64_t b_lo_lo = b.low & 0xFFFFFFFFULL;
    uint64_t b_lo_hi = b.low >> 32;

    /* Compute partial products */
    uint64_t lo_lo = a_lo_lo * b_lo_lo;
    uint64_t hi_lo = a_lo_hi * b_lo_lo;
    uint64_t lo_hi = a_lo_lo * b_lo_hi;
    uint64_t hi_hi = a_lo_hi * b_lo_hi;

    /* Sum with proper carry handling */
    uint64_t cross = (lo_lo >> 32) + (hi_lo & 0xFFFFFFFFULL) + lo_hi;

    result.low = (cross << 32) | (lo_lo & 0xFFFFFFFFULL);
    result.high = (hi_lo >> 32) + (cross >> 32) + hi_hi +
                  a.high * b.low + a.low * b.high;

    return result;
}

/**
 * Add two 128-bit integers with carry.
 */
pcg128_t pcg128_add(pcg128_t a, pcg128_t b) {
    pcg128_t result;
    result.low = a.low + b.low;
    result.high = a.high + b.high + (result.low < a.low ? 1 : 0);
    return result;
}

/* ============ PCG64 Output Function ============ */

/**
 * XSL-RR output function for PCG64.
 * Takes the 128-bit state and produces a 64-bit output.
 *
 * XSL: XOR the high and low 64-bit halves
 * RR: Random rotate based on high bits
 */
static inline uint64_t pcg_output_xsl_rr_128_64(pcg128_t state) {
    uint64_t xorshifted = state.high ^ state.low;
    int rot = (int)(state.high >> 58);  /* Use top 6 bits for rotation */
    return (xorshifted >> rot) | (xorshifted << ((-rot) & 63));
}

/* ============ PCG64 LCG Step ============ */

/**
 * Advance the LCG by one step.
 * new_state = state * multiplier + increment
 */
static inline void pcg_step(pcg_state_t *state) {
    pcg128_t mult = {PCG_DEFAULT_MULTIPLIER_HIGH, PCG_DEFAULT_MULTIPLIER_LOW};
    state->state = pcg128_add(pcg128_mult(state->state, mult), state->inc);
}

/* ============ Public API ============ */

EXPORT pcg64_state* pcg64_create(void) {
    pcg64_state *state = (pcg64_state *)malloc(sizeof(pcg64_state));
    if (state) {
        state->pcg_state = (pcg_state_t *)malloc(sizeof(pcg_state_t));
        state->has_uint32 = 0;
        state->uinteger = 0;
        if (!state->pcg_state) {
            free(state);
            return NULL;
        }
        /* Initialize to zero state - must call pcg64_seed before use */
        memset(state->pcg_state, 0, sizeof(pcg_state_t));
    }
    return state;
}

EXPORT void pcg64_free(pcg64_state *state) {
    if (state) {
        if (state->pcg_state) {
            free(state->pcg_state);
        }
        free(state);
    }
}

EXPORT void pcg64_seed(pcg64_state *state,
                        uint64_t seed_high, uint64_t seed_low,
                        uint64_t inc_high, uint64_t inc_low) {
    if (!state || !state->pcg_state) {
        return;
    }

    /* Ensure increment is odd (required for full period) */
    state->pcg_state->inc.high = inc_high;
    state->pcg_state->inc.low = (inc_low << 1) | 1ULL;

    /* Initialize state to zero, step once */
    state->pcg_state->state.high = 0;
    state->pcg_state->state.low = 0;
    pcg_step(state->pcg_state);

    /* Add seed and step again */
    state->pcg_state->state.high += seed_high;
    state->pcg_state->state.low += seed_low;
    /* Handle carry from low to high */
    if (state->pcg_state->state.low < seed_low) {
        state->pcg_state->state.high++;
    }
    pcg_step(state->pcg_state);

    /* Clear buffered 32-bit value */
    state->has_uint32 = 0;
    state->uinteger = 0;
}

EXPORT uint64_t pcg64_next64(pcg64_state *state) {
    /* Generate output from current state, then advance */
    uint64_t result = pcg_output_xsl_rr_128_64(state->pcg_state->state);
    pcg_step(state->pcg_state);
    return result;
}

EXPORT uint32_t pcg64_next32(pcg64_state *state) {
    /* Use buffered value if available */
    if (state->has_uint32) {
        state->has_uint32 = 0;
        return state->uinteger;
    }

    /* Generate 64-bit value, use lower 32 bits, buffer upper 32 bits */
    uint64_t next = pcg64_next64(state);
    state->has_uint32 = 1;
    state->uinteger = (uint32_t)(next >> 32);
    return (uint32_t)next;
}

EXPORT double pcg64_next_double(pcg64_state *state) {
    return uint64_to_double(pcg64_next64(state));
}

/**
 * Advance the generator by delta steps using O(log n) algorithm.
 *
 * Based on Brown, "Random Number Generation with Arbitrary Stride",
 * Transactions of the American Nuclear Society (Nov. 1994).
 *
 * The algorithm computes state' = state * mult^delta + inc * (mult^delta - 1) / (mult - 1)
 * using repeated squaring.
 */
EXPORT void pcg64_advance(pcg64_state *state, uint64_t delta_high, uint64_t delta_low) {
    pcg128_t delta = {delta_high, delta_low};
    pcg128_t cur_mult = {PCG_DEFAULT_MULTIPLIER_HIGH, PCG_DEFAULT_MULTIPLIER_LOW};
    pcg128_t cur_plus = state->pcg_state->inc;
    pcg128_t acc_mult = {0, 1};  /* Start with identity for multiplication */
    pcg128_t acc_plus = {0, 0};  /* Start with zero for addition */

    /* Binary decomposition of delta */
    while (delta.high > 0 || delta.low > 0) {
        if (delta.low & 1) {
            acc_mult = pcg128_mult(acc_mult, cur_mult);
            acc_plus = pcg128_add(pcg128_mult(acc_plus, cur_mult), cur_plus);
        }

        /* cur_plus = (cur_mult + 1) * cur_plus */
        pcg128_t one = {0, 1};
        cur_plus = pcg128_mult(pcg128_add(cur_mult, one), cur_plus);
        cur_mult = pcg128_mult(cur_mult, cur_mult);

        /* delta >>= 1 */
        delta.low = (delta.low >> 1) | (delta.high << 63);
        delta.high >>= 1;
    }

    /* Apply accumulated transformation */
    state->pcg_state->state = pcg128_add(
        pcg128_mult(acc_mult, state->pcg_state->state),
        acc_plus
    );
}

EXPORT void pcg64_get_state(pcg64_state *state, uint64_t *out) {
    if (!state || !state->pcg_state || !out) {
        return;
    }

    out[0] = state->pcg_state->state.high;
    out[1] = state->pcg_state->state.low;
    out[2] = state->pcg_state->inc.high;
    out[3] = state->pcg_state->inc.low;
    out[4] = (uint64_t)state->has_uint32;
    out[5] = (uint64_t)state->uinteger;
}

EXPORT void pcg64_set_state(pcg64_state *state, uint64_t *in) {
    if (!state || !in) {
        return;
    }

    if (!state->pcg_state) {
        state->pcg_state = (pcg_state_t *)malloc(sizeof(pcg_state_t));
        if (!state->pcg_state) {
            return;
        }
    }

    state->pcg_state->state.high = in[0];
    state->pcg_state->state.low = in[1];
    state->pcg_state->inc.high = in[2];
    state->pcg_state->inc.low = in[3];
    state->has_uint32 = (int)in[4];
    state->uinteger = (uint32_t)in[5];
}

/* ============ BitGenerator Interface ============ */

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
    if (!bitgen || !state) {
        return;
    }

    bitgen->state = state;
    bitgen->next_uint64 = _pcg64_next64_wrapper;
    bitgen->next_uint32 = _pcg64_next32_wrapper;
    bitgen->next_double = _pcg64_next_double_wrapper;
    bitgen->next_raw = _pcg64_next64_wrapper;
}

/* ============ Convenience Functions for WASM ============ */

/**
 * Create and seed a PCG64 generator in one call.
 * Useful for JavaScript bindings.
 */
EXPORT pcg64_state* pcg64_create_seeded(uint64_t seed_high, uint64_t seed_low,
                                         uint64_t inc_high, uint64_t inc_low) {
    pcg64_state *state = pcg64_create();
    if (state) {
        pcg64_seed(state, seed_high, seed_low, inc_high, inc_low);
    }
    return state;
}

/**
 * Seed the PCG64 generator using 32-bit parts.
 * This is useful for JavaScript which doesn't have native 64-bit integers.
 */
EXPORT void pcg64_seed_parts(pcg64_state *state,
                              uint32_t seed_hh, uint32_t seed_hl,
                              uint32_t seed_lh, uint32_t seed_ll,
                              uint32_t inc_hh, uint32_t inc_hl,
                              uint32_t inc_lh, uint32_t inc_ll) {
    uint64_t seed_high = ((uint64_t)seed_hh << 32) | seed_hl;
    uint64_t seed_low = ((uint64_t)seed_lh << 32) | seed_ll;
    uint64_t inc_high = ((uint64_t)inc_hh << 32) | inc_hl;
    uint64_t inc_low = ((uint64_t)inc_lh << 32) | inc_ll;
    pcg64_seed(state, seed_high, seed_low, inc_high, inc_low);
}

/**
 * Generate next 64-bit value and return low 32 bits.
 * High 32 bits returned via out parameter.
 * This is useful for JavaScript which doesn't have native 64-bit integers.
 */
EXPORT uint32_t pcg64_next64_parts(pcg64_state *state, uint32_t *high_out) {
    uint64_t val = pcg64_next64(state);
    if (high_out) {
        *high_out = (uint32_t)(val >> 32);
    }
    return (uint32_t)val;
}

/**
 * Fill an array with random uint64 values.
 */
EXPORT void pcg64_fill_uint64(pcg64_state *state, uint64_t *out, int64_t count) {
    for (int64_t i = 0; i < count; i++) {
        out[i] = pcg64_next64(state);
    }
}

/**
 * Fill an array with random doubles in [0, 1).
 */
EXPORT void pcg64_fill_double(pcg64_state *state, double *out, int64_t count) {
    for (int64_t i = 0; i < count; i++) {
        out[i] = pcg64_next_double(state);
    }
}
