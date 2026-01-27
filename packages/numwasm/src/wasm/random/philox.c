/**
 * NumJS Random - Philox BitGenerator Implementation
 *
 * Philox4x64-10 counter-based pseudo-random number generator.
 * Based on the Random123 library reference implementation.
 */

#include "philox.h"
#include <stdlib.h>
#include <string.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/* ============ 64-bit Multiplication Helpers ============ */

/**
 * Multiply two 64-bit values and return both high and low 64-bit parts.
 * This implements 64x64 -> 128 bit multiplication without __uint128_t.
 *
 * @param a     First multiplicand
 * @param b     Second multiplicand
 * @param hi    Output for upper 64 bits
 * @param lo    Output for lower 64 bits
 */
static inline void mulhilo64(uint64_t a, uint64_t b, uint64_t *hi, uint64_t *lo) {
    /* Split into 32-bit pieces */
    uint64_t a_lo = (uint32_t)a;
    uint64_t a_hi = a >> 32;
    uint64_t b_lo = (uint32_t)b;
    uint64_t b_hi = b >> 32;

    /* Partial products */
    uint64_t p0 = a_lo * b_lo;
    uint64_t p1 = a_lo * b_hi;
    uint64_t p2 = a_hi * b_lo;
    uint64_t p3 = a_hi * b_hi;

    /* Sum middle products with carry */
    uint64_t mid = p1 + p2;
    uint64_t carry = (mid < p1) ? 1ULL : 0ULL;  /* Overflow in middle sum */

    /* Combine results */
    *lo = p0 + (mid << 32);
    uint64_t lo_carry = (*lo < p0) ? 1ULL : 0ULL;

    *hi = p3 + (mid >> 32) + (carry << 32) + lo_carry;
}

/* ============ Philox Core Algorithm ============ */

/**
 * Single Philox round function.
 * Applies multiplication-based mixing to the counter with the current key.
 */
static inline void philox4x64_round(uint64_t *ctr, const uint64_t *key) {
    uint64_t hi0, lo0, hi1, lo1;

    mulhilo64(ctr[0], PHILOX_M0, &hi0, &lo0);
    mulhilo64(ctr[2], PHILOX_M1, &hi1, &lo1);

    /* Feistel-like mixing */
    uint64_t new_ctr0 = hi1 ^ ctr[1] ^ key[0];
    uint64_t new_ctr1 = lo1;
    uint64_t new_ctr2 = hi0 ^ ctr[3] ^ key[1];
    uint64_t new_ctr3 = lo0;

    ctr[0] = new_ctr0;
    ctr[1] = new_ctr1;
    ctr[2] = new_ctr2;
    ctr[3] = new_ctr3;
}

/**
 * Full Philox4x64-10 function.
 * Applies 10 rounds with key bumping.
 *
 * @param ctr  Counter (input, will be modified in place for output)
 * @param key  Key (will not be modified)
 */
static void philox4x64_10(uint64_t *ctr, const uint64_t *key) {
    uint64_t k[2] = {key[0], key[1]};

    for (int round = 0; round < PHILOX_ROUNDS; round++) {
        philox4x64_round(ctr, k);

        /* Bump key for next round */
        if (round < PHILOX_ROUNDS - 1) {
            k[0] += PHILOX_W0;
            k[1] += PHILOX_W1;
        }
    }
}

/**
 * Generate a block of 4 random values.
 */
static void philox_generate_block(philox_state *state) {
    /* Copy counter to buffer (this will be transformed) */
    memcpy(state->buffer, state->ctr.v, sizeof(state->buffer));

    /* Apply Philox function */
    philox4x64_10(state->buffer, state->key.v);

    /* Increment counter (256-bit addition) */
    state->ctr.v[0]++;
    if (state->ctr.v[0] == 0) {
        state->ctr.v[1]++;
        if (state->ctr.v[1] == 0) {
            state->ctr.v[2]++;
            if (state->ctr.v[2] == 0) {
                state->ctr.v[3]++;
            }
        }
    }

    state->buffer_pos = 0;
}

/* ============ Public API ============ */

EXPORT philox_state* philox_create(void) {
    philox_state *state = (philox_state *)malloc(sizeof(philox_state));
    if (state) {
        memset(state, 0, sizeof(philox_state));
        state->buffer_pos = PHILOX_BUFFER_SIZE;  /* Buffer empty */
    }
    return state;
}

EXPORT void philox_free(philox_state *state) {
    if (state) {
        free(state);
    }
}

EXPORT void philox_seed(philox_state *state, uint64_t key0, uint64_t key1) {
    if (!state) {
        return;
    }

    state->key.v[0] = key0;
    state->key.v[1] = key1;

    /* Reset counter to zero */
    memset(state->ctr.v, 0, sizeof(state->ctr.v));

    /* Invalidate buffer */
    state->buffer_pos = PHILOX_BUFFER_SIZE;
    state->has_uint32 = 0;
    state->uinteger = 0;
}

EXPORT void philox_seed_parts(philox_state *state, uint32_t *parts) {
    if (!state || !parts) {
        return;
    }

    uint64_t key0 = ((uint64_t)parts[0] << 32) | parts[1];
    uint64_t key1 = ((uint64_t)parts[2] << 32) | parts[3];

    philox_seed(state, key0, key1);
}

EXPORT void philox_set_counter(philox_state *state, uint64_t *counter, uint64_t *key) {
    if (!state) {
        return;
    }

    if (counter) {
        memcpy(state->ctr.v, counter, sizeof(state->ctr.v));
    }
    if (key) {
        memcpy(state->key.v, key, sizeof(state->key.v));
    }

    /* Invalidate buffer */
    state->buffer_pos = PHILOX_BUFFER_SIZE;
    state->has_uint32 = 0;
    state->uinteger = 0;
}

EXPORT uint64_t philox_next64(philox_state *state) {
    if (state->buffer_pos >= PHILOX_BUFFER_SIZE) {
        philox_generate_block(state);
    }
    return state->buffer[state->buffer_pos++];
}

EXPORT uint32_t philox_next32(philox_state *state) {
    /* Use buffered value if available */
    if (state->has_uint32) {
        state->has_uint32 = 0;
        return state->uinteger;
    }

    /* Generate 64-bit value, use lower 32 bits, buffer upper 32 bits */
    uint64_t next = philox_next64(state);
    state->has_uint32 = 1;
    state->uinteger = (uint32_t)(next >> 32);
    return (uint32_t)next;
}

EXPORT double philox_next_double(philox_state *state) {
    return uint64_to_double(philox_next64(state));
}

EXPORT uint32_t philox_next64_parts(philox_state *state, uint32_t *high_out) {
    uint64_t val = philox_next64(state);
    if (high_out) {
        *high_out = (uint32_t)(val >> 32);
    }
    return (uint32_t)val;
}

EXPORT void philox_jump(philox_state *state) {
    if (!state) {
        return;
    }

    /* Jump by 2^128: increment counter[2] by 1 */
    state->ctr.v[2]++;
    if (state->ctr.v[2] == 0) {
        state->ctr.v[3]++;
    }

    /* Invalidate buffer */
    state->buffer_pos = PHILOX_BUFFER_SIZE;
    state->has_uint32 = 0;
    state->uinteger = 0;
}

EXPORT void philox_advance(philox_state *state, uint64_t *step) {
    if (!state || !step) {
        return;
    }

    /* 256-bit addition: ctr += step */
    uint64_t carry = 0;

    for (int i = 0; i < 4; i++) {
        uint64_t old = state->ctr.v[i];
        state->ctr.v[i] += step[i] + carry;

        /* Check for overflow */
        if (carry) {
            carry = (state->ctr.v[i] <= old) ? 1 : 0;
        } else {
            carry = (state->ctr.v[i] < old) ? 1 : 0;
        }
    }

    /* Invalidate buffer */
    state->buffer_pos = PHILOX_BUFFER_SIZE;
    state->has_uint32 = 0;
    state->uinteger = 0;
}

EXPORT void philox_advance64(philox_state *state, uint64_t step) {
    uint64_t step_arr[4] = {step, 0, 0, 0};
    philox_advance(state, step_arr);
}

EXPORT void philox_get_state(philox_state *state,
                              uint64_t *ctr_out,
                              uint64_t *key_out,
                              int *buffer_pos_out,
                              uint64_t *buffer_out,
                              int *has_uint32_out,
                              uint32_t *uinteger_out) {
    if (!state) {
        return;
    }

    if (ctr_out) {
        memcpy(ctr_out, state->ctr.v, sizeof(state->ctr.v));
    }
    if (key_out) {
        memcpy(key_out, state->key.v, sizeof(state->key.v));
    }
    if (buffer_pos_out) {
        *buffer_pos_out = state->buffer_pos;
    }
    if (buffer_out) {
        memcpy(buffer_out, state->buffer, sizeof(state->buffer));
    }
    if (has_uint32_out) {
        *has_uint32_out = state->has_uint32;
    }
    if (uinteger_out) {
        *uinteger_out = state->uinteger;
    }
}

EXPORT void philox_set_state(philox_state *state,
                              uint64_t *ctr_in,
                              uint64_t *key_in,
                              int buffer_pos,
                              uint64_t *buffer_in,
                              int has_uint32,
                              uint32_t uinteger) {
    if (!state) {
        return;
    }

    if (ctr_in) {
        memcpy(state->ctr.v, ctr_in, sizeof(state->ctr.v));
    }
    if (key_in) {
        memcpy(state->key.v, key_in, sizeof(state->key.v));
    }
    state->buffer_pos = buffer_pos;
    if (buffer_in) {
        memcpy(state->buffer, buffer_in, sizeof(state->buffer));
    }
    state->has_uint32 = has_uint32;
    state->uinteger = uinteger;
}

/* ============ BitGenerator Interface ============ */

static uint64_t _philox_next64_wrapper(void *state) {
    return philox_next64((philox_state *)state);
}

static uint32_t _philox_next32_wrapper(void *state) {
    return philox_next32((philox_state *)state);
}

static double _philox_next_double_wrapper(void *state) {
    return philox_next_double((philox_state *)state);
}

EXPORT void philox_init_bitgen(bitgen_t *bitgen, philox_state *state) {
    if (!bitgen || !state) {
        return;
    }

    bitgen->state = state;
    bitgen->next_uint64 = _philox_next64_wrapper;
    bitgen->next_uint32 = _philox_next32_wrapper;
    bitgen->next_double = _philox_next_double_wrapper;
    bitgen->next_raw = _philox_next64_wrapper;
}

/* ============ Bulk Generation ============ */

EXPORT void philox_fill_uint64(philox_state *state, uint64_t *out, int64_t count) {
    for (int64_t i = 0; i < count; i++) {
        out[i] = philox_next64(state);
    }
}

EXPORT void philox_fill_double(philox_state *state, double *out, int64_t count) {
    for (int64_t i = 0; i < count; i++) {
        out[i] = philox_next_double(state);
    }
}
