/**
 * NumJS Random - MT19937 BitGenerator Implementation
 *
 * MT19937 (Mersenne Twister) pseudo-random number generator.
 * Based on the reference implementation by Makoto Matsumoto and Takuji Nishimura.
 */

#include "mt19937.h"
#include <stdlib.h>
#include <string.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/* ============ MT19937 Core Algorithm ============ */

/**
 * Generate the next N values in the state array (twist operation).
 */
static void mt19937_gen(mt19937_internal_state *s) {
    static const uint32_t mag01[2] = {0x0UL, MT19937_MATRIX_A};
    uint32_t y;
    int kk;

    for (kk = 0; kk < MT19937_N - MT19937_M; kk++) {
        y = (s->key[kk] & MT19937_UPPER_MASK) | (s->key[kk + 1] & MT19937_LOWER_MASK);
        s->key[kk] = s->key[kk + MT19937_M] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    for (; kk < MT19937_N - 1; kk++) {
        y = (s->key[kk] & MT19937_UPPER_MASK) | (s->key[kk + 1] & MT19937_LOWER_MASK);
        s->key[kk] = s->key[kk + (MT19937_M - MT19937_N)] ^ (y >> 1) ^ mag01[y & 0x1UL];
    }
    y = (s->key[MT19937_N - 1] & MT19937_UPPER_MASK) | (s->key[0] & MT19937_LOWER_MASK);
    s->key[MT19937_N - 1] = s->key[MT19937_M - 1] ^ (y >> 1) ^ mag01[y & 0x1UL];

    s->pos = 0;
}

/**
 * Generate a single 32-bit value with tempering.
 */
static inline uint32_t mt19937_genrand(mt19937_internal_state *s) {
    uint32_t y;

    if (s->pos >= MT19937_N) {
        mt19937_gen(s);
    }

    y = s->key[s->pos++];

    /* Tempering */
    y ^= (y >> 11);
    y ^= (y << 7) & 0x9d2c5680UL;
    y ^= (y << 15) & 0xefc60000UL;
    y ^= (y >> 18);

    return y;
}

/* ============ Public API ============ */

EXPORT mt19937_state* mt19937_create(void) {
    mt19937_state *state = (mt19937_state *)malloc(sizeof(mt19937_state));
    if (state) {
        state->state = (mt19937_internal_state *)malloc(sizeof(mt19937_internal_state));
        state->has_uint32 = 0;
        state->uinteger = 0;
        if (!state->state) {
            free(state);
            return NULL;
        }
        /* Initialize to zero state - must call mt19937_seed before use */
        memset(state->state->key, 0, sizeof(state->state->key));
        state->state->pos = MT19937_N + 1;  /* Indicates uninitialized */
    }
    return state;
}

EXPORT void mt19937_free(mt19937_state *state) {
    if (state) {
        if (state->state) {
            free(state->state);
        }
        free(state);
    }
}

EXPORT void mt19937_seed(mt19937_state *state, uint32_t seed) {
    if (!state || !state->state) {
        return;
    }

    mt19937_internal_state *s = state->state;

    s->key[0] = seed;
    for (int i = 1; i < MT19937_N; i++) {
        /* See Knuth TAOCP Vol2. 3rd Ed. P.106 for multiplier */
        s->key[i] = (1812433253UL * (s->key[i - 1] ^ (s->key[i - 1] >> 30)) + i);
    }
    s->pos = MT19937_N;

    /* Clear buffered 32-bit value */
    state->has_uint32 = 0;
    state->uinteger = 0;
}

EXPORT void mt19937_seed_array(mt19937_state *state, uint32_t *init_key, int key_length) {
    if (!state || !state->state || !init_key || key_length <= 0) {
        return;
    }

    mt19937_internal_state *s = state->state;

    /* First initialize with a fixed seed */
    mt19937_seed(state, 19650218UL);

    int i = 1;
    int j = 0;
    int k = (MT19937_N > key_length) ? MT19937_N : key_length;

    for (; k > 0; k--) {
        /* Non-linear mixing */
        s->key[i] = (s->key[i] ^ ((s->key[i - 1] ^ (s->key[i - 1] >> 30)) * 1664525UL))
                    + init_key[j] + j;
        i++;
        j++;
        if (i >= MT19937_N) {
            s->key[0] = s->key[MT19937_N - 1];
            i = 1;
        }
        if (j >= key_length) {
            j = 0;
        }
    }

    for (k = MT19937_N - 1; k > 0; k--) {
        s->key[i] = (s->key[i] ^ ((s->key[i - 1] ^ (s->key[i - 1] >> 30)) * 1566083941UL))
                    - i;
        i++;
        if (i >= MT19937_N) {
            s->key[0] = s->key[MT19937_N - 1];
            i = 1;
        }
    }

    /* MSB is 1; assuring non-zero initial array */
    s->key[0] = 0x80000000UL;

    /* Clear buffered 32-bit value */
    state->has_uint32 = 0;
    state->uinteger = 0;
}

EXPORT uint32_t mt19937_next32(mt19937_state *state) {
    return mt19937_genrand(state->state);
}

EXPORT uint64_t mt19937_next64(mt19937_state *state) {
    /* Combine two 32-bit values into one 64-bit value */
    uint64_t high = (uint64_t)mt19937_next32(state);
    uint64_t low = (uint64_t)mt19937_next32(state);
    return (high << 32) | low;
}

EXPORT double mt19937_next_double(mt19937_state *state) {
    /* Use the MT19937 convention: combines 27 + 26 bits */
    uint32_t a = mt19937_next32(state) >> 5;   /* 27 bits */
    uint32_t b = mt19937_next32(state) >> 6;   /* 26 bits */
    return (a * 67108864.0 + b) / 9007199254740992.0;  /* 2^53 */
}

EXPORT uint32_t mt19937_next64_parts(mt19937_state *state, uint32_t *high_out) {
    uint64_t val = mt19937_next64(state);
    if (high_out) {
        *high_out = (uint32_t)(val >> 32);
    }
    return (uint32_t)val;
}

EXPORT void mt19937_get_state(mt19937_state *state, uint32_t *key_out,
                               int *pos_out, int *has_uint32_out, uint32_t *uinteger_out) {
    if (!state || !state->state) {
        return;
    }

    if (key_out) {
        memcpy(key_out, state->state->key, MT19937_N * sizeof(uint32_t));
    }
    if (pos_out) {
        *pos_out = state->state->pos;
    }
    if (has_uint32_out) {
        *has_uint32_out = state->has_uint32;
    }
    if (uinteger_out) {
        *uinteger_out = state->uinteger;
    }
}

EXPORT void mt19937_set_state(mt19937_state *state, uint32_t *key_in,
                               int pos, int has_uint32, uint32_t uinteger) {
    if (!state || !key_in) {
        return;
    }

    if (!state->state) {
        state->state = (mt19937_internal_state *)malloc(sizeof(mt19937_internal_state));
        if (!state->state) {
            return;
        }
    }

    memcpy(state->state->key, key_in, MT19937_N * sizeof(uint32_t));
    state->state->pos = pos;
    state->has_uint32 = has_uint32;
    state->uinteger = uinteger;
}

/* ============ BitGenerator Interface ============ */

static uint64_t _mt19937_next64_wrapper(void *state) {
    return mt19937_next64((mt19937_state *)state);
}

static uint32_t _mt19937_next32_wrapper(void *state) {
    return mt19937_next32((mt19937_state *)state);
}

static double _mt19937_next_double_wrapper(void *state) {
    return mt19937_next_double((mt19937_state *)state);
}

EXPORT void mt19937_init_bitgen(bitgen_t *bitgen, mt19937_state *state) {
    if (!bitgen || !state) {
        return;
    }

    bitgen->state = state;
    bitgen->next_uint64 = _mt19937_next64_wrapper;
    bitgen->next_uint32 = _mt19937_next32_wrapper;
    bitgen->next_double = _mt19937_next_double_wrapper;
    bitgen->next_raw = _mt19937_next64_wrapper;
}

/* ============ Bulk Generation ============ */

EXPORT void mt19937_fill_uint32(mt19937_state *state, uint32_t *out, int64_t count) {
    for (int64_t i = 0; i < count; i++) {
        out[i] = mt19937_next32(state);
    }
}

EXPORT void mt19937_fill_uint64(mt19937_state *state, uint64_t *out, int64_t count) {
    for (int64_t i = 0; i < count; i++) {
        out[i] = mt19937_next64(state);
    }
}

EXPORT void mt19937_fill_double(mt19937_state *state, double *out, int64_t count) {
    for (int64_t i = 0; i < count; i++) {
        out[i] = mt19937_next_double(state);
    }
}
