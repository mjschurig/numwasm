/**
 * NumJS Random - SFC64 BitGenerator Implementation
 *
 * SFC64 (Small Fast Chaotic 64) pseudo-random number generator.
 * Based on PractRand's sfc64 by Chris Doty-Humphrey.
 */

#include "sfc64.h"
#include <stdlib.h>
#include <string.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/* ============ SFC64 Core Algorithm ============ */

/**
 * Rotate left operation for 64-bit values.
 */
static inline uint64_t rotl(uint64_t value, unsigned int rot) {
    return (value << rot) | (value >> (64 - rot));
}

/**
 * SFC64 core step function.
 * Generates one 64-bit output and advances the state.
 */
static inline uint64_t sfc64_next(sfc64_internal_state *s) {
    /* Combine a + b + counter, then increment counter */
    uint64_t tmp = s->a + s->b + s->w++;

    /* Update state using mix operations */
    s->a = s->b ^ (s->b >> 11);
    s->b = s->c + (s->c << 3);
    s->c = rotl(s->c, 24) + tmp;

    return tmp;
}

/* ============ Public API ============ */

EXPORT sfc64_state* sfc64_create(void) {
    sfc64_state *state = (sfc64_state *)malloc(sizeof(sfc64_state));
    if (state) {
        state->state = (sfc64_internal_state *)malloc(sizeof(sfc64_internal_state));
        state->has_uint32 = 0;
        state->uinteger = 0;
        if (!state->state) {
            free(state);
            return NULL;
        }
        /* Initialize to zero state - must call sfc64_seed before use */
        memset(state->state, 0, sizeof(sfc64_internal_state));
    }
    return state;
}

EXPORT void sfc64_free(sfc64_state *state) {
    if (state) {
        if (state->state) {
            free(state->state);
        }
        free(state);
    }
}

EXPORT void sfc64_seed(sfc64_state *state, uint64_t *seed) {
    if (!state || !state->state || !seed) {
        return;
    }

    /* Initialize state from seed values */
    state->state->a = seed[0];
    state->state->b = seed[1];
    state->state->c = seed[2];
    state->state->w = seed[3];

    /* Clear buffered 32-bit value */
    state->has_uint32 = 0;
    state->uinteger = 0;

    /* Mix state by running 12 iterations (warmup) */
    for (int i = 0; i < 12; i++) {
        (void)sfc64_next(state->state);
    }
}

EXPORT void sfc64_seed_parts(sfc64_state *state, uint32_t *parts) {
    if (!state || !state->state || !parts) {
        return;
    }

    /* Combine 32-bit parts into 64-bit values */
    uint64_t seed[4];
    seed[0] = ((uint64_t)parts[0] << 32) | parts[1];
    seed[1] = ((uint64_t)parts[2] << 32) | parts[3];
    seed[2] = ((uint64_t)parts[4] << 32) | parts[5];
    seed[3] = ((uint64_t)parts[6] << 32) | parts[7];

    sfc64_seed(state, seed);
}

EXPORT uint64_t sfc64_next64(sfc64_state *state) {
    return sfc64_next(state->state);
}

EXPORT uint32_t sfc64_next32(sfc64_state *state) {
    /* Use buffered value if available */
    if (state->has_uint32) {
        state->has_uint32 = 0;
        return state->uinteger;
    }

    /* Generate 64-bit value, use lower 32 bits, buffer upper 32 bits */
    uint64_t next = sfc64_next64(state);
    state->has_uint32 = 1;
    state->uinteger = (uint32_t)(next >> 32);
    return (uint32_t)next;
}

EXPORT double sfc64_next_double(sfc64_state *state) {
    return uint64_to_double(sfc64_next64(state));
}

EXPORT uint32_t sfc64_next64_parts(sfc64_state *state, uint32_t *high_out) {
    uint64_t val = sfc64_next64(state);
    if (high_out) {
        *high_out = (uint32_t)(val >> 32);
    }
    return (uint32_t)val;
}

EXPORT void sfc64_get_state(sfc64_state *state, uint64_t *out) {
    if (!state || !state->state || !out) {
        return;
    }

    out[0] = state->state->a;
    out[1] = state->state->b;
    out[2] = state->state->c;
    out[3] = state->state->w;
    out[4] = (uint64_t)state->has_uint32;
    out[5] = (uint64_t)state->uinteger;
}

EXPORT void sfc64_set_state(sfc64_state *state, uint64_t *in) {
    if (!state || !in) {
        return;
    }

    if (!state->state) {
        state->state = (sfc64_internal_state *)malloc(sizeof(sfc64_internal_state));
        if (!state->state) {
            return;
        }
    }

    state->state->a = in[0];
    state->state->b = in[1];
    state->state->c = in[2];
    state->state->w = in[3];
    state->has_uint32 = (int)in[4];
    state->uinteger = (uint32_t)in[5];
}

/* ============ BitGenerator Interface ============ */

static uint64_t _sfc64_next64_wrapper(void *state) {
    return sfc64_next64((sfc64_state *)state);
}

static uint32_t _sfc64_next32_wrapper(void *state) {
    return sfc64_next32((sfc64_state *)state);
}

static double _sfc64_next_double_wrapper(void *state) {
    return sfc64_next_double((sfc64_state *)state);
}

EXPORT void sfc64_init_bitgen(bitgen_t *bitgen, sfc64_state *state) {
    if (!bitgen || !state) {
        return;
    }

    bitgen->state = state;
    bitgen->next_uint64 = _sfc64_next64_wrapper;
    bitgen->next_uint32 = _sfc64_next32_wrapper;
    bitgen->next_double = _sfc64_next_double_wrapper;
    bitgen->next_raw = _sfc64_next64_wrapper;
}

/* ============ Bulk Generation ============ */

EXPORT void sfc64_fill_uint64(sfc64_state *state, uint64_t *out, int64_t count) {
    for (int64_t i = 0; i < count; i++) {
        out[i] = sfc64_next64(state);
    }
}

EXPORT void sfc64_fill_double(sfc64_state *state, double *out, int64_t count) {
    for (int64_t i = 0; i < count; i++) {
        out[i] = sfc64_next_double(state);
    }
}
