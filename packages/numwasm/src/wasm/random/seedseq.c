/**
 * NumJS Random - SeedSequence Implementation
 *
 * Implements entropy mixing for reproducible seed generation.
 */

#include "seedseq.h"
#include <stdlib.h>
#include <string.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/* ============ Internal Hash Functions ============ */

/**
 * Hash mixing function based on MurmurHash3's finalizer.
 * Mixes a value with the current hash state.
 */
static inline uint32_t hashmix(uint32_t value, uint32_t *hash) {
    value ^= *hash;
    *hash *= SEED_SEQ_MULT_A;
    value *= SEED_SEQ_MULT_A;
    value ^= value >> 16;
    return value;
}

/**
 * Rotate left by n bits.
 */
static inline uint32_t rotl32(uint32_t x, int n) {
    return (x << n) | (x >> (32 - n));
}

/* ============ Public API ============ */

EXPORT void seed_seq_init(seed_sequence_t *seq,
                           const uint32_t *entropy, int32_t entropy_len,
                           const uint32_t *spawn_key, int32_t spawn_key_len,
                           int32_t pool_size) {
    if (!seq) {
        return;
    }

    /* Default pool size is 4 (128 bits) */
    seq->pool_size = pool_size > 0 ? pool_size : 4;
    seq->pool = (uint32_t *)calloc(seq->pool_size, sizeof(uint32_t));
    seq->n_children = 0;

    if (!seq->pool) {
        seq->pool_size = 0;
        return;
    }

    /* Mix entropy into pool */
    if (entropy && entropy_len > 0) {
        uint32_t hash = SEED_SEQ_INIT_A;
        for (int32_t i = 0; i < entropy_len; i++) {
            int32_t idx = i % seq->pool_size;
            seq->pool[idx] ^= hashmix(entropy[i], &hash);
        }
    }

    /* Mix spawn_key into pool with different initial hash */
    if (spawn_key && spawn_key_len > 0) {
        uint32_t hash = SEED_SEQ_INIT_B;
        for (int32_t i = 0; i < spawn_key_len; i++) {
            int32_t idx = i % seq->pool_size;
            seq->pool[idx] ^= hashmix(spawn_key[i], &hash);
        }
    }

    /* Final mixing pass to ensure all bits affect all others */
    /* This is important for short entropy inputs */
    for (int32_t round = 0; round < 3; round++) {
        uint32_t hash = SEED_SEQ_INIT_A + round;
        for (int32_t i = 0; i < seq->pool_size; i++) {
            /* Mix with neighbors */
            uint32_t prev = seq->pool[(i + seq->pool_size - 1) % seq->pool_size];
            uint32_t curr = seq->pool[i];
            uint32_t next = seq->pool[(i + 1) % seq->pool_size];

            curr ^= rotl32(prev, 7);
            curr ^= rotl32(next, 13);
            seq->pool[i] = hashmix(curr, &hash);
        }
    }
}

EXPORT void seed_seq_generate(seed_sequence_t *seq, uint32_t *out, int32_t n_words) {
    if (!seq || !seq->pool || !out || n_words <= 0) {
        return;
    }

    /* Generate output words by mixing pool with index */
    uint32_t hash = SEED_SEQ_INIT_A;

    for (int32_t i = 0; i < n_words; i++) {
        /* Combine pool value with index for unique output per position */
        uint32_t value = seq->pool[i % seq->pool_size];
        value ^= (uint32_t)i;

        /* Additional mixing based on pool state */
        value ^= rotl32(seq->pool[(i + 1) % seq->pool_size], 11);
        value ^= rotl32(seq->pool[(i + 2) % seq->pool_size], 23);

        out[i] = hashmix(value, &hash);
    }

    /* Second pass to improve statistical properties */
    hash = SEED_SEQ_INIT_B;
    for (int32_t i = 0; i < n_words; i++) {
        uint32_t prev = out[(i + n_words - 1) % n_words];
        out[i] ^= rotl32(hashmix(prev, &hash), 17);
    }
}

EXPORT void seed_seq_free(seed_sequence_t *seq) {
    if (seq && seq->pool) {
        free(seq->pool);
        seq->pool = NULL;
        seq->pool_size = 0;
    }
}

EXPORT void seed_seq_generate_words(const uint32_t *entropy, int32_t entropy_len,
                                     const uint32_t *spawn_key, int32_t spawn_key_len,
                                     int32_t pool_size,
                                     uint32_t *out, int32_t n_words) {
    seed_sequence_t seq;
    seed_seq_init(&seq, entropy, entropy_len, spawn_key, spawn_key_len, pool_size);
    seed_seq_generate(&seq, out, n_words);
    seed_seq_free(&seq);
}

/* ============ Additional Utilities ============ */

/**
 * Mix a 64-bit seed into two 32-bit words.
 * Useful for converting a single integer seed to entropy array.
 */
EXPORT void seed_seq_mix64(uint64_t seed, uint32_t *out) {
    if (!out) return;

    /* Split into low and high parts */
    out[0] = (uint32_t)(seed & 0xFFFFFFFFULL);
    out[1] = (uint32_t)(seed >> 32);

    /* Mix using SplitMix-style transformation */
    uint64_t z = seed + 0x9e3779b97f4a7c15ULL;
    z = (z ^ (z >> 30)) * 0xbf58476d1ce4e5b9ULL;
    z = (z ^ (z >> 27)) * 0x94d049bb133111ebULL;
    z = z ^ (z >> 31);

    out[0] ^= (uint32_t)(z & 0xFFFFFFFFULL);
    out[1] ^= (uint32_t)(z >> 32);
}

/**
 * Generate entropy from timestamp (for JavaScript fallback).
 * This is not cryptographically secure but provides reasonable uniqueness.
 */
EXPORT void seed_seq_from_time(uint32_t timestamp_low, uint32_t timestamp_high,
                                uint32_t counter, uint32_t *out, int32_t n_words) {
    /* Create entropy from timestamp components */
    uint32_t entropy[4];
    entropy[0] = timestamp_low;
    entropy[1] = timestamp_high;
    entropy[2] = counter;
    entropy[3] = timestamp_low ^ timestamp_high ^ counter;

    seed_seq_generate_words(entropy, 4, NULL, 0, 4, out, n_words);
}
