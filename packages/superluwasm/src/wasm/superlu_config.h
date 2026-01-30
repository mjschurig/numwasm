/*
 * SuperLU configuration for WebAssembly/Emscripten build
 */

#ifndef SUPERLU_CONFIG_H
#define SUPERLU_CONFIG_H

/* Enable COLAMD ordering */
#define HAVE_COLAMD 1

/* Use 32-bit indices (sufficient for WASM memory limits) */
typedef int int_t;

#endif /* SUPERLU_CONFIG_H */
