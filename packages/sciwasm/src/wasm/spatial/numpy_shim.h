#ifndef NUMPY_SHIM_H
#define NUMPY_SHIM_H

#include <stdint.h>

// Shim for numpy types - we don't use numpy in WASM
#ifdef __EMSCRIPTEN__
typedef int64_t npy_intp;
typedef int32_t npy_int32;
typedef double npy_float64;
#else
typedef long npy_intp;
typedef int npy_int32;
typedef double npy_float64;
#endif

#endif  // NUMPY_SHIM_H
