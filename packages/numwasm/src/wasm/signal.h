#ifndef NUMJS_SIGNAL_H
#define NUMJS_SIGNAL_H

#include "ndarray.h"

/* Convolution/correlation modes */
#define CONVOLVE_MODE_FULL  0
#define CONVOLVE_MODE_SAME  1
#define CONVOLVE_MODE_VALID 2

/**
 * 1D discrete convolution of two arrays.
 *
 * @param a First input array (1D)
 * @param v Second input array (1D)
 * @param mode Output mode: FULL (0), SAME (1), or VALID (2)
 * @return New array containing the discrete linear convolution
 */
NDArray* ndarray_convolve(NDArray* a, NDArray* v, int32_t mode);

/**
 * 1D cross-correlation of two arrays.
 *
 * @param a First input array (1D)
 * @param v Second input array (1D)
 * @param mode Output mode: FULL (0), SAME (1), or VALID (2)
 * @return New array containing the cross-correlation
 */
NDArray* ndarray_correlate(NDArray* a, NDArray* v, int32_t mode);

#endif /* NUMJS_SIGNAL_H */
