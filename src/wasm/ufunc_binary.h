/**
 * NumJS-WASM Binary Universal Functions
 *
 * Declares ~25 binary element-wise operations with broadcasting support.
 * Each function takes two NDArrays and returns a new NDArray with the result.
 */

#ifndef NUMJS_UFUNC_BINARY_H
#define NUMJS_UFUNC_BINARY_H

#include "ndarray.h"

/* ============ Arithmetic Operations ============ */

/** Add arguments element-wise. */
NDArray* ufunc_add(NDArray* x1, NDArray* x2);

/** Subtract arguments, element-wise. */
NDArray* ufunc_subtract(NDArray* x1, NDArray* x2);

/** Multiply arguments element-wise. */
NDArray* ufunc_multiply(NDArray* x1, NDArray* x2);

/** Divide arguments element-wise (true division). */
NDArray* ufunc_divide(NDArray* x1, NDArray* x2);
NDArray* ufunc_true_divide(NDArray* x1, NDArray* x2);  /* Alias */

/** Return the largest integer smaller or equal to the division of the inputs. */
NDArray* ufunc_floor_divide(NDArray* x1, NDArray* x2);

/** Return element-wise remainder of division. */
NDArray* ufunc_remainder(NDArray* x1, NDArray* x2);
NDArray* ufunc_mod(NDArray* x1, NDArray* x2);  /* Alias */

/** Returns element-wise remainder of floor_divide (like Python %). */
NDArray* ufunc_fmod(NDArray* x1, NDArray* x2);

/** First array elements raised to powers from second array, element-wise. */
NDArray* ufunc_power(NDArray* x1, NDArray* x2);

/* ============ Comparison Operations ============ */

/** Return (x1 == x2) element-wise. */
NDArray* ufunc_equal(NDArray* x1, NDArray* x2);

/** Return (x1 != x2) element-wise. */
NDArray* ufunc_not_equal(NDArray* x1, NDArray* x2);

/** Return (x1 < x2) element-wise. */
NDArray* ufunc_less(NDArray* x1, NDArray* x2);

/** Return (x1 <= x2) element-wise. */
NDArray* ufunc_less_equal(NDArray* x1, NDArray* x2);

/** Return (x1 > x2) element-wise. */
NDArray* ufunc_greater(NDArray* x1, NDArray* x2);

/** Return (x1 >= x2) element-wise. */
NDArray* ufunc_greater_equal(NDArray* x1, NDArray* x2);

/* ============ Extrema Operations ============ */

/** Element-wise maximum of array elements. */
NDArray* ufunc_maximum(NDArray* x1, NDArray* x2);

/** Element-wise minimum of array elements. */
NDArray* ufunc_minimum(NDArray* x1, NDArray* x2);

/** Element-wise maximum (propagates NaN). */
NDArray* ufunc_fmax(NDArray* x1, NDArray* x2);

/** Element-wise minimum (propagates NaN). */
NDArray* ufunc_fmin(NDArray* x1, NDArray* x2);

/* ============ Logical Operations ============ */

/** Compute the truth value of x1 AND x2 element-wise. */
NDArray* ufunc_logical_and(NDArray* x1, NDArray* x2);

/** Compute the truth value of x1 OR x2 element-wise. */
NDArray* ufunc_logical_or(NDArray* x1, NDArray* x2);

/** Compute the truth value of x1 XOR x2, element-wise. */
NDArray* ufunc_logical_xor(NDArray* x1, NDArray* x2);

/* ============ Bitwise Operations ============ */

/** Compute the bit-wise AND of two arrays element-wise. */
NDArray* ufunc_bitwise_and(NDArray* x1, NDArray* x2);

/** Compute the bit-wise OR of two arrays element-wise. */
NDArray* ufunc_bitwise_or(NDArray* x1, NDArray* x2);

/** Compute the bit-wise XOR of two arrays element-wise. */
NDArray* ufunc_bitwise_xor(NDArray* x1, NDArray* x2);

/** Shift the bits of an integer to the left. */
NDArray* ufunc_left_shift(NDArray* x1, NDArray* x2);

/** Shift the bits of an integer to the right. */
NDArray* ufunc_right_shift(NDArray* x1, NDArray* x2);

/* ============ Special Mathematical Functions ============ */

/** Element-wise arc tangent of x1/x2 choosing the quadrant correctly. */
NDArray* ufunc_arctan2(NDArray* x1, NDArray* x2);

/** Given the "legs" of a right triangle, return its hypotenuse. */
NDArray* ufunc_hypot(NDArray* x1, NDArray* x2);

/** Change the sign of x1 to that of x2, element-wise. */
NDArray* ufunc_copysign(NDArray* x1, NDArray* x2);

/** Return element-wise True where signbit is set (less than zero). */
NDArray* ufunc_signbit(NDArray* x);

/** Compute the logarithm of the sum of exponentials of the inputs. */
NDArray* ufunc_logaddexp(NDArray* x1, NDArray* x2);

/** Compute log(2**x1 + 2**x2). */
NDArray* ufunc_logaddexp2(NDArray* x1, NDArray* x2);

#endif /* NUMJS_UFUNC_BINARY_H */
