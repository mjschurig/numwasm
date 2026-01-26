/**
 * NumJS-WASM Unary Universal Functions
 *
 * Declares ~30 unary element-wise operations.
 * Each function takes an NDArray and returns a new NDArray with the result.
 */

#ifndef NUMJS_UFUNC_UNARY_H
#define NUMJS_UFUNC_UNARY_H

#include "ndarray.h"

/* ============ Arithmetic Operations ============ */

/** Numerical negative, element-wise. */
NDArray* ufunc_negative(NDArray* input);

/** Numerical positive (returns a copy). */
NDArray* ufunc_positive(NDArray* input);

/** Calculate the absolute value element-wise. */
NDArray* ufunc_absolute(NDArray* input);
NDArray* ufunc_abs(NDArray* input);  /* Alias */

/** Returns element-wise sign: -1, 0, or +1. */
NDArray* ufunc_sign(NDArray* input);

/* ============ Powers and Roots ============ */

/** Return the non-negative square-root of an array, element-wise. */
NDArray* ufunc_sqrt(NDArray* input);

/** Return the element-wise square of the input. */
NDArray* ufunc_square(NDArray* input);

/** Return the cube-root of an array, element-wise. */
NDArray* ufunc_cbrt(NDArray* input);

/** Return the reciprocal of the argument, element-wise. */
NDArray* ufunc_reciprocal(NDArray* input);

/* ============ Exponential Functions ============ */

/** Calculate the exponential of all elements in the input array. */
NDArray* ufunc_exp(NDArray* input);

/** Calculate 2**x for all elements in the array. */
NDArray* ufunc_exp2(NDArray* input);

/** Calculate exp(x) - 1 for all elements in the array. */
NDArray* ufunc_expm1(NDArray* input);

/* ============ Logarithmic Functions ============ */

/** Natural logarithm, element-wise. */
NDArray* ufunc_log(NDArray* input);

/** Base-2 logarithm of x. */
NDArray* ufunc_log2(NDArray* input);

/** Base-10 logarithm, element-wise. */
NDArray* ufunc_log10(NDArray* input);

/** Return the natural logarithm of one plus the input array, element-wise. */
NDArray* ufunc_log1p(NDArray* input);

/* ============ Trigonometric Functions ============ */

/** Trigonometric sine, element-wise. */
NDArray* ufunc_sin(NDArray* input);

/** Cosine element-wise. */
NDArray* ufunc_cos(NDArray* input);

/** Compute tangent element-wise. */
NDArray* ufunc_tan(NDArray* input);

/** Inverse sine, element-wise. */
NDArray* ufunc_arcsin(NDArray* input);

/** Trigonometric inverse cosine, element-wise. */
NDArray* ufunc_arccos(NDArray* input);

/** Trigonometric inverse tangent, element-wise. */
NDArray* ufunc_arctan(NDArray* input);

/* ============ Hyperbolic Functions ============ */

/** Hyperbolic sine, element-wise. */
NDArray* ufunc_sinh(NDArray* input);

/** Hyperbolic cosine, element-wise. */
NDArray* ufunc_cosh(NDArray* input);

/** Compute hyperbolic tangent element-wise. */
NDArray* ufunc_tanh(NDArray* input);

/** Inverse hyperbolic sine element-wise. */
NDArray* ufunc_arcsinh(NDArray* input);

/** Inverse hyperbolic cosine, element-wise. */
NDArray* ufunc_arccosh(NDArray* input);

/** Inverse hyperbolic tangent element-wise. */
NDArray* ufunc_arctanh(NDArray* input);

/* ============ Rounding Functions ============ */

/** Return the floor of the input, element-wise. */
NDArray* ufunc_floor(NDArray* input);

/** Return the ceiling of the input, element-wise. */
NDArray* ufunc_ceil(NDArray* input);

/** Return the truncated value of the input, element-wise. */
NDArray* ufunc_trunc(NDArray* input);

/** Round elements of the array to the nearest integer. */
NDArray* ufunc_rint(NDArray* input);

/** Round to nearest even value (banker's rounding). */
NDArray* ufunc_round(NDArray* input);

/* ============ Angle Conversion ============ */

/** Convert angles from radians to degrees. */
NDArray* ufunc_degrees(NDArray* input);
NDArray* ufunc_rad2deg(NDArray* input);  /* Alias */

/** Convert angles from degrees to radians. */
NDArray* ufunc_radians(NDArray* input);
NDArray* ufunc_deg2rad(NDArray* input);  /* Alias */

/* ============ Logical Operations ============ */

/** Compute the truth value of NOT x element-wise. */
NDArray* ufunc_logical_not(NDArray* input);

/* ============ Bitwise Operations (Integer) ============ */

/** Compute bit-wise inversion, element-wise. */
NDArray* ufunc_invert(NDArray* input);
NDArray* ufunc_bitwise_not(NDArray* input);  /* Alias */

#endif /* NUMJS_UFUNC_UNARY_H */
