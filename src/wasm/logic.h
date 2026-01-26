/**
 * NumJS-WASM Logic & Comparison Functions
 *
 * Provides truth testing, element-wise predicates, and array comparison operations.
 * Follows NumPy's logic and comparison API.
 *
 * Reference: numpy/_core/fromnumeric.py, numpy/_core/numeric.py
 */

#ifndef NUMJS_LOGIC_H
#define NUMJS_LOGIC_H

#include "ndarray.h"

/* ============ Element-wise Predicates ============ */

/**
 * Test element-wise for finiteness (not infinity and not NaN).
 *
 * @param arr Input array
 * @return    Boolean array where True indicates finite value
 */
NDArray* ndarray_isfinite(NDArray* arr);

/**
 * Test element-wise for positive or negative infinity.
 *
 * @param arr Input array
 * @return    Boolean array where True indicates infinity
 */
NDArray* ndarray_isinf(NDArray* arr);

/**
 * Test element-wise for NaN.
 *
 * @param arr Input array
 * @return    Boolean array where True indicates NaN
 */
NDArray* ndarray_isnan(NDArray* arr);

/**
 * Test element-wise for negative infinity.
 *
 * @param arr Input array
 * @return    Boolean array where True indicates negative infinity
 */
NDArray* ndarray_isneginf(NDArray* arr);

/**
 * Test element-wise for positive infinity.
 *
 * @param arr Input array
 * @return    Boolean array where True indicates positive infinity
 */
NDArray* ndarray_isposinf(NDArray* arr);

/**
 * Test element-wise if imaginary part is non-zero.
 * For non-complex types, returns all False.
 *
 * @param arr Input array
 * @return    Boolean array
 */
NDArray* ndarray_iscomplex_elem(NDArray* arr);

/**
 * Test element-wise if imaginary part is zero.
 * For non-complex types, returns all True.
 *
 * @param arr Input array
 * @return    Boolean array
 */
NDArray* ndarray_isreal_elem(NDArray* arr);

/* ============ Reductions ============ */

/**
 * Test if all elements evaluate to True.
 *
 * @param arr Input array
 * @return    0-d boolean array (True if all non-zero)
 */
NDArray* ndarray_all(NDArray* arr);

/**
 * Test if all elements along an axis evaluate to True.
 *
 * @param arr      Input array
 * @param axis     Axis to reduce along
 * @param keepdims If non-zero, keep reduced axis as size-1 dimension
 * @return         Reduced boolean array
 */
NDArray* ndarray_all_axis(NDArray* arr, int32_t axis, int keepdims);

/**
 * Test if any element evaluates to True.
 *
 * @param arr Input array
 * @return    0-d boolean array (True if any non-zero)
 */
NDArray* ndarray_any(NDArray* arr);

/**
 * Test if any element along an axis evaluates to True.
 *
 * @param arr      Input array
 * @param axis     Axis to reduce along
 * @param keepdims If non-zero, keep reduced axis as size-1 dimension
 * @return         Reduced boolean array
 */
NDArray* ndarray_any_axis(NDArray* arr, int32_t axis, int keepdims);

/* ============ Comparison Functions ============ */

/**
 * Element-wise comparison with tolerance.
 *
 * Formula: |a - b| <= (atol + rtol * |b|)
 *
 * Broadcasts inputs to common shape.
 *
 * @param a         First array
 * @param b         Second array (reference)
 * @param rtol      Relative tolerance
 * @param atol      Absolute tolerance
 * @param equal_nan If non-zero, NaN == NaN
 * @return          Boolean array of element-wise close status
 */
NDArray* ndarray_isclose(NDArray* a, NDArray* b, double rtol, double atol, int equal_nan);

/**
 * Check if all elements are close within tolerance.
 *
 * @param a         First array
 * @param b         Second array
 * @param rtol      Relative tolerance
 * @param atol      Absolute tolerance
 * @param equal_nan If non-zero, NaN == NaN
 * @return          1 if all close, 0 otherwise
 */
int ndarray_allclose(NDArray* a, NDArray* b, double rtol, double atol, int equal_nan);

/**
 * Check if arrays have same shape and all elements equal.
 *
 * @param a1        First array
 * @param a2        Second array
 * @param equal_nan If non-zero, NaN == NaN
 * @return          1 if equal, 0 otherwise
 */
int ndarray_array_equal(NDArray* a1, NDArray* a2, int equal_nan);

/**
 * Check if arrays are equivalent (allowing broadcasting).
 *
 * @param a1 First array
 * @param a2 Second array
 * @return   1 if equivalent, 0 otherwise
 */
int ndarray_array_equiv(NDArray* a1, NDArray* a2);

#endif /* NUMJS_LOGIC_H */
