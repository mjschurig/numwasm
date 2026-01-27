/**
 * NumJS-WASM Index Functions Header
 *
 * Adapted from NumPy's item_selection.c and multiarraymodule.c.
 * Provides index-based array manipulation functions.
 */

#ifndef NUMJS_INDEXING_H
#define NUMJS_INDEXING_H

#include "ndarray.h"

/* ============ Clip Modes ============ */
/* Adapted from NumPy's NPY_CLIPMODE */
#define CLIP_RAISE  0  /* Raise error on out-of-bounds */
#define CLIP_WRAP   1  /* Wrap around (modulo) */
#define CLIP_CLIP   2  /* Clip to valid range */

/* ============ Take ============ */

/**
 * Take elements from an array along an axis.
 *
 * Adapted from NumPy's PyArray_TakeFrom() in item_selection.c line 232.
 *
 * @param arr       Source array
 * @param indices   Array of indices (must be integer type)
 * @param axis      Axis along which to take (negative allowed)
 * @param clipmode  How to handle out-of-bounds indices (CLIP_*)
 * @return          New array with taken elements, or NULL on error
 */
NDArray* ndarray_take(NDArray* arr, NDArray* indices, int32_t axis, int32_t clipmode);

/**
 * Take elements from flattened array.
 *
 * @param arr       Source array (will be flattened conceptually)
 * @param indices   Array of indices
 * @param clipmode  How to handle out-of-bounds indices
 * @return          New 1D array with taken elements
 */
NDArray* ndarray_take_flat(NDArray* arr, NDArray* indices, int32_t clipmode);

/* ============ Put ============ */

/**
 * Put values into an array at specified flat indices.
 *
 * Adapted from NumPy's PyArray_PutTo() in item_selection.c line 369.
 *
 * @param arr       Target array (modified in place)
 * @param indices   Array of flat indices
 * @param values    Values to put (broadcast if smaller)
 * @param clipmode  How to handle out-of-bounds indices
 * @return          0 on success, -1 on error
 */
int ndarray_put(NDArray* arr, NDArray* indices, NDArray* values, int32_t clipmode);

/* ============ Nonzero ============ */

/**
 * Count nonzero elements in array.
 *
 * @param arr  Input array
 * @return     Number of nonzero elements
 */
size_t ndarray_count_nonzero(NDArray* arr);

/**
 * Find indices of nonzero elements.
 *
 * Adapted from NumPy's PyArray_Nonzero() in item_selection.c line 2806.
 *
 * Returns a 2D array of shape (num_nonzero, ndim) where each row
 * contains the multi-dimensional index of a nonzero element.
 *
 * @param arr   Input array
 * @return      2D index array, or NULL on error
 */
NDArray* ndarray_nonzero(NDArray* arr);

/**
 * Return indices of nonzero elements in a flattened array.
 *
 * @param arr   Input array
 * @return      1D array of flat indices where arr != 0
 */
NDArray* ndarray_flatnonzero(NDArray* arr);

/* ============ Where ============ */

/**
 * Return elements chosen from x or y depending on condition.
 *
 * Adapted from NumPy's PyArray_Where() in multiarraymodule.c.
 *
 * @param condition Boolean condition array
 * @param x         Values where condition is true
 * @param y         Values where condition is false
 * @return          New array with selected values, or NULL on error
 *
 * Note: Arrays are broadcast together. If x and y are NULL,
 * returns same as nonzero(condition).
 */
NDArray* ndarray_where(NDArray* condition, NDArray* x, NDArray* y);

/* ============ Compress ============ */

/**
 * Return selected slices of an array along given axis.
 *
 * @param condition Boolean 1D array selecting which slices to keep
 * @param arr       Input array
 * @param axis      Axis along which to select
 * @return          Compressed array
 */
NDArray* ndarray_compress(NDArray* condition, NDArray* arr, int32_t axis);

/* ============ Extract ============ */

/**
 * Return elements of an array that satisfy some condition.
 *
 * @param condition Boolean array, same shape as arr
 * @param arr       Input array
 * @return          1D array of elements where condition is true
 */
NDArray* ndarray_extract(NDArray* condition, NDArray* arr);

/* ============ Choose ============ */

/**
 * Construct an array from an index array and a set of arrays to choose from.
 *
 * @param indices   Array of indices into choices
 * @param choices   Array of choice arrays (num_choices length)
 * @param num_choices Number of choice arrays
 * @param clipmode  How to handle out-of-bounds indices
 * @return          New array with chosen values
 */
NDArray* ndarray_choose(NDArray* indices, NDArray** choices, int32_t num_choices,
                        int32_t clipmode);

/* ============ Diagonal Indexing ============ */

/**
 * Return specified diagonals.
 *
 * @param arr    Input array (2D or higher)
 * @param offset Offset from main diagonal
 * @param axis1  First axis of 2D sub-arrays
 * @param axis2  Second axis of 2D sub-arrays
 * @return       Array of diagonals
 */
NDArray* ndarray_diagonal(NDArray* arr, int32_t offset, int32_t axis1, int32_t axis2);

#endif /* NUMJS_INDEXING_H */
