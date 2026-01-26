#ifndef NUMJS_SEARCHING_H
#define NUMJS_SEARCHING_H

#include "ndarray.h"

#define SEARCH_LEFT   0
#define SEARCH_RIGHT  1

/**
 * Return indices of maximum values along axis.
 *
 * @param arr      Source array
 * @param axis     Axis to reduce (-1 for last, INT32_MIN to flatten)
 * @param keepdims If true, reduced axis is left with size 1
 * @return         Array of indices (Int32) or NULL on error
 */
NDArray* ndarray_argmax(NDArray* arr, int32_t axis, bool keepdims);

/**
 * Return indices of minimum values along axis.
 *
 * @param arr      Source array
 * @param axis     Axis to reduce (-1 for last, INT32_MIN to flatten)
 * @param keepdims If true, reduced axis is left with size 1
 * @return         Array of indices (Int32) or NULL on error
 */
NDArray* ndarray_argmin(NDArray* arr, int32_t axis, bool keepdims);

/**
 * Find indices where elements should be inserted to maintain order.
 *
 * @param sorted_arr Sorted 1D array
 * @param values     Values to insert
 * @param side       SEARCH_LEFT or SEARCH_RIGHT
 * @param sorter     Optional indices that sort the array (can be NULL)
 * @return           Array of indices or NULL on error
 */
NDArray* ndarray_searchsorted(NDArray* sorted_arr, NDArray* values,
                               int32_t side, NDArray* sorter);

#endif /* NUMJS_SEARCHING_H */
