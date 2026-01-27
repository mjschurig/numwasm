#ifndef NUMJS_SORTING_H
#define NUMJS_SORTING_H

#include "ndarray.h"

/* Sort Kind Constants */
#define SORT_QUICKSORT  0  /* Introsort (quicksort + heapsort fallback) */
#define SORT_MERGESORT  1  /* Stable merge sort */
#define SORT_HEAPSORT   2  /* Heapsort */
#define SORT_STABLE     3  /* Stable sort (alias for mergesort) */

/**
 * Sort array in-place along specified axis.
 *
 * @param arr   Array to sort (modified in place)
 * @param axis  Axis along which to sort (-1 for last axis)
 * @param kind  Sorting algorithm (SORT_* constant)
 * @return      0 on success, -1 on error
 */
int ndarray_sort(NDArray* arr, int32_t axis, int32_t kind);

/**
 * Return a sorted copy of array.
 *
 * @param arr   Source array
 * @param axis  Axis along which to sort (-1 for last, INT32_MIN to flatten)
 * @param kind  Sorting algorithm
 * @return      New sorted array or NULL on error
 */
NDArray* ndarray_sort_copy(NDArray* arr, int32_t axis, int32_t kind);

/**
 * Return indices that would sort array along axis.
 *
 * @param arr   Source array
 * @param axis  Axis along which to sort (-1 for last, INT32_MIN to flatten)
 * @param kind  Sorting algorithm
 * @return      Array of indices (Int32) or NULL on error
 */
NDArray* ndarray_argsort(NDArray* arr, int32_t axis, int32_t kind);

/**
 * Return a partitioned copy of array.
 * Elements before kth are smaller, elements after are larger.
 *
 * @param arr   Source array
 * @param kth   Index to partition around
 * @param axis  Axis along which to partition
 * @return      Partitioned array or NULL on error
 */
NDArray* ndarray_partition(NDArray* arr, int32_t kth, int32_t axis);

/**
 * Return indices that would partition array.
 */
NDArray* ndarray_argpartition(NDArray* arr, int32_t kth, int32_t axis);

#endif /* NUMJS_SORTING_H */
