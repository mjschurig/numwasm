/**
 * NumJS-WASM Set Operations Header
 *
 * Adapted from NumPy's _arraysetops_impl.py.
 * Provides set operations on arrays: unique, union, intersection, etc.
 */

#ifndef NUMJS_SETOPS_H
#define NUMJS_SETOPS_H

#include "ndarray.h"

/* ============ Isin Algorithm Kind ============ */
#define ISIN_AUTO   0  /* Automatically choose algorithm */
#define ISIN_SORT   1  /* Use sorting-based algorithm */
#define ISIN_TABLE  2  /* Use table lookup (integers only) */

/* ============ UniqueResult Structure ============ */

/**
 * Result structure for unique() with optional return values.
 *
 * Allows returning multiple arrays from a single WASM call:
 * - values: the unique elements (always present)
 * - indices: indices of first occurrence (if return_index)
 * - inverse: indices to reconstruct original (if return_inverse)
 * - counts: count of each unique element (if return_counts)
 */
typedef struct {
    NDArray* values;    /* Sorted unique values (always present) */
    NDArray* indices;   /* First occurrence indices, or NULL */
    NDArray* inverse;   /* Reconstruction indices, or NULL */
    NDArray* counts;    /* Element counts, or NULL */
} UniqueResult;

/* ============ UniqueResult Memory Management ============ */

/**
 * Free a UniqueResult structure and all its arrays.
 *
 * @param result  UniqueResult to free (may be NULL)
 */
void unique_result_free(UniqueResult* result);

/* ============ Unique ============ */

/**
 * Find unique elements of an array.
 *
 * Adapted from NumPy's unique() in _arraysetops_impl.py line 350.
 *
 * @param arr            Input array (will be flattened)
 * @param return_index   If true, return indices of first occurrences
 * @param return_inverse If true, return indices to reconstruct input
 * @param return_counts  If true, return count of each unique element
 * @param equal_nan      If true, treat NaN values as equal
 * @return               UniqueResult structure, or NULL on error
 */
UniqueResult* ndarray_unique(NDArray* arr, bool return_index, bool return_inverse,
                             bool return_counts, bool equal_nan);

/**
 * Return only the unique values (simplified interface).
 *
 * Equivalent to unique(arr)[0] in NumPy.
 *
 * @param arr  Input array (will be flattened)
 * @return     Sorted array of unique values, or NULL on error
 */
NDArray* ndarray_unique_values(NDArray* arr);

/* ============ Set Combination Operations ============ */

/**
 * Find the union of two arrays.
 *
 * Adapted from NumPy's union1d() in _arraysetops_impl.py line 809.
 *
 * Returns sorted unique values that are in either of the input arrays.
 *
 * @param ar1  First input array
 * @param ar2  Second input array
 * @return     Sorted array of unique union values, or NULL on error
 */
NDArray* ndarray_union1d(NDArray* ar1, NDArray* ar2);

/**
 * Find the intersection of two arrays.
 *
 * Adapted from NumPy's intersect1d() in _arraysetops_impl.py line 667.
 *
 * Returns sorted unique values that are in both input arrays.
 *
 * @param ar1             First input array
 * @param ar2             Second input array
 * @param assume_unique   If true, assume inputs are already unique
 * @param return_indices  If true, return indices via out params
 * @param indices1_out    Output for ar1 indices (if return_indices)
 * @param indices2_out    Output for ar2 indices (if return_indices)
 * @return                Sorted array of common values, or NULL on error
 */
NDArray* ndarray_intersect1d(NDArray* ar1, NDArray* ar2, bool assume_unique,
                             bool return_indices, NDArray** indices1_out,
                             NDArray** indices2_out);

/**
 * Find the set difference of two arrays.
 *
 * Adapted from NumPy's setdiff1d() in _arraysetops_impl.py line 843.
 *
 * Returns sorted unique values in ar1 that are not in ar2.
 *
 * @param ar1            First input array
 * @param ar2            Second input array
 * @param assume_unique  If true, assume inputs are already unique
 * @return               Sorted array of difference values, or NULL on error
 */
NDArray* ndarray_setdiff1d(NDArray* ar1, NDArray* ar2, bool assume_unique);

/**
 * Find the set exclusive-or of two arrays.
 *
 * Adapted from NumPy's setxor1d() in _arraysetops_impl.py line 763.
 *
 * Returns sorted unique values in exactly one (not both) of the input arrays.
 *
 * @param ar1            First input array
 * @param ar2            Second input array
 * @param assume_unique  If true, assume inputs are already unique
 * @return               Sorted array of symmetric difference, or NULL on error
 */
NDArray* ndarray_setxor1d(NDArray* ar1, NDArray* ar2, bool assume_unique);

/* ============ Membership Testing ============ */

/**
 * Test whether each element of ar1 is also present in ar2.
 *
 * Adapted from NumPy's isin() in _arraysetops_impl.py line 890.
 *
 * Returns a boolean array of the same shape as ar1.
 *
 * @param ar1     Input array
 * @param ar2     Values to test against
 * @param invert  If true, return ~(ar1 in ar2)
 * @param assume_unique  If true, assume inputs are unique
 * @param kind    Algorithm hint (ISIN_AUTO, ISIN_SORT, ISIN_TABLE)
 * @return        Boolean array of membership results, or NULL on error
 */
NDArray* ndarray_isin(NDArray* ar1, NDArray* ar2, bool invert,
                      bool assume_unique, int32_t kind);

/**
 * Test whether each element of ar1 is also present in ar2 (flat version).
 *
 * Adapted from NumPy's in1d() in _arraysetops_impl.py line 936.
 *
 * Returns a flattened boolean array.
 *
 * @param ar1     Input array (will be flattened)
 * @param ar2     Values to test against
 * @param invert  If true, return ~(ar1 in ar2)
 * @param assume_unique  If true, assume inputs are unique
 * @param kind    Algorithm hint (ISIN_AUTO, ISIN_SORT, ISIN_TABLE)
 * @return        1D boolean array of membership results, or NULL on error
 */
NDArray* ndarray_in1d(NDArray* ar1, NDArray* ar2, bool invert,
                      bool assume_unique, int32_t kind);

#endif /* NUMJS_SETOPS_H */
