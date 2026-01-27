#ifndef NUMJS_MANIPULATION_H
#define NUMJS_MANIPULATION_H

#include "ndarray.h"

/**
 * Concatenate arrays along an existing axis.
 *
 * @param arrays    Array of NDArray pointers
 * @param n_arrays  Number of arrays
 * @param axis      Axis along which to concatenate (supports negative)
 * @return          New concatenated array or NULL on error
 */
NDArray* ndarray_concatenate(NDArray** arrays, int32_t n_arrays, int32_t axis);

#endif /* NUMJS_MANIPULATION_H */
