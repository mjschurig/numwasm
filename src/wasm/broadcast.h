/**
 * NumJS-WASM Broadcasting Functions
 *
 * Adapted from NumPy's broadcasting implementation in:
 * numpy/_core/src/multiarray/nditer_constr.c
 *
 * Broadcasting allows arrays with different shapes to be used together
 * in element-wise operations by virtually expanding smaller arrays.
 */

#ifndef NUMJS_BROADCAST_H
#define NUMJS_BROADCAST_H

#include "ndarray.h"

/**
 * Maximum number of dimensions supported
 */
#define NPY_MAXDIMS 32

/**
 * Compute the broadcast shape for two shapes.
 *
 * Broadcasting rules (from NumPy):
 * - Shapes are compared element-wise from right to left
 * - Dimensions are compatible if they are equal or one of them is 1
 * - The result dimension is the larger of the two
 *
 * Adapted from NumPy nditer_constr.c lines 1511-1521
 *
 * @param shape1    First shape array
 * @param ndim1     Number of dimensions in shape1
 * @param shape2    Second shape array
 * @param ndim2     Number of dimensions in shape2
 * @param out_shape Output shape array (caller allocates, size >= max(ndim1, ndim2))
 * @param out_ndim  Output: resulting number of dimensions
 * @return          0 on success, -1 on incompatible shapes
 */
int broadcast_shapes(const int32_t* shape1, int32_t ndim1,
                     const int32_t* shape2, int32_t ndim2,
                     int32_t* out_shape, int32_t* out_ndim);

/**
 * Compute the broadcast shape for multiple shapes.
 *
 * @param shapes     Array of shape pointers
 * @param ndims      Array of ndim values
 * @param num_arrays Number of shapes to broadcast
 * @param out_shape  Output shape array
 * @param out_ndim   Output: resulting number of dimensions
 * @return           0 on success, -1 on incompatible shapes
 */
int broadcast_shapes_multi(const int32_t** shapes, const int32_t* ndims,
                           int32_t num_arrays,
                           int32_t* out_shape, int32_t* out_ndim);

/**
 * Compute strides for broadcasting an array to a target shape.
 *
 * The key broadcasting trick: set stride to 0 for dimensions where
 * the array's shape is 1 but the broadcast shape is larger.
 * This makes the pointer not advance, effectively repeating values.
 *
 * Adapted from NumPy nditer_constr.c lines 1594-1615
 *
 * @param arr           Source array
 * @param target_shape  Target broadcast shape
 * @param target_ndim   Target number of dimensions
 * @param out_strides   Output strides array (caller allocates)
 * @return              0 on success, -1 on incompatible shapes
 */
int broadcast_strides(NDArray* arr, const int32_t* target_shape, int32_t target_ndim,
                      int32_t* out_strides);

/**
 * Create a broadcast view of an array.
 *
 * The returned array:
 * - Has the target shape
 * - Has stride=0 for dimensions being broadcast
 * - Is marked as non-writeable (writes would be repeated)
 * - Shares data with the source array
 *
 * @param arr           Source array
 * @param target_shape  Target shape to broadcast to
 * @param target_ndim   Target number of dimensions
 * @return              Broadcast view or NULL on error
 */
NDArray* ndarray_broadcast_to(NDArray* arr, const int32_t* target_shape, int32_t target_ndim);

/**
 * Check if two shapes are broadcast-compatible.
 *
 * @param shape1  First shape
 * @param ndim1   First ndim
 * @param shape2  Second shape
 * @param ndim2   Second ndim
 * @return        1 if compatible, 0 if not
 */
int shapes_are_broadcastable(const int32_t* shape1, int32_t ndim1,
                             const int32_t* shape2, int32_t ndim2);

#endif /* NUMJS_BROADCAST_H */
