#ifndef NUMJS_NDARRAY_H
#define NUMJS_NDARRAY_H

#include <stdint.h>
#include <stddef.h>
#include <stdbool.h>
#include "dtype.h"

/*
 * NumJS NDArray - A minimal n-dimensional array implementation
 *
 * Inspired by NumPy's PyArrayObject from:
 * numpy/_core/include/numpy/ndarraytypes.h lines 772-819
 *
 * This is a simplified version for WebAssembly that removes
 * Python C API dependencies.
 */

/* ============ Memory Ownership Flags ============ */

#define NDARRAY_OWNDATA       0x0001  /* Array owns its data buffer */
#define NDARRAY_WRITEABLE     0x0002  /* Array data is writeable */
#define NDARRAY_C_CONTIGUOUS  0x0004  /* Data is C-contiguous (row-major) */
#define NDARRAY_F_CONTIGUOUS  0x0008  /* Data is Fortran-contiguous (column-major) */
#define NDARRAY_ALIGNED       0x0010  /* Data is properly aligned */

/* ============ NDArray Structure ============ */

/*
 * NDArray structure - simplified from PyArrayObject_fields
 *
 * Original NumPy struct has:
 *   - PyObject_HEAD (Python object header)
 *   - data, nd, dimensions, strides
 *   - base (for views)
 *   - descr (dtype descriptor object)
 *   - flags, weakreflist
 *
 * We keep only the essential fields for array operations.
 */
typedef struct NDArray_s {
    void* data;              /* Pointer to raw data buffer */
    int32_t ndim;            /* Number of dimensions */
    int32_t* shape;          /* Size in each dimension */
    int32_t* strides;        /* Bytes to jump in each dimension */
    DType dtype;             /* Data type */
    int32_t flags;           /* Memory ownership and layout flags */
    size_t size;             /* Total number of elements */
    struct NDArray_s* base;  /* Base array if this is a view (NULL if owns data) */
} NDArray;

/* ============ Array Creation ============ */

/*
 * Create a new NDArray with given shape, initialized to zeros.
 * The array owns its data buffer.
 *
 * @param ndim   Number of dimensions
 * @param shape  Array of dimension sizes (will be copied)
 * @param dtype  Data type for elements
 * @return       New NDArray or NULL on allocation failure
 */
NDArray* ndarray_create(int32_t ndim, int32_t* shape, DType dtype);

/*
 * Create an NDArray from existing data.
 * Data is copied into a new buffer owned by the array.
 *
 * @param data   Source data buffer (will be copied)
 * @param ndim   Number of dimensions
 * @param shape  Array of dimension sizes
 * @param dtype  Data type for elements
 * @return       New NDArray or NULL on allocation failure
 */
NDArray* ndarray_from_data(void* data, int32_t ndim, int32_t* shape, DType dtype);

/*
 * Create a new NDArray with given shape, uninitialized.
 * Faster than create() as it skips zeroing memory.
 *
 * @param ndim   Number of dimensions
 * @param shape  Array of dimension sizes (will be copied)
 * @param dtype  Data type for elements
 * @return       New NDArray or NULL on allocation failure
 */
NDArray* ndarray_empty(int32_t ndim, int32_t* shape, DType dtype);

/*
 * Create a new NDArray filled with a constant value.
 *
 * @param ndim   Number of dimensions
 * @param shape  Array of dimension sizes (will be copied)
 * @param dtype  Data type for elements
 * @param value  Value to fill with
 * @return       New NDArray or NULL on allocation failure
 */
NDArray* ndarray_full(int32_t ndim, int32_t* shape, DType dtype, double value);

/*
 * Create a 0-dimensional (scalar) array.
 *
 * @param value  Scalar value
 * @param dtype  Data type
 * @return       New NDArray or NULL on allocation failure
 */
NDArray* ndarray_scalar(double value, DType dtype);

/*
 * Free an NDArray and its associated memory.
 *
 * @param arr    Array to free (can be NULL)
 */
void ndarray_free(NDArray* arr);

/*
 * Create a deep copy of an array.
 *
 * @param arr    Array to copy
 * @return       New NDArray or NULL on allocation failure
 */
NDArray* ndarray_copy(NDArray* arr);

/*
 * Create a copy with a different dtype.
 *
 * @param arr    Array to convert
 * @param dtype  Target data type
 * @return       New NDArray or NULL on allocation failure
 */
NDArray* ndarray_astype(NDArray* arr, DType dtype);

/* ============ Views ============ */

/*
 * Create a view of an array with different shape/strides.
 * The view shares data with the source array.
 *
 * @param src     Source array
 * @param ndim    Number of dimensions for the view
 * @param shape   Shape for the view
 * @param strides Strides for the view (in bytes)
 * @return        New view or NULL on failure
 */
NDArray* ndarray_view(NDArray* src, int32_t ndim, int32_t* shape, int32_t* strides);

/*
 * Create a view with a byte offset into the data.
 *
 * @param src         Source array
 * @param ndim        Number of dimensions
 * @param shape       Shape for the view
 * @param strides     Strides for the view
 * @param byte_offset Offset into source data
 * @return            New view or NULL on failure
 */
NDArray* ndarray_view_with_offset(NDArray* src, int32_t ndim, int32_t* shape,
                                   int32_t* strides, size_t byte_offset);

/* ============ Shape Manipulation ============ */

/*
 * Reshape array to new shape.
 * Returns a view if the array is contiguous.
 *
 * @param arr       Source array
 * @param new_ndim  Number of dimensions
 * @param new_shape New shape (can include one -1 for auto-calculation)
 * @return          View with new shape or NULL if not possible
 */
NDArray* ndarray_reshape(NDArray* arr, int32_t new_ndim, int32_t* new_shape);

/*
 * Transpose array (reverse axes or custom permutation).
 * Always returns a view.
 *
 * @param arr  Source array
 * @param axes Permutation of axes (NULL for reverse)
 * @return     Transposed view
 */
NDArray* ndarray_transpose(NDArray* arr, int32_t* axes);

/*
 * Flatten array to 1D.
 * Returns view if contiguous, NULL otherwise (caller should use flatten).
 *
 * @param arr  Source array
 * @return     1D view or NULL if copy needed
 */
NDArray* ndarray_ravel(NDArray* arr);

/*
 * Flatten array to 1D (always creates a copy).
 *
 * @param arr  Source array
 * @return     New 1D array
 */
NDArray* ndarray_flatten(NDArray* arr);

/*
 * Remove size-1 dimensions.
 *
 * @param arr  Source array
 * @param axis Specific axis to squeeze (-1 for all size-1 axes)
 * @return     Squeezed view
 */
NDArray* ndarray_squeeze(NDArray* arr, int32_t axis);

/*
 * Add a size-1 dimension at the specified position.
 *
 * @param arr  Source array
 * @param axis Position for new axis (supports negative indexing)
 * @return     Expanded view
 */
NDArray* ndarray_expand_dims(NDArray* arr, int32_t axis);

/*
 * Swap two axes.
 *
 * @param arr   Source array
 * @param axis1 First axis
 * @param axis2 Second axis
 * @return      View with swapped axes
 */
NDArray* ndarray_swapaxes(NDArray* arr, int32_t axis1, int32_t axis2);

/* ============ Element Access ============ */

/*
 * Compute flat index from multi-dimensional indices.
 * Returns SIZE_MAX on error (out of bounds or invalid args).
 *
 * @param arr     Array
 * @param indices Multi-dimensional indices
 * @param ndim    Number of indices
 * @return        Flat index or SIZE_MAX on error
 */
size_t ndarray_flat_index(NDArray* arr, int32_t* indices, int32_t ndim);

/*
 * Check if indices are within bounds.
 *
 * @param arr     Array
 * @param indices Multi-dimensional indices
 * @param ndim    Number of indices
 * @return        true if valid, false otherwise
 */
bool ndarray_check_bounds(NDArray* arr, int32_t* indices, int32_t ndim);

/*
 * Get element value at multi-dimensional indices.
 * Returns value converted to double.
 *
 * @param arr     Array
 * @param indices Multi-dimensional indices
 * @param ndim    Number of indices
 * @return        Element value as double
 */
double ndarray_get_item(NDArray* arr, int32_t* indices, int32_t ndim);

/*
 * Set element value at multi-dimensional indices.
 *
 * @param arr     Array
 * @param indices Multi-dimensional indices
 * @param ndim    Number of indices
 * @param value   Value to set (converted to array's dtype)
 */
void ndarray_set_item(NDArray* arr, int32_t* indices, int32_t ndim, double value);

/*
 * Get element value at flat index.
 *
 * @param arr       Array
 * @param flat_idx  Flat index
 * @return          Element value as double
 */
double ndarray_get_flat(NDArray* arr, size_t flat_idx);

/*
 * Set element value at flat index.
 *
 * @param arr       Array
 * @param flat_idx  Flat index
 * @param value     Value to set
 */
void ndarray_set_flat(NDArray* arr, size_t flat_idx, double value);

/*
 * Get complex element (real part) at flat index.
 */
double ndarray_get_complex_real(NDArray* arr, size_t flat_idx);

/*
 * Get complex element (imaginary part) at flat index.
 */
double ndarray_get_complex_imag(NDArray* arr, size_t flat_idx);

/*
 * Set complex element at flat index.
 */
void ndarray_set_complex(NDArray* arr, size_t flat_idx, double real, double imag);

/* ============ Array Operations ============ */

/*
 * Compute the sum of all elements using pairwise summation.
 * Uses NumPy's pairwise summation algorithm for O(lg n) rounding error.
 *
 * @param arr    Input array
 * @return       Sum of all elements as double
 */
double ndarray_sum(NDArray* arr);

/*
 * Fill array with a constant value.
 *
 * @param arr    Array to fill
 * @param value  Value to fill with (converted to array's dtype)
 */
void ndarray_fill(NDArray* arr, double value);

/* ============ Array Properties ============ */

/*
 * Get number of dimensions.
 */
int32_t ndarray_get_ndim(NDArray* arr);

/*
 * Get pointer to shape array.
 */
int32_t* ndarray_get_shape(NDArray* arr);

/*
 * Get pointer to strides array.
 */
int32_t* ndarray_get_strides(NDArray* arr);

/*
 * Get pointer to raw data buffer.
 */
void* ndarray_get_data(NDArray* arr);

/*
 * Get total number of elements.
 */
size_t ndarray_get_size(NDArray* arr);

/*
 * Get data type.
 */
int32_t ndarray_get_dtype(NDArray* arr);

/*
 * Get flags.
 */
int32_t ndarray_get_flags(NDArray* arr);

/*
 * Get base array (NULL if this array owns its data).
 */
NDArray* ndarray_get_base(NDArray* arr);

/* ============ Contiguity ============ */

/*
 * Check if array is C-contiguous (row-major).
 */
bool ndarray_is_c_contiguous(NDArray* arr);

/*
 * Check if array is F-contiguous (column-major).
 */
bool ndarray_is_f_contiguous(NDArray* arr);

/*
 * Update contiguity flags based on current shape/strides.
 */
void ndarray_update_flags(NDArray* arr);

/* ============ Slicing ============ */

/*
 * Index specification types for slicing
 * Adapted from NumPy's mapping.c HAS_* constants
 */
#define INDEX_TYPE_INTEGER  0  /* Integer index (removes dimension) */
#define INDEX_TYPE_SLICE    1  /* Slice (keeps dimension) */
#define INDEX_TYPE_NEWAXIS  2  /* New axis (adds dimension of size 1) */
#define INDEX_TYPE_ELLIPSIS 3  /* Ellipsis (expands to full slices) */

/*
 * Index specification for multi-dimensional slicing.
 * Adapted from NumPy's npy_index_info structure.
 */
typedef struct {
    int32_t type;   /* INDEX_TYPE_* constant */
    int32_t value;  /* Integer index value (if type=INTEGER), or ellipsis length */
    int32_t start;  /* Slice start (if type=SLICE) */
    int32_t stop;   /* Slice stop (if type=SLICE) */
    int32_t step;   /* Slice step (if type=SLICE) */
} IndexSpec;

/*
 * Create a sliced view of an array.
 *
 * Adapted from NumPy's get_view_from_index() in mapping.c lines 865-946.
 *
 * @param arr         Source array
 * @param indices     Array of IndexSpec structs
 * @param num_indices Number of indices
 * @return            New view or NULL on error
 */
NDArray* ndarray_slice(NDArray* arr, IndexSpec* indices, int32_t num_indices);

/*
 * Get a sub-array by integer index along axis 0.
 * Returns a view with ndim-1 dimensions.
 *
 * @param arr   Source array
 * @param index Index along first axis (supports negative indices)
 * @return      View into arr or NULL on error
 */
NDArray* ndarray_get_subarray(NDArray* arr, int32_t index);

/* ============ View Extensions ============ */

/*
 * Create a view with different dtype interpretation.
 * The array must be C-contiguous.
 * Last dimension is adjusted for size difference.
 *
 * @param arr   Source array
 * @param dtype New dtype
 * @return      View with new dtype or NULL on error
 */
NDArray* ndarray_view_dtype(NDArray* arr, DType dtype);

/*
 * Return array as C-contiguous.
 * Returns view if already contiguous, otherwise copies.
 *
 * @param arr   Source array
 * @return      Contiguous array
 */
NDArray* ndarray_ascontiguousarray(NDArray* arr);

/*
 * Return array as Fortran-contiguous.
 * Returns view if already F-contiguous, otherwise copies.
 *
 * @param arr   Source array
 * @return      F-contiguous array
 */
NDArray* ndarray_asfortranarray(NDArray* arr);

/* ============ WASM Memory Helpers ============ */

/*
 * Allocate memory (for use from JavaScript).
 */
void* wasm_malloc(size_t size);

/*
 * Free memory (for use from JavaScript).
 */
void wasm_free(void* ptr);

#endif /* NUMJS_NDARRAY_H */
