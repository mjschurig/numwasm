# Level 2: Views, Slicing & Broadcasting Implementation Plan

Level 2 extends the NumJS-WASM foundation with three major components: View System Extensions, Slicing (basic and advanced), and Broadcasting. This layer enables array indexing, shape-compatible operations, and efficient memory views.

---

## Current State (Level 1 Complete)

```
src/wasm/
├── ndarray.h          # NDArray struct with base pointer for views
├── ndarray.c          # Views, element access, shape manipulation
├── dtype.h/c          # DType utilities
└── pairwise_sum.h/c   # Pairwise summation

src/ts/
├── types.ts           # DType enum, WasmModule interface, flags
├── NDArray.ts         # Factory methods, get/set, reshape, transpose, etc.
├── iterators.ts       # FlatIterator, nditer, ndenumerate, ndindex
├── dtype.ts           # Type promotion utilities
├── wasm-loader.ts     # WASM module loading
└── index.ts           # Public exports
```

**Existing Infrastructure:**
- View creation: `ndarray_view()`, `ndarray_view_with_offset()`
- Element access: `get()`, `set()`, `getFlat()`, `setFlat()`, `item()`
- Shape ops: `reshape()`, `transpose()`, `ravel()`, `flatten()`, `squeeze()`, `expandDims()`, `swapaxes()`
- Contiguity checks: `ndarray_is_c_contiguous()`, `ndarray_is_f_contiguous()`
- Flags: OWNDATA, WRITEABLE, C_CONTIGUOUS, F_CONTIGUOUS, ALIGNED

---

## Level 2 Implementation Tree

```
LEVEL 2: VIEWS, SLICING & BROADCASTING
│
├── 2.1 Slice Object & Parsing (TypeScript)
│   ├── 2.1.1 Slice class (start, stop, step)
│   ├── 2.1.2 slice() factory function
│   ├── 2.1.3 ellipsis and newaxis constants
│   ├── 2.1.4 IndexElement type union
│   ├── 2.1.5 expandEllipsis() utility
│   └── 2.1.6 ParsedIndex interface
│
├── 2.2 Basic Slicing (C + TypeScript)
│   ├── 2.2.1 C: SliceSpec and IndexSpec structs
│   ├── 2.2.2 C: ndarray_slice() → creates sliced view
│   ├── 2.2.3 C: ndarray_get_subarray() → integer index along axis 0
│   ├── 2.2.4 TS: NDArray.slice() method
│   ├── 2.2.5 TS: NDArray.at() method
│   └── 2.2.6 TS: _buildIndexSpecs() helper
│
│   Dependencies: 2.1.*
│
├── 2.3 Broadcasting (C + TypeScript)
│   ├── 2.3.1 C: broadcast_shapes() → compute output shape
│   ├── 2.3.2 C: broadcast_shapes_multi() → multiple arrays
│   ├── 2.3.3 C: broadcast_strides() → compute strides (0 for broadcast dims)
│   ├── 2.3.4 C: ndarray_broadcast_to() → create broadcast view
│   ├── 2.3.5 TS: broadcastShapes() function
│   ├── 2.3.6 TS: broadcastTo() function
│   └── 2.3.7 TS: broadcastArrays() function
│
│   Dependencies: Level 1 views
│
├── 2.4 View System Extensions (C + TypeScript)
│   ├── 2.4.1 C: ndarray_view_dtype() → view with different dtype
│   ├── 2.4.2 C: ndarray_ascontiguousarray() → ensure C-contiguous
│   ├── 2.4.3 C: ndarray_asfortranarray() → ensure F-contiguous
│   ├── 2.4.4 TS: NDArray.view() method
│   ├── 2.4.5 TS: ascontiguousarray() function
│   └── 2.4.6 TS: asfortranarray() function
│
│   Dependencies: Level 1 views
│
├── 2.5 Index Functions (C + TypeScript)
│   ├── 2.5.1 C: ndarray_take() → select by indices along axis
│   ├── 2.5.2 C: ndarray_put() → assign by flat indices
│   ├── 2.5.3 C: ndarray_nonzero() → find non-zero indices
│   ├── 2.5.4 C: ndarray_where() → conditional selection
│   ├── 2.5.5 TS: take(), put(), nonzero(), where() wrappers
│   └── 2.5.6 TS: compress(), extract() functions
│
│   Dependencies: 2.3.* (broadcasting for where)
│
└── 2.6 Index Generation (TypeScript)
    ├── 2.6.1 indices(dimensions) → grid indices
    ├── 2.6.2 diag_indices(n, ndim) → diagonal indices
    ├── 2.6.3 tril_indices(n, k, m) → lower triangle indices
    ├── 2.6.4 triu_indices(n, k, m) → upper triangle indices
    ├── 2.6.5 ravel_multi_index(multi_index, dims)
    └── 2.6.6 unravel_index(indices, shape)

    Dependencies: 2.2.*
```

---

## Detailed Implementation Specifications

### 2.1 Slice Object & Parsing (TypeScript)

**File:** `src/ts/slice.ts` (new file)

```typescript
/**
 * Represents a slice with start:stop:step semantics.
 * Mirrors Python's slice object behavior.
 */
export class Slice {
  readonly start: number | null;
  readonly stop: number | null;
  readonly step: number | null;

  constructor(
    start: number | null = null,
    stop: number | null = null,
    step: number | null = null
  ) {
    this.start = start;
    this.stop = stop;
    this.step = step;

    if (step === 0) {
      throw new Error('slice step cannot be zero');
    }
  }

  /**
   * Compute concrete indices for a dimension of given length.
   * @returns [start, stop, step, length]
   */
  indices(length: number): [number, number, number, number] {
    const step = this.step ?? 1;
    let start: number;
    let stop: number;

    if (step > 0) {
      start = this.start ?? 0;
      stop = this.stop ?? length;

      // Handle negative indices
      if (start < 0) start = Math.max(0, length + start);
      if (stop < 0) stop = Math.max(0, length + stop);

      // Clamp to bounds
      start = Math.min(start, length);
      stop = Math.min(stop, length);
    } else {
      start = this.start ?? length - 1;
      stop = this.stop ?? -length - 1;

      if (start < 0) start = Math.max(-1, length + start);
      if (stop < 0) stop = length + stop;

      start = Math.min(start, length - 1);
      stop = Math.max(stop, -1);
    }

    // Calculate slice length
    let sliceLength: number;
    if (step > 0) {
      sliceLength = Math.max(0, Math.ceil((stop - start) / step));
    } else {
      sliceLength = Math.max(0, Math.ceil((stop - start) / step));
    }

    return [start, stop, step, sliceLength];
  }
}

/**
 * Factory function for cleaner slice syntax.
 * @example
 * slice(1, 5)        // 1:5
 * slice(null, 5)     // :5
 * slice(1, null, 2)  // 1::2
 * slice(null, null, -1) // ::-1 (reverse)
 */
export function slice(
  start: number | null = null,
  stop: number | null = null,
  step: number | null = null
): Slice {
  return new Slice(start, stop, step);
}

/** Sentinel for ellipsis (...) - expands to fill remaining dimensions */
export const ellipsis = Symbol('ellipsis');

/** Sentinel for newaxis - inserts a new dimension of size 1 */
export const newaxis = Symbol('newaxis');

/** Valid index elements in a slice operation */
export type IndexElement =
  | number           // Integer index
  | Slice            // Slice object
  | typeof ellipsis  // Ellipsis (...)
  | typeof newaxis;  // New axis (np.newaxis)

/**
 * Expand ellipsis and validate index dimensions.
 */
export function expandEllipsis(
  indices: IndexElement[],
  ndim: number
): IndexElement[] {
  let numIndices = 0;
  let hasEllipsis = false;

  for (const idx of indices) {
    if (idx === ellipsis) {
      if (hasEllipsis) {
        throw new Error('an index can only have a single ellipsis');
      }
      hasEllipsis = true;
    } else if (idx !== newaxis) {
      numIndices++;
    }
  }

  if (numIndices > ndim) {
    throw new Error(`too many indices for array: array is ${ndim}-dimensional`);
  }

  if (!hasEllipsis) {
    return indices;
  }

  // Expand ellipsis to full slices
  const expanded: IndexElement[] = [];
  const ellipsisLength = ndim - numIndices;

  for (const idx of indices) {
    if (idx === ellipsis) {
      for (let i = 0; i < ellipsisLength; i++) {
        expanded.push(new Slice());
      }
    } else {
      expanded.push(idx);
    }
  }

  return expanded;
}
```

---

### 2.2 Basic Slicing (C)

**File:** `src/wasm/ndarray.h` (additions)

```c
/* ============ Slicing ============ */

/** Slice specification for a single dimension */
typedef struct {
    int32_t start;
    int32_t stop;
    int32_t step;
} SliceSpec;

/** Index specification types */
#define INDEX_TYPE_INTEGER 0
#define INDEX_TYPE_SLICE   1
#define INDEX_TYPE_NEWAXIS 2

/** Index specification for multi-dimensional slicing */
typedef struct {
    int type;        /* INDEX_TYPE_* constant */
    int32_t value;   /* integer index (if type=INTEGER) */
    SliceSpec slice; /* slice parameters (if type=SLICE) */
} IndexSpec;

/**
 * Create a sliced view of an array.
 *
 * @param arr         Source array
 * @param indices     Array of IndexSpec structs
 * @param num_indices Number of indices
 * @param out_shape   Output: resulting shape (caller allocates, size >= 32)
 * @param out_ndim    Output: resulting ndim
 * @return            New view or NULL on error
 */
NDArray* ndarray_slice(NDArray* arr, IndexSpec* indices, int32_t num_indices,
                       int32_t* out_shape, int32_t* out_ndim);

/**
 * Get a sub-array by integer index along axis 0.
 * Returns a view with ndim-1 dimensions.
 *
 * @param arr   Source array
 * @param index Index along first axis (supports negative)
 * @return      View into arr or NULL on error
 */
NDArray* ndarray_get_subarray(NDArray* arr, int32_t index);
```

**File:** `src/wasm/ndarray.c` (additions)

```c
EXPORT NDArray* ndarray_slice(NDArray* arr, IndexSpec* indices, int32_t num_indices,
                               int32_t* out_shape, int32_t* out_ndim)
{
    if (!arr || !indices) return NULL;

    int32_t new_shape[32];
    int32_t new_strides[32];
    int new_ndim = 0;
    size_t byte_offset = 0;

    int arr_axis = 0;  /* Current axis in source array */

    for (int i = 0; i < num_indices; i++) {
        IndexSpec* idx = &indices[i];

        switch (idx->type) {
            case INDEX_TYPE_INTEGER: {
                /* Integer index: reduce dimension, add offset */
                if (arr_axis >= arr->ndim) return NULL;

                int32_t index = idx->value;
                int32_t dim_size = arr->shape[arr_axis];

                /* Handle negative indices */
                if (index < 0) index += dim_size;
                if (index < 0 || index >= dim_size) return NULL;

                byte_offset += (size_t)index * (size_t)arr->strides[arr_axis];
                arr_axis++;
                break;
            }

            case INDEX_TYPE_SLICE: {
                /* Slice: compute new shape/stride for this dimension */
                if (arr_axis >= arr->ndim) return NULL;

                int32_t dim_size = arr->shape[arr_axis];
                int32_t start = idx->slice.start;
                int32_t stop = idx->slice.stop;
                int32_t step = idx->slice.step;

                /* Resolve defaults */
                if (step == 0) step = 1;

                /* Resolve negative indices and clamp to bounds */
                if (step > 0) {
                    if (start < 0) start = (start < -dim_size) ? 0 : dim_size + start;
                    else if (start > dim_size) start = dim_size;

                    if (stop < 0) stop = (stop < -dim_size) ? 0 : dim_size + stop;
                    else if (stop > dim_size) stop = dim_size;
                } else {
                    if (start < 0) start = (start < -dim_size) ? -1 : dim_size + start;
                    else if (start >= dim_size) start = dim_size - 1;

                    if (stop < 0) stop = (stop < -dim_size) ? -1 : dim_size + stop;
                    else if (stop >= dim_size) stop = dim_size;
                }

                /* Calculate slice length */
                int32_t length;
                if (step > 0) {
                    length = (stop > start) ? (stop - start + step - 1) / step : 0;
                } else {
                    length = (start > stop) ? (start - stop - step - 1) / (-step) : 0;
                }

                new_shape[new_ndim] = length;
                new_strides[new_ndim] = arr->strides[arr_axis] * step;
                byte_offset += (size_t)start * (size_t)arr->strides[arr_axis];

                new_ndim++;
                arr_axis++;
                break;
            }

            case INDEX_TYPE_NEWAXIS: {
                /* Insert new dimension of size 1 */
                new_shape[new_ndim] = 1;
                new_strides[new_ndim] = (new_ndim > 0) ?
                    new_strides[new_ndim - 1] : dtype_size(arr->dtype);
                new_ndim++;
                /* Don't advance arr_axis */
                break;
            }

            default:
                return NULL;
        }
    }

    /* Fill remaining dimensions with full slices */
    while (arr_axis < arr->ndim) {
        new_shape[new_ndim] = arr->shape[arr_axis];
        new_strides[new_ndim] = arr->strides[arr_axis];
        new_ndim++;
        arr_axis++;
    }

    /* Return results via output parameters */
    if (out_shape) memcpy(out_shape, new_shape, new_ndim * sizeof(int32_t));
    if (out_ndim) *out_ndim = new_ndim;

    return ndarray_view_with_offset(arr, new_ndim, new_shape, new_strides, byte_offset);
}

EXPORT NDArray* ndarray_get_subarray(NDArray* arr, int32_t index)
{
    if (!arr || arr->ndim == 0) return NULL;

    int32_t idx = index;
    if (idx < 0) idx += arr->shape[0];
    if (idx < 0 || idx >= arr->shape[0]) return NULL;

    size_t offset = (size_t)idx * (size_t)arr->strides[0];

    if (arr->ndim == 1) {
        /* Return 0-d array (scalar view) */
        return ndarray_scalar_view(arr, offset);
    }

    /* Return (ndim-1) dimensional view */
    return ndarray_view_with_offset(arr, arr->ndim - 1,
                                     arr->shape + 1,
                                     arr->strides + 1,
                                     offset);
}
```

---

### 2.3 Broadcasting (C)

**File:** `src/wasm/broadcast.h` (new file)

```c
#ifndef NUMJS_BROADCAST_H
#define NUMJS_BROADCAST_H

#include "ndarray.h"

/**
 * Compute the broadcast shape for two shapes.
 * Returns 0 on success, -1 on incompatible shapes.
 */
int broadcast_shapes(const int32_t* shape1, int32_t ndim1,
                     const int32_t* shape2, int32_t ndim2,
                     int32_t* out_shape, int32_t* out_ndim);

/**
 * Compute the broadcast shape for multiple shapes.
 */
int broadcast_shapes_multi(const int32_t** shapes, const int32_t* ndims,
                            int32_t num_arrays,
                            int32_t* out_shape, int32_t* out_ndim);

/**
 * Compute strides for broadcasting an array to a target shape.
 * Sets stride to 0 for dimensions being broadcast.
 */
int broadcast_strides(NDArray* arr, const int32_t* target_shape, int32_t target_ndim,
                      int32_t* out_strides);

/**
 * Create a broadcast view of an array.
 * The returned array has stride=0 for broadcast dimensions.
 * The view is marked as non-writeable.
 */
NDArray* ndarray_broadcast_to(NDArray* arr, const int32_t* target_shape, int32_t target_ndim);

#endif /* NUMJS_BROADCAST_H */
```

**File:** `src/wasm/broadcast.c` (new file)

```c
#include "broadcast.h"
#include <stdlib.h>
#include <string.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

EXPORT int broadcast_shapes(const int32_t* shape1, int32_t ndim1,
                             const int32_t* shape2, int32_t ndim2,
                             int32_t* out_shape, int32_t* out_ndim)
{
    int32_t result_ndim = (ndim1 > ndim2) ? ndim1 : ndim2;

    /* Work backwards from the trailing dimensions */
    for (int i = 0; i < result_ndim; i++) {
        int idx1 = ndim1 - 1 - i;
        int idx2 = ndim2 - 1 - i;

        int32_t dim1 = (idx1 >= 0) ? shape1[idx1] : 1;
        int32_t dim2 = (idx2 >= 0) ? shape2[idx2] : 1;

        if (dim1 == dim2) {
            out_shape[result_ndim - 1 - i] = dim1;
        } else if (dim1 == 1) {
            out_shape[result_ndim - 1 - i] = dim2;
        } else if (dim2 == 1) {
            out_shape[result_ndim - 1 - i] = dim1;
        } else {
            /* Incompatible shapes */
            return -1;
        }
    }

    *out_ndim = result_ndim;
    return 0;
}

EXPORT int broadcast_strides(NDArray* arr, const int32_t* target_shape, int32_t target_ndim,
                              int32_t* out_strides)
{
    if (!arr || !target_shape || !out_strides) return -1;

    int diff = target_ndim - arr->ndim;

    for (int i = 0; i < target_ndim; i++) {
        int arr_idx = i - diff;

        if (arr_idx < 0) {
            /* Prepended dimension: stride = 0 (broadcast) */
            out_strides[i] = 0;
        } else if (arr->shape[arr_idx] == 1) {
            /* Size-1 dimension being broadcast: stride = 0 */
            out_strides[i] = 0;
        } else if (arr->shape[arr_idx] == target_shape[i]) {
            /* Matching dimension: use original stride */
            out_strides[i] = arr->strides[arr_idx];
        } else {
            /* Incompatible shapes */
            return -1;
        }
    }

    return 0;
}

EXPORT NDArray* ndarray_broadcast_to(NDArray* arr, const int32_t* target_shape, int32_t target_ndim)
{
    if (!arr || !target_shape) return NULL;

    /* Validate broadcast compatibility */
    int diff = target_ndim - arr->ndim;
    if (diff < 0) return NULL;

    for (int i = 0; i < arr->ndim; i++) {
        int target_idx = i + diff;
        if (arr->shape[i] != 1 && arr->shape[i] != target_shape[target_idx]) {
            return NULL;
        }
    }

    /* Compute broadcast strides */
    int32_t new_strides[32];
    if (broadcast_strides(arr, target_shape, target_ndim, new_strides) != 0) {
        return NULL;
    }

    /* Create view with broadcast shape and strides */
    int32_t* shape_copy = malloc(target_ndim * sizeof(int32_t));
    if (!shape_copy) return NULL;
    memcpy(shape_copy, target_shape, target_ndim * sizeof(int32_t));

    NDArray* view = ndarray_view(arr, target_ndim, shape_copy, new_strides);
    free(shape_copy);

    /* Mark as not writeable since writes would be repeated */
    if (view) {
        view->flags &= ~NDARRAY_WRITEABLE;
    }

    return view;
}
```

---

### 2.4 View System Extensions (C)

**File:** `src/wasm/ndarray.h` (additions)

```c
/**
 * Create a view with different dtype interpretation.
 * Must be C-contiguous. Last dimension is adjusted for size difference.
 */
NDArray* ndarray_view_dtype(NDArray* arr, DType dtype);

/**
 * Return array as C-contiguous.
 * Returns view if already contiguous, otherwise copies.
 */
NDArray* ndarray_ascontiguousarray(NDArray* arr);

/**
 * Return array as Fortran-contiguous.
 * Returns view if already F-contiguous, otherwise copies.
 */
NDArray* ndarray_asfortranarray(NDArray* arr);
```

---

### 2.5 Index Functions (C)

**File:** `src/wasm/indexing.h` (new file)

```c
#ifndef NUMJS_INDEXING_H
#define NUMJS_INDEXING_H

#include "ndarray.h"

/**
 * Take elements from an array along an axis.
 *
 * @param arr       Source array
 * @param indices   Array of indices to take
 * @param n_indices Number of indices
 * @param axis      Axis along which to take (negative for from end)
 * @return          New array with selected elements
 */
NDArray* ndarray_take(NDArray* arr, const int32_t* indices, size_t n_indices, int32_t axis);

/**
 * Replace elements at flat indices with values.
 *
 * @param arr     Array to modify (in place)
 * @param indices Flat indices
 * @param values  Values to place
 * @param n       Number of indices/values
 */
void ndarray_put(NDArray* arr, const size_t* indices, const double* values, size_t n);

/**
 * Find indices of non-zero elements.
 *
 * @param arr         Source array
 * @param out_indices Output: array of indices per dimension
 * @param out_count   Output: number of non-zero elements
 * @return            0 on success, -1 on error
 */
int ndarray_nonzero(NDArray* arr, int32_t*** out_indices, size_t* out_count);

/**
 * Return elements chosen from x or y depending on condition.
 * All arrays are broadcast together.
 */
NDArray* ndarray_where(NDArray* condition, NDArray* x, NDArray* y);

#endif /* NUMJS_INDEXING_H */
```

---

### 2.2-2.6 TypeScript Implementations

**File:** `src/ts/NDArray.ts` (additions to class)

```typescript
import { Slice, slice, ellipsis, newaxis, IndexElement, expandEllipsis } from './slice.js';

// Inside NDArray class:

/**
 * Create a sliced view of the array.
 *
 * @param indices - Array of index specifications
 * @returns Sliced view
 *
 * @example
 * arr.slice([1])                    // Second row
 * arr.slice([slice(1, 3)])          // Rows 1-2
 * arr.slice([slice(), 2])           // Third column
 * arr.slice([slice(null, null, 2)]) // Every other row
 * arr.slice([newaxis, slice()])     // Add dimension
 * arr.slice([ellipsis, 0])          // Last axis, first element
 */
slice(indices: IndexElement[]): NDArray {
  this.ensureNotDisposed();

  const expanded = expandEllipsis(indices, this.ndim);
  const specs = this._buildIndexSpecs(expanded);

  // Allocate and call C function...
  // (implementation details as in spec above)
}

/**
 * Get element or sub-array at index along first axis.
 * Returns scalar for 1D arrays, view for higher dimensions.
 */
at(index: number): number | NDArray {
  this.ensureNotDisposed();

  if (this.ndim === 0) {
    throw new Error('Cannot index 0-dimensional array');
  }

  const resultPtr = this._module._ndarray_get_subarray(this._ptr, index);
  if (resultPtr === 0) {
    throw new RangeError(
      `Index ${index} out of bounds for axis 0 with size ${this.shape[0]}`
    );
  }

  const view = new NDArray(resultPtr, this._module);

  // If result is 0-d (scalar), return the value
  if (view.ndim === 0) {
    const val = view.item();
    view.dispose();
    return val;
  }

  return view;
}

/**
 * View array with different dtype interpretation.
 */
view(dtype: DType): NDArray {
  this.ensureNotDisposed();

  const resultPtr = this._module._ndarray_view_dtype(this._ptr, dtype);
  if (resultPtr === 0) {
    throw new Error('Cannot create dtype view: array must be contiguous');
  }

  return new NDArray(resultPtr, this._module);
}
```

**File:** `src/ts/broadcast.ts` (new file)

```typescript
import { NDArray } from './NDArray.js';

/**
 * Compute broadcast shape for multiple shapes.
 * @throws If shapes are not broadcastable
 */
export function broadcastShapes(...shapes: number[][]): number[] {
  if (shapes.length === 0) return [];

  let result = shapes[0].slice();
  for (let i = 1; i < shapes.length; i++) {
    result = _broadcastTwo(result, shapes[i]);
  }
  return result;
}

function _broadcastTwo(shape1: number[], shape2: number[]): number[] {
  const ndim = Math.max(shape1.length, shape2.length);
  const result: number[] = new Array(ndim);

  for (let i = 0; i < ndim; i++) {
    const dim1 = shape1[shape1.length - 1 - i] ?? 1;
    const dim2 = shape2[shape2.length - 1 - i] ?? 1;

    if (dim1 === dim2) {
      result[ndim - 1 - i] = dim1;
    } else if (dim1 === 1) {
      result[ndim - 1 - i] = dim2;
    } else if (dim2 === 1) {
      result[ndim - 1 - i] = dim1;
    } else {
      throw new Error(
        `operands could not be broadcast together with shapes ` +
        `(${shape1.join(',')}) (${shape2.join(',')})`
      );
    }
  }

  return result;
}

/**
 * Broadcast an array to a new shape.
 * Returns a read-only view with stride=0 for broadcast dimensions.
 */
export function broadcastTo(arr: NDArray, shape: number[]): NDArray {
  // ... implementation calling _ndarray_broadcast_to
}

/**
 * Broadcast multiple arrays to a common shape.
 * Returns array of views all having the same shape.
 */
export function broadcastArrays(...arrays: NDArray[]): NDArray[] {
  if (arrays.length === 0) return [];

  const shapes = arrays.map(arr => arr.shape);
  const targetShape = broadcastShapes(...shapes);

  return arrays.map(arr => {
    if (arraysEqual(arr.shape, targetShape)) {
      return arr;
    }
    return broadcastTo(arr, targetShape);
  });
}
```

**File:** `src/ts/indexing.ts` (new file)

```typescript
import { NDArray } from './NDArray.js';

/**
 * Take elements from an array along an axis.
 */
export function take(arr: NDArray, indices: number[], axis: number = 0): NDArray {
  // ... implementation calling _ndarray_take
}

/**
 * Replace elements at flat indices.
 */
export function put(arr: NDArray, indices: number[], values: number[]): void {
  // ... implementation calling _ndarray_put
}

/**
 * Return indices of non-zero elements.
 * Returns tuple of arrays, one per dimension.
 */
export function nonzero(arr: NDArray): NDArray[] {
  // ... implementation calling _ndarray_nonzero
}

/**
 * Return elements chosen from x or y depending on condition.
 */
export function where(condition: NDArray, x: NDArray, y: NDArray): NDArray {
  // ... implementation calling _ndarray_where
}

/**
 * Return selected slices of array along axis where condition is True.
 */
export function compress(condition: NDArray, arr: NDArray, axis?: number): NDArray {
  // ... implementation
}

/**
 * Return elements of array that satisfy condition.
 */
export function extract(condition: NDArray, arr: NDArray): NDArray {
  // ... implementation
}
```

---

## File Changes Summary

### New Files to Create

```
src/wasm/
├── broadcast.h        # Broadcasting declarations
├── broadcast.c        # Broadcasting implementation
├── indexing.h         # Index functions declarations
└── indexing.c         # Index functions implementation

src/ts/
├── slice.ts           # Slice class, ellipsis, newaxis, parsing
├── broadcast.ts       # broadcastTo, broadcastArrays, broadcastShapes
└── indexing.ts        # take, put, nonzero, where, compress, extract
```

### Files to Modify

```
src/wasm/ndarray.h
├── Add SliceSpec, IndexSpec structs
├── Add INDEX_TYPE_* constants
├── Declare ndarray_slice()
├── Declare ndarray_get_subarray()
├── Declare ndarray_view_dtype()
├── Declare ndarray_ascontiguousarray()
└── Declare ndarray_asfortranarray()

src/wasm/ndarray.c
├── Implement ndarray_slice()
├── Implement ndarray_get_subarray()
├── Implement ndarray_view_dtype()
├── Implement ndarray_ascontiguousarray()
└── Implement ndarray_asfortranarray()

src/ts/NDArray.ts
├── Import from slice.ts
├── Add slice() method
├── Add at() method
└── Add view() method

src/ts/types.ts
├── Add IndexSpec interface
└── Add new WASM function declarations

src/ts/index.ts
├── Export Slice, slice, ellipsis, newaxis
├── Export broadcastTo, broadcastArrays, broadcastShapes
└── Export take, put, nonzero, where, compress, extract

scripts/build-wasm.sh
├── Add broadcast.c to compilation
├── Add indexing.c to compilation
└── Add new EXPORTED_FUNCTIONS
```

---

## Build Script Updates

Add to `scripts/build-wasm.sh` EXPORTED_FUNCTIONS:

```bash
# Slicing
"_ndarray_slice",
"_ndarray_get_subarray",

# Broadcasting
"_broadcast_shapes",
"_broadcast_shapes_multi",
"_broadcast_strides",
"_ndarray_broadcast_to",

# View extensions
"_ndarray_view_dtype",
"_ndarray_ascontiguousarray",
"_ndarray_asfortranarray",

# Index functions
"_ndarray_take",
"_ndarray_put",
"_ndarray_nonzero",
"_ndarray_where"
```

Add new source files:

```bash
emcc \
    "$SRC_DIR/ndarray.c" \
    "$SRC_DIR/dtype.c" \
    "$SRC_DIR/pairwise_sum.c" \
    "$SRC_DIR/broadcast.c" \
    "$SRC_DIR/indexing.c" \
    -o "$OUT_DIR/numjs.cjs" \
    ...
```

---

## WasmModule Interface Updates

Add to `src/ts/types.ts`:

```typescript
// Slicing
_ndarray_slice(ptr: number, indicesPtr: number, numIndices: number,
               outShapePtr: number, outNdimPtr: number): number;
_ndarray_get_subarray(ptr: number, index: number): number;

// Broadcasting
_broadcast_shapes(shape1Ptr: number, ndim1: number,
                  shape2Ptr: number, ndim2: number,
                  outShapePtr: number, outNdimPtr: number): number;
_broadcast_shapes_multi(shapesPtr: number, ndimsPtr: number, numArrays: number,
                        outShapePtr: number, outNdimPtr: number): number;
_broadcast_strides(ptr: number, targetShapePtr: number, targetNdim: number,
                   outStridesPtr: number): number;
_ndarray_broadcast_to(ptr: number, targetShapePtr: number, targetNdim: number): number;

// View extensions
_ndarray_view_dtype(ptr: number, dtype: number): number;
_ndarray_ascontiguousarray(ptr: number): number;
_ndarray_asfortranarray(ptr: number): number;

// Index functions
_ndarray_take(ptr: number, indicesPtr: number, nIndices: number, axis: number): number;
_ndarray_put(ptr: number, indicesPtr: number, valuesPtr: number, n: number): void;
_ndarray_nonzero(ptr: number, outIndicesPtr: number, outCountPtr: number): number;
_ndarray_where(condPtr: number, xPtr: number, yPtr: number): number;
```

---

## Implementation Order

```
Week 1: Slice Infrastructure
├── Day 1: Slice class and TypeScript utilities (slice.ts)
├── Day 2: Index parsing (expandEllipsis, _buildIndexSpecs)
├── Day 3: C: IndexSpec struct, ndarray_slice()
├── Day 4: C: ndarray_get_subarray()
└── Day 5: TypeScript slice() and at() methods + tests

Week 2: Broadcasting
├── Day 1: C: broadcast_shapes(), broadcast_shapes_multi()
├── Day 2: C: broadcast_strides(), ndarray_broadcast_to()
├── Day 3: TypeScript broadcast.ts functions
├── Day 4: Integration tests
└── Day 5: Edge cases and documentation

Week 3: View Extensions & Index Functions
├── Day 1: C: ndarray_view_dtype()
├── Day 2: C: ndarray_ascontiguousarray(), ndarray_asfortranarray()
├── Day 3: C: ndarray_take(), ndarray_put()
├── Day 4: C: ndarray_nonzero(), ndarray_where()
└── Day 5: TypeScript wrappers (indexing.ts)

Week 4: Index Generation & Polish
├── Day 1: TS: indices(), diag_indices()
├── Day 2: TS: tril_indices(), triu_indices()
├── Day 3: TS: ravel_multi_index(), unravel_index()
├── Day 4: Comprehensive test suite
└── Day 5: Documentation, examples, cleanup
```

---

## Verification Plan

After Level 2 completion, verify:

```bash
# Build
npm run build

# Run all tests
npm test

# Specific Level 2 tests should pass:

# Slicing
✓ slice() creates views with correct shape/strides
✓ slice() handles negative indices
✓ slice() handles negative step (reverse)
✓ slice() with ellipsis expands correctly
✓ slice() with newaxis adds dimensions
✓ at() returns sub-array or scalar
✓ Sliced views share data with base

# Broadcasting
✓ broadcastShapes() computes correct output shape
✓ broadcastShapes() throws on incompatible shapes
✓ broadcastTo() creates view with stride=0 for broadcast dims
✓ broadcastTo() result is read-only
✓ broadcastArrays() broadcasts multiple arrays to common shape

# View Extensions
✓ view(dtype) reinterprets data correctly
✓ ascontiguousarray() returns view if already contiguous
✓ ascontiguousarray() returns copy if non-contiguous
✓ asfortranarray() returns F-ordered array

# Index Functions
✓ take() selects elements by indices
✓ take() works along specified axis
✓ put() assigns elements at flat indices
✓ nonzero() returns indices of non-zero elements
✓ where() selects from x or y based on condition
✓ where() broadcasts inputs correctly
```

Generate NumPy comparison vectors:

```python
# tests/python/generate_level2_tests.py
import numpy as np
import json

tests = {
    "slicing": [
        {"shape": [5], "slice": "1:4", "expected_shape": [3]},
        {"shape": [5], "slice": "::2", "expected_shape": [3]},
        {"shape": [5], "slice": "::-1", "expected_shape": [5]},
        {"shape": [3, 4], "slice": "1", "expected_shape": [4]},
        {"shape": [3, 4], "slice": ":, 2", "expected_shape": [3]},
    ],
    "broadcasting": [
        {"shapes": [[3, 1], [1, 4]], "expected": [3, 4]},
        {"shapes": [[2, 3], [3]], "expected": [2, 3]},
        {"shapes": [[1], [5, 4]], "expected": [5, 4]},
    ],
    "take": [
        {"shape": [5], "indices": [0, 2, 4], "axis": 0, "expected_shape": [3]},
        {"shape": [3, 4], "indices": [0, 2], "axis": 1, "expected_shape": [3, 2]},
    ],
    "nonzero": [
        {"data": [0, 1, 0, 2, 0], "expected_indices": [[1, 3]]},
        {"data": [[1, 0], [0, 2]], "expected_indices": [[0, 1], [0, 1]]},
    ],
}

with open("tests/fixtures/level2_vectors.json", "w") as f:
    json.dump(tests, f, indent=2)
```

---

## Critical Dependencies for Level 3

Level 2 completion enables:

- **Level 3.1 Ufunc Infrastructure**: Broadcasting is essential for binary operations
- **Level 3.2 Element-wise Operations**: Requires broadcasting for mismatched shapes
- **Level 3.3 Reductions with Axis**: Can leverage slicing for axis iteration

Level 2 MUST be complete before proceeding to Level 3.
