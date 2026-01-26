# Phase 5: Array Manipulation Implementation Plan

Phase 5 extends NumJS-WASM with array manipulation functions: joining, splitting, tiling, and rearranging arrays. This enables building complex data structures from simpler arrays.

---

## Current State (Phases 1-4 Complete)

```
src/wasm/
├── ndarray.h/c        # Core array, views, slicing, shape manipulation
├── dtype.h/c          # DType utilities
├── broadcast.h/c      # Broadcasting
├── indexing.h/c       # Index operations (take, put, nonzero, where)
└── pairwise_sum.h/c   # Accurate summation

src/ts/
├── NDArray.ts         # Core class with reshape, transpose, slice, etc.
├── types.ts           # DType enum, WasmModule interface
├── broadcast.ts       # broadcastTo, broadcastArrays, broadcastShapes
├── indexing.ts        # take, put, nonzero, where, atleast_Nd, indices, etc.
├── slice.ts           # Slice class, ellipsis, newaxis
├── iterators.ts       # FlatIterator, nditer, ndenumerate
├── dtype.ts           # Type promotion utilities
└── index.ts           # Public exports
```

**Existing Utilities:**
- Views: `ndarray_view()`, `ndarray_view_with_offset()`
- Shape ops: `reshape()`, `transpose()`, `squeeze()`, `expandDims()`, `swapaxes()`
- Memory: `compute_strides()`, `compute_size()`, `flat_to_byte_offset()`
- Contiguity: `ndarray_is_c_contiguous()`, `ndarray_is_f_contiguous()`
- Broadcasting: `broadcastShapes()`, `broadcastTo()`, `broadcastArrays()`
- atleast_Nd: `atleast_1d()`, `atleast_2d()`, `atleast_3d()` (in indexing.ts)

---

## Phase 5 Implementation Tree

```
PHASE 5: ARRAY MANIPULATION
│
├── 5.1 Joining Arrays (C + TypeScript)
│   ├── 5.1.1 C: ndarray_concatenate() → join along axis
│   ├── 5.1.2 TS: concatenate(arrays, axis, dtype)
│   ├── 5.1.3 TS: stack(arrays, axis) → insert new axis
│   ├── 5.1.4 TS: vstack(tup) → vertical stack (row_stack alias)
│   ├── 5.1.5 TS: hstack(tup) → horizontal stack
│   ├── 5.1.6 TS: dstack(tup) → depth stack (3rd axis)
│   ├── 5.1.7 TS: column_stack(tup) → stack as columns
│   ├── 5.1.8 TS: block(nested_arrays) → assemble from blocks
│   └── 5.1.9 TS: append(arr, values, axis)
│
│   Dependencies: atleast_Nd, broadcasting
│
├── 5.2 Splitting Arrays (C + TypeScript)
│   ├── 5.2.1 C: ndarray_split() → split into sections
│   ├── 5.2.2 TS: split(arr, indices_or_sections, axis)
│   ├── 5.2.3 TS: array_split(arr, indices_or_sections, axis)
│   ├── 5.2.4 TS: vsplit(arr, indices_or_sections)
│   ├── 5.2.5 TS: hsplit(arr, indices_or_sections)
│   ├── 5.2.6 TS: dsplit(arr, indices_or_sections)
│   └── 5.2.7 TS: unstack(arr, axis)
│
│   Dependencies: slicing, swapaxes
│
├── 5.3 Tiling & Repeating (C + TypeScript)
│   ├── 5.3.1 C: ndarray_tile() → repeat array by reps
│   ├── 5.3.2 C: ndarray_repeat() → repeat elements
│   ├── 5.3.3 TS: tile(arr, reps)
│   ├── 5.3.4 TS: repeat(arr, repeats, axis)
│   └── 5.3.5 TS: pad(arr, pad_width, mode, ...)
│
│   Dependencies: concatenate, reshape
│
├── 5.4 Rearranging (C + TypeScript)
│   ├── 5.4.1 C: ndarray_flip() → reverse along axis
│   ├── 5.4.2 C: ndarray_roll() → shift elements
│   ├── 5.4.3 TS: flip(arr, axis)
│   ├── 5.4.4 TS: fliplr(arr), flipud(arr)
│   ├── 5.4.5 TS: roll(arr, shift, axis)
│   ├── 5.4.6 TS: rot90(arr, k, axes)
│   ├── 5.4.7 TS: resize(arr, new_shape)
│   ├── 5.4.8 TS: trim_zeros(arr, trim)
│   ├── 5.4.9 TS: insert(arr, obj, values, axis)
│   └── 5.4.10 TS: delete(arr, obj, axis)
│
│   Dependencies: slicing, concatenate
│
└── 5.5 Copying Extensions (TypeScript)
    ├── 5.5.1 TS: copyto(dst, src, casting, where)
    └── 5.5.2 TS: asarray(a, dtype, order)

    Dependencies: broadcasting, dtype conversion
```

---

## Detailed Implementation Specifications

### 5.1 Joining Arrays

#### 5.1.1 C: ndarray_concatenate()

**File:** `src/wasm/manipulation.h` (new file)

```c
#ifndef NUMJS_MANIPULATION_H
#define NUMJS_MANIPULATION_H

#include "ndarray.h"

/**
 * Concatenate arrays along an existing axis.
 *
 * @param arrays    Array of NDArray pointers
 * @param n_arrays  Number of arrays
 * @param axis      Axis along which to concatenate
 * @return          New concatenated array or NULL on error
 */
NDArray* ndarray_concatenate(NDArray** arrays, int32_t n_arrays, int32_t axis);

/**
 * Tile an array by repeating it along each axis.
 *
 * @param arr   Source array
 * @param reps  Number of repetitions along each axis
 * @param ndim  Number of repetition dimensions
 * @return      Tiled array
 */
NDArray* ndarray_tile(NDArray* arr, int32_t* reps, int32_t ndim);

/**
 * Repeat elements of an array.
 *
 * @param arr      Source array
 * @param repeats  Number of repetitions for each element (or single value)
 * @param n_reps   Length of repeats array (1 for broadcast)
 * @param axis     Axis along which to repeat (-1 for flattened)
 * @return         Array with repeated elements
 */
NDArray* ndarray_repeat(NDArray* arr, int32_t* repeats, int32_t n_reps, int32_t axis);

/**
 * Flip array along specified axis.
 *
 * @param arr   Source array
 * @param axis  Axis to flip (-1 for all axes)
 * @return      View with flipped axis (uses negative stride)
 */
NDArray* ndarray_flip(NDArray* arr, int32_t axis);

/**
 * Roll array elements along axis.
 *
 * @param arr    Source array
 * @param shift  Number of positions to shift
 * @param axis   Axis along which to roll (-1 for flattened)
 * @return       Rolled array (copy)
 */
NDArray* ndarray_roll(NDArray* arr, int32_t shift, int32_t axis);

/**
 * Split array into multiple sub-arrays.
 *
 * @param arr           Source array
 * @param indices       Split points or number of sections
 * @param n_indices     Length of indices array
 * @param axis          Axis along which to split
 * @param out_arrays    Output: array of NDArray pointers (caller allocates)
 * @param out_count     Output: number of resulting arrays
 * @return              0 on success, -1 on error
 */
int ndarray_split(NDArray* arr, int32_t* indices, int32_t n_indices,
                  int32_t axis, NDArray*** out_arrays, int32_t* out_count);

#endif /* NUMJS_MANIPULATION_H */
```

**File:** `src/wasm/manipulation.c` (new file)

```c
#include "manipulation.h"
#include <stdlib.h>
#include <string.h>

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

EXPORT NDArray* ndarray_concatenate(NDArray** arrays, int32_t n_arrays, int32_t axis)
{
    if (!arrays || n_arrays == 0) return NULL;

    NDArray* first = arrays[0];
    if (!first) return NULL;

    /* Normalize axis */
    int32_t ndim = first->ndim;
    if (axis < 0) axis += ndim;
    if (axis < 0 || axis >= ndim) return NULL;

    /* Validate all arrays have same shape except along axis */
    int32_t total_axis_size = first->shape[axis];
    for (int i = 1; i < n_arrays; i++) {
        if (!arrays[i] || arrays[i]->ndim != ndim) return NULL;

        for (int d = 0; d < ndim; d++) {
            if (d != axis && arrays[i]->shape[d] != first->shape[d]) {
                return NULL;  /* Shape mismatch */
            }
        }
        total_axis_size += arrays[i]->shape[axis];
    }

    /* Create output shape */
    int32_t* out_shape = malloc(ndim * sizeof(int32_t));
    if (!out_shape) return NULL;
    memcpy(out_shape, first->shape, ndim * sizeof(int32_t));
    out_shape[axis] = total_axis_size;

    /* Create output array */
    NDArray* result = ndarray_create(ndim, out_shape, first->dtype);
    free(out_shape);
    if (!result) return NULL;

    /* Copy data from each input array */
    size_t axis_offset = 0;
    for (int i = 0; i < n_arrays; i++) {
        NDArray* arr = arrays[i];

        /* For each position in the array, copy to result */
        size_t arr_size = arr->size;
        for (size_t flat = 0; flat < arr_size; flat++) {
            /* Convert flat index to multi-index */
            int32_t multi_idx[32];
            size_t remainder = flat;
            for (int d = ndim - 1; d >= 0; d--) {
                multi_idx[d] = remainder % arr->shape[d];
                remainder /= arr->shape[d];
            }

            /* Adjust axis index for destination */
            multi_idx[axis] += axis_offset;

            /* Get value and set in result */
            double val = ndarray_get_flat(arr, flat);

            /* Convert multi-index to flat index in result */
            size_t result_flat = 0;
            size_t multiplier = 1;
            for (int d = ndim - 1; d >= 0; d--) {
                result_flat += multi_idx[d] * multiplier;
                multiplier *= result->shape[d];
            }

            ndarray_set_flat(result, result_flat, val);
        }

        axis_offset += arr->shape[axis];
    }

    return result;
}

EXPORT NDArray* ndarray_flip(NDArray* arr, int32_t axis)
{
    if (!arr) return NULL;

    /* Normalize axis */
    int32_t ndim = arr->ndim;
    if (axis < -ndim || axis >= ndim) return NULL;
    if (axis < 0) axis += ndim;

    /* Create new strides with negative stride for flipped axis */
    int32_t* new_strides = malloc(ndim * sizeof(int32_t));
    if (!new_strides) return NULL;
    memcpy(new_strides, arr->strides, ndim * sizeof(int32_t));
    new_strides[axis] = -arr->strides[axis];

    /* Compute byte offset to start at end of axis */
    size_t byte_offset = (arr->shape[axis] - 1) * arr->strides[axis];

    NDArray* view = ndarray_view_with_offset(arr, ndim, arr->shape, new_strides, byte_offset);
    free(new_strides);

    return view;
}

EXPORT NDArray* ndarray_roll(NDArray* arr, int32_t shift, int32_t axis)
{
    if (!arr) return NULL;

    int32_t ndim = arr->ndim;

    /* Handle flattened case */
    if (axis < -ndim || axis >= ndim) {
        /* Flatten, roll, reshape */
        NDArray* flat = ndarray_flatten(arr);
        if (!flat) return NULL;

        int32_t n = flat->size;
        shift = ((shift % n) + n) % n;  /* Normalize */

        NDArray* result = ndarray_empty(1, &n, flat->dtype);
        if (!result) {
            ndarray_free(flat);
            return NULL;
        }

        for (int32_t i = 0; i < n; i++) {
            int32_t src_idx = (i - shift + n) % n;
            ndarray_set_flat(result, i, ndarray_get_flat(flat, src_idx));
        }

        ndarray_free(flat);

        /* Reshape back to original shape */
        NDArray* final = ndarray_reshape(result, arr->ndim, arr->shape);
        ndarray_free(result);
        return final;
    }

    if (axis < 0) axis += ndim;

    int32_t axis_size = arr->shape[axis];
    shift = ((shift % axis_size) + axis_size) % axis_size;

    if (shift == 0) return ndarray_copy(arr);

    /* Create result array */
    NDArray* result = ndarray_empty(ndim, arr->shape, arr->dtype);
    if (!result) return NULL;

    /* Copy with shift */
    for (size_t flat = 0; flat < arr->size; flat++) {
        /* Convert to multi-index */
        int32_t multi_idx[32];
        size_t remainder = flat;
        for (int d = ndim - 1; d >= 0; d--) {
            multi_idx[d] = remainder % arr->shape[d];
            remainder /= arr->shape[d];
        }

        /* Compute source index with shift */
        int32_t src_axis_idx = (multi_idx[axis] - shift + axis_size) % axis_size;
        multi_idx[axis] = src_axis_idx;

        /* Convert back to flat */
        size_t src_flat = 0;
        size_t multiplier = 1;
        for (int d = ndim - 1; d >= 0; d--) {
            src_flat += multi_idx[d] * multiplier;
            multiplier *= arr->shape[d];
        }

        ndarray_set_flat(result, flat, ndarray_get_flat(arr, src_flat));
    }

    return result;
}
```

---

### 5.1.2-5.1.9 TypeScript: Stack Functions

**File:** `src/ts/manipulation.ts` (new file)

```typescript
import { NDArray } from './NDArray.js';
import { DType } from './types.js';
import { atleast_1d, atleast_2d, atleast_3d } from './indexing.js';
import { broadcastShapes } from './broadcast.js';

/**
 * Join arrays along an existing axis.
 */
export async function concatenate(
  arrays: NDArray[],
  axis: number = 0,
  dtype?: DType
): Promise<NDArray> {
  if (arrays.length === 0) {
    throw new Error('need at least one array to concatenate');
  }

  const module = arrays[0]._wasmModule;

  // Allocate array of pointers
  const ptrsPtr = module._malloc(arrays.length * 4);
  for (let i = 0; i < arrays.length; i++) {
    module.setValue(ptrsPtr + i * 4, arrays[i]._wasmPtr, 'i32');
  }

  const resultPtr = module._ndarray_concatenate(ptrsPtr, arrays.length, axis);
  module._free(ptrsPtr);

  if (resultPtr === 0) {
    throw new Error('concatenate failed: incompatible shapes or invalid axis');
  }

  let result = NDArray._fromPtr(resultPtr, module);

  if (dtype !== undefined && result.dtype !== dtype) {
    const converted = result.astype(dtype);
    result.dispose();
    result = converted;
  }

  return result;
}

/**
 * Join arrays along a new axis.
 */
export async function stack(
  arrays: NDArray[],
  axis: number = 0
): Promise<NDArray> {
  if (arrays.length === 0) {
    throw new Error('need at least one array to stack');
  }

  // Validate all arrays have same shape
  const shape = arrays[0].shape;
  for (let i = 1; i < arrays.length; i++) {
    if (!arraysEqual(arrays[i].shape, shape)) {
      throw new Error('all input arrays must have the same shape');
    }
  }

  // Normalize axis
  const ndim = arrays[0].ndim + 1;
  axis = axis < 0 ? axis + ndim : axis;
  if (axis < 0 || axis > arrays[0].ndim) {
    throw new Error(`axis ${axis} out of bounds`);
  }

  // Expand each array along the new axis
  const expanded = arrays.map(arr => arr.expandDims(axis));

  // Concatenate along the new axis
  const result = await concatenate(expanded, axis);

  // Clean up expanded views
  for (const exp of expanded) {
    exp.dispose();
  }

  return result;
}

/**
 * Stack arrays vertically (row-wise).
 */
export async function vstack(tup: NDArray[]): Promise<NDArray> {
  // atleast_2d ensures 1D arrays become (1, N)
  const arrs = atleast_2d(...tup);
  const arrays = Array.isArray(arrs) ? arrs : [arrs];
  return concatenate(arrays, 0);
}

// Alias
export const row_stack = vstack;

/**
 * Stack arrays horizontally (column-wise).
 */
export async function hstack(tup: NDArray[]): Promise<NDArray> {
  const arrs = atleast_1d(...tup);
  const arrays = Array.isArray(arrs) ? arrs : [arrs];

  // For 1D arrays, concatenate along axis 0
  // For 2D+ arrays, concatenate along axis 1
  if (arrays[0].ndim === 1) {
    return concatenate(arrays, 0);
  }
  return concatenate(arrays, 1);
}

/**
 * Stack arrays along the third axis (depth).
 */
export async function dstack(tup: NDArray[]): Promise<NDArray> {
  const arrs = atleast_3d(...tup);
  const arrays = Array.isArray(arrs) ? arrs : [arrs];
  return concatenate(arrays, 2);
}

/**
 * Stack 1D arrays as columns into a 2D array.
 */
export async function column_stack(tup: NDArray[]): Promise<NDArray> {
  const arrays: NDArray[] = [];

  for (const arr of tup) {
    if (arr.ndim < 2) {
      // 1D to column: [N] -> [N, 1]
      arrays.push(arr.reshape([arr.size, 1]));
    } else {
      arrays.push(arr);
    }
  }

  return concatenate(arrays, 1);
}

/**
 * Assemble arrays from nested lists of blocks.
 */
export async function block(arrays: (NDArray | NDArray[])[]): Promise<NDArray> {
  // Handle simple 1D list case
  if (arrays.every(arr => arr instanceof NDArray)) {
    return hstack(arrays as NDArray[]);
  }

  // Handle 2D block case: list of lists
  const rows: NDArray[] = [];
  for (const row of arrays) {
    if (Array.isArray(row)) {
      rows.push(await hstack(row));
    } else {
      rows.push(row);
    }
  }

  return vstack(rows);
}

/**
 * Append values to the end of an array.
 */
export async function append(
  arr: NDArray,
  values: NDArray,
  axis?: number
): Promise<NDArray> {
  if (axis === undefined) {
    // Flatten both and concatenate
    const flat1 = arr.ravel();
    const flat2 = values.ravel();
    const result = await concatenate([flat1, flat2], 0);
    flat1.dispose();
    flat2.dispose();
    return result;
  }

  return concatenate([arr, values], axis);
}

function arraysEqual(a: number[], b: number[]): boolean {
  if (a.length !== b.length) return false;
  for (let i = 0; i < a.length; i++) {
    if (a[i] !== b[i]) return false;
  }
  return true;
}
```

---

### 5.2 Splitting Arrays

**Add to `src/ts/manipulation.ts`:**

```typescript
/**
 * Split array into multiple sub-arrays.
 * Raises error if array cannot be split into equal sections.
 */
export async function split(
  arr: NDArray,
  indices_or_sections: number | number[],
  axis: number = 0
): Promise<NDArray[]> {
  // Normalize axis
  axis = axis < 0 ? axis + arr.ndim : axis;
  if (axis < 0 || axis >= arr.ndim) {
    throw new Error(`axis ${axis} out of bounds for array with ${arr.ndim} dimensions`);
  }

  const axisSize = arr.shape[axis];

  if (typeof indices_or_sections === 'number') {
    // Number of sections
    const sections = indices_or_sections;
    if (axisSize % sections !== 0) {
      throw new Error('array split does not result in an equal division');
    }
    return array_split(arr, indices_or_sections, axis);
  }

  return array_split(arr, indices_or_sections, axis);
}

/**
 * Split array into multiple sub-arrays.
 * Allows unequal division (unlike split).
 */
export async function array_split(
  arr: NDArray,
  indices_or_sections: number | number[],
  axis: number = 0
): Promise<NDArray[]> {
  // Normalize axis
  axis = axis < 0 ? axis + arr.ndim : axis;

  const axisSize = arr.shape[axis];
  let divPoints: number[];

  if (typeof indices_or_sections === 'number') {
    // Number of sections - compute division points
    const sections = indices_or_sections;
    const [base, extras] = [Math.floor(axisSize / sections), axisSize % sections];

    // First 'extras' sections get base+1 elements
    divPoints = [0];
    for (let i = 0; i < sections; i++) {
      divPoints.push(divPoints[i] + base + (i < extras ? 1 : 0));
    }
  } else {
    // Explicit indices
    divPoints = [0, ...indices_or_sections, axisSize];
  }

  // Create sub-arrays using slicing
  const result: NDArray[] = [];
  const { slice } = await import('./slice.js');

  for (let i = 0; i < divPoints.length - 1; i++) {
    const start = divPoints[i];
    const stop = divPoints[i + 1];

    // Build slice indices
    const indices: any[] = [];
    for (let d = 0; d < arr.ndim; d++) {
      if (d === axis) {
        indices.push(slice(start, stop));
      } else {
        indices.push(slice(null, null));  // Full slice
      }
    }

    const subArr = arr.slice(indices);
    result.push(subArr.copy());  // Copy to own data
    subArr.dispose();
  }

  return result;
}

/**
 * Split array into sub-arrays vertically (row-wise).
 */
export function vsplit(
  arr: NDArray,
  indices_or_sections: number | number[]
): Promise<NDArray[]> {
  if (arr.ndim < 2) {
    throw new Error('vsplit only works on arrays of 2 or more dimensions');
  }
  return split(arr, indices_or_sections, 0);
}

/**
 * Split array into sub-arrays horizontally (column-wise).
 */
export function hsplit(
  arr: NDArray,
  indices_or_sections: number | number[]
): Promise<NDArray[]> {
  if (arr.ndim < 1) {
    throw new Error('hsplit only works on arrays of 1 or more dimensions');
  }
  if (arr.ndim > 1) {
    return split(arr, indices_or_sections, 1);
  }
  return split(arr, indices_or_sections, 0);
}

/**
 * Split array into sub-arrays along the 3rd axis (depth).
 */
export function dsplit(
  arr: NDArray,
  indices_or_sections: number | number[]
): Promise<NDArray[]> {
  if (arr.ndim < 3) {
    throw new Error('dsplit only works on arrays of 3 or more dimensions');
  }
  return split(arr, indices_or_sections, 2);
}

/**
 * Unpack array along specified axis (inverse of stack).
 */
export async function unstack(arr: NDArray, axis: number = 0): Promise<NDArray[]> {
  // Normalize axis
  axis = axis < 0 ? axis + arr.ndim : axis;

  if (arr.ndim === 0) {
    throw new Error('cannot unstack 0-d array');
  }

  const result: NDArray[] = [];
  for (let i = 0; i < arr.shape[axis]; i++) {
    // Use at() for axis 0, otherwise slice
    if (axis === 0) {
      const sub = arr.at(i);
      result.push(sub.copy());
      sub.dispose();
    } else {
      const { slice } = await import('./slice.js');
      const indices: any[] = [];
      for (let d = 0; d < arr.ndim; d++) {
        indices.push(d === axis ? i : slice(null, null));
      }
      const sub = arr.slice(indices);
      result.push(sub.copy());
      sub.dispose();
    }
  }

  return result;
}
```

---

### 5.3 Tiling & Repeating

**Add to `src/ts/manipulation.ts`:**

```typescript
/**
 * Construct array by repeating arr the number of times given by reps.
 */
export async function tile(arr: NDArray, reps: number | number[]): Promise<NDArray> {
  const repsArr = Array.isArray(reps) ? reps : [reps];

  // Extend reps to match or exceed arr.ndim
  const d = Math.max(arr.ndim, repsArr.length);

  // Pad reps with 1s at the beginning if needed
  const paddedReps = new Array(d - repsArr.length).fill(1).concat(repsArr);

  // Pad array shape with 1s at the beginning if needed
  let tiled = arr;
  if (arr.ndim < d) {
    const newShape = new Array(d - arr.ndim).fill(1).concat(arr.shape);
    tiled = arr.reshape(newShape);
  }

  // Compute output shape
  const outShape = tiled.shape.map((s, i) => s * paddedReps[i]);

  // Create result
  const result = await NDArray.empty(outShape, { dtype: arr.dtype });

  // Fill by iterating over repetition grid
  for (let outFlat = 0; outFlat < result.size; outFlat++) {
    // Convert to multi-index
    const outIdx = flatToMulti(outFlat, outShape);

    // Map to source index (modulo original shape)
    const srcIdx = outIdx.map((idx, d) => idx % tiled.shape[d]);
    const srcFlat = multiToFlat(srcIdx, tiled.shape);

    result.setFlat(outFlat, tiled.getFlat(srcFlat));
  }

  if (tiled !== arr) {
    tiled.dispose();
  }

  return result;
}

/**
 * Repeat elements of an array.
 */
export async function repeat(
  arr: NDArray,
  repeats: number | number[],
  axis?: number
): Promise<NDArray> {
  if (axis === undefined) {
    // Flatten and repeat
    const flat = arr.ravel();
    const result = await repeatAlongAxis(flat, repeats, 0);
    flat.dispose();
    return result;
  }

  return repeatAlongAxis(arr, repeats, axis);
}

async function repeatAlongAxis(
  arr: NDArray,
  repeats: number | number[],
  axis: number
): Promise<NDArray> {
  axis = axis < 0 ? axis + arr.ndim : axis;

  const axisSize = arr.shape[axis];
  const repsArr = typeof repeats === 'number'
    ? new Array(axisSize).fill(repeats)
    : repeats;

  if (repsArr.length !== axisSize) {
    throw new Error('repeats must have same length as axis');
  }

  const totalRepeats = repsArr.reduce((a, b) => a + b, 0);

  // Compute output shape
  const outShape = [...arr.shape];
  outShape[axis] = totalRepeats;

  const result = await NDArray.empty(outShape, { dtype: arr.dtype });

  // Use swapaxes trick: move axis to front
  const swapped = arr.swapaxes(axis, 0);
  const outSwapped = result.swapaxes(axis, 0);

  let outIdx = 0;
  for (let i = 0; i < axisSize; i++) {
    const slice = swapped.at(i);
    const numReps = repsArr[i];

    for (let r = 0; r < numReps; r++) {
      // Copy slice to output
      const outSlice = outSwapped.at(outIdx);
      copyArrayData(slice, outSlice);
      outSlice.dispose();
      outIdx++;
    }
    slice.dispose();
  }

  swapped.dispose();
  outSwapped.dispose();

  return result;
}

/**
 * Pad an array.
 * Modes: 'constant', 'edge', 'reflect', 'symmetric', 'wrap'
 */
export async function pad(
  arr: NDArray,
  pad_width: number | [number, number] | [number, number][],
  mode: string = 'constant',
  constant_values: number = 0
): Promise<NDArray> {
  // Normalize pad_width to [[before, after], ...] format
  let padSpec: [number, number][];

  if (typeof pad_width === 'number') {
    padSpec = arr.shape.map(() => [pad_width, pad_width]);
  } else if (Array.isArray(pad_width) && typeof pad_width[0] === 'number') {
    const [before, after] = pad_width as [number, number];
    padSpec = arr.shape.map(() => [before, after]);
  } else {
    padSpec = pad_width as [number, number][];
  }

  if (padSpec.length !== arr.ndim) {
    throw new Error('pad_width must have same length as array dimensions');
  }

  // Compute output shape
  const outShape = arr.shape.map((s, i) => s + padSpec[i][0] + padSpec[i][1]);

  // Create padded array with constant fill
  const result = await NDArray.full(outShape, constant_values, { dtype: arr.dtype });

  // Build slice to place original data
  const { slice } = await import('./slice.js');
  const indices: any[] = padSpec.map(([before, _], i) =>
    slice(before, before + arr.shape[i])
  );

  // Copy original data into the center
  const target = result.slice(indices);
  copyArrayData(arr, target);
  target.dispose();

  // TODO: Handle edge modes (edge, reflect, symmetric, wrap)

  return result;
}

function copyArrayData(src: NDArray, dst: NDArray): void {
  for (let i = 0; i < src.size; i++) {
    dst.setFlat(i, src.getFlat(i));
  }
}
```

---

### 5.4 Rearranging

**Add to `src/ts/manipulation.ts`:**

```typescript
/**
 * Reverse the order of elements along the given axis.
 */
export function flip(arr: NDArray, axis?: number | number[]): NDArray {
  if (axis === undefined) {
    // Flip all axes
    const axes = Array.from({ length: arr.ndim }, (_, i) => i);
    return flipMulti(arr, axes);
  }

  if (typeof axis === 'number') {
    return flipSingle(arr, axis);
  }

  return flipMulti(arr, axis);
}

function flipSingle(arr: NDArray, axis: number): NDArray {
  axis = axis < 0 ? axis + arr.ndim : axis;

  // Use slicing with step -1
  const { slice } = require('./slice.js');
  const indices: any[] = [];
  for (let d = 0; d < arr.ndim; d++) {
    if (d === axis) {
      indices.push(slice(null, null, -1));  // Reverse this axis
    } else {
      indices.push(slice(null, null));  // Full slice
    }
  }

  return arr.slice(indices);
}

function flipMulti(arr: NDArray, axes: number[]): NDArray {
  let result = arr;
  for (const axis of axes) {
    const flipped = flipSingle(result, axis);
    if (result !== arr) result.dispose();
    result = flipped;
  }
  return result;
}

/**
 * Flip array left to right (axis 1).
 */
export function fliplr(arr: NDArray): NDArray {
  if (arr.ndim < 2) {
    throw new Error('fliplr requires array with at least 2 dimensions');
  }
  return flip(arr, 1);
}

/**
 * Flip array up to down (axis 0).
 */
export function flipud(arr: NDArray): NDArray {
  if (arr.ndim < 1) {
    throw new Error('flipud requires array with at least 1 dimension');
  }
  return flip(arr, 0);
}

/**
 * Roll array elements along axis.
 */
export async function roll(
  arr: NDArray,
  shift: number | number[],
  axis?: number | number[]
): Promise<NDArray> {
  if (axis === undefined) {
    // Roll flattened array
    const flat = arr.flatten();
    const rolled = await rollSingle(flat, shift as number, 0);
    const result = rolled.reshape(arr.shape);
    flat.dispose();
    rolled.dispose();
    return result;
  }

  if (typeof axis === 'number' && typeof shift === 'number') {
    return rollSingle(arr, shift, axis);
  }

  // Multiple axes
  const shifts = Array.isArray(shift) ? shift : [shift];
  const axes = Array.isArray(axis) ? axis : [axis];

  if (shifts.length !== axes.length) {
    throw new Error('shift and axis must have same length');
  }

  let result = arr;
  for (let i = 0; i < axes.length; i++) {
    const rolled = await rollSingle(result, shifts[i], axes[i]);
    if (result !== arr) result.dispose();
    result = rolled;
  }

  return result;
}

async function rollSingle(arr: NDArray, shift: number, axis: number): Promise<NDArray> {
  axis = axis < 0 ? axis + arr.ndim : axis;
  const axisSize = arr.shape[axis];

  // Normalize shift
  shift = ((shift % axisSize) + axisSize) % axisSize;
  if (shift === 0) return arr.copy();

  // Split and concatenate
  const { slice } = await import('./slice.js');

  const indices1: any[] = [];
  const indices2: any[] = [];
  for (let d = 0; d < arr.ndim; d++) {
    if (d === axis) {
      indices1.push(slice(axisSize - shift, null));
      indices2.push(slice(null, axisSize - shift));
    } else {
      indices1.push(slice(null, null));
      indices2.push(slice(null, null));
    }
  }

  const part1 = arr.slice(indices1);
  const part2 = arr.slice(indices2);
  const result = await concatenate([part1, part2], axis);
  part1.dispose();
  part2.dispose();

  return result;
}

/**
 * Rotate array by 90 degrees in the plane specified by axes.
 */
export function rot90(
  arr: NDArray,
  k: number = 1,
  axes: [number, number] = [0, 1]
): NDArray {
  if (arr.ndim < 2) {
    throw new Error('rot90 requires array with at least 2 dimensions');
  }

  k = ((k % 4) + 4) % 4;  // Normalize to 0-3

  const [ax1, ax2] = axes;

  if (k === 0) return arr.copy();

  if (k === 1) {
    // Transpose then flip
    return flipSingle(arr.swapaxes(ax1, ax2), ax1);
  }

  if (k === 2) {
    // Flip both axes
    return flipSingle(flipSingle(arr, ax1), ax2);
  }

  // k === 3
  return flipSingle(arr.swapaxes(ax1, ax2), ax2);
}

/**
 * Return array with new shape, repeating or truncating as needed.
 */
export async function resize(arr: NDArray, new_shape: number[]): Promise<NDArray> {
  const newSize = new_shape.reduce((a, b) => a * b, 1);
  const result = await NDArray.empty(new_shape, { dtype: arr.dtype });

  // Fill by cycling through original data
  for (let i = 0; i < newSize; i++) {
    result.setFlat(i, arr.getFlat(i % arr.size));
  }

  return result;
}

/**
 * Trim leading and/or trailing zeros from a 1D array.
 */
export function trim_zeros(arr: NDArray, trim: string = 'fb'): NDArray {
  if (arr.ndim !== 1) {
    throw new Error('trim_zeros only works on 1D arrays');
  }

  const data = arr.toArray();
  let start = 0;
  let end = data.length;

  if (trim.includes('f')) {
    while (start < end && data[start] === 0) start++;
  }

  if (trim.includes('b')) {
    while (end > start && data[end - 1] === 0) end--;
  }

  const { slice } = require('./slice.js');
  return arr.slice([slice(start, end)]);
}

// Helper functions
function flatToMulti(flatIdx: number, shape: number[]): number[] {
  const result = new Array(shape.length);
  let remainder = flatIdx;
  for (let i = shape.length - 1; i >= 0; i--) {
    result[i] = remainder % shape[i];
    remainder = Math.floor(remainder / shape[i]);
  }
  return result;
}

function multiToFlat(multiIdx: number[], shape: number[]): number {
  let flat = 0;
  let multiplier = 1;
  for (let i = shape.length - 1; i >= 0; i--) {
    flat += multiIdx[i] * multiplier;
    multiplier *= shape[i];
  }
  return flat;
}
```

---

### 5.5 Copying Extensions

**Add to `src/ts/manipulation.ts`:**

```typescript
/**
 * Copy values from src to dst.
 */
export async function copyto(
  dst: NDArray,
  src: NDArray,
  casting: string = 'same_kind',
  where?: NDArray
): Promise<void> {
  if (where === undefined) {
    // Simple copy
    if (dst.size !== src.size) {
      throw new Error('dst and src must have same size');
    }
    for (let i = 0; i < dst.size; i++) {
      dst.setFlat(i, src.getFlat(i));
    }
  } else {
    // Conditional copy
    if (dst.size !== src.size || dst.size !== where.size) {
      throw new Error('dst, src, and where must have same size');
    }
    for (let i = 0; i < dst.size; i++) {
      if (where.getFlat(i) !== 0) {
        dst.setFlat(i, src.getFlat(i));
      }
    }
  }
}

/**
 * Convert input to array.
 */
export async function asarray(
  a: NDArray | number[] | number,
  dtype?: DType
): Promise<NDArray> {
  if (a instanceof NDArray) {
    if (dtype !== undefined && a.dtype !== dtype) {
      return a.astype(dtype);
    }
    return a;  // Return as-is (not a copy)
  }

  return NDArray.fromArray(a as number[], undefined, { dtype });
}
```

---

## File Changes Summary

### New Files to Create

```
src/wasm/
└── manipulation.h/c   # C implementations for concatenate, tile, repeat, flip, roll, split

src/ts/
└── manipulation.ts    # All manipulation functions
```

### Files to Modify

```
src/wasm/CMakeLists.txt or build script
├── Add manipulation.c to compilation

src/ts/types.ts
├── Add WasmModule function declarations:
│   _ndarray_concatenate
│   _ndarray_tile
│   _ndarray_repeat
│   _ndarray_flip
│   _ndarray_roll
│   _ndarray_split

src/ts/index.ts
├── Export all manipulation functions:
│   concatenate, stack, vstack, hstack, dstack, column_stack, block, append
│   split, array_split, vsplit, hsplit, dsplit, unstack
│   tile, repeat, pad
│   flip, fliplr, flipud, roll, rot90, resize, trim_zeros, insert, delete
│   copyto, asarray

scripts/build-wasm.sh
├── Add manipulation.c to compilation
└── Add new EXPORTED_FUNCTIONS
```

---

## Build Script Updates

Add to `scripts/build-wasm.sh` EXPORTED_FUNCTIONS:

```bash
# Manipulation
"_ndarray_concatenate",
"_ndarray_tile",
"_ndarray_repeat",
"_ndarray_flip",
"_ndarray_roll",
"_ndarray_split"
```

Add new source file:

```bash
emcc \
    "$SRC_DIR/ndarray.c" \
    "$SRC_DIR/dtype.c" \
    "$SRC_DIR/pairwise_sum.c" \
    "$SRC_DIR/broadcast.c" \
    "$SRC_DIR/indexing.c" \
    "$SRC_DIR/manipulation.c" \
    -o "$OUT_DIR/numjs.cjs" \
    ...
```

---

## WasmModule Interface Updates

Add to `src/ts/types.ts`:

```typescript
// Manipulation functions
_ndarray_concatenate(arraysPtr: number, nArrays: number, axis: number): number;
_ndarray_tile(ptr: number, repsPtr: number, ndim: number): number;
_ndarray_repeat(ptr: number, repeatsPtr: number, nReps: number, axis: number): number;
_ndarray_flip(ptr: number, axis: number): number;
_ndarray_roll(ptr: number, shift: number, axis: number): number;
_ndarray_split(ptr: number, indicesPtr: number, nIndices: number, axis: number,
               outArraysPtr: number, outCountPtr: number): number;
```

---

## Implementation Order

```
Week 1: Joining Functions
├── Day 1: C: ndarray_concatenate()
├── Day 2: TS: concatenate() wrapper + tests
├── Day 3: TS: stack(), vstack(), hstack()
├── Day 4: TS: dstack(), column_stack()
└── Day 5: TS: block(), append() + tests

Week 2: Splitting Functions
├── Day 1: C: ndarray_split() (optional, can be TS-only)
├── Day 2: TS: split(), array_split()
├── Day 3: TS: vsplit(), hsplit(), dsplit()
├── Day 4: TS: unstack()
└── Day 5: Integration tests

Week 3: Tiling & Repeating
├── Day 1: C: ndarray_tile(), ndarray_repeat()
├── Day 2: TS: tile() wrapper + tests
├── Day 3: TS: repeat() wrapper + tests
├── Day 4: TS: pad() (constant mode)
└── Day 5: TS: pad() (edge, reflect, symmetric, wrap modes)

Week 4: Rearranging & Polish
├── Day 1: TS: flip(), fliplr(), flipud()
├── Day 2: TS: roll(), rot90()
├── Day 3: TS: resize(), trim_zeros()
├── Day 4: TS: insert(), delete()
└── Day 5: copyto(), asarray(), comprehensive tests
```

---

## Verification Plan

After Phase 5 completion, verify:

```bash
# Build
npm run build

# Run all tests
npm test

# Specific Phase 5 tests should pass:

# Joining
✓ concatenate() joins arrays along axis
✓ concatenate() validates shape compatibility
✓ stack() inserts new axis and joins
✓ vstack() stacks vertically (handles 1D as row)
✓ hstack() stacks horizontally (1D along axis 0, 2D+ along axis 1)
✓ dstack() stacks along third axis
✓ column_stack() makes columns from 1D arrays
✓ block() assembles nested blocks
✓ append() adds values to end

# Splitting
✓ split() divides equally, raises on unequal
✓ array_split() allows unequal division
✓ vsplit(), hsplit(), dsplit() split along correct axes
✓ unstack() is inverse of stack()

# Tiling
✓ tile() repeats array according to reps
✓ repeat() repeats elements
✓ pad() adds padding with various modes

# Rearranging
✓ flip() reverses along axis
✓ fliplr(), flipud() flip horizontally/vertically
✓ roll() shifts elements cyclically
✓ rot90() rotates by 90 degrees
✓ resize() reshapes with repeat/truncate
✓ trim_zeros() removes leading/trailing zeros
✓ insert() adds elements
✓ delete() removes elements
```

Generate NumPy comparison vectors:

```python
# tests/python/generate_phase5_tests.py
import numpy as np
import json

tests = {
    "concatenate": [
        {"arrays": [[[1, 2], [3, 4]], [[5, 6]]], "axis": 0, "expected": [[1, 2], [3, 4], [5, 6]]},
        {"arrays": [[[1], [2]], [[3], [4]]], "axis": 1, "expected": [[1, 3], [2, 4]]},
    ],
    "stack": [
        {"arrays": [[1, 2], [3, 4]], "axis": 0, "expected_shape": [2, 2]},
        {"arrays": [[1, 2], [3, 4]], "axis": 1, "expected_shape": [2, 2]},
    ],
    "split": [
        {"shape": [6], "sections": 3, "axis": 0, "expected_shapes": [[2], [2], [2]]},
        {"shape": [6], "indices": [2, 4], "axis": 0, "expected_shapes": [[2], [2], [2]]},
    ],
    "tile": [
        {"shape": [2, 3], "reps": [2, 1], "expected_shape": [4, 3]},
        {"shape": [2], "reps": 3, "expected_shape": [6]},
    ],
    "flip": [
        {"data": [[1, 2], [3, 4]], "axis": 0, "expected": [[3, 4], [1, 2]]},
        {"data": [[1, 2], [3, 4]], "axis": 1, "expected": [[2, 1], [4, 3]]},
    ],
    "roll": [
        {"data": [1, 2, 3, 4, 5], "shift": 2, "expected": [4, 5, 1, 2, 3]},
        {"data": [1, 2, 3, 4, 5], "shift": -2, "expected": [3, 4, 5, 1, 2]},
    ],
}

with open("tests/fixtures/phase5_vectors.json", "w") as f:
    json.dump(tests, f, indent=2)
```

---

## Critical Dependencies for Phase 6

Phase 5 completion enables:

- **Phase 6.1 Sorting**: Uses manipulation for rearranging sorted results
- **Phase 6.2 Statistics**: Uses concatenate for axis-based operations
- **Phase 7 Logic**: Uses manipulation for array comparison results

Phase 5 should be complete before proceeding to Phase 6 (Sorting, Searching & Statistics).
