/**
 * NumJS Index Functions
 *
 * TypeScript wrappers for WASM index operations.
 * Adapted from NumPy's item_selection.c
 */

import { NDArray } from './NDArray.js';
import { DType, CLIP_RAISE, CLIP_WRAP, CLIP_CLIP } from './types.js';

export { CLIP_RAISE, CLIP_WRAP, CLIP_CLIP };

/**
 * Clip mode type
 */
export type ClipMode = typeof CLIP_RAISE | typeof CLIP_WRAP | typeof CLIP_CLIP;

/**
 * Take elements from an array along an axis.
 *
 * @param arr - Source array
 * @param indices - Array of indices to take
 * @param axis - Axis along which to take (default: 0)
 * @param mode - How to handle out-of-bounds indices (default: 'raise')
 * @returns New array with taken elements
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([[1, 2], [3, 4], [5, 6]]);
 * const taken = await take(arr, await NDArray.fromArray([0, 2]));
 * // [[1, 2], [5, 6]]
 * ```
 */
export async function take(
  arr: NDArray,
  indices: NDArray,
  axis: number = 0,
  mode: ClipMode = CLIP_RAISE
): Promise<NDArray> {
  const module = arr._wasmModule;

  const resultPtr = module._ndarray_take(
    arr._wasmPtr,
    indices._wasmPtr,
    axis,
    mode
  );

  if (resultPtr === 0) {
    throw new Error('take failed: index out of bounds or invalid axis');
  }

  return NDArray._fromPtr(resultPtr, module);
}

/**
 * Take elements from a flattened array.
 *
 * @param arr - Source array (will be flattened conceptually)
 * @param indices - Array of flat indices
 * @param mode - How to handle out-of-bounds indices (default: 'raise')
 * @returns New array with taken elements (same shape as indices)
 */
export async function takeFlat(
  arr: NDArray,
  indices: NDArray,
  mode: ClipMode = CLIP_RAISE
): Promise<NDArray> {
  const module = arr._wasmModule;

  const resultPtr = module._ndarray_take_flat(
    arr._wasmPtr,
    indices._wasmPtr,
    mode
  );

  if (resultPtr === 0) {
    throw new Error('takeFlat failed: index out of bounds');
  }

  return NDArray._fromPtr(resultPtr, module);
}

/**
 * Put values into an array at specified flat indices.
 *
 * @param arr - Target array (modified in place)
 * @param indices - Array of flat indices
 * @param values - Values to put (broadcast if smaller)
 * @param mode - How to handle out-of-bounds indices (default: 'raise')
 *
 * @example
 * ```typescript
 * const arr = await NDArray.zeros([3, 3]);
 * const indices = await NDArray.fromArray([0, 4, 8]);
 * const values = await NDArray.fromArray([1, 1, 1]);
 * await put(arr, indices, values);  // Sets diagonal to 1
 * ```
 */
export async function put(
  arr: NDArray,
  indices: NDArray,
  values: NDArray,
  mode: ClipMode = CLIP_RAISE
): Promise<void> {
  const module = arr._wasmModule;

  const result = module._ndarray_put(
    arr._wasmPtr,
    indices._wasmPtr,
    values._wasmPtr,
    mode
  );

  if (result !== 0) {
    throw new Error('put failed: index out of bounds or array not writeable');
  }
}

/**
 * Count the number of nonzero elements in an array.
 *
 * @param arr - Input array
 * @returns Number of nonzero elements
 */
export async function countNonzero(arr: NDArray): Promise<number> {
  const module = arr._wasmModule;
  return module._ndarray_count_nonzero(arr._wasmPtr);
}

/**
 * Find indices of nonzero elements.
 *
 * Returns a 2D array of shape (num_nonzero, ndim) where each row
 * contains the multi-dimensional index of a nonzero element.
 *
 * @param arr - Input array
 * @returns 2D index array
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([[1, 0], [0, 2]]);
 * const indices = await nonzero(arr);
 * // [[0, 0], [1, 1]] - indices of nonzero elements
 * ```
 */
export async function nonzero(arr: NDArray): Promise<NDArray> {
  const module = arr._wasmModule;

  const resultPtr = module._ndarray_nonzero(arr._wasmPtr);

  if (resultPtr === 0) {
    throw new Error('nonzero failed: cannot call on 0-d array');
  }

  return NDArray._fromPtr(resultPtr, module);
}

/**
 * Return flat indices of nonzero elements.
 *
 * @param arr - Input array
 * @returns 1D array of flat indices where arr != 0
 */
export async function flatnonzero(arr: NDArray): Promise<NDArray> {
  const module = arr._wasmModule;

  const resultPtr = module._ndarray_flatnonzero(arr._wasmPtr);

  if (resultPtr === 0) {
    throw new Error('flatnonzero failed');
  }

  return NDArray._fromPtr(resultPtr, module);
}

/**
 * Return elements chosen from x or y depending on condition.
 *
 * If x and y are provided, returns an array where elements are taken
 * from x where condition is true, and from y where condition is false.
 *
 * If only condition is provided, returns the indices where condition is true
 * (equivalent to nonzero).
 *
 * @param condition - Boolean condition array
 * @param x - Values where condition is true (optional)
 * @param y - Values where condition is false (optional)
 * @returns Selected values or indices
 *
 * @example
 * ```typescript
 * const cond = await NDArray.fromArray([true, false, true]);
 * const x = await NDArray.fromArray([1, 2, 3]);
 * const y = await NDArray.fromArray([10, 20, 30]);
 * const result = await where(cond, x, y);  // [1, 20, 3]
 *
 * // Without x and y, returns indices where true
 * const indices = await where(cond);  // [[0], [2]]
 * ```
 */
export async function where(
  condition: NDArray,
  x?: NDArray,
  y?: NDArray
): Promise<NDArray> {
  const module = condition._wasmModule;

  const xPtr = x ? x._wasmPtr : 0;
  const yPtr = y ? y._wasmPtr : 0;

  const resultPtr = module._ndarray_where(condition._wasmPtr, xPtr, yPtr);

  if (resultPtr === 0) {
    if (x !== undefined && y === undefined) {
      throw new Error('where requires both x and y, or neither');
    }
    throw new Error('where failed: shapes not broadcastable');
  }

  return NDArray._fromPtr(resultPtr, module);
}

/**
 * Return selected slices of an array along given axis.
 *
 * @param condition - Boolean 1D array selecting which slices to keep
 * @param arr - Input array
 * @param axis - Axis along which to select (default: 0)
 * @returns Compressed array
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([[1, 2], [3, 4], [5, 6]]);
 * const cond = await NDArray.fromArray([true, false, true]);
 * const result = await compress(cond, arr);  // [[1, 2], [5, 6]]
 * ```
 */
export async function compress(
  condition: NDArray,
  arr: NDArray,
  axis: number = 0
): Promise<NDArray> {
  const module = arr._wasmModule;

  const resultPtr = module._ndarray_compress(
    condition._wasmPtr,
    arr._wasmPtr,
    axis
  );

  if (resultPtr === 0) {
    throw new Error('compress failed: condition must be 1D or axis not supported');
  }

  return NDArray._fromPtr(resultPtr, module);
}

/**
 * Return elements of an array that satisfy some condition.
 *
 * @param condition - Boolean array, same shape as arr
 * @param arr - Input array
 * @returns 1D array of elements where condition is true
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([[1, 2], [3, 4]]);
 * const cond = await NDArray.fromArray([[true, false], [false, true]]);
 * const result = await extract(cond, arr);  // [1, 4]
 * ```
 */
export async function extract(
  condition: NDArray,
  arr: NDArray
): Promise<NDArray> {
  const module = arr._wasmModule;

  const resultPtr = module._ndarray_extract(condition._wasmPtr, arr._wasmPtr);

  if (resultPtr === 0) {
    throw new Error('extract failed: condition and arr must have same size');
  }

  return NDArray._fromPtr(resultPtr, module);
}

/**
 * Construct an array from an index array and a set of arrays to choose from.
 *
 * @param indices - Array of indices into choices
 * @param choices - Array of choice arrays
 * @param mode - How to handle out-of-bounds indices (default: 'raise')
 * @returns New array with chosen values
 *
 * @example
 * ```typescript
 * const indices = await NDArray.fromArray([0, 1, 0, 1]);
 * const choices = [
 *   await NDArray.fromArray([1, 1, 1, 1]),
 *   await NDArray.fromArray([2, 2, 2, 2])
 * ];
 * const result = await choose(indices, choices);  // [1, 2, 1, 2]
 * ```
 */
export async function choose(
  indices: NDArray,
  choices: NDArray[],
  mode: ClipMode = CLIP_RAISE
): Promise<NDArray> {
  if (choices.length === 0) {
    throw new Error('choices must not be empty');
  }

  const module = indices._wasmModule;

  // Allocate array of pointers to choice arrays
  const choicesPtr = module._malloc(choices.length * 4);
  for (let i = 0; i < choices.length; i++) {
    module.setValue(choicesPtr + i * 4, choices[i]._wasmPtr, 'i32');
  }

  const resultPtr = module._ndarray_choose(
    indices._wasmPtr,
    choicesPtr,
    choices.length,
    mode
  );

  module._free(choicesPtr);

  if (resultPtr === 0) {
    throw new Error('choose failed: index out of bounds or incompatible shapes');
  }

  return NDArray._fromPtr(resultPtr, module);
}

/**
 * Return specified diagonals from a 2D array.
 *
 * @param arr - Input array (must be 2D)
 * @param offset - Offset from main diagonal (default: 0)
 * @returns 1D array of diagonal elements
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([[1, 2, 3], [4, 5, 6], [7, 8, 9]]);
 * const main = await diagonal(arr);      // [1, 5, 9]
 * const upper = await diagonal(arr, 1);  // [2, 6]
 * const lower = await diagonal(arr, -1); // [4, 8]
 * ```
 */
export async function diagonal(
  arr: NDArray,
  offset: number = 0,
  axis1: number = 0,
  axis2: number = 1
): Promise<NDArray> {
  const module = arr._wasmModule;

  const resultPtr = module._ndarray_diagonal(
    arr._wasmPtr,
    offset,
    axis1,
    axis2
  );

  if (resultPtr === 0) {
    throw new Error('diagonal failed: array must be at least 2D');
  }

  return NDArray._fromPtr(resultPtr, module);
}

/* ============ Index Generation Functions ============ */
/* Note: argsort, argmax, argmin moved to statistics.ts (WASM-accelerated) */

/**
 * Convert a flat index into a tuple of coordinate arrays.
 *
 * @param indices - Array of flat indices
 * @param shape - Shape of the array
 * @returns Array of coordinate arrays, one per dimension
 */
export async function unravelIndex(
  indices: number[] | NDArray,
  shape: number[]
): Promise<number[][]> {
  const flatIndices =
    indices instanceof NDArray ? (indices.toArray() as number[]) : indices;

  const result: number[][] = shape.map(() => []);

  for (const idx of flatIndices) {
    let remainder = idx;
    for (let d = shape.length - 1; d >= 0; d--) {
      result[d].push(remainder % shape[d]);
      remainder = Math.floor(remainder / shape[d]);
    }
  }

  return result;
}

/**
 * Convert a tuple of coordinate arrays into flat indices.
 *
 * @param coords - Array of coordinate arrays, one per dimension
 * @param shape - Shape of the array
 * @returns Array of flat indices
 */
export async function ravelMultiIndex(
  coords: number[][],
  shape: number[]
): Promise<number[]> {
  if (coords.length !== shape.length) {
    throw new Error('Number of coordinate arrays must match number of dimensions');
  }

  const n = coords[0].length;
  const result: number[] = [];

  for (let i = 0; i < n; i++) {
    let flatIdx = 0;
    let multiplier = 1;
    for (let d = shape.length - 1; d >= 0; d--) {
      flatIdx += coords[d][i] * multiplier;
      multiplier *= shape[d];
    }
    result.push(flatIdx);
  }

  return result;
}

/**
 * Return coordinate matrices from coordinate vectors.
 *
 * @param xi - 1D coordinate arrays
 * @param indexing - Cartesian ('xy') or matrix ('ij') indexing (default: 'xy')
 * @returns Array of coordinate matrices
 */
export async function meshgrid(
  ...xi: (number[] | NDArray)[]
): Promise<NDArray[]> {
  const arrays = await Promise.all(
    xi.map(async (x) =>
      x instanceof NDArray ? (x.toArray() as number[]) : x
    )
  );

  const ndim = arrays.length;
  if (ndim === 0) return [];

  // Compute output shape
  const shape = arrays.map((arr) => arr.length);

  const results: NDArray[] = [];

  for (let dim = 0; dim < ndim; dim++) {
    // Create shape for this result (all 1s except for dim)
    const resultShape = shape.slice();

    // Build the output data
    const totalSize = resultShape.reduce((a, b) => a * b, 1);
    const data: number[] = new Array(totalSize);

    // Calculate strides
    const strides: number[] = new Array(ndim);
    strides[ndim - 1] = 1;
    for (let i = ndim - 2; i >= 0; i--) {
      strides[i] = strides[i + 1] * shape[i + 1];
    }

    // Fill in the data
    for (let i = 0; i < totalSize; i++) {
      // Convert flat index to multi-index
      let remainder = i;
      const multiIdx: number[] = new Array(ndim);
      for (let d = 0; d < ndim; d++) {
        multiIdx[d] = Math.floor(remainder / strides[d]);
        remainder = remainder % strides[d];
      }

      // For this dimension, use the coordinate from arrays[dim]
      data[i] = arrays[dim][multiIdx[dim]];
    }

    results.push(await NDArray.fromArray(data, resultShape));
  }

  return results;
}

/* ============ Shape Manipulation Functions ============ */

/**
 * Convert inputs to arrays with at least one dimension.
 * Scalar inputs are converted to 1D arrays. Higher-dimensional inputs are unchanged.
 *
 * @param arrs - Input arrays
 * @returns Single array if one input, tuple of arrays if multiple inputs
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray(5);  // 0-d scalar
 * const result = await atleast_1d(arr);    // shape: [1]
 * ```
 */
export function atleast_1d(...arrs: NDArray[]): NDArray | NDArray[] {
  const results = arrs.map(arr => {
    if (arr.ndim === 0) {
      return arr.reshape([1]);
    }
    return arr;  // Already 1D+, return as-is (view)
  });
  return arrs.length === 1 ? results[0] : results;
}

/**
 * View inputs as arrays with at least two dimensions.
 *
 * @param arrs - Input arrays
 * @returns Single array if one input, tuple of arrays if multiple inputs
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2, 3]);  // 1D
 * const result = await atleast_2d(arr);            // shape: [1, 3]
 * ```
 */
export function atleast_2d(...arrs: NDArray[]): NDArray | NDArray[] {
  const results = arrs.map(arr => {
    if (arr.ndim === 0) {
      return arr.reshape([1, 1]);
    } else if (arr.ndim === 1) {
      return arr.expandDims(0);  // [N] -> [1, N]
    }
    return arr;  // Already 2D+, return as-is
  });
  return arrs.length === 1 ? results[0] : results;
}

/**
 * View inputs as arrays with at least three dimensions.
 *
 * @param arrs - Input arrays
 * @returns Single array if one input, tuple of arrays if multiple inputs
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([[1, 2], [3, 4]]);  // 2D
 * const result = await atleast_3d(arr);                   // shape: [2, 2, 1]
 * ```
 */
export function atleast_3d(...arrs: NDArray[]): NDArray | NDArray[] {
  const results = arrs.map(arr => {
    if (arr.ndim === 0) {
      return arr.reshape([1, 1, 1]);
    } else if (arr.ndim === 1) {
      // [N] -> [1, N, 1]
      return arr.expandDims(0).expandDims(2);
    } else if (arr.ndim === 2) {
      // [M, N] -> [M, N, 1]
      return arr.expandDims(2);
    }
    return arr;  // Already 3D+, return as-is
  });
  return arrs.length === 1 ? results[0] : results;
}

/* ============ Index Generation Functions ============ */

/**
 * Find the indices of array elements that are non-zero, grouped by element.
 *
 * @param arr - Input array
 * @returns 2D array of shape (N, arr.ndim) where N is number of nonzero elements
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([[1, 0], [0, 2]]);
 * const result = await argwhere(arr);
 * // [[0, 0], [1, 1]] - indices of nonzero elements
 * ```
 */
export async function argwhere(arr: NDArray): Promise<NDArray> {
  // For 0-d arrays, return empty with 0 columns
  if (arr.ndim === 0) {
    if (arr.item() !== 0) {
      // Single nonzero element at "index" []
      return NDArray.empty([1, 0], { dtype: DType.Int32 });
    } else {
      return NDArray.empty([0, 0], { dtype: DType.Int32 });
    }
  }

  // Our nonzero() returns shape (num_nonzero, ndim), which is argwhere format
  return nonzero(arr);
}

/**
 * Return an array representing the indices of a grid.
 *
 * @param dimensions - The shape of the grid
 * @param dtype - Data type of the result (default: Int32)
 * @param sparse - Return sparse representation (default: false)
 * @returns If sparse is false, array of shape (N, d1, d2, ..., dN).
 *          If sparse is true, tuple of N arrays.
 *
 * @example
 * ```typescript
 * // Dense: returns array of shape (2, 3, 4)
 * const grid = await indices([3, 4]);
 * // grid[0] contains row indices, grid[1] contains column indices
 *
 * // Sparse: returns tuple of arrays
 * const [rows, cols] = await indices([3, 4], DType.Int32, true);
 * // rows.shape = [3, 1], cols.shape = [1, 4]
 * ```
 */
export async function indices(
  dimensions: number[],
  dtype: DType = DType.Int32,
  sparse: boolean = false
): Promise<NDArray | NDArray[]> {
  const N = dimensions.length;

  if (N === 0) {
    if (sparse) {
      return [];
    }
    return NDArray.empty([0], { dtype });
  }

  if (sparse) {
    // Return tuple of N arrays, each with shape (1,...,dim_i,...,1)
    const result: NDArray[] = [];
    for (let i = 0; i < N; i++) {
      const shape = new Array(N).fill(1);
      shape[i] = dimensions[i];
      const idx = await NDArray.arange(0, dimensions[i]);
      const converted = idx.astype(dtype);
      idx.dispose();
      result.push(converted.reshape(shape));
    }
    return result;
  } else {
    // Return single array of shape (N, ...dimensions)
    const totalShape = [N, ...dimensions];
    const result = await NDArray.empty(totalShape, { dtype });

    // Calculate strides for the dimensions part
    const dimStrides: number[] = new Array(N);
    dimStrides[N - 1] = 1;
    for (let i = N - 2; i >= 0; i--) {
      dimStrides[i] = dimStrides[i + 1] * dimensions[i + 1];
    }

    const subSize = dimensions.reduce((a, b) => a * b, 1);

    // Fill each dimension slice
    for (let dim = 0; dim < N; dim++) {
      for (let flat = 0; flat < subSize; flat++) {
        // Convert flat index to multi-index
        let remainder = flat;
        const multiIdx: number[] = new Array(N);
        for (let d = 0; d < N; d++) {
          multiIdx[d] = Math.floor(remainder / dimStrides[d]);
          remainder = remainder % dimStrides[d];
        }
        // Set the index value for this dimension
        result.set(multiIdx[dim], dim, ...multiIdx);
      }
    }

    return result;
  }
}

/**
 * Construct an open mesh from multiple sequences.
 * This function takes N 1-D sequences and returns N outputs with N dimensions each,
 * such that the shape in the ith dimension is the length of the ith input.
 *
 * @param args - 1-D sequences (arrays or NDArrays)
 * @returns Tuple of N-D arrays
 *
 * @example
 * ```typescript
 * const [rows, cols] = await ix_([0, 1], [2, 3, 4]);
 * // rows.shape = [2, 1], rows = [[0], [1]]
 * // cols.shape = [1, 3], cols = [[2, 3, 4]]
 * ```
 */
export async function ix_(
  ...args: (number[] | NDArray)[]
): Promise<NDArray[]> {
  const nd = args.length;
  if (nd === 0) return [];

  const result: NDArray[] = [];

  for (let k = 0; k < nd; k++) {
    let arr: NDArray;
    if (args[k] instanceof NDArray) {
      arr = args[k] as NDArray;
    } else {
      arr = await NDArray.fromArray(args[k] as number[]);
    }

    if (arr.ndim !== 1) {
      throw new Error('Cross index must be 1 dimensional');
    }

    // Handle boolean arrays (convert via nonzero)
    if (arr.dtype === DType.Bool) {
      const nz = await flatnonzero(arr);
      if (!(args[k] instanceof NDArray)) {
        arr.dispose();
      }
      arr = nz;
    }

    // Reshape: (1,)*k + (size,) + (1,)*(nd-k-1)
    const shape = new Array(nd).fill(1);
    shape[k] = arr.size;

    // If we created the array, we can reshape in place
    if (args[k] instanceof NDArray) {
      result.push(arr.reshape(shape));
    } else {
      const reshaped = arr.reshape(shape);
      result.push(reshaped);
    }
  }

  return result;
}

/**
 * Return the indices to access the main diagonal of an array.
 *
 * @param n - Size of the array (along each dimension)
 * @param ndim - Number of dimensions (default: 2)
 * @returns Tuple of ndim index arrays
 *
 * @example
 * ```typescript
 * const [i, j] = await diag_indices(3);
 * // i = [0, 1, 2], j = [0, 1, 2]
 * // Use to access arr[i, j] which gives diagonal elements
 * ```
 */
export async function diag_indices(
  n: number,
  ndim: number = 2
): Promise<NDArray[]> {
  const idx = await NDArray.arange(0, n);
  const converted = idx.astype(DType.Int32);
  idx.dispose();

  const result: NDArray[] = [];
  for (let i = 0; i < ndim; i++) {
    result.push(i === 0 ? converted : converted.copy());
  }
  return result;
}

/**
 * Return the indices for the lower-triangle of an (n, m) array.
 *
 * @param n - Number of rows
 * @param k - Diagonal offset (default: 0, main diagonal)
 * @param m - Number of columns (default: n)
 * @returns Tuple of two arrays (row indices, column indices)
 *
 * @example
 * ```typescript
 * const [rows, cols] = await tril_indices(3);
 * // rows = [0, 1, 1, 2, 2, 2]
 * // cols = [0, 0, 1, 0, 1, 2]
 * ```
 */
export async function tril_indices(
  n: number,
  k: number = 0,
  m?: number
): Promise<[NDArray, NDArray]> {
  m = m ?? n;

  const rows: number[] = [];
  const cols: number[] = [];

  for (let i = 0; i < n; i++) {
    for (let j = 0; j < m; j++) {
      // Lower triangle: j <= i + k
      if (j <= i + k) {
        rows.push(i);
        cols.push(j);
      }
    }
  }

  return [
    await NDArray.fromArray(rows, [rows.length], { dtype: DType.Int32 }),
    await NDArray.fromArray(cols, [cols.length], { dtype: DType.Int32 })
  ];
}

/**
 * Return the indices for the upper-triangle of an (n, m) array.
 *
 * @param n - Number of rows
 * @param k - Diagonal offset (default: 0, main diagonal)
 * @param m - Number of columns (default: n)
 * @returns Tuple of two arrays (row indices, column indices)
 *
 * @example
 * ```typescript
 * const [rows, cols] = await triu_indices(3);
 * // rows = [0, 0, 0, 1, 1, 2]
 * // cols = [0, 1, 2, 1, 2, 2]
 * ```
 */
export async function triu_indices(
  n: number,
  k: number = 0,
  m?: number
): Promise<[NDArray, NDArray]> {
  m = m ?? n;

  const rows: number[] = [];
  const cols: number[] = [];

  for (let i = 0; i < n; i++) {
    for (let j = 0; j < m; j++) {
      // Upper triangle: j >= i + k
      if (j >= i + k) {
        rows.push(i);
        cols.push(j);
      }
    }
  }

  return [
    await NDArray.fromArray(rows, [rows.length], { dtype: DType.Int32 }),
    await NDArray.fromArray(cols, [cols.length], { dtype: DType.Int32 })
  ];
}

/* ============ Advanced Indexing Functions ============ */

/**
 * Helper to convert flat index to multi-dimensional index.
 */
function unravelIndexSingle(flatIndex: number, shape: number[]): number[] {
  const result: number[] = new Array(shape.length);
  let remainder = flatIndex;

  for (let i = shape.length - 1; i >= 0; i--) {
    result[i] = remainder % shape[i];
    remainder = Math.floor(remainder / shape[i]);
  }

  return result;
}

/**
 * Take values from the input array by matching 1d index and data slices.
 *
 * @param arr - Source array
 * @param indices - Indices to take along axis
 * @param axis - The axis to take along
 * @returns Array with values taken along the specified axis
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([[10, 20, 30], [40, 50, 60]]);
 * const idx = await NDArray.fromArray([[0, 2], [1, 0]], { dtype: DType.Int32 });
 * const result = await take_along_axis(arr, idx, 1);
 * // [[10, 30], [50, 40]]
 * ```
 */
export async function take_along_axis(
  arr: NDArray,
  indicesArr: NDArray,
  axis: number
): Promise<NDArray> {
  // Normalize axis
  const ndim = arr.ndim;
  axis = axis < 0 ? axis + ndim : axis;

  if (axis < 0 || axis >= ndim) {
    throw new Error(`axis ${axis} is out of bounds for array with ${ndim} dimensions`);
  }

  // Validate that indices has same ndim
  if (indicesArr.ndim !== ndim) {
    throw new Error('indices must have same number of dimensions as arr');
  }

  // Create result array with same shape as indices
  const result = await NDArray.empty(indicesArr.shape, { dtype: arr.dtype });

  // Iterate over all positions in indices
  const arrShape = arr.shape;
  for (let i = 0; i < indicesArr.size; i++) {
    // Get multi-index in indices array
    const multiIdx = unravelIndexSingle(i, indicesArr.shape);

    // Get the index value at this position
    let idxVal = Math.floor(indicesArr.getFlat(i));

    // Handle negative indices
    if (idxVal < 0) {
      idxVal += arrShape[axis];
    }

    // Build source index: same as multiIdx but with idxVal at axis
    const srcIdx = [...multiIdx];
    srcIdx[axis] = idxVal;

    // Copy value
    result.setFlat(i, arr.get(...srcIdx));
  }

  return result;
}

/**
 * Put values into the destination array by matching 1d index and data slices.
 *
 * @param arr - Destination array (modified in place)
 * @param indices - Indices to put along axis
 * @param values - Values to put
 * @param axis - The axis to put along
 *
 * @example
 * ```typescript
 * const arr = await NDArray.zeros([2, 3]);
 * const idx = await NDArray.fromArray([[0, 2], [1, 0]], { dtype: DType.Int32 });
 * const vals = await NDArray.fromArray([[1, 2], [3, 4]]);
 * await put_along_axis(arr, idx, vals, 1);
 * // arr = [[1, 0, 2], [0, 3, 4]]
 * ```
 */
export async function put_along_axis(
  arr: NDArray,
  indicesArr: NDArray,
  values: NDArray,
  axis: number
): Promise<void> {
  // Normalize axis
  const ndim = arr.ndim;
  axis = axis < 0 ? axis + ndim : axis;

  if (axis < 0 || axis >= ndim) {
    throw new Error(`axis ${axis} is out of bounds for array with ${ndim} dimensions`);
  }

  // Validate dimensions match
  if (indicesArr.ndim !== ndim) {
    throw new Error('indices must have same number of dimensions as arr');
  }
  if (values.ndim !== ndim) {
    throw new Error('values must have same number of dimensions as arr');
  }

  // Iterate over all positions in indices
  const arrShape = arr.shape;
  for (let i = 0; i < indicesArr.size; i++) {
    // Get multi-index in indices array
    const multiIdx = unravelIndexSingle(i, indicesArr.shape);

    // Get the index value at this position
    let idxVal = Math.floor(indicesArr.getFlat(i));

    // Handle negative indices
    if (idxVal < 0) {
      idxVal += arrShape[axis];
    }

    // Build destination index: same as multiIdx but with idxVal at axis
    const dstIdx = [...multiIdx];
    dstIdx[axis] = idxVal;

    // Get value (broadcast if values is smaller)
    const val = values.getFlat(i % values.size);

    // Set value
    arr.set(val, ...dstIdx);
  }
}

/**
 * Change elements of an array based on conditional and input values.
 * Sets a.flat[n] = values[n % len(values)] for each n where mask.flat[n]==True.
 *
 * @param arr - Array to modify (in place)
 * @param mask - Boolean mask array
 * @param values - Values to put (cycled if smaller than number of True in mask)
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2, 3, 4]);
 * const mask = await NDArray.fromArray([1, 0, 1, 0], { dtype: DType.Bool });
 * const vals = await NDArray.fromArray([10, 20]);
 * await putmask(arr, mask, vals);
 * // arr = [10, 2, 20, 4]
 * ```
 */
export async function putmask(
  arr: NDArray,
  mask: NDArray,
  values: NDArray
): Promise<void> {
  const n = arr.size;
  const valSize = values.size;

  if (mask.size !== n) {
    throw new Error('mask and arr must have the same size');
  }

  let valIdx = 0;
  for (let i = 0; i < n; i++) {
    if (mask.getFlat(i) !== 0) {
      arr.setFlat(i, values.getFlat(valIdx % valSize));
      valIdx++;
    }
  }
}

/**
 * Change elements of an array based on conditional and input values.
 * Similar to putmask, but places values sequentially.
 *
 * @param arr - Array to modify (in place)
 * @param mask - Boolean mask array
 * @param vals - Values to place (used sequentially, cycled if needed)
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2, 3, 4]);
 * const mask = await NDArray.fromArray([1, 0, 1, 0], { dtype: DType.Bool });
 * const vals = await NDArray.fromArray([10, 20]);
 * await place(arr, mask, vals);
 * // arr = [10, 2, 20, 4]
 * ```
 */
export async function place(
  arr: NDArray,
  mask: NDArray,
  vals: NDArray
): Promise<void> {
  // For our implementation, place and putmask behave the same way
  // In NumPy, the difference is subtle and relates to how values are cycled
  return putmask(arr, mask, vals);
}

/**
 * Return an array drawn from elements in choicelist, depending on conditions.
 *
 * @param condlist - List of boolean arrays (conditions)
 * @param choicelist - List of arrays from which output elements are taken
 * @param defaultValue - Value used when all conditions are False (default: 0)
 * @returns Array with elements from choicelist based on conditions
 *
 * @example
 * ```typescript
 * const x = await NDArray.fromArray([0, 1, 2, 3, 4]);
 * const cond1 = await NDArray.fromArray([0, 0, 1, 1, 0], { dtype: DType.Bool });
 * const cond2 = await NDArray.fromArray([1, 0, 0, 0, 1], { dtype: DType.Bool });
 * const choice1 = await NDArray.fromArray([10, 10, 10, 10, 10]);
 * const choice2 = await NDArray.fromArray([20, 20, 20, 20, 20]);
 * const result = await select([cond1, cond2], [choice1, choice2], 0);
 * // [20, 0, 10, 10, 20]
 * ```
 */
export async function select(
  condlist: NDArray[],
  choicelist: NDArray[],
  defaultValue: number = 0
): Promise<NDArray> {
  if (condlist.length !== choicelist.length) {
    throw new Error('condlist and choicelist must have same length');
  }

  if (condlist.length === 0) {
    throw new Error('select requires at least one condition');
  }

  // Get the common shape from the first array
  const shape = condlist[0].shape;
  const size = condlist[0].size;

  // Verify all arrays have compatible sizes
  for (const cond of condlist) {
    if (cond.size !== size) {
      throw new Error('all conditions must have the same shape');
    }
  }
  for (const choice of choicelist) {
    if (choice.size !== size) {
      throw new Error('all choices must have the same shape as conditions');
    }
  }

  // Initialize with default
  const result = await NDArray.full(shape, defaultValue);

  // Apply in reverse order (first condition takes precedence)
  for (let i = condlist.length - 1; i >= 0; i--) {
    const cond = condlist[i];
    const choice = choicelist[i];
    for (let j = 0; j < size; j++) {
      if (cond.getFlat(j) !== 0) {
        result.setFlat(j, choice.getFlat(j));
      }
    }
  }

  return result;
}
