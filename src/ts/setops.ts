/**
 * NumJS Set Operations
 *
 * TypeScript wrappers for WASM set operations.
 * Adapted from NumPy's _arraysetops_impl.py
 */

import { NDArray } from './NDArray.js';

/* ============ Isin Algorithm Kind ============ */
export const ISIN_AUTO = 0;
export const ISIN_SORT = 1;
export const ISIN_TABLE = 2;

export type IsinKind = typeof ISIN_AUTO | typeof ISIN_SORT | typeof ISIN_TABLE;

/* ============ UniqueResult Interface ============ */

/**
 * Result from unique() when returning additional arrays.
 */
export interface UniqueResult {
  /** Sorted unique values */
  values: NDArray;
  /** Indices of first occurrence in original array (if return_index) */
  indices?: NDArray;
  /** Indices to reconstruct original from unique (if return_inverse) */
  inverse?: NDArray;
  /** Count of each unique element (if return_counts) */
  counts?: NDArray;
}

/* ============ Unique Options ============ */

export interface UniqueOptions {
  /** Return indices of first occurrences */
  returnIndex?: boolean;
  /** Return indices to reconstruct original */
  returnInverse?: boolean;
  /** Return counts of each unique element */
  returnCounts?: boolean;
  /** Treat NaN values as equal (default: true) */
  equalNan?: boolean;
}

export interface IsinOptions {
  /** Invert the result (return elements NOT in ar2) */
  invert?: boolean;
  /** Assume inputs are already unique */
  assumeUnique?: boolean;
  /** Algorithm hint */
  kind?: IsinKind;
}

export interface IntersectOptions {
  /** Assume inputs are already unique */
  assumeUnique?: boolean;
  /** Return indices in original arrays */
  returnIndices?: boolean;
}

export interface SetDiffOptions {
  /** Assume inputs are already unique */
  assumeUnique?: boolean;
}

export interface SetXorOptions {
  /** Assume inputs are already unique */
  assumeUnique?: boolean;
}

/* ============ Unique Functions ============ */

/**
 * Find the unique elements of an array.
 *
 * @param arr - Input array (will be flattened)
 * @param options - Configuration options
 * @returns Unique values, or UniqueResult if any return_* option is true
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2, 2, 3, 1]);
 * const result = await unique(arr);
 * // NDArray [1, 2, 3]
 *
 * const full = await unique(arr, { returnIndex: true, returnCounts: true });
 * // full.values: [1, 2, 3]
 * // full.indices: [0, 1, 3]
 * // full.counts: [2, 2, 1]
 * ```
 */
export async function unique(
  arr: NDArray,
  options: UniqueOptions = {}
): Promise<NDArray | UniqueResult> {
  const module = arr._wasmModule;
  const {
    returnIndex = false,
    returnInverse = false,
    returnCounts = false,
    equalNan = true,
  } = options;

  const resultPtr = module._ndarray_unique(
    arr._wasmPtr,
    returnIndex,
    returnInverse,
    returnCounts,
    equalNan
  );

  if (resultPtr === 0) {
    throw new Error('unique failed');
  }

  // Read UniqueResult struct from WASM memory
  // UniqueResult { values, indices, inverse, counts } - 4 pointers (4 bytes each on wasm32)
  const valuesPtr = module.getValue(resultPtr, 'i32');
  const indicesPtr = module.getValue(resultPtr + 4, 'i32');
  const inversePtr = module.getValue(resultPtr + 8, 'i32');
  const countsPtr = module.getValue(resultPtr + 12, 'i32');

  const values = NDArray._fromPtr(valuesPtr, module);

  // If no additional returns requested, free struct and return just values
  if (!returnIndex && !returnInverse && !returnCounts) {
    // Clear pointers in struct to prevent double-free, then free struct
    module.setValue(resultPtr, 0, 'i32');
    module._unique_result_free(resultPtr);
    return values;
  }

  const result: UniqueResult = { values };

  if (returnIndex && indicesPtr !== 0) {
    result.indices = NDArray._fromPtr(indicesPtr, module);
    module.setValue(resultPtr + 4, 0, 'i32'); // Prevent double-free
  }

  if (returnInverse && inversePtr !== 0) {
    result.inverse = NDArray._fromPtr(inversePtr, module);
    module.setValue(resultPtr + 8, 0, 'i32');
  }

  if (returnCounts && countsPtr !== 0) {
    result.counts = NDArray._fromPtr(countsPtr, module);
    module.setValue(resultPtr + 12, 0, 'i32');
  }

  // Clear values pointer and free struct
  module.setValue(resultPtr, 0, 'i32');
  module._unique_result_free(resultPtr);

  return result;
}

/**
 * Return only the unique values (simplified interface).
 *
 * @param arr - Input array
 * @returns Sorted array of unique values
 */
export async function uniqueValues(arr: NDArray): Promise<NDArray> {
  const module = arr._wasmModule;
  const resultPtr = module._ndarray_unique_values(arr._wasmPtr);

  if (resultPtr === 0) {
    throw new Error('uniqueValues failed');
  }

  return NDArray._fromPtr(resultPtr, module);
}

/**
 * Return unique values and indices of first occurrence.
 *
 * @param arr - Input array
 * @returns Object with values and indices arrays
 */
export async function uniqueIndex(
  arr: NDArray
): Promise<{ values: NDArray; indices: NDArray }> {
  const result = (await unique(arr, { returnIndex: true })) as UniqueResult;
  return { values: result.values, indices: result.indices! };
}

/**
 * Return unique values and indices to reconstruct original.
 *
 * @param arr - Input array
 * @returns Object with values and inverse arrays
 */
export async function uniqueInverse(
  arr: NDArray
): Promise<{ values: NDArray; inverse: NDArray }> {
  const result = (await unique(arr, { returnInverse: true })) as UniqueResult;
  return { values: result.values, inverse: result.inverse! };
}

/**
 * Return unique values and their counts.
 *
 * @param arr - Input array
 * @returns Object with values and counts arrays
 */
export async function uniqueCounts(
  arr: NDArray
): Promise<{ values: NDArray; counts: NDArray }> {
  const result = (await unique(arr, { returnCounts: true })) as UniqueResult;
  return { values: result.values, counts: result.counts! };
}

/**
 * Return unique values, indices, inverse, and counts.
 *
 * @param arr - Input array
 * @returns Full UniqueResult with all arrays
 */
export async function uniqueAll(arr: NDArray): Promise<UniqueResult> {
  return (await unique(arr, {
    returnIndex: true,
    returnInverse: true,
    returnCounts: true,
  })) as UniqueResult;
}

/* ============ Set Combination Operations ============ */

/**
 * Find the union of two arrays.
 *
 * Returns the sorted unique values that are in either of the input arrays.
 *
 * @param ar1 - First input array
 * @param ar2 - Second input array
 * @returns Sorted array of unique union values
 *
 * @example
 * ```typescript
 * const a = await NDArray.fromArray([1, 2, 3]);
 * const b = await NDArray.fromArray([2, 3, 4]);
 * const result = await union1d(a, b);
 * // NDArray [1, 2, 3, 4]
 * ```
 */
export async function union1d(ar1: NDArray, ar2: NDArray): Promise<NDArray> {
  const module = ar1._wasmModule;
  const resultPtr = module._ndarray_union1d(ar1._wasmPtr, ar2._wasmPtr);

  if (resultPtr === 0) {
    throw new Error('union1d failed');
  }

  return NDArray._fromPtr(resultPtr, module);
}

/**
 * Find the intersection of two arrays.
 *
 * Returns the sorted unique values that are in both input arrays.
 *
 * @param ar1 - First input array
 * @param ar2 - Second input array
 * @param options - Configuration options
 * @returns Sorted array of common values, or object with indices if returnIndices
 *
 * @example
 * ```typescript
 * const a = await NDArray.fromArray([1, 3, 4, 3]);
 * const b = await NDArray.fromArray([3, 1, 2, 1]);
 * const result = await intersect1d(a, b);
 * // NDArray [1, 3]
 * ```
 */
export async function intersect1d(
  ar1: NDArray,
  ar2: NDArray,
  options: IntersectOptions = {}
): Promise<NDArray | { values: NDArray; indices1: NDArray; indices2: NDArray }> {
  const module = ar1._wasmModule;
  const { assumeUnique = false, returnIndices = false } = options;

  if (returnIndices) {
    // Allocate pointers for output indices
    const idx1PtrPtr = module._malloc(4);
    const idx2PtrPtr = module._malloc(4);

    const resultPtr = module._ndarray_intersect1d(
      ar1._wasmPtr,
      ar2._wasmPtr,
      assumeUnique,
      true,
      idx1PtrPtr,
      idx2PtrPtr
    );

    if (resultPtr === 0) {
      module._free(idx1PtrPtr);
      module._free(idx2PtrPtr);
      throw new Error('intersect1d failed');
    }

    const idx1Ptr = module.getValue(idx1PtrPtr, 'i32');
    const idx2Ptr = module.getValue(idx2PtrPtr, 'i32');

    module._free(idx1PtrPtr);
    module._free(idx2PtrPtr);

    return {
      values: NDArray._fromPtr(resultPtr, module),
      indices1: NDArray._fromPtr(idx1Ptr, module),
      indices2: NDArray._fromPtr(idx2Ptr, module),
    };
  }

  const resultPtr = module._ndarray_intersect1d(
    ar1._wasmPtr,
    ar2._wasmPtr,
    assumeUnique,
    false,
    0,
    0
  );

  if (resultPtr === 0) {
    throw new Error('intersect1d failed');
  }

  return NDArray._fromPtr(resultPtr, module);
}

/**
 * Find the set exclusive-or of two arrays.
 *
 * Returns the sorted unique values in exactly one (not both) of the input arrays.
 *
 * @param ar1 - First input array
 * @param ar2 - Second input array
 * @param options - Configuration options
 * @returns Sorted array of symmetric difference values
 *
 * @example
 * ```typescript
 * const a = await NDArray.fromArray([1, 2, 3]);
 * const b = await NDArray.fromArray([2, 3, 4]);
 * const result = await setxor1d(a, b);
 * // NDArray [1, 4]
 * ```
 */
export async function setxor1d(
  ar1: NDArray,
  ar2: NDArray,
  options: SetXorOptions = {}
): Promise<NDArray> {
  const module = ar1._wasmModule;
  const { assumeUnique = false } = options;

  const resultPtr = module._ndarray_setxor1d(
    ar1._wasmPtr,
    ar2._wasmPtr,
    assumeUnique
  );

  if (resultPtr === 0) {
    throw new Error('setxor1d failed');
  }

  return NDArray._fromPtr(resultPtr, module);
}

/**
 * Find the set difference of two arrays.
 *
 * Returns the sorted unique values in ar1 that are not in ar2.
 *
 * @param ar1 - First input array
 * @param ar2 - Second input array
 * @param options - Configuration options
 * @returns Sorted array of difference values
 *
 * @example
 * ```typescript
 * const a = await NDArray.fromArray([1, 2, 3, 4]);
 * const b = await NDArray.fromArray([2, 4]);
 * const result = await setdiff1d(a, b);
 * // NDArray [1, 3]
 * ```
 */
export async function setdiff1d(
  ar1: NDArray,
  ar2: NDArray,
  options: SetDiffOptions = {}
): Promise<NDArray> {
  const module = ar1._wasmModule;
  const { assumeUnique = false } = options;

  const resultPtr = module._ndarray_setdiff1d(
    ar1._wasmPtr,
    ar2._wasmPtr,
    assumeUnique
  );

  if (resultPtr === 0) {
    throw new Error('setdiff1d failed');
  }

  return NDArray._fromPtr(resultPtr, module);
}

/* ============ Membership Testing ============ */

/**
 * Test whether each element of ar1 is also present in ar2.
 *
 * Returns a boolean array of the same shape as ar1.
 *
 * @param ar1 - Input array
 * @param ar2 - Values to test against
 * @param options - Configuration options
 * @returns Boolean array of membership results
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([[0, 2], [4, 6]]);
 * const test = await NDArray.fromArray([1, 2, 4, 8]);
 * const result = await isin(arr, test);
 * // NDArray [[false, true], [true, false]]
 * ```
 */
export async function isin(
  ar1: NDArray,
  ar2: NDArray,
  options: IsinOptions = {}
): Promise<NDArray> {
  const module = ar1._wasmModule;
  const {
    invert = false,
    assumeUnique = false,
    kind = ISIN_AUTO,
  } = options;

  const resultPtr = module._ndarray_isin(
    ar1._wasmPtr,
    ar2._wasmPtr,
    invert,
    assumeUnique,
    kind
  );

  if (resultPtr === 0) {
    throw new Error('isin failed');
  }

  return NDArray._fromPtr(resultPtr, module);
}

/**
 * Test whether each element of ar1 is also present in ar2 (flat version).
 *
 * @param ar1 - Input array (will be flattened)
 * @param ar2 - Values to test against
 * @param options - Configuration options
 * @returns 1D boolean array of membership results
 */
export async function in1d(
  ar1: NDArray,
  ar2: NDArray,
  options: IsinOptions = {}
): Promise<NDArray> {
  const module = ar1._wasmModule;
  const {
    invert = false,
    assumeUnique = false,
    kind = ISIN_AUTO,
  } = options;

  const resultPtr = module._ndarray_in1d(
    ar1._wasmPtr,
    ar2._wasmPtr,
    invert,
    assumeUnique,
    kind
  );

  if (resultPtr === 0) {
    throw new Error('in1d failed');
  }

  return NDArray._fromPtr(resultPtr, module);
}

/* ============ ediff1d (TypeScript-only) ============ */

export interface Ediff1dOptions {
  /** Values to prepend to the differences */
  toPrepend?: NDArray | number[];
  /** Values to append to the differences */
  toAppend?: NDArray | number[];
}

/**
 * Compute the differences between consecutive elements of an array.
 *
 * This is implemented in TypeScript as it's straightforward and doesn't
 * require WASM for performance.
 *
 * @param arr - Input array
 * @param options - Configuration options
 * @returns 1D array of differences with optional prepend/append
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2, 4, 7]);
 * const result = await ediff1d(arr);
 * // NDArray [1, 2, 3]
 *
 * const withEnds = await ediff1d(arr, { toPrepend: [0], toAppend: [10] });
 * // NDArray [0, 1, 2, 3, 10]
 * ```
 */
export async function ediff1d(
  arr: NDArray,
  options: Ediff1dOptions = {}
): Promise<NDArray> {
  const { toPrepend, toAppend } = options;

  // Flatten input
  const flat = await arr.flatten();
  const n = flat.size;

  if (n === 0) {
    // Handle empty array
    const parts: number[] = [];
    if (toPrepend) {
      const prepend = Array.isArray(toPrepend)
        ? toPrepend
        : await toPrepend.toArray();
      parts.push(...(prepend as number[]));
    }
    if (toAppend) {
      const append = Array.isArray(toAppend)
        ? toAppend
        : await toAppend.toArray();
      parts.push(...(append as number[]));
    }
    if (parts.length === 0) {
      flat.dispose();
      return NDArray.fromArray([]);
    }
    flat.dispose();
    return NDArray.fromArray(parts);
  }

  // Compute differences
  const diffs: number[] = [];

  // Prepend if specified
  if (toPrepend) {
    const prepend = Array.isArray(toPrepend)
      ? toPrepend
      : await toPrepend.toArray();
    diffs.push(...(prepend as number[]));
  }

  // Compute element-wise differences
  for (let i = 1; i < n; i++) {
    const curr = flat.getFlat(i);
    const prev = flat.getFlat(i - 1);
    diffs.push(curr - prev);
  }

  // Append if specified
  if (toAppend) {
    const append = Array.isArray(toAppend)
      ? toAppend
      : await toAppend.toArray();
    diffs.push(...(append as number[]));
  }

  flat.dispose();

  return NDArray.fromArray(diffs);
}
