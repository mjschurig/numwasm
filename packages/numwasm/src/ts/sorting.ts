/**
 * Sorting functions for NDArray
 *
 * Provides WASM-accelerated sorting operations with full NumPy-compatible
 * axis support.
 */

import { NDArray } from "./_core/NDArray.js";
import { DType } from "./types.js";

export type SortKind = "quicksort" | "mergesort" | "heapsort" | "stable";

const KIND_MAP: Record<SortKind, number> = {
  quicksort: 0,
  mergesort: 1,
  heapsort: 2,
  stable: 3,
};

/**
 * Return a sorted copy of an array.
 *
 * @param a - Array to sort
 * @param axis - Axis along which to sort. null flattens the array.
 *               Default is -1 (last axis).
 * @param kind - Sorting algorithm: 'quicksort', 'mergesort', 'heapsort', or 'stable'
 * @returns Sorted copy of the array
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([[3, 1], [2, 4]]);
 * const sorted = await sort(arr);       // Sort along last axis
 * const flat = await sort(arr, null);   // Sort flattened
 * ```
 */
export async function sort(
  a: NDArray,
  axis: number | null = -1,
  kind: SortKind = "quicksort",
): Promise<NDArray> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_sort_copy(
    a._wasmPtr,
    axisVal,
    KIND_MAP[kind],
  );
  if (resultPtr === 0) throw new Error("sort failed");
  return NDArray._fromPtr(resultPtr, a._wasmModule);
}

/**
 * Returns the indices that would sort an array.
 *
 * @param a - Array to sort
 * @param axis - Axis along which to sort. null flattens the array.
 *               Default is -1 (last axis).
 * @param kind - Sorting algorithm
 * @returns Array of indices that sort the array
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([3, 1, 2]);
 * const idx = await argsort(arr);  // [1, 2, 0]
 * ```
 */
export async function argsort(
  a: NDArray,
  axis: number | null = -1,
  kind: SortKind = "quicksort",
): Promise<NDArray> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_argsort(
    a._wasmPtr,
    axisVal,
    KIND_MAP[kind],
  );
  if (resultPtr === 0) throw new Error("argsort failed");
  return NDArray._fromPtr(resultPtr, a._wasmModule);
}

/**
 * Return a partitioned copy of an array.
 *
 * Creates a copy of the array with its elements rearranged so that
 * elements smaller than the kth element are moved before it and
 * elements greater are moved after it.
 *
 * @param a - Array to partition
 * @param kth - Index to partition around
 * @param axis - Axis along which to partition (default: -1)
 * @returns Partitioned copy of the array
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([3, 1, 4, 1, 5, 9, 2]);
 * const part = await partition(arr, 3);  // Elements < arr[3] before, > after
 * ```
 */
export async function partition(
  a: NDArray,
  kth: number,
  axis: number = -1,
): Promise<NDArray> {
  const resultPtr = a._wasmModule._ndarray_partition(a._wasmPtr, kth, axis);
  if (resultPtr === 0) throw new Error("partition failed");
  return NDArray._fromPtr(resultPtr, a._wasmModule);
}

/**
 * Returns the indices that would partition an array.
 *
 * @param a - Array to partition
 * @param kth - Index to partition around
 * @param axis - Axis along which to partition (default: -1)
 * @returns Array of indices
 */
export async function argpartition(
  a: NDArray,
  kth: number,
  axis: number = -1,
): Promise<NDArray> {
  const resultPtr = a._wasmModule._ndarray_argpartition(a._wasmPtr, kth, axis);
  if (resultPtr === 0) throw new Error("argpartition failed");
  return NDArray._fromPtr(resultPtr, a._wasmModule);
}

/**
 * Sort a complex array using the real part first, then the imaginary part.
 *
 * The input array is flattened before sorting. Non-complex arrays are
 * converted to complex in the output.
 *
 * @param a - Input array (will be flattened)
 * @returns Sorted 1D complex array
 *
 * @example
 * ```typescript
 * const c = await NDArray.fromArray([1+2j, 2-1j, 3-2j, 3+5j], { dtype: DType.Complex128 });
 * const sorted = sort_complex(c);
 * // [1+2j, 2-1j, 3-2j, 3+5j] (sorted by real, then imaginary)
 *
 * // Real arrays are converted to complex:
 * const r = await NDArray.fromArray([5, 3, 6, 2, 1]);
 * const sorted2 = sort_complex(r);
 * // [1+0j, 2+0j, 3+0j, 5+0j, 6+0j]
 * ```
 */
export function sort_complex(a: NDArray): NDArray {
  // Flatten the array
  const flat = a.flatten();
  const size = flat.size;
  const module = a._wasmModule;

  // Determine output dtype based on input
  let outputDtype: DType;
  if (flat.dtype === DType.Complex64 || flat.dtype === DType.Complex128) {
    outputDtype = flat.dtype;
  } else if (
    flat.dtype === DType.Int8 ||
    flat.dtype === DType.Uint8 ||
    flat.dtype === DType.Int16 ||
    flat.dtype === DType.Uint16
  ) {
    outputDtype = DType.Complex64;
  } else {
    outputDtype = DType.Complex128;
  }

  // Convert to complex if needed
  let source: NDArray;
  if (flat.dtype !== outputDtype) {
    source = flat.astype(outputDtype);
  } else {
    source = flat;
  }

  // Extract values for sorting
  const values: Array<{ real: number; imag: number }> = [];
  for (let i = 0; i < size; i++) {
    const real = module._ndarray_get_complex_real(source._wasmPtr, i);
    const imag = module._ndarray_get_complex_imag(source._wasmPtr, i);
    values.push({ real, imag });
  }

  // Sort by real first, then by imaginary
  values.sort((a, b) => {
    if (a.real !== b.real) return a.real - b.real;
    return a.imag - b.imag;
  });

  // Create result array
  const shapePtr = module._malloc(4);
  module.setValue(shapePtr, size, "i32");
  const resultPtr = module._ndarray_empty(1, shapePtr, outputDtype);
  module._free(shapePtr);

  if (resultPtr === 0) {
    throw new Error("Failed to create result array for sort_complex");
  }

  // Fill result array with sorted values
  for (let i = 0; i < size; i++) {
    module._ndarray_set_complex(resultPtr, i, values[i].real, values[i].imag);
  }

  return NDArray._fromPtr(resultPtr, module);
}

/**
 * Perform an indirect stable sort using a sequence of keys.
 *
 * Given multiple sorting keys, lexsort returns an array of integer indices
 * that describes the sort order by multiple columns. The last key in the
 * sequence is used for the primary sort order, the second-to-last key for
 * the secondary sort order, and so on.
 *
 * @param keys - Array of keys to sort by. Each key should be 1-D with the same length.
 *               The last key is the primary sort key.
 * @returns Array of indices that sort the keys
 *
 * @example
 * ```typescript
 * // Sort by last name, then by first name
 * const surnames = await NDArray.fromArray(['Smith', 'Jones', 'Smith', 'Adams']);
 * const firstNames = await NDArray.fromArray(['John', 'Mary', 'Alice', 'Bob']);
 * const idx = await lexsort([firstNames, surnames]);
 * // idx sorts first by surname, then by first name within each surname
 * ```
 *
 * @example
 * ```typescript
 * // Sort by column 2, then by column 1
 * const a = await NDArray.fromArray([1, 5, 1, 4, 3, 4, 4]);
 * const b = await NDArray.fromArray([9, 4, 0, 4, 0, 2, 1]);
 * const idx = await lexsort([b, a]);
 * // First sorts by 'a', then by 'b' within equal 'a' values
 * ```
 */
export async function lexsort(keys: NDArray[]): Promise<NDArray> {
  if (keys.length === 0) {
    throw new Error("lexsort requires at least one key array");
  }

  // Validate all keys have the same size
  const n = keys[0].size;
  for (let i = 1; i < keys.length; i++) {
    if (keys[i].size !== n) {
      throw new Error("all keys must have the same size");
    }
  }

  // Create array of indices [0, 1, 2, ..., n-1]
  const indices = new Int32Array(n);
  for (let i = 0; i < n; i++) {
    indices[i] = i;
  }

  // Extract values from all keys for comparison
  const keyValues: number[][] = [];
  for (const key of keys) {
    const vals: number[] = [];
    for (let i = 0; i < n; i++) {
      vals.push(key.getFlat(i) as number);
    }
    keyValues.push(vals);
  }

  // Sort indices using a stable sort, comparing from last key to first
  // (last key is primary sort key)
  const indicesArray = Array.from(indices);
  indicesArray.sort((i, j) => {
    // Compare from last key (primary) to first key (least significant)
    for (let k = keyValues.length - 1; k >= 0; k--) {
      const a = keyValues[k][i];
      const b = keyValues[k][j];
      if (a < b) return -1;
      if (a > b) return 1;
    }
    // Stable sort: preserve original order for equal elements
    return i - j;
  });

  return NDArray.fromArray(indicesArray, undefined, { dtype: DType.Int32 });
}
