/**
 * Sorting functions for NDArray
 *
 * Provides WASM-accelerated sorting operations with full NumPy-compatible
 * axis support.
 */

import { NDArray } from './NDArray.js';

export type SortKind = 'quicksort' | 'mergesort' | 'heapsort' | 'stable';

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
  kind: SortKind = 'quicksort'
): Promise<NDArray> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_sort_copy(
    a._wasmPtr,
    axisVal,
    KIND_MAP[kind]
  );
  if (resultPtr === 0) throw new Error('sort failed');
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
  kind: SortKind = 'quicksort'
): Promise<NDArray> {
  const axisVal = axis === null ? -2147483648 : axis;
  const resultPtr = a._wasmModule._ndarray_argsort(
    a._wasmPtr,
    axisVal,
    KIND_MAP[kind]
  );
  if (resultPtr === 0) throw new Error('argsort failed');
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
  axis: number = -1
): Promise<NDArray> {
  const resultPtr = a._wasmModule._ndarray_partition(a._wasmPtr, kth, axis);
  if (resultPtr === 0) throw new Error('partition failed');
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
  axis: number = -1
): Promise<NDArray> {
  const resultPtr = a._wasmModule._ndarray_argpartition(a._wasmPtr, kth, axis);
  if (resultPtr === 0) throw new Error('argpartition failed');
  return NDArray._fromPtr(resultPtr, a._wasmModule);
}
