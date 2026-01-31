/**
 * Array insert and delete functions for NumJS-WASM
 *
 * Functions for inserting and deleting elements from arrays.
 */

import { NDArray } from "../_core/NDArray.js";
import { DType } from "../types.js";
import { Slice, slice } from "../slice.js";
import { take } from "../indexing.js";
import { concatenate } from "./join.js";

/**
 * Insert values along the given axis before the given indices.
 *
 * @param arr - Input array
 * @param obj - Index or indices before which values are inserted
 * @param values - Values to insert
 * @param axis - Axis along which to insert (default: flatten first)
 * @returns Array with inserted values
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2, 3, 4, 5]);
 * const result = await insert(arr, 2, [10, 11]);
 * // [1, 2, 10, 11, 3, 4, 5]
 * ```
 */
export function insert(
  arr: NDArray,
  obj: number | number[],
  values: NDArray,
  axis?: number,
): NDArray {
  if (axis === undefined) {
    // Insert into flattened array
    const flat = arr.ravel();
    const result = insertAlongAxis(flat, obj, values, 0);
    flat.dispose();
    return result;
  }

  return insertAlongAxis(arr, obj, values, axis);
}

function insertAlongAxis(
  arr: NDArray,
  obj: number | number[],
  values: NDArray,
  axis: number,
): NDArray {
  axis = axis < 0 ? axis + arr.ndim : axis;

  const indices = typeof obj === "number" ? [obj] : [...obj];
  const sortedIndices = indices.sort((a, b) => a - b);

  const parts: NDArray[] = [];
  const toDispose: NDArray[] = [];
  let lastIdx = 0;

  for (let i = 0; i < sortedIndices.length; i++) {
    const idx = sortedIndices[i];

    // Add part before insertion point
    if (idx > lastIdx) {
      const sliceIndices: Slice[] = [];
      for (let d = 0; d < arr.ndim; d++) {
        sliceIndices.push(d === axis ? slice(lastIdx, idx) : slice(null, null));
      }
      const part = arr.slice(sliceIndices);
      parts.push(part);
      toDispose.push(part);
    }

    // Add inserted values
    parts.push(values);
    lastIdx = idx;
  }

  // Add remaining part
  if (lastIdx < arr.shape[axis]) {
    const sliceIndices: Slice[] = [];
    for (let d = 0; d < arr.ndim; d++) {
      sliceIndices.push(d === axis ? slice(lastIdx, null) : slice(null, null));
    }
    const part = arr.slice(sliceIndices);
    parts.push(part);
    toDispose.push(part);
  }

  const result = concatenate(parts, axis);

  for (const p of toDispose) {
    p.dispose();
  }

  return result;
}

/**
 * Return a new array with sub-arrays along an axis deleted.
 *
 * @param arr - Input array
 * @param obj - Index or indices to delete
 * @param axis - Axis along which to delete (default: flatten first)
 * @returns Array with elements deleted
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2, 3, 4, 5]);
 * const result = await deleteArr(arr, [1, 3]);
 * // [1, 3, 5]
 * ```
 */
export async function deleteArr(
  arr: NDArray,
  obj: number | number[],
  axis?: number,
): Promise<NDArray> {
  if (axis === undefined) {
    const flat = arr.ravel();
    const result = await deleteAlongAxis(flat, obj, 0);
    flat.dispose();
    return result;
  }

  return deleteAlongAxis(arr, obj, axis);
}

async function deleteAlongAxis(
  arr: NDArray,
  obj: number | number[],
  axis: number,
): Promise<NDArray> {
  axis = axis < 0 ? axis + arr.ndim : axis;

  const indices = typeof obj === "number" ? [obj] : obj;
  const deleteSet = new Set(
    indices.map((i) => (i < 0 ? i + arr.shape[axis] : i)),
  );

  // Keep indices not in deleteSet
  const keepIndices: number[] = [];
  for (let i = 0; i < arr.shape[axis]; i++) {
    if (!deleteSet.has(i)) keepIndices.push(i);
  }

  // Use take to select kept indices
  const indicesArr = await NDArray.fromArray(
    keepIndices,
    [keepIndices.length],
    {
      dtype: DType.Int32,
    },
  );
  const result = take(arr, indicesArr, axis);
  indicesArr.dispose();

  return result;
}

/* ============ Copying Functions ============ */

/**
 * Copy values from one array to another, broadcasting as necessary.
 *
 * @param dst - Destination array
 * @param src - Source array
 * @param casting - Type casting rule (default: 'same_kind')
 * @param where - Optional boolean mask array
 */
export function copyto(
  dst: NDArray,
  src: NDArray,
  _casting: string = "same_kind",
  where?: NDArray,
): void {
  if (dst.size !== src.size) {
    throw new Error("dst and src must have same size");
  }
  if (where !== undefined && where.size !== dst.size) {
    throw new Error("where mask must have same size as dst");
  }

  const module = dst._wasmModule;
  const wherePtr = where ? where._wasmPtr : 0;
  const result = module._ndarray_copyto(dst._wasmPtr, src._wasmPtr, wherePtr);

  if (result !== 0) {
    throw new Error("copyto failed");
  }
}

/**
 * Convert the input to an array.
 *
 * @param a - Input data
 * @param dtype - Optional data type
 * @returns Array from input
 */
export async function asarray(
  a: NDArray | number[] | number,
  dtype?: DType,
): Promise<NDArray> {
  if (a instanceof NDArray) {
    if (dtype !== undefined && a.dtype !== dtype) {
      return a.astype(dtype);
    }
    return a; // Return as-is (not a copy)
  }

  const data = typeof a === "number" ? [a] : a;
  const shape = typeof a === "number" ? [] : [data.length];
  return NDArray.fromArray(data, shape, { dtype });
}
