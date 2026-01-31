/**
 * Array tiling and repetition functions for NumJS-WASM
 *
 * Functions for tiling, repeating, and padding arrays.
 */

import { NDArray } from "../_core/NDArray.js";
import { slice, Slice } from "../slice.js";
import { concatenate } from "./join.js";

/* ============ Helper Functions ============ */

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

function copyArrayData(src: NDArray, dst: NDArray): void {
  // Use WASM copyto for efficient copying
  const module = dst._wasmModule;
  module._ndarray_copyto(dst._wasmPtr, src._wasmPtr, 0);
}

/* ============ Tiling Functions ============ */

/**
 * Construct an array by repeating arr the number of times given by reps.
 *
 * @param arr - Input array
 * @param reps - Number of repetitions along each axis
 * @returns Tiled array
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2]);
 * const result = await tile(arr, 3);
 * // [1, 2, 1, 2, 1, 2]
 *
 * const arr2 = await NDArray.fromArray([[1, 2], [3, 4]]);
 * const result2 = await tile(arr2, [2, 1]);
 * // [[1, 2], [3, 4], [1, 2], [3, 4]]
 * ```
 */
export async function tile(
  arr: NDArray,
  reps: number | number[],
): Promise<NDArray> {
  const repsArr = Array.isArray(reps) ? reps : [reps];

  // Extend reps to match or exceed arr.ndim
  const d = Math.max(arr.ndim, repsArr.length);

  // Pad reps with 1s at the beginning if needed
  const paddedReps = new Array(d - repsArr.length).fill(1).concat(repsArr);

  // Pad array shape with 1s at the beginning if needed
  let current = arr;
  let needsDispose = false;
  if (arr.ndim < d) {
    const newShape = new Array(d - arr.ndim).fill(1).concat([...arr.shape]);
    current = arr.reshape(newShape);
    needsDispose = true;
  }

  // Tile along each axis using concatenate (like NumPy)
  for (let axis = 0; axis < d; axis++) {
    const nrep = paddedReps[axis];
    if (nrep === 1) continue;

    // Create array of nrep copies and concatenate along axis
    const copies: NDArray[] = [];
    for (let i = 0; i < nrep; i++) {
      copies.push(current);
    }
    const tiled = concatenate(copies, axis);

    // Dispose previous if it was created by us
    if (needsDispose) {
      current.dispose();
    }
    current = tiled;
    needsDispose = true;
  }

  return current;
}

/**
 * Repeat elements of an array.
 *
 * @param arr - Input array
 * @param repeats - Number of repetitions for each element
 * @param axis - Axis along which to repeat (default: flatten first)
 * @returns Array with repeated elements
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2, 3]);
 * const result = await repeat(arr, 2);
 * // [1, 1, 2, 2, 3, 3]
 *
 * const arr2 = await NDArray.fromArray([[1, 2], [3, 4]]);
 * const result2 = await repeat(arr2, 2, 0);
 * // [[1, 2], [1, 2], [3, 4], [3, 4]]
 * ```
 */
export async function repeat(
  arr: NDArray,
  repeats: number | number[],
  axis?: number,
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
  axis: number,
): Promise<NDArray> {
  axis = axis < 0 ? axis + arr.ndim : axis;

  const axisSize = arr.shape[axis];
  const repsArr =
    typeof repeats === "number" ? new Array(axisSize).fill(repeats) : repeats;

  if (repsArr.length !== axisSize) {
    throw new Error("repeats must have same length as axis");
  }

  // Check if all repeats are the same (common case)
  const allSame = repsArr.every((r) => r === repsArr[0]);

  if (allSame && repsArr[0] > 0) {
    // Fast path: use expandDims + tile-like approach with concatenate
    // For repeat([1,2,3], 2) -> [1,1,2,2,3,3]
    // We expand to shape [..., axisSize, 1, ...], tile along new axis, then reshape
    const nrep = repsArr[0];

    // Add a new axis after the repeat axis
    const expanded = arr.expandDims(axis + 1);

    // Tile along the new axis using concatenate
    const copies: NDArray[] = [];
    for (let i = 0; i < nrep; i++) {
      copies.push(expanded);
    }
    const tiled = concatenate(copies, axis + 1);
    expanded.dispose();

    // Reshape to merge the repeated axis
    const newShape = [...arr.shape];
    newShape[axis] = arr.shape[axis] * nrep;
    const result = tiled.reshape(newShape);
    tiled.dispose();

    return result;
  }

  // General case: different repeat counts per element
  // Build arrays for each element and concatenate
  const slices: NDArray[] = [];
  const toDispose: NDArray[] = [];

  for (let i = 0; i < axisSize; i++) {
    const numReps = repsArr[i];
    if (numReps === 0) continue;

    // Get slice at index i along axis
    const sliceSpec = arr.shape.map((_, dim) =>
      dim === axis ? slice(i, i + 1) : slice(null),
    );
    const srcSlice = arr.slice(sliceSpec);
    toDispose.push(srcSlice);

    // Add numReps copies
    for (let r = 0; r < numReps; r++) {
      slices.push(srcSlice);
    }
  }

  const result =
    slices.length > 0
      ? concatenate(slices, axis)
      : await NDArray.empty(
          [...arr.shape.slice(0, axis), 0, ...arr.shape.slice(axis + 1)],
          { dtype: arr.dtype },
        );

  // Dispose temporary slices
  for (const s of toDispose) {
    s.dispose();
  }

  return result;
}

/**
 * Pad an array.
 *
 * @param arr - Input array
 * @param pad_width - Number of values padded to edges of each axis
 * @param mode - Padding mode (default: 'constant')
 * @param constant_values - Value to use for constant padding (default: 0)
 * @returns Padded array
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2, 3]);
 * const result = await pad(arr, 2, 'constant', 0);
 * // [0, 0, 1, 2, 3, 0, 0]
 *
 * const arr2 = await NDArray.fromArray([[1, 2], [3, 4]]);
 * const result2 = await pad(arr2, [[1, 1], [2, 2]], 'constant', 0);
 * // 0 padding around the original array
 * ```
 */
export async function pad(
  arr: NDArray,
  pad_width: number | [number, number] | [number, number][],
  mode: string = "constant",
  constant_values: number = 0,
): Promise<NDArray> {
  // Normalize pad_width to [[before, after], ...] format
  let padSpec: [number, number][];

  if (typeof pad_width === "number") {
    padSpec = arr.shape.map(() => [pad_width, pad_width]);
  } else if (Array.isArray(pad_width) && typeof pad_width[0] === "number") {
    const [before, after] = pad_width as [number, number];
    padSpec = arr.shape.map(() => [before, after]);
  } else {
    padSpec = pad_width as [number, number][];
  }

  if (padSpec.length !== arr.ndim) {
    throw new Error("pad_width must have same length as array dimensions");
  }

  // Compute output shape
  const outShape = arr.shape.map((s, i) => s + padSpec[i][0] + padSpec[i][1]);

  // Create padded array with constant fill
  const result = await NDArray.full(outShape, constant_values, {
    dtype: arr.dtype,
  });

  // Build slice to place original data
  const indices: Slice[] = padSpec.map(([before], i) =>
    slice(before, before + arr.shape[i]),
  );

  // Copy original data into the center
  const target = result.slice(indices);
  copyArrayData(arr, target);
  target.dispose();

  // Handle other modes
  if (mode !== "constant") {
    // For edge mode, reflect, etc., we'd need additional implementation
    // For now, only constant mode is fully supported
    if (mode === "edge") {
      padEdge(result, arr.shape, padSpec);
    }
    // Other modes could be added: 'reflect', 'symmetric', 'wrap'
  }

  return result;
}

function padEdge(
  result: NDArray,
  originalShape: number[],
  padSpec: [number, number][],
): void {
  // Edge padding - replicate edge values
  const outShape = result.shape;

  for (let i = 0; i < result.size; i++) {
    const outIdx = flatToMulti(i, outShape);

    // Map to source index, clamping to original bounds
    const srcIdx = outIdx.map((idx, d) => {
      const [before] = padSpec[d];
      const shifted = idx - before;
      return Math.max(0, Math.min(shifted, originalShape[d] - 1));
    });

    // Check if this is a padded position
    const isPadded = outIdx.some((idx, d) => {
      const [before] = padSpec[d];
      return idx < before || idx >= before + originalShape[d];
    });

    if (isPadded) {
      // Get source value and set in result
      const srcFlat = multiToFlat(
        srcIdx.map((idx, d) => idx + padSpec[d][0]),
        outShape,
      );
      result.setFlat(i, result.getFlat(srcFlat));
    }
  }
}
