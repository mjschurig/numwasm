/**
 * Array splitting functions for NumJS-WASM
 *
 * Functions for splitting and unpacking arrays.
 */

import { NDArray } from "../_core/NDArray.js";
import { Slice, slice } from "../slice.js";

/**
 * Split an array into multiple sub-arrays as views.
 * Allows unequal division (unlike split).
 *
 * @param arr - Array to split
 * @param indices_or_sections - Number of sections or split indices
 * @param axis - Axis along which to split (default: 0)
 * @returns List of sub-arrays
 */
export function array_split(
  arr: NDArray,
  indices_or_sections: number | number[],
  axis: number = 0,
): NDArray[] {
  // Normalize axis
  axis = axis < 0 ? axis + arr.ndim : axis;
  if (axis < 0 || axis >= arr.ndim) {
    throw new Error(
      `axis ${axis} out of bounds for array with ${arr.ndim} dimensions`,
    );
  }

  const axisSize = arr.shape[axis];
  let divPoints: number[];

  if (typeof indices_or_sections === "number") {
    // Number of sections - compute division points
    const sections = indices_or_sections;
    const base = Math.floor(axisSize / sections);
    const extras = axisSize % sections;

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

  for (let i = 0; i < divPoints.length - 1; i++) {
    const start = divPoints[i];
    const stop = divPoints[i + 1];

    // Build slice indices
    const indices: (Slice | number)[] = [];
    for (let d = 0; d < arr.ndim; d++) {
      if (d === axis) {
        indices.push(slice(start, stop));
      } else {
        indices.push(slice(null, null)); // Full slice
      }
    }

    const subArr = arr.slice(indices);
    result.push(subArr.copy()); // Copy to own data
    subArr.dispose();
  }

  return result;
}

/**
 * Split an array into multiple sub-arrays of equal size.
 * Raises error if array cannot be split into equal sections.
 *
 * @param arr - Array to split
 * @param indices_or_sections - Number of equal sections or split indices
 * @param axis - Axis along which to split (default: 0)
 * @returns List of sub-arrays
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2, 3, 4, 5, 6]);
 * const result = await split(arr, 3);
 * // [[1, 2], [3, 4], [5, 6]]
 * ```
 */
export function split(
  arr: NDArray,
  indices_or_sections: number | number[],
  axis: number = 0,
): NDArray[] {
  // Normalize axis
  const normalizedAxis = axis < 0 ? axis + arr.ndim : axis;

  if (typeof indices_or_sections === "number") {
    const axisSize = arr.shape[normalizedAxis];
    if (axisSize % indices_or_sections !== 0) {
      throw new Error("array split does not result in an equal division");
    }
  }

  return array_split(arr, indices_or_sections, axis);
}

/**
 * Split an array into multiple sub-arrays vertically (row-wise).
 *
 * @param arr - Array with at least 2 dimensions
 * @param indices_or_sections - Number of sections or split indices
 * @returns List of sub-arrays
 */
export function vsplit(
  arr: NDArray,
  indices_or_sections: number | number[],
): NDArray[] {
  if (arr.ndim < 2) {
    throw new Error("vsplit only works on arrays of 2 or more dimensions");
  }
  return split(arr, indices_or_sections, 0);
}

/**
 * Split an array into multiple sub-arrays horizontally (column-wise).
 *
 * @param arr - Array
 * @param indices_or_sections - Number of sections or split indices
 * @returns List of sub-arrays
 */
export function hsplit(
  arr: NDArray,
  indices_or_sections: number | number[],
): NDArray[] {
  if (arr.ndim < 1) {
    throw new Error("hsplit only works on arrays of 1 or more dimensions");
  }
  if (arr.ndim > 1) {
    return split(arr, indices_or_sections, 1);
  }
  return split(arr, indices_or_sections, 0);
}

/**
 * Split array into multiple sub-arrays along the 3rd axis (depth).
 *
 * @param arr - Array with at least 3 dimensions
 * @param indices_or_sections - Number of sections or split indices
 * @returns List of sub-arrays
 */
export function dsplit(
  arr: NDArray,
  indices_or_sections: number | number[],
): NDArray[] {
  if (arr.ndim < 3) {
    throw new Error("dsplit only works on arrays of 3 or more dimensions");
  }
  return split(arr, indices_or_sections, 2);
}

/**
 * Unpack the array along a given axis.
 * The inverse of stack.
 *
 * @param arr - Array to unpack
 * @param axis - Axis along which to unpack (default: 0)
 * @returns List of arrays, each with one fewer dimension
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([[1, 2], [3, 4]]);
 * const [a, b] = await unstack(arr, 0);
 * // a: [1, 2], b: [3, 4]
 * ```
 */
export function unstack(arr: NDArray, axis: number = 0): NDArray[] {
  // Normalize axis
  axis = axis < 0 ? axis + arr.ndim : axis;

  if (arr.ndim === 0) {
    throw new Error("cannot unstack 0-d array");
  }

  const numArrays = arr.shape[axis];
  return split(arr, numArrays, axis).map((sub) => {
    // Remove the singleton dimension created by split
    // squeeze() returns a view, so we need to copy it before disposing sub
    const squeezed = sub.squeeze(axis);
    const copied = squeezed.copy();
    squeezed.dispose();
    sub.dispose();
    return copied;
  });
}
