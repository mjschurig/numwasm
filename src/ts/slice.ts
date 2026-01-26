/**
 * Slice and indexing utilities for NumJS-WASM
 *
 * Provides Python-style slice semantics for array indexing.
 * Adapted from NumPy's indexing behavior.
 */

/**
 * Represents a slice with start:stop:step semantics.
 * Mirrors Python's slice object behavior.
 *
 * @example
 * new Slice(1, 5)        // 1:5
 * new Slice(null, 5)     // :5
 * new Slice(1, null, 2)  // 1::2
 * new Slice(null, null, -1) // ::-1 (reverse)
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
   * Implements the same logic as Python's slice.indices().
   *
   * @param length - Length of the dimension
   * @returns [start, stop, step, sliceLength] where sliceLength is the number of elements
   */
  indices(length: number): [number, number, number, number] {
    const step = this.step ?? 1;
    let start: number;
    let stop: number;

    if (step > 0) {
      // Forward iteration
      start = this.start ?? 0;
      stop = this.stop ?? length;

      // Handle negative indices
      if (start < 0) {
        start = Math.max(0, length + start);
      }
      if (stop < 0) {
        stop = Math.max(0, length + stop);
      }

      // Clamp to bounds
      start = Math.min(start, length);
      stop = Math.min(stop, length);
    } else {
      // Backward iteration
      start = this.start ?? length - 1;
      stop = this.stop ?? -length - 1;

      // Handle negative indices
      if (start < 0) {
        start = Math.max(-1, length + start);
      } else if (start >= length) {
        start = length - 1;
      }

      if (stop < -1) {
        stop = -1;
      } else if (stop >= length) {
        stop = length;
      }
    }

    // Calculate slice length
    let sliceLength: number;
    if (step > 0) {
      sliceLength = stop > start ? Math.ceil((stop - start) / step) : 0;
    } else {
      sliceLength = start > stop ? Math.ceil((start - stop) / -step) : 0;
    }

    return [start, stop, step, sliceLength];
  }

  /**
   * String representation for debugging
   */
  toString(): string {
    const startStr = this.start !== null ? String(this.start) : '';
    const stopStr = this.stop !== null ? String(this.stop) : '';
    const stepStr = this.step !== null && this.step !== 1 ? `:${this.step}` : '';
    return `${startStr}:${stopStr}${stepStr}`;
  }
}

/**
 * Factory function for creating slices with cleaner syntax.
 *
 * @param start - Start index (null for default)
 * @param stop - Stop index (null for default)
 * @param step - Step size (null for 1)
 * @returns Slice object
 *
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

/**
 * Sentinel for ellipsis (...) in indexing.
 * Expands to fill remaining dimensions with full slices.
 *
 * @example
 * arr.slice([ellipsis, 0])  // arr[..., 0] in NumPy
 */
export const ellipsis: unique symbol = Symbol('ellipsis');
export type Ellipsis = typeof ellipsis;

/**
 * Sentinel for newaxis in indexing.
 * Inserts a new dimension of size 1.
 *
 * @example
 * arr.slice([newaxis, slice()])  // arr[np.newaxis, :] in NumPy
 */
export const newaxis: unique symbol = Symbol('newaxis');
export type Newaxis = typeof newaxis;

/**
 * Type for valid index elements in a slice operation.
 */
export type IndexElement =
  | number // Integer index (selects single element, removes dimension)
  | Slice // Slice object (selects range, keeps dimension)
  | Ellipsis // Ellipsis (expands to fill dimensions)
  | Newaxis; // New axis (inserts dimension of size 1)

/**
 * Index specification types for C interop
 */
export const INDEX_TYPE_INTEGER = 0;
export const INDEX_TYPE_SLICE = 1;
export const INDEX_TYPE_NEWAXIS = 2;
export const INDEX_TYPE_ELLIPSIS = 3;

/**
 * Parsed index specification for a single index element.
 * Used for communication with C layer.
 */
export interface IndexSpec {
  type: number;
  value: number; // Integer index value (if type=INTEGER), or ellipsis length
  start: number; // Slice start (if type=SLICE)
  stop: number; // Slice stop (if type=SLICE)
  step: number; // Slice step (if type=SLICE)
}

/**
 * Expand ellipsis and validate index dimensions.
 *
 * @param indices - Array of index elements (may contain ellipsis)
 * @param ndim - Number of dimensions in the array being indexed
 * @returns Expanded index array with ellipsis replaced by full slices
 * @throws Error if multiple ellipses or too many indices
 */
export function expandEllipsis(
  indices: IndexElement[],
  ndim: number
): IndexElement[] {
  // Count non-ellipsis, non-newaxis indices (these consume dimensions)
  let numDimIndices = 0;
  let ellipsisIndex = -1;

  for (let i = 0; i < indices.length; i++) {
    const idx = indices[i];
    if (idx === ellipsis) {
      if (ellipsisIndex !== -1) {
        throw new Error('an index can only have a single ellipsis (...)');
      }
      ellipsisIndex = i;
    } else if (idx !== newaxis) {
      numDimIndices++;
    }
  }

  if (numDimIndices > ndim) {
    throw new Error(
      `too many indices for array: array is ${ndim}-dimensional, ` +
        `but ${numDimIndices} were indexed`
    );
  }

  // If no ellipsis, return as-is
  if (ellipsisIndex === -1) {
    return indices;
  }

  // Expand ellipsis to full slices
  const ellipsisLength = ndim - numDimIndices;
  const expanded: IndexElement[] = [];

  for (let i = 0; i < indices.length; i++) {
    if (i === ellipsisIndex) {
      // Replace ellipsis with full slices
      for (let j = 0; j < ellipsisLength; j++) {
        expanded.push(new Slice());
      }
    } else {
      expanded.push(indices[i]);
    }
  }

  return expanded;
}

/**
 * Build IndexSpec array from parsed indices for C interop.
 *
 * @param indices - Expanded index elements (no ellipsis)
 * @param shape - Shape of the array being indexed
 * @returns Array of IndexSpec objects
 */
export function buildIndexSpecs(
  indices: IndexElement[],
  shape: number[]
): IndexSpec[] {
  const specs: IndexSpec[] = [];
  let dimIndex = 0;

  for (const idx of indices) {
    if (typeof idx === 'number') {
      // Integer index
      let value = idx;
      const dimSize = shape[dimIndex];

      // Handle negative indices
      if (value < 0) {
        value += dimSize;
      }

      // Bounds check
      if (value < 0 || value >= dimSize) {
        throw new RangeError(
          `index ${idx} is out of bounds for axis ${dimIndex} with size ${dimSize}`
        );
      }

      specs.push({
        type: INDEX_TYPE_INTEGER,
        value,
        start: 0,
        stop: 0,
        step: 1,
      });
      dimIndex++;
    } else if (idx instanceof Slice) {
      // Slice
      const dimSize = shape[dimIndex];
      const [start, stop, step, _length] = idx.indices(dimSize);

      specs.push({
        type: INDEX_TYPE_SLICE,
        value: 0,
        start,
        stop,
        step,
      });
      dimIndex++;
    } else if (idx === newaxis) {
      // Newaxis - doesn't consume a dimension
      specs.push({
        type: INDEX_TYPE_NEWAXIS,
        value: 0,
        start: 0,
        stop: 0,
        step: 1,
      });
    }
    // Note: ellipsis should already be expanded
  }

  return specs;
}

/**
 * Compute the result shape from index specifications.
 *
 * @param specs - Array of IndexSpec
 * @param shape - Original array shape
 * @returns Result shape after indexing
 */
export function computeResultShape(specs: IndexSpec[], shape: number[]): number[] {
  const resultShape: number[] = [];
  let dimIndex = 0;

  for (const spec of specs) {
    switch (spec.type) {
      case INDEX_TYPE_INTEGER:
        // Integer index removes dimension
        dimIndex++;
        break;

      case INDEX_TYPE_SLICE: {
        // Slice keeps dimension with computed length
        const step = spec.step || 1;
        let length: number;

        if (step > 0) {
          length = Math.max(0, Math.ceil((spec.stop - spec.start) / step));
        } else {
          length = Math.max(0, Math.ceil((spec.start - spec.stop) / -step));
        }

        resultShape.push(length);
        dimIndex++;
        break;
      }

      case INDEX_TYPE_NEWAXIS:
        // Newaxis adds dimension of size 1
        resultShape.push(1);
        break;
    }
  }

  // Add remaining dimensions (implicit full slices)
  while (dimIndex < shape.length) {
    resultShape.push(shape[dimIndex]);
    dimIndex++;
  }

  return resultShape;
}
