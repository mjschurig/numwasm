/**
 * NumJS Iterators
 *
 * Provides iteration utilities for NDArray, following NumPy conventions.
 */

import type { NDArray } from "../_core/NDArray.js";

/**
 * Iterator over array elements in row-major (C) order.
 * Implements JavaScript Iterator protocol.
 *
 * @example
 * ```typescript
 * const arr = await NDArray.fromArray([1, 2, 3, 4, 5, 6], [2, 3]);
 * for (const value of arr.flat) {
 *     console.log(value);  // 1, 2, 3, 4, 5, 6
 * }
 * ```
 */
export class FlatIterator implements IterableIterator<number> {
  private _arr: NDArray;
  private _index: number = 0;
  private _size: number;
  private _shape: number[];

  constructor(arr: NDArray) {
    this._arr = arr;
    this._size = arr.size;
    this._shape = arr.shape;
  }

  [Symbol.iterator](): IterableIterator<number> {
    return this;
  }

  next(): IteratorResult<number> {
    if (this._index >= this._size) {
      return { done: true, value: undefined };
    }

    const indices = this._unravelIndex(this._index);
    const value = this._arr.get(...indices);
    this._index++;

    return { done: false, value };
  }

  /**
   * Convert flat index to N-dimensional indices (row-major order).
   */
  private _unravelIndex(flatIndex: number): number[] {
    const indices: number[] = new Array(this._shape.length);
    let remaining = flatIndex;

    for (let i = this._shape.length - 1; i >= 0; i--) {
      indices[i] = remaining % this._shape[i];
      remaining = Math.floor(remaining / this._shape[i]);
    }

    return indices;
  }

  /**
   * Current flat index position.
   */
  get index(): number {
    return this._index;
  }

  /**
   * Current N-dimensional coordinates.
   */
  get coords(): number[] {
    return this._unravelIndex(this._index);
  }

  /**
   * Reset iterator to beginning.
   */
  reset(): void {
    this._index = 0;
  }
}

/**
 * Iterate over array elements.
 *
 * @param arr - Array to iterate
 * @yields Element values in row-major order
 *
 * @example
 * ```typescript
 * for (const x of nditer(arr)) {
 *     console.log(x);
 * }
 * ```
 */
export function* nditer(arr: NDArray): Generator<number> {
  const size = arr.size;
  const shape = arr.shape;
  const ndim = shape.length;
  const indices = new Array(ndim).fill(0);

  for (let flat = 0; flat < size; flat++) {
    yield arr.get(...indices);

    // Increment indices (row-major order)
    for (let d = ndim - 1; d >= 0; d--) {
      indices[d]++;
      if (indices[d] < shape[d]) break;
      indices[d] = 0;
    }
  }
}

/**
 * Iterate over array with indices and values.
 *
 * @param arr - Array to iterate
 * @yields [indices, value] pairs
 *
 * @example
 * ```typescript
 * for (const [idx, val] of ndenumerate(arr)) {
 *     console.log(`arr[${idx.join(', ')}] = ${val}`);
 * }
 * // Output for 2x3 array:
 * // arr[0, 0] = 1
 * // arr[0, 1] = 2
 * // arr[0, 2] = 3
 * // arr[1, 0] = 4
 * // arr[1, 1] = 5
 * // arr[1, 2] = 6
 * ```
 */
export function* ndenumerate(arr: NDArray): Generator<[number[], number]> {
  const size = arr.size;
  const shape = arr.shape;
  const ndim = shape.length;
  const indices = new Array(ndim).fill(0);

  for (let flat = 0; flat < size; flat++) {
    yield [indices.slice(), arr.get(...indices)];

    // Increment indices (row-major order)
    for (let d = ndim - 1; d >= 0; d--) {
      indices[d]++;
      if (indices[d] < shape[d]) break;
      indices[d] = 0;
    }
  }
}

/**
 * Generate all index combinations for a shape.
 *
 * @param shape - Dimensions
 * @yields Index arrays
 *
 * @example
 * ```typescript
 * for (const idx of ndindex(2, 3)) {
 *     console.log(idx);
 * }
 * // Output:
 * // [0, 0]
 * // [0, 1]
 * // [0, 2]
 * // [1, 0]
 * // [1, 1]
 * // [1, 2]
 * ```
 */
export function* ndindex(...shape: number[]): Generator<number[]> {
  const ndim = shape.length;
  if (ndim === 0) return;

  const size = shape.reduce((a, b) => a * b, 1);
  if (size === 0) return;

  const indices = new Array(ndim).fill(0);

  for (let flat = 0; flat < size; flat++) {
    yield indices.slice();

    // Increment indices (row-major order)
    for (let d = ndim - 1; d >= 0; d--) {
      indices[d]++;
      if (indices[d] < shape[d]) break;
      indices[d] = 0;
    }
  }
}
