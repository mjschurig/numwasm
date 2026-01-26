/**
 * NumJS Masked Arrays - Creation Functions
 *
 * Factory functions for creating masked arrays.
 * Compatible with NumPy's numpy.ma module.
 */

import { NDArray } from '../NDArray.js';
import { DType } from '../types.js';
import { MaskedArray } from './core.js';
import { nomask, MaskType, getDefaultFillValue } from './types.js';

/**
 * Create a masked array (async).
 *
 * @param data - Input data (array, NDArray, or MaskedArray)
 * @param mask - Boolean mask (true = masked)
 * @param dtype - Data type
 * @param copy - If true, copy data
 * @param fill_value - Fill value for masked elements
 * @returns Promise resolving to MaskedArray
 *
 * @example
 * ```typescript
 * const ma = await masked_array([1, 2, 3], [false, true, false]);
 * // ma.data = [1, 2, 3], ma.mask = [false, true, false]
 * ```
 */
export async function masked_array(
  data: NDArray | number[] | MaskedArray,
  mask: MaskType | boolean[] | boolean = nomask,
  dtype: DType = DType.Float64,
  _copy: boolean = false,
  fill_value: number | null = null
): Promise<MaskedArray> {
  return MaskedArray.create(data, mask as boolean[] | MaskType, dtype, fill_value, false, true);
}

/**
 * Alias for masked_array.
 */
export const array = masked_array;

/**
 * Mask elements equal to a value.
 */
export async function masked_equal(x: NDArray | number[], value: number): Promise<MaskedArray> {
  const arr = Array.isArray(x) ? await NDArray.fromArray(x) : x;
  const mask = await NDArray.empty(arr.shape, { dtype: DType.Bool });

  for (let i = 0; i < arr.size; i++) {
    mask.setFlat(i, arr.getFlat(i) === value ? 1 : 0);
  }

  return new MaskedArray(arr, mask);
}

/**
 * Mask elements not equal to a value.
 */
export async function masked_not_equal(x: NDArray | number[], value: number): Promise<MaskedArray> {
  const arr = Array.isArray(x) ? await NDArray.fromArray(x) : x;
  const mask = await NDArray.empty(arr.shape, { dtype: DType.Bool });

  for (let i = 0; i < arr.size; i++) {
    mask.setFlat(i, arr.getFlat(i) !== value ? 1 : 0);
  }

  return new MaskedArray(arr, mask);
}

/**
 * Mask elements greater than a value.
 */
export async function masked_greater(x: NDArray | number[], value: number): Promise<MaskedArray> {
  const arr = Array.isArray(x) ? await NDArray.fromArray(x) : x;
  const mask = await NDArray.empty(arr.shape, { dtype: DType.Bool });

  for (let i = 0; i < arr.size; i++) {
    mask.setFlat(i, arr.getFlat(i) > value ? 1 : 0);
  }

  return new MaskedArray(arr, mask);
}

/**
 * Mask elements greater than or equal to a value.
 */
export async function masked_greater_equal(
  x: NDArray | number[],
  value: number
): Promise<MaskedArray> {
  const arr = Array.isArray(x) ? await NDArray.fromArray(x) : x;
  const mask = await NDArray.empty(arr.shape, { dtype: DType.Bool });

  for (let i = 0; i < arr.size; i++) {
    mask.setFlat(i, arr.getFlat(i) >= value ? 1 : 0);
  }

  return new MaskedArray(arr, mask);
}

/**
 * Mask elements less than a value.
 */
export async function masked_less(x: NDArray | number[], value: number): Promise<MaskedArray> {
  const arr = Array.isArray(x) ? await NDArray.fromArray(x) : x;
  const mask = await NDArray.empty(arr.shape, { dtype: DType.Bool });

  for (let i = 0; i < arr.size; i++) {
    mask.setFlat(i, arr.getFlat(i) < value ? 1 : 0);
  }

  return new MaskedArray(arr, mask);
}

/**
 * Mask elements less than or equal to a value.
 */
export async function masked_less_equal(x: NDArray | number[], value: number): Promise<MaskedArray> {
  const arr = Array.isArray(x) ? await NDArray.fromArray(x) : x;
  const mask = await NDArray.empty(arr.shape, { dtype: DType.Bool });

  for (let i = 0; i < arr.size; i++) {
    mask.setFlat(i, arr.getFlat(i) <= value ? 1 : 0);
  }

  return new MaskedArray(arr, mask);
}

/**
 * Mask elements inside an interval [v1, v2].
 */
export async function masked_inside(
  x: NDArray | number[],
  v1: number,
  v2: number
): Promise<MaskedArray> {
  const arr = Array.isArray(x) ? await NDArray.fromArray(x) : x;
  const mask = await NDArray.empty(arr.shape, { dtype: DType.Bool });
  const lo = Math.min(v1, v2);
  const hi = Math.max(v1, v2);

  for (let i = 0; i < arr.size; i++) {
    const val = arr.getFlat(i);
    mask.setFlat(i, val >= lo && val <= hi ? 1 : 0);
  }

  return new MaskedArray(arr, mask);
}

/**
 * Mask elements outside an interval (v1, v2).
 */
export async function masked_outside(
  x: NDArray | number[],
  v1: number,
  v2: number
): Promise<MaskedArray> {
  const arr = Array.isArray(x) ? await NDArray.fromArray(x) : x;
  const mask = await NDArray.empty(arr.shape, { dtype: DType.Bool });
  const lo = Math.min(v1, v2);
  const hi = Math.max(v1, v2);

  for (let i = 0; i < arr.size; i++) {
    const val = arr.getFlat(i);
    mask.setFlat(i, val < lo || val > hi ? 1 : 0);
  }

  return new MaskedArray(arr, mask);
}

/**
 * Mask elements where condition is true.
 */
export async function masked_where(
  condition: NDArray | boolean[],
  x: NDArray | number[]
): Promise<MaskedArray> {
  const arr = Array.isArray(x) ? await NDArray.fromArray(x) : x;
  let cond: NDArray;

  if (Array.isArray(condition)) {
    cond = await NDArray.fromArray(
      condition.map((c) => (c ? 1 : 0)),
      undefined,
      { dtype: DType.Bool }
    );
  } else {
    cond = condition.dtype === DType.Bool ? condition : condition.astype(DType.Bool);
  }

  return new MaskedArray(arr, cond);
}

/**
 * Mask NaN and Inf values.
 */
export async function masked_invalid(x: NDArray | number[]): Promise<MaskedArray> {
  const arr = Array.isArray(x) ? await NDArray.fromArray(x) : x;
  const mask = await NDArray.empty(arr.shape, { dtype: DType.Bool });

  for (let i = 0; i < arr.size; i++) {
    const val = arr.getFlat(i);
    mask.setFlat(i, !Number.isFinite(val) ? 1 : 0);
  }

  return new MaskedArray(arr, mask);
}

/**
 * Mask values close to a given value.
 */
export async function masked_values(
  x: NDArray | number[],
  value: number,
  rtol: number = 1e-5,
  atol: number = 1e-8
): Promise<MaskedArray> {
  const arr = Array.isArray(x) ? await NDArray.fromArray(x) : x;
  const mask = await NDArray.empty(arr.shape, { dtype: DType.Bool });

  for (let i = 0; i < arr.size; i++) {
    const val = arr.getFlat(i);
    const close = Math.abs(val - value) <= atol + rtol * Math.abs(value);
    mask.setFlat(i, close ? 1 : 0);
  }

  return new MaskedArray(arr, mask);
}

/**
 * Create masked array of zeros.
 */
export async function zeros(
  shape: number | number[],
  dtype: DType = DType.Float64
): Promise<MaskedArray> {
  const shapeArr = typeof shape === 'number' ? [shape] : shape;
  const data = await NDArray.zeros(shapeArr, { dtype });
  return new MaskedArray(data);
}

/**
 * Create masked array of ones.
 */
export async function ones(
  shape: number | number[],
  dtype: DType = DType.Float64
): Promise<MaskedArray> {
  const shapeArr = typeof shape === 'number' ? [shape] : shape;
  const data = await NDArray.ones(shapeArr, { dtype });
  return new MaskedArray(data);
}

/**
 * Create empty masked array.
 */
export async function empty(
  shape: number | number[],
  dtype: DType = DType.Float64
): Promise<MaskedArray> {
  const shapeArr = typeof shape === 'number' ? [shape] : shape;
  const data = await NDArray.empty(shapeArr, { dtype });
  return new MaskedArray(data);
}

/**
 * Create fully masked array.
 */
export async function masked_all(
  shape: number | number[],
  dtype: DType = DType.Float64
): Promise<MaskedArray> {
  const shapeArr = typeof shape === 'number' ? [shape] : shape;
  const data = await NDArray.zeros(shapeArr, { dtype });
  const mask = await NDArray.full(shapeArr, 1, { dtype: DType.Bool });
  return new MaskedArray(data, mask);
}

/**
 * Create fully masked array like another.
 */
export async function masked_all_like(prototype: MaskedArray | NDArray): Promise<MaskedArray> {
  const shape = prototype.shape;
  const dtype = prototype.dtype;
  const data = await NDArray.zeros(shape, { dtype });
  const mask = await NDArray.full(shape, 1, { dtype: DType.Bool });
  return new MaskedArray(data, mask);
}

/**
 * Create zeros_like masked array.
 */
export async function zeros_like(prototype: MaskedArray | NDArray): Promise<MaskedArray> {
  const shape = prototype.shape;
  const dtype = prototype.dtype;
  const data = await NDArray.zeros(shape, { dtype });
  return new MaskedArray(data);
}

/**
 * Create ones_like masked array.
 */
export async function ones_like(prototype: MaskedArray | NDArray): Promise<MaskedArray> {
  const shape = prototype.shape;
  const dtype = prototype.dtype;
  const data = await NDArray.ones(shape, { dtype });
  return new MaskedArray(data);
}

/**
 * Create empty_like masked array.
 */
export async function empty_like(prototype: MaskedArray | NDArray): Promise<MaskedArray> {
  const shape = prototype.shape;
  const dtype = prototype.dtype;
  const data = await NDArray.empty(shape, { dtype });
  return new MaskedArray(data);
}

/**
 * Create masked array from function.
 */
export async function fromfunction(
  shape: number[],
  fn: (...indices: number[]) => [number, boolean],
  dtype: DType = DType.Float64
): Promise<MaskedArray> {
  const data = await NDArray.empty(shape, { dtype });
  const mask = await NDArray.zeros(shape, { dtype: DType.Bool });

  const size = shape.reduce((a, b) => a * b, 1);

  for (let flatIdx = 0; flatIdx < size; flatIdx++) {
    // Convert flat index to multi-index
    const indices: number[] = [];
    let remaining = flatIdx;
    for (let i = shape.length - 1; i >= 0; i--) {
      indices.unshift(remaining % shape[i]);
      remaining = Math.floor(remaining / shape[i]);
    }

    const [value, isMasked] = fn(...indices);
    data.setFlat(flatIdx, value);
    mask.setFlat(flatIdx, isMasked ? 1 : 0);
  }

  return new MaskedArray(data, mask, getDefaultFillValue(dtype));
}
