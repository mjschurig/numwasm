/**
 * NumJS Masked Arrays - Fill Value Operations
 *
 * Utilities for managing fill values in masked arrays.
 * Compatible with NumPy's numpy.ma module.
 */

import { NDArray } from '../NDArray.js';
import { DType } from '../types.js';
import {
  getDefaultFillValue,
  MaskType,
} from './types.js';

// Forward declaration to avoid circular dependency
type MaskedArrayLike = {
  _data: NDArray;
  _mask: MaskType;
  _fill_value: number | string | boolean;
  fill_value: number | string | boolean;
};

/**
 * Get the default fill value for an object based on its dtype.
 *
 * @param obj - NDArray, MaskedArray, or dtype value
 * @returns Default fill value for the dtype
 *
 * @example
 * ```typescript
 * default_fill_value(arr);  // Returns default for arr's dtype
 * default_fill_value(DType.Float64);  // Returns 1e20
 * default_fill_value(DType.Int32);  // Returns 999999
 * ```
 */
export function default_fill_value(
  obj: NDArray | { dtype: DType } | DType
): number | string | boolean {
  let dtype: DType;

  if (typeof obj === 'number') {
    // Direct DType value
    dtype = obj;
  } else if (obj instanceof NDArray) {
    dtype = obj.dtype;
  } else if (typeof obj === 'object' && 'dtype' in obj) {
    dtype = obj.dtype;
  } else {
    return 1e20;
  }

  return getDefaultFillValue(dtype);
}

/**
 * Find a common fill value for two masked arrays.
 *
 * Returns the common fill value if both arrays have the same fill value,
 * otherwise returns null.
 *
 * @param a - First array
 * @param b - Second array
 * @returns Common fill value, or null if different
 *
 * @example
 * ```typescript
 * const fv = common_fill_value(ma1, ma2);
 * if (fv !== null) {
 *   console.log('Arrays have common fill value:', fv);
 * }
 * ```
 */
export function common_fill_value(
  a: unknown,
  b: unknown
): number | string | boolean | null {
  // Get fill values from both arrays
  const fv1 = getFillValue(a);
  const fv2 = getFillValue(b);

  if (fv1 === null || fv2 === null) {
    return null;
  }

  // Check if they're equal
  if (fv1 === fv2) {
    return fv1;
  }

  // For numbers, check with tolerance for floating point
  if (typeof fv1 === 'number' && typeof fv2 === 'number') {
    if (Math.abs(fv1 - fv2) < 1e-10 * Math.max(Math.abs(fv1), Math.abs(fv2))) {
      return fv1;
    }
  }

  return null;
}

/**
 * Set the fill value of a masked array.
 *
 * @param a - MaskedArray to modify
 * @param fill_value - New fill value
 *
 * @example
 * ```typescript
 * set_fill_value(ma, 0);
 * // ma.fill_value is now 0
 * ```
 */
export function set_fill_value(
  a: unknown,
  fill_value: number | string | boolean
): void {
  if (
    a !== null &&
    typeof a === 'object' &&
    '_fill_value' in a
  ) {
    (a as MaskedArrayLike)._fill_value = fill_value;
  } else {
    throw new Error('Cannot set fill_value on non-MaskedArray');
  }
}

/**
 * Get the fill value from an object.
 *
 * @param obj - MaskedArray or other object
 * @returns Fill value, or null if not a MaskedArray
 */
function getFillValue(obj: unknown): number | string | boolean | null {
  if (
    obj !== null &&
    typeof obj === 'object' &&
    '_fill_value' in obj
  ) {
    return (obj as MaskedArrayLike)._fill_value;
  }

  if (obj instanceof NDArray) {
    return getDefaultFillValue(obj.dtype);
  }

  return null;
}

/**
 * Get the appropriate fill value for a reduction operation.
 *
 * @param dtype - Data type
 * @param operation - Reduction operation ('sum', 'prod', 'min', 'max')
 * @returns Fill value that acts as identity for the operation
 *
 * @example
 * ```typescript
 * getReductionFillValue(DType.Float64, 'sum');  // 0
 * getReductionFillValue(DType.Float64, 'prod');  // 1
 * getReductionFillValue(DType.Float64, 'min');  // +Infinity
 * getReductionFillValue(DType.Float64, 'max');  // -Infinity
 * ```
 */
export function getReductionFillValue(
  dtype: DType,
  operation: 'sum' | 'prod' | 'min' | 'max'
): number {
  switch (operation) {
    case 'sum':
      return 0;
    case 'prod':
      return 1;
    case 'min':
      // Use maximum value so masked elements don't affect min
      switch (dtype) {
        case DType.Int8:
          return 127;
        case DType.Int16:
          return 32767;
        case DType.Int32:
          return 2147483647;
        case DType.Uint8:
          return 255;
        case DType.Uint16:
          return 65535;
        case DType.Uint32:
          return 4294967295;
        default:
          return Number.POSITIVE_INFINITY;
      }
    case 'max':
      // Use minimum value so masked elements don't affect max
      switch (dtype) {
        case DType.Int8:
          return -128;
        case DType.Int16:
          return -32768;
        case DType.Int32:
          return -2147483648;
        case DType.Uint8:
        case DType.Uint16:
        case DType.Uint32:
          return 0;
        default:
          return Number.NEGATIVE_INFINITY;
      }
  }
}
