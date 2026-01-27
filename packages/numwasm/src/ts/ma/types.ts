/**
 * NumJS Masked Arrays - Types and Constants
 *
 * Core type definitions, constants, and error classes for the masked arrays module.
 * Compatible with NumPy's numpy.ma module.
 */

import { NDArray } from '../NDArray.js';
import { DType } from '../types.js';

/**
 * Sentinel value indicating no mask.
 * When mask is nomask, all elements are valid.
 */
export const nomask: unique symbol = Symbol('nomask');

/**
 * Type for mask: either a boolean NDArray or nomask sentinel.
 */
export type MaskType = NDArray | typeof nomask;

/**
 * Masked scalar constant.
 * Represents a single masked (invalid) value.
 * Used when accessing a masked element.
 */
export const masked = Object.freeze({
  __masked__: true as const,
  toString(): string {
    return '--';
  },
  valueOf(): number {
    return NaN;
  },
});

/**
 * Type for the masked constant.
 */
export type MaskedConstant = typeof masked;

/**
 * Check if a value is the masked constant.
 *
 * @param x - Value to check
 * @returns true if x is the masked constant
 *
 * @example
 * ```typescript
 * const val = ma.getFlat(0);
 * if (isMaskedConstant(val)) {
 *   console.log('Element is masked');
 * }
 * ```
 */
export function isMaskedConstant(x: unknown): x is MaskedConstant {
  return (
    x !== null &&
    typeof x === 'object' &&
    '__masked__' in x &&
    (x as Record<string, unknown>).__masked__ === true
  );
}

/**
 * Error for masked array operations.
 */
export class MaskedArrayError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'MaskedArrayError';
  }
}

/**
 * Error for mask operations.
 */
export class MaskError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'MaskError';
  }
}

/**
 * Default fill values for each dtype.
 * Used when masked values need to be filled for operations.
 */
export const defaultFillValues: Partial<Record<DType, number | boolean | string>> = {
  [DType.Bool]: true,
  [DType.Int8]: 999999,
  [DType.Int16]: 999999,
  [DType.Int32]: 999999,
  [DType.Int64]: 999999,
  [DType.Uint8]: 999999,
  [DType.Uint16]: 999999,
  [DType.Uint32]: 999999,
  [DType.Uint64]: 999999,
  [DType.Float16]: 1e20,
  [DType.Float32]: 1e20,
  [DType.Float64]: 1e20,
  [DType.Complex64]: 1e20,
  [DType.Complex128]: 1e20,
  [DType.String]: 'N/A',
};

/**
 * Get the default fill value for a dtype.
 *
 * @param dtype - Data type
 * @returns Default fill value for the dtype
 */
export function getDefaultFillValue(dtype: DType): number | boolean | string {
  return defaultFillValues[dtype] ?? 1e20;
}

/**
 * Maximum fill values for dtype (used for min operations).
 * When computing min, masked values should be replaced with +Infinity
 * so they don't affect the result.
 *
 * @param dtype - Data type
 * @returns Maximum value for the dtype
 */
export function maximumFillValue(dtype: DType): number {
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
    case DType.Float16:
    case DType.Float32:
    case DType.Float64:
      return Number.POSITIVE_INFINITY;
    default:
      return 1e20;
  }
}

/**
 * Minimum fill values for dtype (used for max operations).
 * When computing max, masked values should be replaced with -Infinity
 * so they don't affect the result.
 *
 * @param dtype - Data type
 * @returns Minimum value for the dtype
 */
export function minimumFillValue(dtype: DType): number {
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
    case DType.Float16:
    case DType.Float32:
    case DType.Float64:
      return Number.NEGATIVE_INFINITY;
    default:
      return -1e20;
  }
}
