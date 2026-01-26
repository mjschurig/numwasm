/**
 * NumJS Masked Arrays - Mask Operations
 *
 * Low-level mask manipulation functions for working with boolean masks.
 * Compatible with NumPy's numpy.ma module.
 */

import { NDArray } from '../NDArray.js';
import { DType } from '../types.js';
import { nomask, MaskType, isMaskedConstant } from './types.js';

// Forward declaration - MaskedArray is imported dynamically to avoid circular deps
type MaskedArrayLike = {
  _data: NDArray;
  _mask: MaskType;
};

/**
 * Create a boolean mask from an array (async version).
 *
 * @param m - Input data (array, boolean, or mask)
 * @param copy - If true, return a copy of the mask
 * @param shrink - If true, return nomask if all values are False
 * @returns Promise resolving to boolean mask array or nomask
 *
 * @example
 * ```typescript
 * const mask = await make_mask([false, true, false]);
 * // Returns boolean NDArray
 *
 * const mask2 = await make_mask([false, false, false], false, true);
 * // Returns nomask (all False, shrink=true)
 * ```
 */
export async function make_mask(
  m: MaskType | boolean[] | boolean | NDArray | number[],
  copy: boolean = false,
  shrink: boolean = true
): Promise<MaskType> {
  // nomask passes through
  if (m === nomask) {
    return nomask;
  }

  // Boolean scalar
  if (typeof m === 'boolean') {
    if (!m && shrink) {
      return nomask;
    }
    // Can't create a scalar mask - caller should handle this
    return nomask;
  }

  // Convert to NDArray if needed
  let maskArr: NDArray;
  if (m instanceof NDArray) {
    if (m.dtype === DType.Bool) {
      maskArr = copy ? m.copy() : m;
    } else {
      // Convert to boolean
      maskArr = await NDArray.empty(m.shape, { dtype: DType.Bool });
      for (let i = 0; i < m.size; i++) {
        maskArr.setFlat(i, m.getFlat(i) ? 1 : 0);
      }
    }
  } else if (Array.isArray(m)) {
    // Convert boolean/number array to NDArray
    const boolData = m.map((v) => (v ? 1 : 0));
    maskArr = await NDArray.fromArray(boolData, undefined, { dtype: DType.Bool });
  } else {
    return nomask;
  }

  // Shrink to nomask if all False
  if (shrink && !anyTrue(maskArr)) {
    return nomask;
  }

  return maskArr;
}

/**
 * Create a boolean mask from an array (sync version, requires existing NDArray).
 * Use this when you already have an NDArray mask.
 */
export function make_mask_sync(
  m: MaskType | NDArray,
  copy: boolean = false,
  shrink: boolean = true
): MaskType {
  // nomask passes through
  if (m === nomask) {
    return nomask;
  }

  if (!(m instanceof NDArray)) {
    throw new Error('make_mask_sync requires NDArray input, use make_mask for arrays');
  }

  let maskArr: NDArray;
  if (m.dtype === DType.Bool) {
    maskArr = copy ? m.copy() : m;
  } else {
    // Need async to create new array - can't do this sync
    throw new Error('Converting non-bool array to mask requires async make_mask');
  }

  // Shrink to nomask if all False
  if (shrink && !anyTrue(maskArr)) {
    return nomask;
  }

  return maskArr;
}

/**
 * Create a mask with all values False (no masking).
 *
 * @param shape - Shape of the mask to create
 * @returns nomask (since all values would be False)
 *
 * @example
 * ```typescript
 * const mask = make_mask_none([3, 4]);
 * // Returns nomask
 * ```
 */
export function make_mask_none(_shape: number[]): MaskType {
  // An all-False mask is equivalent to nomask
  return nomask;
}

/**
 * Get the mask from an array.
 *
 * @param a - Input array (MaskedArray, NDArray, or other)
 * @returns Mask of the array, or nomask if no mask
 *
 * @example
 * ```typescript
 * const mask = getmask(maskedArray);
 * if (mask !== nomask) {
 *   console.log('Array has a mask');
 * }
 * ```
 */
export function getmask(a: unknown): MaskType {
  // Check if it's a MaskedArray-like object
  if (
    a !== null &&
    typeof a === 'object' &&
    '_mask' in a &&
    ('_data' in a || 'data' in a)
  ) {
    return (a as MaskedArrayLike)._mask;
  }

  // Regular arrays have no mask
  return nomask;
}

/**
 * Get the mask from an array as an NDArray (never nomask).
 *
 * @param a - Input array
 * @returns Promise resolving to mask as boolean NDArray (all False if no mask)
 *
 * @example
 * ```typescript
 * const maskArr = await getmaskarray(maskedArray);
 * // Always returns an NDArray, even if all False
 * ```
 */
export async function getmaskarray(a: unknown): Promise<NDArray> {
  const mask = getmask(a);

  if (mask === nomask) {
    // Create an all-False mask with the same shape
    let shape: number[] = [1];
    if (a instanceof NDArray) {
      shape = a.shape;
    } else if (
      a !== null &&
      typeof a === 'object' &&
      '_data' in a &&
      (a as MaskedArrayLike)._data instanceof NDArray
    ) {
      shape = (a as MaskedArrayLike)._data.shape;
    }
    return await NDArray.zeros(shape, { dtype: DType.Bool });
  }

  return mask;
}

/**
 * Get the data from an array.
 *
 * @param a - Input array (MaskedArray, NDArray, or other)
 * @returns Underlying data array
 *
 * @example
 * ```typescript
 * const data = getdata(maskedArray);
 * // Returns the underlying NDArray
 * ```
 */
export function getdata(a: unknown): NDArray {
  // Check if it's a MaskedArray-like object
  if (
    a !== null &&
    typeof a === 'object' &&
    '_data' in a &&
    (a as MaskedArrayLike)._data instanceof NDArray
  ) {
    return (a as MaskedArrayLike)._data;
  }

  // If it's already an NDArray, return it
  if (a instanceof NDArray) {
    return a;
  }

  throw new Error('Cannot get data from non-array object');
}

/**
 * Check if an object is a valid mask.
 *
 * @param m - Object to check
 * @returns true if m is a valid mask (boolean NDArray or nomask)
 *
 * @example
 * ```typescript
 * is_mask(nomask);  // true
 * is_mask(boolArray);  // true
 * is_mask(floatArray);  // false
 * ```
 */
export function is_mask(m: unknown): boolean {
  if (m === nomask) {
    return true;
  }

  if (m instanceof NDArray) {
    return m.dtype === DType.Bool;
  }

  return false;
}

/**
 * Check if an array has any masked elements.
 *
 * @param x - Input to check
 * @returns true if any element is masked
 *
 * @example
 * ```typescript
 * is_masked(maskedArray);  // true if any element is masked
 * is_masked(regularArray);  // false (no mask)
 * ```
 */
export function is_masked(x: unknown): boolean {
  // Check for masked constant
  if (isMaskedConstant(x)) {
    return true;
  }

  // Get the mask
  const mask = getmask(x);

  if (mask === nomask) {
    return false;
  }

  // Check if any element is True (masked)
  return anyTrue(mask);
}

/**
 * Combine two masks with logical OR (async version).
 *
 * @param m1 - First mask
 * @param m2 - Second mask
 * @param shrink - If true, return nomask if result is all False
 * @returns Promise resolving to combined mask (m1 OR m2)
 *
 * @example
 * ```typescript
 * const combined = await mask_or(mask1, mask2);
 * // Element is masked if masked in either mask
 * ```
 */
export async function mask_or(
  m1: MaskType,
  m2: MaskType,
  shrink: boolean = true
): Promise<MaskType> {
  // Handle nomask cases
  if (m1 === nomask && m2 === nomask) {
    return nomask;
  }
  if (m1 === nomask) {
    return shrink && !anyTrue(m2 as NDArray) ? nomask : m2;
  }
  if (m2 === nomask) {
    return shrink && !anyTrue(m1 as NDArray) ? nomask : m1;
  }

  // Both are NDArrays - combine with OR
  const arr1 = m1 as NDArray;
  const arr2 = m2 as NDArray;

  // Shapes must be compatible (broadcasting would be needed for different shapes)
  // For simplicity, we require same shape
  if (!shapesEqual(arr1.shape, arr2.shape)) {
    throw new Error(
      `Mask shapes must match: ${arr1.shape} vs ${arr2.shape}`
    );
  }

  const result = await NDArray.zeros(arr1.shape, { dtype: DType.Bool });
  for (let i = 0; i < result.size; i++) {
    result.setFlat(i, arr1.getFlat(i) || arr2.getFlat(i) ? 1 : 0);
  }

  // Shrink to nomask if all False
  if (shrink && !anyTrue(result)) {
    return nomask;
  }

  return result;
}

/**
 * Combine two masks with logical OR (sync version for NDArray masks).
 * Use when you already have NDArray masks and want sync behavior.
 */
export function mask_or_sync(
  m1: MaskType,
  m2: MaskType,
  result: NDArray,
  shrink: boolean = true
): MaskType {
  // Handle nomask cases
  if (m1 === nomask && m2 === nomask) {
    return nomask;
  }
  if (m1 === nomask) {
    return shrink && !anyTrue(m2 as NDArray) ? nomask : m2;
  }
  if (m2 === nomask) {
    return shrink && !anyTrue(m1 as NDArray) ? nomask : m1;
  }

  // Both are NDArrays - combine with OR
  const arr1 = m1 as NDArray;
  const arr2 = m2 as NDArray;

  if (!shapesEqual(arr1.shape, arr2.shape)) {
    throw new Error(`Mask shapes must match: ${arr1.shape} vs ${arr2.shape}`);
  }

  for (let i = 0; i < result.size; i++) {
    result.setFlat(i, arr1.getFlat(i) || arr2.getFlat(i) ? 1 : 0);
  }

  if (shrink && !anyTrue(result)) {
    return nomask;
  }

  return result;
}

/**
 * Flatten a mask to 1D.
 *
 * @param mask - Input mask
 * @returns Flattened 1D mask
 */
export function flatten_mask(mask: MaskType): MaskType {
  if (mask === nomask) {
    return nomask;
  }

  return (mask as NDArray).flatten();
}

/**
 * Reshape a mask to a new shape.
 *
 * @param mask - Input mask
 * @param shape - New shape
 * @returns Reshaped mask
 */
export function reshape_mask(mask: MaskType, shape: number[]): MaskType {
  if (mask === nomask) {
    return nomask;
  }

  return (mask as NDArray).reshape(shape);
}

/**
 * Broadcast a mask to a target shape (async version).
 *
 * @param mask - Input mask
 * @param targetShape - Shape to broadcast to
 * @returns Promise resolving to broadcast mask
 */
export async function broadcast_mask(mask: MaskType, targetShape: number[]): Promise<MaskType> {
  if (mask === nomask) {
    return nomask;
  }

  const maskArr = mask as NDArray;

  // If already the right shape, return as-is
  if (shapesEqual(maskArr.shape, targetShape)) {
    return maskArr;
  }

  // Manual broadcasting for simple cases
  const result = await NDArray.zeros(targetShape, { dtype: DType.Bool });

  // Calculate how to broadcast
  const maskShape = maskArr.shape;
  const maskNdim = maskShape.length;
  const targetNdim = targetShape.length;

  // Pad mask shape with leading 1s
  const paddedMaskShape = new Array(targetNdim).fill(1);
  for (let i = 0; i < maskNdim; i++) {
    paddedMaskShape[targetNdim - maskNdim + i] = maskShape[i];
  }

  // Check broadcast compatibility
  for (let i = 0; i < targetNdim; i++) {
    if (paddedMaskShape[i] !== 1 && paddedMaskShape[i] !== targetShape[i]) {
      throw new Error(
        `Cannot broadcast mask shape ${maskShape} to ${targetShape}`
      );
    }
  }

  // Copy values with broadcasting
  for (let i = 0; i < result.size; i++) {
    // Calculate indices in target shape
    const indices = flatToMultiIndex(i, targetShape);

    // Map to mask indices (use modulo for broadcast dimensions)
    const maskIndices = new Array(maskNdim);
    for (let j = 0; j < maskNdim; j++) {
      const targetIdx = targetNdim - maskNdim + j;
      maskIndices[j] =
        maskShape[j] === 1 ? 0 : indices[targetIdx];
    }

    const maskFlatIdx = multiToFlatIndex(maskIndices, maskShape);
    result.setFlat(i, maskArr.getFlat(maskFlatIdx));
  }

  return result;
}

/**
 * Broadcast a mask to a target shape (sync version, result pre-allocated).
 */
export function broadcast_mask_sync(mask: MaskType, targetShape: number[], result: NDArray): MaskType {
  if (mask === nomask) {
    return nomask;
  }

  const maskArr = mask as NDArray;

  if (shapesEqual(maskArr.shape, targetShape)) {
    return maskArr;
  }

  const maskShape = maskArr.shape;
  const maskNdim = maskShape.length;
  const targetNdim = targetShape.length;

  const paddedMaskShape = new Array(targetNdim).fill(1);
  for (let i = 0; i < maskNdim; i++) {
    paddedMaskShape[targetNdim - maskNdim + i] = maskShape[i];
  }

  for (let i = 0; i < targetNdim; i++) {
    if (paddedMaskShape[i] !== 1 && paddedMaskShape[i] !== targetShape[i]) {
      throw new Error(`Cannot broadcast mask shape ${maskShape} to ${targetShape}`);
    }
  }

  for (let i = 0; i < result.size; i++) {
    const indices = flatToMultiIndex(i, targetShape);
    const maskIndices = new Array(maskNdim);
    for (let j = 0; j < maskNdim; j++) {
      const targetIdx = targetNdim - maskNdim + j;
      maskIndices[j] = maskShape[j] === 1 ? 0 : indices[targetIdx];
    }
    const maskFlatIdx = multiToFlatIndex(maskIndices, maskShape);
    result.setFlat(i, maskArr.getFlat(maskFlatIdx));
  }

  return result;
}

// ============ Helper Functions ============

/**
 * Check if any element in a boolean array is true.
 */
function anyTrue(arr: NDArray): boolean {
  for (let i = 0; i < arr.size; i++) {
    if (arr.getFlat(i)) {
      return true;
    }
  }
  return false;
}

/**
 * Check if all elements in a boolean array are true.
 */
export function allTrue(arr: NDArray): boolean {
  for (let i = 0; i < arr.size; i++) {
    if (!arr.getFlat(i)) {
      return false;
    }
  }
  return true;
}

/**
 * Check if two shapes are equal.
 */
function shapesEqual(s1: number[], s2: number[]): boolean {
  if (s1.length !== s2.length) {
    return false;
  }
  for (let i = 0; i < s1.length; i++) {
    if (s1[i] !== s2[i]) {
      return false;
    }
  }
  return true;
}

/**
 * Convert flat index to multi-dimensional indices.
 */
function flatToMultiIndex(flatIdx: number, shape: number[]): number[] {
  const indices = new Array(shape.length);
  let remaining = flatIdx;

  for (let i = shape.length - 1; i >= 0; i--) {
    indices[i] = remaining % shape[i];
    remaining = Math.floor(remaining / shape[i]);
  }

  return indices;
}

/**
 * Convert multi-dimensional indices to flat index.
 */
function multiToFlatIndex(indices: number[], shape: number[]): number {
  let flatIdx = 0;
  let multiplier = 1;

  for (let i = shape.length - 1; i >= 0; i--) {
    flatIdx += indices[i] * multiplier;
    multiplier *= shape[i];
  }

  return flatIdx;
}
