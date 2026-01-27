/**
 * NumJS Masked Arrays Module
 *
 * Provides masked array functionality compatible with NumPy's numpy.ma module.
 * Masked arrays support missing or invalid data by maintaining a boolean mask
 * alongside the data array.
 *
 * @example
 * ```typescript
 * import { ma } from 'numjs';
 *
 * // Create a masked array
 * const data = NDArray.fromArray([1, 2, 3, 4, 5]);
 * const arr = new ma.MaskedArray(data, [false, true, false, false, true]);
 *
 * // Operations exclude masked elements
 * arr.mean();  // (1 + 3 + 4) / 3 = 2.667
 * arr.filled(0);  // [1, 0, 3, 4, 0]
 * arr.compressed();  // [1, 3, 4]
 *
 * // Create masked arrays with convenience functions
 * const invalid = ma.masked_invalid([1, NaN, 3, Infinity, 5]);
 * const threshold = ma.masked_greater([1, 2, 3, 4, 5], 3);
 * ```
 */

// Types and constants
export {
  nomask,
  masked,
  MaskedArrayError,
  MaskError,
  defaultFillValues,
  getDefaultFillValue,
  maximumFillValue,
  minimumFillValue,
  isMaskedConstant,
} from './types.js';

export type { MaskType, MaskedConstant } from './types.js';

// Core class
export { MaskedArray } from './core.js';

// Mask operations
export {
  make_mask,
  make_mask_sync,
  make_mask_none,
  getmask,
  getmaskarray,
  getdata,
  is_mask,
  is_masked,
  mask_or,
  mask_or_sync,
  flatten_mask,
  reshape_mask,
  broadcast_mask,
  broadcast_mask_sync,
  allTrue,
} from './mask_ops.js';

// Fill value operations
export {
  default_fill_value,
  common_fill_value,
  set_fill_value,
  getReductionFillValue,
} from './fill_value.js';

// Creation functions
export {
  masked_array,
  array,
  masked_equal,
  masked_not_equal,
  masked_greater,
  masked_greater_equal,
  masked_less,
  masked_less_equal,
  masked_inside,
  masked_outside,
  masked_where,
  masked_invalid,
  masked_values,
  zeros,
  ones,
  empty,
  masked_all,
  masked_all_like,
  zeros_like,
  ones_like,
  empty_like,
  fromfunction,
} from './creation.js';

// Extras (statistics and utilities)
export {
  average,
  median,
  cov,
  corrcoef,
  notmasked_edges,
  notmasked_contiguous,
  flatnotmasked_edges,
  flatnotmasked_contiguous,
  clump_masked,
  clump_unmasked,
  apply_along_axis,
} from './extras.js';

export type { SliceInfo } from './extras.js';

// Re-import for namespace object
import {
  nomask,
  masked,
  MaskedArrayError,
  MaskError,
  defaultFillValues,
  getDefaultFillValue,
  maximumFillValue,
  minimumFillValue,
  isMaskedConstant,
} from './types.js';

import { MaskedArray } from './core.js';

import {
  make_mask,
  make_mask_sync,
  make_mask_none,
  getmask,
  getmaskarray,
  getdata,
  is_mask,
  is_masked,
  mask_or,
  mask_or_sync,
  flatten_mask,
  reshape_mask,
  broadcast_mask,
  broadcast_mask_sync,
  allTrue,
} from './mask_ops.js';

import {
  default_fill_value,
  common_fill_value,
  set_fill_value,
  getReductionFillValue,
} from './fill_value.js';

import {
  masked_array,
  array,
  masked_equal,
  masked_not_equal,
  masked_greater,
  masked_greater_equal,
  masked_less,
  masked_less_equal,
  masked_inside,
  masked_outside,
  masked_where,
  masked_invalid,
  masked_values,
  zeros,
  ones,
  empty,
  masked_all,
  masked_all_like,
  zeros_like,
  ones_like,
  empty_like,
  fromfunction,
} from './creation.js';

import {
  average,
  median,
  cov,
  corrcoef,
  notmasked_edges,
  notmasked_contiguous,
  flatnotmasked_edges,
  flatnotmasked_contiguous,
  clump_masked,
  clump_unmasked,
  apply_along_axis,
} from './extras.js';

/**
 * Masked arrays module namespace object for convenient grouped import.
 *
 * @example
 * ```typescript
 * import { ma } from 'numjs';
 *
 * const arr = ma.masked_invalid([1, NaN, 3]);
 * const mean = arr.mean();
 * ```
 */
export const ma = {
  // Constants
  nomask,
  masked,

  // Error classes
  MaskedArrayError,
  MaskError,

  // Type utilities
  defaultFillValues,
  getDefaultFillValue,
  maximumFillValue,
  minimumFillValue,
  isMaskedConstant,

  // Core class
  MaskedArray,

  // Mask operations
  make_mask,
  make_mask_sync,
  make_mask_none,
  getmask,
  getmaskarray,
  getdata,
  is_mask,
  is_masked,
  mask_or,
  mask_or_sync,
  flatten_mask,
  reshape_mask,
  broadcast_mask,
  broadcast_mask_sync,
  allTrue,

  // Fill value operations
  default_fill_value,
  common_fill_value,
  set_fill_value,
  getReductionFillValue,

  // Creation functions
  masked_array,
  array,
  masked_equal,
  masked_not_equal,
  masked_greater,
  masked_greater_equal,
  masked_less,
  masked_less_equal,
  masked_inside,
  masked_outside,
  masked_where,
  masked_invalid,
  masked_values,
  zeros,
  ones,
  empty,
  masked_all,
  masked_all_like,
  zeros_like,
  ones_like,
  empty_like,
  fromfunction,

  // Extras (statistics and utilities)
  average,
  median,
  cov,
  corrcoef,
  notmasked_edges,
  notmasked_contiguous,
  flatnotmasked_edges,
  flatnotmasked_contiguous,
  clump_masked,
  clump_unmasked,
  apply_along_axis,
};
