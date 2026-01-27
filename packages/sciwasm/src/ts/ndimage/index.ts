/**
 * N-dimensional image processing.
 * @module ndimage
 */

import { NotImplementedError } from '../errors.js';

/**
 * Multi-dimensional convolution.
 * Mirrors scipy.ndimage.convolve.
 */
export function convolve(
  _input: number[][],
  _weights: number[][],
  _options?: { mode?: 'reflect' | 'constant' | 'nearest' | 'wrap'; cval?: number }
): number[][] {
  throw new NotImplementedError('sciwasm.ndimage.convolve');
}

/**
 * Multi-dimensional Gaussian filter.
 * Mirrors scipy.ndimage.gaussian_filter.
 */
export function gaussian_filter(
  _input: number[][],
  _sigma: number,
  _options?: { mode?: 'reflect' | 'constant' | 'nearest' | 'wrap' }
): number[][] {
  throw new NotImplementedError('sciwasm.ndimage.gaussian_filter');
}

/**
 * Label features in an array.
 * Mirrors scipy.ndimage.label.
 */
export function label(
  _input: number[][],
  _structure?: number[][]
): { labeled: number[][]; num_features: number } {
  throw new NotImplementedError('sciwasm.ndimage.label');
}

/**
 * Binary erosion on a binary array.
 * Mirrors scipy.ndimage.binary_erosion.
 */
export function binary_erosion(
  _input: boolean[][],
  _structure?: boolean[][],
  _iterations?: number
): boolean[][] {
  throw new NotImplementedError('sciwasm.ndimage.binary_erosion');
}

/**
 * Binary dilation on a binary array.
 * Mirrors scipy.ndimage.binary_dilation.
 */
export function binary_dilation(
  _input: boolean[][],
  _structure?: boolean[][],
  _iterations?: number
): boolean[][] {
  throw new NotImplementedError('sciwasm.ndimage.binary_dilation');
}
