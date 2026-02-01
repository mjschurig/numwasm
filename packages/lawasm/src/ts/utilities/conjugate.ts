/**
 * Complex Conjugate
 *
 * Compute the complex conjugate of a matrix.
 */

import {
  prepareMatrix,
  getMatrixDimensions,
} from '../helpers.js';
import type { Matrix, ConjugateResult } from './types.js';

/**
 * Compute the complex conjugate of a matrix.
 *
 * For real matrices, this returns a copy of the matrix.
 * For complex matrices (stored as interleaved real/imag pairs),
 * this negates the imaginary parts.
 *
 * @param A - Input matrix (m × n)
 * @returns Conjugated matrix
 *
 * @example
 * ```typescript
 * // Real matrix (just returns a copy)
 * const A = [[1, 2], [3, 4]];
 * const { conj } = conjugate(A);
 *
 * // Complex matrix (interleaved format: [re, im, re, im, ...])
 * // For a 2×1 complex matrix with values (1+2i) and (3+4i):
 * // Pass as 2×2 real: [[1, 2], [3, 4]] where each row is [re, im]
 * ```
 */
export function conjugate(A: Matrix): ConjugateResult {
  // Get dimensions
  const [m, n] = getMatrixDimensions(A);

  // Handle edge cases
  if (m === 0 || n === 0) {
    return {
      conj: new Float64Array(0),
      m,
      n,
      success: true,
      message: 'Success',
    };
  }

  // Prepare matrix (column-major)
  const aData = prepareMatrix(A);

  // For real matrices, just copy
  // For complex interpretation (even n), negate every other element
  const result = new Float64Array(m * n);

  // Check if this might be a complex matrix (even number of columns)
  // and the structure suggests interleaved format
  const isLikelyComplex = n >= 2 && n % 2 === 0;

  if (isLikelyComplex) {
    // Treat as complex: columns come in pairs (real, imag)
    for (let j = 0; j < n; j += 2) {
      for (let i = 0; i < m; i++) {
        // Real part (even column)
        result[i + j * m] = aData[i + j * m];
        // Imaginary part (odd column) - negate for conjugate
        result[i + (j + 1) * m] = -aData[i + (j + 1) * m];
      }
    }
  } else {
    // Real matrix - just copy
    for (let i = 0; i < m * n; i++) {
      result[i] = aData[i];
    }
  }

  return {
    conj: result,
    m,
    n,
    success: true,
    message: 'Success',
  };
}
