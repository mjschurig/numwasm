/**
 * Hermitian Rank-k Update (Pure TypeScript)
 *
 * Compute C = alpha*A*A^H + beta*C or C = alpha*A^H*A + beta*C
 *
 * Note: ZHERK is for complex matrices. This implementation provides
 * a real-valued approximation for symmetric positive semi-definite updates.
 */

import {
  prepareMatrix,
  getMatrixDimensions,
} from '../helpers.js';
import type { Matrix, HerkOptions, HerkResult } from './types.js';

/**
 * Compute Hermitian rank-k update C = alpha*A*A^H + beta*C or C = alpha*A^H*A + beta*C.
 *
 * For real matrices, this is equivalent to symmetric rank-k update (SYRK).
 * The result C is Hermitian (symmetric for real matrices).
 *
 * @param A - Input matrix
 * @param options - Update options
 * @returns Hermitian result matrix C
 *
 * @example
 * ```typescript
 * const A = [[1, 2, 3], [4, 5, 6]];
 *
 * // C = A*A^H (2×2 result since A is 2×3)
 * const { C } = herk(A);
 *
 * // C = A^H*A (3×3 result)
 * const { C: C2 } = herk(A, { trans: 'C' });
 * ```
 */
export function herk(A: Matrix, options: HerkOptions = {}): HerkResult {
  const {
    uplo = 'lower',
    trans = 'N',
    alpha = 1.0,
    beta = 0.0,
    C,
  } = options;

  // Get dimensions
  const [rowsA, colsA] = getMatrixDimensions(A);

  // Determine n (dimension of result) and k (common dimension)
  const n = trans === 'N' ? rowsA : colsA;
  const k = trans === 'N' ? colsA : rowsA;

  // Prepare A matrix (column-major)
  const aData = prepareMatrix(A);

  // Prepare or create C
  let cData: Float64Array;
  if (C !== undefined) {
    const [rowsC, colsC] = getMatrixDimensions(C);
    if (rowsC !== n || colsC !== n) {
      throw new Error(
        `C dimensions ${rowsC}×${colsC} incompatible with result ${n}×${n}`
      );
    }
    cData = prepareMatrix(C);
  } else {
    cData = new Float64Array(n * n);
  }

  // For real matrices, Hermitian rank-k update is same as symmetric
  // C = alpha * A * A^T + beta * C  (when trans = 'N')
  // C = alpha * A^T * A + beta * C  (when trans = 'C' or 'T')

  const result = new Float64Array(n * n);

  // Apply beta to existing C
  if (beta !== 0) {
    for (let i = 0; i < n * n; i++) {
      result[i] = beta * cData[i];
    }
  }

  // Compute alpha * A * A^T or alpha * A^T * A
  if (trans === 'N') {
    // C = alpha * A * A^T
    // A is rowsA × colsA, result is rowsA × rowsA
    for (let i = 0; i < n; i++) {
      for (let j = 0; j <= i; j++) {
        let sum = 0;
        for (let l = 0; l < k; l++) {
          // A(i,l) * A(j,l) in column-major: aData[i + l*rowsA] * aData[j + l*rowsA]
          sum += aData[i + l * rowsA] * aData[j + l * rowsA];
        }
        result[i + j * n] += alpha * sum;
        if (i !== j) {
          result[j + i * n] += alpha * sum;
        }
      }
    }
  } else {
    // C = alpha * A^T * A (or A^H * A)
    // A is rowsA × colsA, result is colsA × colsA
    for (let i = 0; i < n; i++) {
      for (let j = 0; j <= i; j++) {
        let sum = 0;
        for (let l = 0; l < k; l++) {
          // A(l,i) * A(l,j) in column-major: aData[l + i*rowsA] * aData[l + j*rowsA]
          sum += aData[l + i * rowsA] * aData[l + j * rowsA];
        }
        result[i + j * n] += alpha * sum;
        if (i !== j) {
          result[j + i * n] += alpha * sum;
        }
      }
    }
  }

  // If uplo is 'upper', we already have the full matrix, but the convention
  // is to have the specified triangle filled. Since we computed the full
  // symmetric matrix, we're fine.

  // Ensure symmetry based on uplo preference (optional cleanup)
  if (uplo === 'upper') {
    for (let j = 0; j < n; j++) {
      for (let i = j + 1; i < n; i++) {
        result[j + i * n] = result[i + j * n];
      }
    }
  } else {
    for (let j = 0; j < n; j++) {
      for (let i = j + 1; i < n; i++) {
        result[i + j * n] = result[j + i * n];
      }
    }
  }

  return {
    C: result,
    n,
    success: true,
    message: 'Success',
  };
}
