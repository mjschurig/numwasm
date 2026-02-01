/**
 * Hermitian Rank-1 Update (Pure TypeScript)
 *
 * Compute A = alpha*x*x^H + A where A is Hermitian
 *
 * Note: ZHER is for complex matrices. This implementation provides
 * a real-valued version (equivalent to symmetric rank-1 update).
 */

import {
  prepareMatrix,
  prepareVector,
  getMatrixDimensions,
} from '../helpers.js';
import type { Vector, HerOptions, HerResult } from './types.js';

/**
 * Compute Hermitian rank-1 update A = alpha*x*x^H + A.
 *
 * For real matrices, this is equivalent to symmetric rank-1 update.
 * The result A is Hermitian (symmetric for real matrices).
 *
 * @param x - Vector (length n)
 * @param options - Update options
 * @returns Hermitian result matrix A (n × n)
 *
 * @example
 * ```typescript
 * const x = [1, 2, 3];
 *
 * // A = x*x^H (3×3 Hermitian result)
 * const { A } = her(x);
 *
 * // A = 2*x*x^H + A0
 * const A0 = [[1, 0, 0], [0, 1, 0], [0, 0, 1]];
 * const { A: A2 } = her(x, { alpha: 2, A: A0 });
 * ```
 */
export function her(x: Vector, options: HerOptions = {}): HerResult {
  const {
    uplo = 'lower',
    alpha = 1.0,
    A,
  } = options;

  // Get vector length
  const xData = prepareVector(x);
  const n = xData.length;

  // Prepare or create A
  let aData: Float64Array;
  if (A !== undefined) {
    const [rowsA, colsA] = getMatrixDimensions(A);
    if (rowsA !== n || colsA !== n) {
      throw new Error(
        `A dimensions ${rowsA}×${colsA} incompatible with result ${n}×${n}`
      );
    }
    aData = prepareMatrix(A);
  } else {
    aData = new Float64Array(n * n);
  }

  // For real matrices, Hermitian rank-1 update is same as symmetric
  // A = alpha * x * x^T + A

  const result = new Float64Array(n * n);

  // Copy existing A
  for (let i = 0; i < n * n; i++) {
    result[i] = aData[i];
  }

  // Add alpha * x * x^T
  for (let i = 0; i < n; i++) {
    for (let j = 0; j <= i; j++) {
      const update = alpha * xData[i] * xData[j];
      result[i + j * n] += update;
      if (i !== j) {
        result[j + i * n] += update;
      }
    }
  }

  // Ensure symmetry based on uplo preference
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
    A: result,
    n,
    success: true,
    message: 'Success',
  };
}
