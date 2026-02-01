/**
 * Matrix Trace
 *
 * Compute the sum of diagonal elements.
 */

import {
  prepareMatrix,
  getMatrixDimensions,
} from '../helpers.js';
import type { Matrix, TraceResult } from './types.js';

/**
 * Compute the trace of a matrix (sum of diagonal elements).
 *
 * The trace is only defined for square matrices.
 *
 * @param A - Input square matrix (n × n)
 * @returns The trace (sum of diagonal elements)
 *
 * @example
 * ```typescript
 * const A = [[1, 2, 3], [4, 5, 6], [7, 8, 9]];
 * const { trace: tr } = trace(A);
 * // tr = 1 + 5 + 9 = 15
 *
 * const I = [[1, 0], [0, 1]];
 * const { trace: tr2 } = trace(I);
 * // tr2 = 2
 * ```
 */
export function trace(A: Matrix): TraceResult {
  // Get dimensions
  const [m, n] = getMatrixDimensions(A);

  // Check squareness
  if (m !== n) {
    throw new Error(`Matrix must be square, got ${m}×${n}`);
  }

  // Handle edge case
  if (n === 0) {
    return {
      trace: 0,
      n: 0,
      success: true,
      message: 'Success',
    };
  }

  // Prepare matrix (column-major)
  const aData = prepareMatrix(A);

  // Sum diagonal elements
  let sum = 0;
  for (let i = 0; i < n; i++) {
    sum += aData[i + i * n]; // A[i,i] in column-major
  }

  return {
    trace: sum,
    n,
    success: true,
    message: 'Success',
  };
}
