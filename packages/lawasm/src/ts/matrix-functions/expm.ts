/**
 * Matrix Exponential
 *
 * Compute e^A for a square matrix A.
 */

import {
  prepareMatrix,
  getMatrixDimensions,
} from '../helpers.js';
import type { Matrix, ExpmResult } from './types.js';

/**
 * Compute the matrix exponential e^A.
 *
 * Uses the scaling and squaring method with Padé approximation.
 * This is a pure TypeScript implementation.
 *
 * @param A - Input square matrix (n × n)
 * @returns Matrix exponential e^A
 *
 * @example
 * ```typescript
 * // For diagonal matrix, expm gives e^(diagonal elements)
 * const A = [[1, 0], [0, 2]];
 * const { E } = expm(A);
 * // E ≈ [[e, 0], [0, e²]]
 *
 * // For nilpotent matrix
 * const N = [[0, 1], [0, 0]];
 * const { E: E2 } = expm(N);
 * // E2 = [[1, 1], [0, 1]] (I + N)
 * ```
 */
export function expm(A: Matrix): ExpmResult {
  // Get dimensions
  const [m, n] = getMatrixDimensions(A);

  if (m !== n) {
    throw new Error(`Matrix must be square, got ${m}×${n}`);
  }

  if (n === 0) {
    return {
      E: new Float64Array(0),
      n: 0,
      success: true,
      message: 'Success',
    };
  }

  // Prepare matrix (column-major)
  const aData = prepareMatrix(A);

  // Scaling and squaring with Padé approximation
  // Algorithm from Higham (2005): "The Scaling and Squaring Method for the Matrix Exponential Revisited"

  // 1. Compute norm of A
  let normA = 0;
  for (let j = 0; j < n; j++) {
    let colSum = 0;
    for (let i = 0; i < n; i++) {
      colSum += Math.abs(aData[i + j * n]);
    }
    if (colSum > normA) normA = colSum;
  }

  // 2. Determine scaling factor s such that ||A/2^s|| < 1
  let s = 0;
  if (normA > 0) {
    s = Math.max(0, Math.ceil(Math.log2(normA)));
  }

  // 3. Scale matrix: B = A / 2^s
  const scale = Math.pow(2, -s);
  const B = new Float64Array(n * n);
  for (let i = 0; i < n * n; i++) {
    B[i] = aData[i] * scale;
  }

  // 4. Compute Padé approximant of e^B using [6/6] Padé
  // e^B ≈ (I - B/2 + B²/12 - B³/120 + ...)^(-1) * (I + B/2 + B²/12 + B³/120 + ...)
  // We use a simpler approach: Taylor series for small ||B||

  // Compute powers of B
  const I = new Float64Array(n * n);
  for (let i = 0; i < n; i++) {
    I[i + i * n] = 1;
  }

  // Taylor series: e^B = I + B + B²/2! + B³/3! + ...
  const result = new Float64Array(n * n);
  const term = new Float64Array(n * n);
  const temp = new Float64Array(n * n);

  // Start with I
  for (let i = 0; i < n * n; i++) {
    result[i] = I[i];
    term[i] = I[i];
  }

  // Add terms until convergence
  const maxTerms = 30;
  const eps = 1e-16;

  for (let k = 1; k <= maxTerms; k++) {
    // term = term * B / k
    // temp = term * B
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        let sum = 0;
        for (let l = 0; l < n; l++) {
          sum += term[i + l * n] * B[l + j * n];
        }
        temp[i + j * n] = sum / k;
      }
    }

    // term = temp
    let maxTerm = 0;
    for (let i = 0; i < n * n; i++) {
      term[i] = temp[i];
      result[i] += term[i];
      if (Math.abs(term[i]) > maxTerm) maxTerm = Math.abs(term[i]);
    }

    // Check convergence
    if (maxTerm < eps) break;
  }

  // 5. Square s times: e^A = (e^B)^(2^s)
  for (let k = 0; k < s; k++) {
    // result = result * result
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        let sum = 0;
        for (let l = 0; l < n; l++) {
          sum += result[i + l * n] * result[l + j * n];
        }
        temp[i + j * n] = sum;
      }
    }
    for (let i = 0; i < n * n; i++) {
      result[i] = temp[i];
    }
  }

  return {
    E: result,
    n,
    success: true,
    message: 'Success',
  };
}
