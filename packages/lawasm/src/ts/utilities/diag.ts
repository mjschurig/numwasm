/**
 * Diagonal Operations
 *
 * Extract diagonal from a matrix or create a diagonal matrix.
 */

import {
  prepareMatrix,
  prepareVector,
  getMatrixDimensions,
} from '../helpers.js';
import type { Matrix, DiagResult } from './types.js';
import type { RealArray } from '../helpers.js';

/**
 * Vector type for diag function.
 */
type Vector = number[] | RealArray;

/**
 * Extract the diagonal from a matrix, or create a diagonal matrix from a vector.
 *
 * If input is a matrix (2D), extracts the diagonal elements.
 * If input is a vector (1D), creates a diagonal matrix.
 *
 * @param A - Input matrix or vector
 * @param k - Diagonal offset (0 = main, positive = above, negative = below)
 * @returns Diagonal elements or diagonal matrix
 *
 * @example
 * ```typescript
 * // Extract diagonal from matrix
 * const A = [[1, 2, 3], [4, 5, 6], [7, 8, 9]];
 * const { diag: d } = diag(A);
 * // d = [1, 5, 9]
 *
 * // Extract superdiagonal
 * const { diag: d1 } = diag(A, 1);
 * // d1 = [2, 6]
 *
 * // Create diagonal matrix from vector
 * const v = [1, 2, 3];
 * const { diag: D } = diag(v);
 * // D is 3×3 diagonal matrix with [1, 2, 3] on diagonal
 * ```
 */
export function diag(A: Matrix | Vector, k: number = 0): DiagResult {
  // Check if input is a vector (1D) or matrix (2D)
  const isVector =
    ArrayBuffer.isView(A) ||
    (Array.isArray(A) && (A.length === 0 || typeof A[0] === 'number'));

  if (isVector) {
    // Create diagonal matrix from vector
    const vData = prepareVector(A as Vector);
    const n = vData.length;

    if (n === 0) {
      return {
        diag: new Float64Array(0),
        n: 0,
        success: true,
        message: 'Success',
      };
    }

    // Size of resulting matrix depends on k
    const size = n + Math.abs(k);
    const result = new Float64Array(size * size);

    // Place vector on k-th diagonal
    for (let i = 0; i < n; i++) {
      let row: number, col: number;
      if (k >= 0) {
        row = i;
        col = i + k;
      } else {
        row = i - k;
        col = i;
      }
      result[row + col * size] = vData[i];
    }

    return {
      diag: result,
      n: size,
      success: true,
      message: 'Success',
    };
  } else {
    // Extract diagonal from matrix
    const [m, n] = getMatrixDimensions(A as Matrix);

    if (m === 0 || n === 0) {
      return {
        diag: new Float64Array(0),
        n: 0,
        success: true,
        message: 'Success',
      };
    }

    const aData = prepareMatrix(A as Matrix);

    // Determine length of diagonal
    let diagLen: number;
    if (k >= 0) {
      diagLen = Math.min(m, n - k);
    } else {
      diagLen = Math.min(m + k, n);
    }

    if (diagLen <= 0) {
      return {
        diag: new Float64Array(0),
        n: 0,
        success: true,
        message: 'Diagonal outside matrix bounds',
      };
    }

    const result = new Float64Array(diagLen);

    // Extract k-th diagonal
    for (let i = 0; i < diagLen; i++) {
      let row: number, col: number;
      if (k >= 0) {
        row = i;
        col = i + k;
      } else {
        row = i - k;
        col = i;
      }
      result[i] = aData[row + col * m];
    }

    return {
      diag: result,
      n: diagLen,
      success: true,
      message: 'Success',
    };
  }
}
