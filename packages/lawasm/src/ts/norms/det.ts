/**
 * Matrix Determinant
 *
 * Compute the determinant of a square matrix using LU factorization.
 */

import { getLAPACKModule } from '../loader.js';
import {
  allocateDoubles,
  allocateInts,
  readDoubles,
  readInt,
  freeAll,
  prepareMatrix,
  getMatrixDimensions,
  getLapackErrorMessage,
} from '../helpers.js';
import type { Matrix, DetOptions, DetResult } from './types.js';

/**
 * Compute the determinant of a square matrix.
 *
 * Uses LU factorization: A = P*L*U, then det(A) = det(P) * det(L) * det(U)
 * Since L has unit diagonal, det(L) = 1.
 * det(U) = product of diagonal elements.
 * det(P) = (-1)^(number of row interchanges).
 *
 * WARNING: For large matrices, the determinant can overflow or underflow.
 * Use logdet() or slogdet() for numerically stable computation.
 *
 * @param A - Input square matrix (n × n)
 * @param options - Computation options
 * @returns The determinant
 *
 * @example
 * ```typescript
 * const A = [[1, 2], [3, 4]];
 * const { det: d } = det(A);
 * // d = -2
 *
 * const B = [[1, 0, 0], [0, 2, 0], [0, 0, 3]];
 * const { det: d2 } = det(B);
 * // d2 = 6
 * ```
 */
export function det(A: Matrix, options: DetOptions = {}): DetResult {
  const Module = getLAPACKModule();

  const { overwriteA: _overwriteA = false } = options;

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);
  if (m !== n) {
    throw new Error(`Matrix A must be square, got ${m}×${n}`);
  }

  // Handle edge case
  if (n === 0) {
    return {
      det: 1, // Empty matrix has determinant 1 by convention
      n: 0,
      info: 0,
      success: true,
      message: 'Success',
    };
  }

  // Prepare matrix
  const aData = prepareMatrix(A);

  // Allocate WASM memory
  const aPtr = allocateDoubles(Module, aData, n * n);
  const ipivPtr = allocateInts(Module, null, n);
  const infoPtr = allocateInts(Module, [0]);
  const mPtr = allocateInts(Module, [n]);
  const nPtr = allocateInts(Module, [n]);
  const ldaPtr = allocateInts(Module, [n]);

  try {
    // LU factorization using DGETRF
    Module._dgetrf_(mPtr, nPtr, aPtr, ldaPtr, ipivPtr, infoPtr);

    const info = readInt(Module, infoPtr);
    if (info < 0) {
      return {
        det: NaN,
        n,
        info,
        success: false,
        message: getLapackErrorMessage('DGETRF', info),
      };
    }

    if (info > 0) {
      // Matrix is singular (U has a zero on diagonal)
      return {
        det: 0,
        n,
        info,
        success: true,
        message: 'Matrix is singular',
      };
    }

    // Read the LU factorization and pivot indices
    const luData = readDoubles(Module, aPtr, n * n);
    const ipiv = new Int32Array(n);
    const baseIdx = ipivPtr >> 2;
    for (let i = 0; i < n; i++) {
      ipiv[i] = Module.HEAP32[baseIdx + i];
    }

    // Compute determinant as product of diagonal elements of U
    // multiplied by sign from permutation
    let detValue = 1.0;
    let sign = 1;

    for (let i = 0; i < n; i++) {
      // Diagonal of U (stored in LU)
      detValue *= luData[i * n + i]; // Column-major: element (i,i)

      // Check pivot - LAPACK uses 1-based indexing
      if (ipiv[i] !== i + 1) {
        sign = -sign;
      }
    }

    detValue *= sign;

    return {
      det: detValue,
      n,
      info: 0,
      success: true,
      message: 'Success',
    };
  } finally {
    freeAll(Module, [aPtr, ipivPtr, infoPtr, mPtr, nPtr, ldaPtr]);
  }
}
