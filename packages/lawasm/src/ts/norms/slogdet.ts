/**
 * Sign and Log-Determinant
 *
 * Compute the sign and natural logarithm of the absolute value of the determinant.
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
import type { Matrix, SlogDetOptions, SlogDetResult } from './types.js';

/**
 * Compute the sign and natural logarithm of the absolute value of the determinant.
 *
 * Returns (sign, logabsdet) such that det(A) = sign * exp(logabsdet).
 *
 * This is the most numerically stable way to work with determinants,
 * especially for large matrices or matrices with very large/small determinants.
 *
 * @param A - Input square matrix (n × n)
 * @param options - Computation options
 * @returns Sign and log of absolute determinant
 *
 * @example
 * ```typescript
 * const A = [[1, 2], [3, 4]];
 * const { sign, logabsdet } = slogdet(A);
 * // sign = -1, logabsdet = log(2) ≈ 0.693
 * // det(A) = -1 * exp(0.693) = -2
 *
 * const B = [[1, 0], [0, 2]];
 * const { sign: s2, logabsdet: l2 } = slogdet(B);
 * // s2 = 1, l2 = log(2) ≈ 0.693
 * // det(B) = 1 * exp(0.693) = 2
 * ```
 */
export function slogdet(A: Matrix, options: SlogDetOptions = {}): SlogDetResult {
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
      sign: 1,
      logabsdet: 0, // log(1) = 0
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
        sign: 0,
        logabsdet: NaN,
        n,
        info,
        success: false,
        message: getLapackErrorMessage('DGETRF', info),
      };
    }

    if (info > 0) {
      // Matrix is singular (U has a zero on diagonal)
      return {
        sign: 0,
        logabsdet: -Infinity,
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

    // Compute sign and log|det|
    let logAbsDet = 0.0;
    let sign = 1;

    for (let i = 0; i < n; i++) {
      const diag = luData[i * n + i]; // Column-major: element (i,i)

      if (diag === 0) {
        return {
          sign: 0,
          logabsdet: -Infinity,
          n,
          info: 0,
          success: true,
          message: 'Matrix is singular',
        };
      }

      // Accumulate log of absolute value
      logAbsDet += Math.log(Math.abs(diag));

      // Track sign from diagonal elements
      if (diag < 0) {
        sign = -sign;
      }

      // Track sign from permutation - LAPACK uses 1-based indexing
      if (ipiv[i] !== i + 1) {
        sign = -sign;
      }
    }

    return {
      sign,
      logabsdet: logAbsDet,
      n,
      info: 0,
      success: true,
      message: 'Success',
    };
  } finally {
    freeAll(Module, [aPtr, ipivPtr, infoPtr, mPtr, nPtr, ldaPtr]);
  }
}
