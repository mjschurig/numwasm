/**
 * Symmetric Positive Definite Matrix Inverse
 *
 * Compute the inverse of a symmetric positive definite matrix using Cholesky factorization.
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
  CHAR,
} from '../helpers.js';
import type { Matrix, InvSymmetricOptions, InvSymmetricResult } from './types.js';

/**
 * Compute the inverse of a symmetric positive definite n×n matrix A.
 *
 * Uses Cholesky factorization: A = L*L^T (or U^T*U),
 * then computes A^(-1) from this factorization.
 *
 * This is more efficient than general matrix inversion when A is known
 * to be symmetric positive definite.
 *
 * @param A - Input symmetric positive definite matrix (n × n)
 * @param options - Computation options
 * @returns Matrix inverse
 *
 * @example
 * ```typescript
 * const A = [[4, 2], [2, 5]];
 * const { inv } = invSymmetric(A);
 * // inv is A^(-1), also symmetric positive definite
 * ```
 */
export function invSymmetric(A: Matrix, options: InvSymmetricOptions = {}): InvSymmetricResult {
  const Module = getLAPACKModule();

  const { uplo = 'lower', overwriteA = false } = options;

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);
  if (m !== n) {
    throw new Error(`Matrix A must be square, got ${m}×${n}`);
  }

  // Prepare matrix
  const aData = prepareMatrix(A);

  // Allocate WASM memory
  const aPtr = allocateDoubles(Module, overwriteA ? null : aData, n * n);
  const infoPtr = allocateInts(Module, [0]);

  // Write data if not overwriting
  if (!overwriteA) {
    const baseIdx = aPtr >> 3;
    for (let i = 0; i < aData.length; i++) {
      Module.HEAPF64[baseIdx + i] = aData[i];
    }
  }

  // Parameter pointers
  const uploPtr = allocateInts(Module, [uplo === 'upper' ? CHAR.U : CHAR.L]);
  const nPtr = allocateInts(Module, [n]);
  const ldaPtr = allocateInts(Module, [n]);

  try {
    // Step 1: Cholesky factorization using DPOTRF
    Module._dpotrf_(uploPtr, nPtr, aPtr, ldaPtr, infoPtr);

    let info = readInt(Module, infoPtr);
    if (info !== 0) {
      if (info > 0) {
        return {
          inv: new Float64Array(0),
          n,
          uplo,
          info,
          success: false,
          message: `Matrix is not positive definite: leading minor of order ${info} is not positive`,
        };
      }
      return {
        inv: new Float64Array(0),
        n,
        uplo,
        info,
        success: false,
        message: getLapackErrorMessage('DPOTRF', info),
      };
    }

    // Step 2: Compute inverse using DPOTRI
    Module._dpotri_(uploPtr, nPtr, aPtr, ldaPtr, infoPtr);

    info = readInt(Module, infoPtr);
    if (info !== 0) {
      return {
        inv: new Float64Array(0),
        n,
        uplo,
        info,
        success: false,
        message: getLapackErrorMessage('DPOTRI', info),
      };
    }

    // Read result
    const invMatrix = readDoubles(Module, aPtr, n * n);

    // Symmetrize the result (DPOTRI only fills one triangle)
    for (let j = 0; j < n; j++) {
      for (let i = j + 1; i < n; i++) {
        if (uplo === 'lower') {
          // Lower triangle is filled, copy to upper
          invMatrix[i * n + j] = invMatrix[j * n + i];
        } else {
          // Upper triangle is filled, copy to lower
          invMatrix[j * n + i] = invMatrix[i * n + j];
        }
      }
    }

    return {
      inv: invMatrix,
      n,
      uplo,
      info: 0,
      success: true,
      message: 'Success',
    };
  } finally {
    freeAll(Module, [aPtr, infoPtr, uploPtr, nPtr, ldaPtr]);
  }
}
