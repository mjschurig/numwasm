/**
 * Cholesky Factorization
 *
 * Compute the Cholesky factorization of a symmetric positive definite matrix.
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
import type { Matrix, CholeskyOptions, CholeskyResult } from './types.js';

/**
 * Compute the Cholesky factorization of a symmetric positive definite matrix A.
 *
 * If upper=false (default): A = L * L^T where L is lower triangular
 * If upper=true: A = U^T * U where U is upper triangular
 *
 * @param A - Symmetric positive definite matrix (n × n).
 *            Only the specified triangle (upper or lower) is referenced.
 * @param options - Factorization options
 * @returns Cholesky factorization result
 *
 * @throws Error if the matrix is not positive definite (info > 0)
 *
 * @example
 * ```typescript
 * const A = [[4, 2], [2, 5]]; // Symmetric positive definite
 * const { factor, success } = cholesky(A);
 * // factor is L such that A = L * L^T
 * ```
 *
 * @example
 * ```typescript
 * // Compute upper triangular factor
 * const { factor } = cholesky(A, { upper: true });
 * // factor is U such that A = U^T * U
 * ```
 */
export function cholesky(A: Matrix, options: CholeskyOptions = {}): CholeskyResult {
  const Module = getLAPACKModule();

  const { upper = false, overwriteA = false } = options;

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);
  if (m !== n) {
    throw new Error(`Matrix A must be square, got ${m}×${n}`);
  }

  // Prepare matrix
  const aData = prepareMatrix(A);

  // Allocate WASM memory
  const uplo = upper ? CHAR.U : CHAR.L;
  const uploPtr = allocateInts(Module, [uplo]);
  const nPtr = allocateInts(Module, [n]);
  const aPtr = allocateDoubles(Module, overwriteA ? null : aData, n * n);
  const ldaPtr = allocateInts(Module, [n]);
  const infoPtr = allocateInts(Module, [0]);

  // Write data if not overwriting
  if (!overwriteA) {
    const baseIdx = aPtr >> 3;
    for (let i = 0; i < aData.length; i++) {
      Module.HEAPF64[baseIdx + i] = aData[i];
    }
  }

  try {
    // Call DPOTRF
    Module._dpotrf_(uploPtr, nPtr, aPtr, ldaPtr, infoPtr);

    const info = readInt(Module, infoPtr);
    const success = info === 0;
    const message = info === 0
      ? 'Success'
      : info > 0
        ? `The leading minor of order ${info} is not positive definite, and the factorization could not be completed.`
        : getLapackErrorMessage('DPOTRF', info);

    // Read the result
    const rawFactor = readDoubles(Module, aPtr, n * n);

    // Zero out the non-factor part
    const factor = new Float64Array(n * n);
    for (let j = 0; j < n; j++) {
      for (let i = 0; i < n; i++) {
        if (upper) {
          // Upper triangular: copy on and above diagonal
          if (i <= j) {
            factor[j * n + i] = rawFactor[j * n + i];
          }
        } else {
          // Lower triangular: copy on and below diagonal
          if (i >= j) {
            factor[j * n + i] = rawFactor[j * n + i];
          }
        }
      }
    }

    return {
      factor,
      upper,
      n,
      info,
      success,
      message,
    };
  } finally {
    freeAll(Module, [uploPtr, nPtr, aPtr, ldaPtr, infoPtr]);
  }
}
