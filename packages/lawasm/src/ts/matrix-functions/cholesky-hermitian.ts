/**
 * Cholesky Factorization for Hermitian Matrices
 *
 * Compute the Cholesky factorization of a Hermitian positive definite matrix.
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
import type { Matrix, CholeskyHermitianResult } from './types.js';

/**
 * Options for Hermitian Cholesky factorization.
 */
export interface CholeskyHermitianOptions {
  /**
   * Which triangle to compute/use.
   * @default 'lower'
   */
  uplo?: 'upper' | 'lower';
}

/**
 * Compute Cholesky factorization of a Hermitian positive definite matrix.
 *
 * For real matrices, this is equivalent to standard Cholesky factorization.
 * For complex Hermitian matrices, A = L*L^H or A = U^H*U.
 *
 * Uses DPOTRF from LAPACK (for real matrices).
 *
 * @param A - Input Hermitian positive definite matrix (n × n)
 * @param options - Factorization options
 * @returns Cholesky factor L or U
 *
 * @example
 * ```typescript
 * const A = [[4, 2], [2, 5]]; // Symmetric positive definite
 * const { L } = choleskyHermitian(A);
 * // L * L' = A
 * ```
 */
export function choleskyHermitian(
  A: Matrix,
  options: CholeskyHermitianOptions = {}
): CholeskyHermitianResult {
  const Module = getLAPACKModule();

  const { uplo = 'lower' } = options;

  const [m, n] = getMatrixDimensions(A);

  if (m !== n) {
    throw new Error(`Matrix must be square, got ${m}×${n}`);
  }

  if (n === 0) {
    return {
      L: new Float64Array(0),
      n: 0,
      uplo,
      info: 0,
      success: true,
      message: 'Success',
    };
  }

  // Prepare matrix
  const aData = prepareMatrix(A);

  // Allocate WASM memory
  const aPtr = allocateDoubles(Module, aData, n * n);
  const infoPtr = allocateInts(Module, [0]);

  // Parameter pointers
  const uploPtr = allocateInts(Module, [uplo === 'upper' ? CHAR.U : CHAR.L]);
  const nPtr = allocateInts(Module, [n]);
  const ldaPtr = allocateInts(Module, [n]);

  try {
    // Compute Cholesky factorization using DPOTRF
    Module._dpotrf_(uploPtr, nPtr, aPtr, ldaPtr, infoPtr);

    const info = readInt(Module, infoPtr);

    if (info > 0) {
      return {
        L: new Float64Array(0),
        n,
        uplo,
        info,
        success: false,
        message: `Matrix is not positive definite: leading minor of order ${info} is not positive`,
      };
    }

    if (info < 0) {
      return {
        L: new Float64Array(0),
        n,
        uplo,
        info,
        success: false,
        message: getLapackErrorMessage('DPOTRF', info),
      };
    }

    // Read result
    const lData = readDoubles(Module, aPtr, n * n);

    // Zero out the non-factor part
    if (uplo === 'lower') {
      // Zero upper triangle
      for (let j = 0; j < n; j++) {
        for (let i = 0; i < j; i++) {
          lData[i + j * n] = 0;
        }
      }
    } else {
      // Zero lower triangle
      for (let j = 0; j < n; j++) {
        for (let i = j + 1; i < n; i++) {
          lData[i + j * n] = 0;
        }
      }
    }

    return {
      L: lData,
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
