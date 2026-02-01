/**
 * General Matrix Inverse
 *
 * Compute the inverse of a general square matrix using LU factorization.
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
import type { Matrix, InvOptions, InvResult } from './types.js';

/**
 * Compute the inverse of a general n×n matrix A.
 *
 * Uses LU factorization with partial pivoting: A = P*L*U,
 * then computes A^(-1) from this factorization.
 *
 * @param A - Input square matrix (n × n)
 * @param options - Computation options
 * @returns Matrix inverse
 *
 * @example
 * ```typescript
 * const A = [[4, 7], [2, 6]];
 * const { inv, success } = inv(A);
 * // inv is A^(-1) such that A * inv = I
 * ```
 *
 * @throws Error if matrix is singular
 */
export function inv(A: Matrix, options: InvOptions = {}): InvResult {
  const Module = getLAPACKModule();

  const { overwriteA = false } = options;

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);
  if (m !== n) {
    throw new Error(`Matrix A must be square, got ${m}×${n}`);
  }

  // Prepare matrix
  const aData = prepareMatrix(A);

  // Allocate WASM memory
  const aPtr = allocateDoubles(Module, overwriteA ? null : aData, n * n);
  const ipivPtr = allocateInts(Module, null, n);
  const infoPtr = allocateInts(Module, [0]);
  const lworkPtr = allocateInts(Module, [-1]);
  const workQueryPtr = allocateDoubles(Module, null, 1);

  // Write data if not overwriting
  if (!overwriteA) {
    const baseIdx = aPtr >> 3;
    for (let i = 0; i < aData.length; i++) {
      Module.HEAPF64[baseIdx + i] = aData[i];
    }
  }

  // Parameter pointers
  const mPtr = allocateInts(Module, [n]);
  const nPtr = allocateInts(Module, [n]);
  const ldaPtr = allocateInts(Module, [n]);

  try {
    // Step 1: LU factorization using DGETRF
    Module._dgetrf_(mPtr, nPtr, aPtr, ldaPtr, ipivPtr, infoPtr);

    let info = readInt(Module, infoPtr);
    if (info !== 0) {
      if (info > 0) {
        return {
          inv: new Float64Array(0),
          n,
          info,
          success: false,
          message: `Matrix is singular: U(${info},${info}) is exactly zero`,
        };
      }
      return {
        inv: new Float64Array(0),
        n,
        info,
        success: false,
        message: getLapackErrorMessage('DGETRF', info),
      };
    }

    // Step 2: Compute inverse using DGETRI
    // Workspace query
    Module._dgetri_(nPtr, aPtr, ldaPtr, ipivPtr, workQueryPtr, lworkPtr, infoPtr);

    info = readInt(Module, infoPtr);
    if (info !== 0) {
      return {
        inv: new Float64Array(0),
        n,
        info,
        success: false,
        message: getLapackErrorMessage('DGETRI', info),
      };
    }

    const lwork = Math.ceil(Module.HEAPF64[workQueryPtr >> 3]);
    const workPtr = allocateDoubles(Module, null, lwork);
    Module.HEAP32[lworkPtr >> 2] = lwork;

    // Compute inverse
    Module._dgetri_(nPtr, aPtr, ldaPtr, ipivPtr, workPtr, lworkPtr, infoPtr);

    info = readInt(Module, infoPtr);
    freeAll(Module, [workPtr]);

    if (info !== 0) {
      return {
        inv: new Float64Array(0),
        n,
        info,
        success: false,
        message: getLapackErrorMessage('DGETRI', info),
      };
    }

    // Read result
    const invMatrix = readDoubles(Module, aPtr, n * n);

    return {
      inv: invMatrix,
      n,
      info: 0,
      success: true,
      message: 'Success',
    };
  } finally {
    freeAll(Module, [
      aPtr, ipivPtr, infoPtr, lworkPtr, workQueryPtr,
      mPtr, nPtr, ldaPtr
    ]);
  }
}
