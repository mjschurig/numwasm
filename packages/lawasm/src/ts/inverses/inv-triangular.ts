/**
 * Triangular Matrix Inverse
 *
 * Compute the inverse of a triangular matrix.
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
import type { Matrix, InvTriangularOptions, InvTriangularResult } from './types.js';

/**
 * Compute the inverse of a triangular n×n matrix A.
 *
 * For upper triangular A, the inverse is also upper triangular.
 * For lower triangular A, the inverse is also lower triangular.
 *
 * @param A - Input triangular matrix (n × n)
 * @param options - Computation options
 * @returns Matrix inverse
 *
 * @example
 * ```typescript
 * // Upper triangular
 * const A = [[2, 3], [0, 4]];
 * const { inv } = invTriangular(A, { upper: true });
 * ```
 *
 * @example
 * ```typescript
 * // Lower triangular
 * const A = [[2, 0], [3, 4]];
 * const { inv } = invTriangular(A, { upper: false });
 * ```
 */
export function invTriangular(A: Matrix, options: InvTriangularOptions = {}): InvTriangularResult {
  const Module = getLAPACKModule();

  const {
    upper = true,
    unitDiagonal = false,
    overwriteA = false,
  } = options;

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
  const uploPtr = allocateInts(Module, [upper ? CHAR.U : CHAR.L]);
  const diagPtr = allocateInts(Module, [unitDiagonal ? CHAR.UNIT : CHAR.N]);
  const nPtr = allocateInts(Module, [n]);
  const ldaPtr = allocateInts(Module, [n]);

  try {
    // Compute inverse using DTRTRI
    Module._dtrtri_(uploPtr, diagPtr, nPtr, aPtr, ldaPtr, infoPtr);

    const info = readInt(Module, infoPtr);
    if (info !== 0) {
      if (info > 0) {
        return {
          inv: new Float64Array(0),
          n,
          upper,
          info,
          success: false,
          message: `Matrix is singular: diagonal element ${info} is zero`,
        };
      }
      return {
        inv: new Float64Array(0),
        n,
        upper,
        info,
        success: false,
        message: getLapackErrorMessage('DTRTRI', info),
      };
    }

    // Read result
    const invMatrix = readDoubles(Module, aPtr, n * n);

    return {
      inv: invMatrix,
      n,
      upper,
      info: 0,
      success: true,
      message: 'Success',
    };
  } finally {
    freeAll(Module, [aPtr, infoPtr, uploPtr, diagPtr, nPtr, ldaPtr]);
  }
}
