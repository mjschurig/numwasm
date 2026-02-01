/**
 * Matrix-Vector Multiplication (DGEMV)
 *
 * Compute y = alpha*op(A)*x + beta*y
 */

import { getLAPACKModule } from '../loader.js';
import {
  allocateDoubles,
  allocateInts,
  readDoubles,
  freeAll,
  prepareMatrix,
  prepareVector,
  getMatrixDimensions,
  CHAR,
} from '../helpers.js';
import type { Matrix, Vector, MatvecOptions, MatvecResult, TransposeOp } from './types.js';

/**
 * Get CHAR code for transpose operation.
 */
function getTransChar(trans: TransposeOp): number {
  switch (trans) {
    case 'T':
      return CHAR.T;
    case 'C':
      return CHAR.C;
    default:
      return CHAR.N;
  }
}

/**
 * Compute matrix-vector multiplication y = alpha*op(A)*x + beta*y.
 *
 * Uses DGEMV from BLAS.
 *
 * @param A - Input matrix (m × n)
 * @param x - Input vector
 * @param options - Multiplication options
 * @returns Result vector y
 *
 * @example
 * ```typescript
 * const A = [[1, 2, 3], [4, 5, 6]];
 * const x = [1, 2, 3];
 *
 * // y = A*x
 * const { y } = matvec(A, x);
 *
 * // y = A^T*x (need x to have length 2 now)
 * const x2 = [1, 2];
 * const { y: y2 } = matvec(A, x2, { trans: 'T' });
 *
 * // y = 2*A*x + 3*y0
 * const y0 = [1, 1];
 * const { y: y3 } = matvec(A, x, { alpha: 2, beta: 3, y: y0 });
 * ```
 */
export function matvec(A: Matrix, x: Vector, options: MatvecOptions = {}): MatvecResult {
  const Module = getLAPACKModule();

  const {
    trans = 'N',
    alpha = 1.0,
    beta = 0.0,
    y,
  } = options;

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);
  const xData = prepareVector(x);

  // Determine expected dimensions
  const xLen = trans === 'N' ? n : m;
  const yLen = trans === 'N' ? m : n;

  // Check x dimension
  if (xData.length !== xLen) {
    throw new Error(
      `Vector x length ${xData.length} incompatible with op(A) requiring length ${xLen}`
    );
  }

  // Prepare or create y
  let yData: Float64Array;
  if (y !== undefined) {
    yData = prepareVector(y);
    if (yData.length !== yLen) {
      throw new Error(
        `Vector y length ${yData.length} incompatible with result length ${yLen}`
      );
    }
  } else {
    yData = new Float64Array(yLen);
  }

  // Prepare A matrix
  const aData = prepareMatrix(A);

  // Leading dimension and increments
  const lda = m;
  const incx = 1;
  const incy = 1;

  // Allocate WASM memory
  const aPtr = allocateDoubles(Module, aData, m * n);
  const xPtr = allocateDoubles(Module, xData, xLen);
  const yPtr = allocateDoubles(Module, yData, yLen);

  // Parameter pointers
  const transPtr = allocateInts(Module, [getTransChar(trans)]);
  const mPtr = allocateInts(Module, [m]);
  const nPtr = allocateInts(Module, [n]);
  const alphaPtr = allocateDoubles(Module, [alpha], 1);
  const betaPtr = allocateDoubles(Module, [beta], 1);
  const ldaPtr = allocateInts(Module, [lda]);
  const incxPtr = allocateInts(Module, [incx]);
  const incyPtr = allocateInts(Module, [incy]);

  try {
    // Call DGEMV
    Module._dgemv_(
      transPtr,
      mPtr,
      nPtr,
      alphaPtr,
      aPtr,
      ldaPtr,
      xPtr,
      incxPtr,
      betaPtr,
      yPtr,
      incyPtr
    );

    // Read result
    const result = readDoubles(Module, yPtr, yLen);

    return {
      y: result,
      n: yLen,
      success: true,
      message: 'Success',
    };
  } finally {
    freeAll(Module, [
      aPtr,
      xPtr,
      yPtr,
      transPtr,
      mPtr,
      nPtr,
      alphaPtr,
      betaPtr,
      ldaPtr,
      incxPtr,
      incyPtr,
    ]);
  }
}
