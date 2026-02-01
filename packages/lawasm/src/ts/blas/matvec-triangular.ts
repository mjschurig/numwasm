/**
 * Triangular Matrix-Vector Multiplication (DTRMV)
 *
 * Compute x = op(A)*x where A is triangular
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
import type { Matrix, Vector, MatvecTriangularOptions, MatvecTriangularResult, TransposeOp } from './types.js';

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
 * Compute triangular matrix-vector multiplication x = op(A)*x.
 *
 * Uses DTRMV from BLAS. The input vector x is overwritten with the result.
 *
 * @param A - Triangular matrix (n × n)
 * @param x - Input/output vector (length n)
 * @param options - Multiplication options
 * @returns Result vector x
 *
 * @example
 * ```typescript
 * // Upper triangular A
 * const A = [[2, 3], [0, 4]];
 * const x = [1, 2];
 *
 * // x = A*x
 * const { x: result } = matvecTriangular(A, x, { upper: true });
 *
 * // x = A^T*x
 * const { x: result2 } = matvecTriangular(A, [1, 2], { upper: true, trans: 'T' });
 * ```
 */
export function matvecTriangular(
  A: Matrix,
  x: Vector,
  options: MatvecTriangularOptions = {}
): MatvecTriangularResult {
  const Module = getLAPACKModule();

  const {
    upper = true,
    trans = 'N',
    unitDiagonal = false,
  } = options;

  // Get dimensions
  const [rowsA, colsA] = getMatrixDimensions(A);
  if (rowsA !== colsA) {
    throw new Error(`Matrix A must be square, got ${rowsA}×${colsA}`);
  }
  const n = rowsA;

  const xData = prepareVector(x);
  if (xData.length !== n) {
    throw new Error(`Vector x length ${xData.length} must equal matrix dimension ${n}`);
  }

  // Prepare matrix
  const aData = prepareMatrix(A);

  // Leading dimension and increment
  const lda = n;
  const incx = 1;

  // Allocate WASM memory
  const aPtr = allocateDoubles(Module, aData, n * n);
  const xPtr = allocateDoubles(Module, xData, n);

  // Parameter pointers
  const uploPtr = allocateInts(Module, [upper ? CHAR.U : CHAR.L]);
  const transPtr = allocateInts(Module, [getTransChar(trans)]);
  const diagPtr = allocateInts(Module, [unitDiagonal ? CHAR.UNIT : CHAR.N]);
  const nPtr = allocateInts(Module, [n]);
  const ldaPtr = allocateInts(Module, [lda]);
  const incxPtr = allocateInts(Module, [incx]);

  try {
    // Call DTRMV
    Module._dtrmv_(
      uploPtr,
      transPtr,
      diagPtr,
      nPtr,
      aPtr,
      ldaPtr,
      xPtr,
      incxPtr
    );

    // Read result
    const result = readDoubles(Module, xPtr, n);

    return {
      x: result,
      n,
      success: true,
      message: 'Success',
    };
  } finally {
    freeAll(Module, [aPtr, xPtr, uploPtr, transPtr, diagPtr, nPtr, ldaPtr, incxPtr]);
  }
}
