/**
 * Triangular Matrix Solve (DTRSM)
 *
 * Solve op(A)*X = alpha*B or X*op(A) = alpha*B
 * where A is triangular.
 */

import { getLAPACKModule } from '../loader.js';
import {
  allocateDoubles,
  allocateInts,
  readDoubles,
  freeAll,
  prepareMatrix,
  getMatrixDimensions,
  CHAR,
} from '../helpers.js';
import type { Matrix, SolveMatrixTriangularOptions, SolveMatrixTriangularResult, TransposeOp } from './types.js';

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
 * Solve a triangular matrix system op(A)*X = alpha*B or X*op(A) = alpha*B.
 *
 * Uses DTRSM from BLAS.
 *
 * @param A - Triangular matrix (n × n)
 * @param B - Right-hand side matrix
 * @param options - Solve options
 * @returns Solution matrix X
 *
 * @example
 * ```typescript
 * // Upper triangular A
 * const A = [[2, 3], [0, 4]];
 * const B = [[1, 2], [3, 4]];
 *
 * // Solve A*X = B (left side)
 * const { X } = solveMatrixTriangular(A, B, { upper: true });
 *
 * // Solve X*A = B (right side)
 * const { X: X2 } = solveMatrixTriangular(A, B, { side: 'R', upper: true });
 *
 * // Solve A^T*X = 2*B
 * const { X: X3 } = solveMatrixTriangular(A, B, { transA: 'T', alpha: 2 });
 * ```
 */
export function solveMatrixTriangular(
  A: Matrix,
  B: Matrix,
  options: SolveMatrixTriangularOptions = {}
): SolveMatrixTriangularResult {
  const Module = getLAPACKModule();

  const {
    side = 'L',
    upper = true,
    transA = 'N',
    unitDiagonal = false,
    alpha = 1.0,
  } = options;

  // Get dimensions
  const [rowsA, colsA] = getMatrixDimensions(A);
  const [m, n] = getMatrixDimensions(B);

  // Check that A is square
  if (rowsA !== colsA) {
    throw new Error(`Matrix A must be square, got ${rowsA}×${colsA}`);
  }

  const nA = rowsA;

  // Check dimension compatibility
  if (side === 'L' && nA !== m) {
    throw new Error(`For left side, A (${nA}×${nA}) and B (${m}×${n}) dimensions incompatible`);
  }
  if (side === 'R' && nA !== n) {
    throw new Error(`For right side, A (${nA}×${nA}) and B (${m}×${n}) dimensions incompatible`);
  }

  // Prepare matrices
  const aData = prepareMatrix(A);
  const bData = prepareMatrix(B);

  // Leading dimensions
  const lda = nA;
  const ldb = m;

  // Allocate WASM memory
  const aPtr = allocateDoubles(Module, aData, nA * nA);
  const bPtr = allocateDoubles(Module, bData, m * n);

  // Parameter pointers
  const sidePtr = allocateInts(Module, [side === 'L' ? CHAR.LEFT : CHAR.RIGHT]);
  const uploPtr = allocateInts(Module, [upper ? CHAR.U : CHAR.L]);
  const transAPtr = allocateInts(Module, [getTransChar(transA)]);
  const diagPtr = allocateInts(Module, [unitDiagonal ? CHAR.UNIT : CHAR.N]);
  const mPtr = allocateInts(Module, [m]);
  const nPtr = allocateInts(Module, [n]);
  const alphaPtr = allocateDoubles(Module, [alpha], 1);
  const ldaPtr = allocateInts(Module, [lda]);
  const ldbPtr = allocateInts(Module, [ldb]);

  try {
    // Call DTRSM
    Module._dtrsm_(
      sidePtr,
      uploPtr,
      transAPtr,
      diagPtr,
      mPtr,
      nPtr,
      alphaPtr,
      aPtr,
      ldaPtr,
      bPtr,
      ldbPtr
    );

    // Read result (B is overwritten with X)
    const result = readDoubles(Module, bPtr, m * n);

    return {
      X: result,
      m,
      n,
      success: true,
      message: 'Success',
    };
  } finally {
    freeAll(Module, [
      aPtr,
      bPtr,
      sidePtr,
      uploPtr,
      transAPtr,
      diagPtr,
      mPtr,
      nPtr,
      alphaPtr,
      ldaPtr,
      ldbPtr,
    ]);
  }
}
