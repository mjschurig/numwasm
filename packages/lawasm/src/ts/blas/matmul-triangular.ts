/**
 * Triangular Matrix Multiplication (DTRMM)
 *
 * Compute B = alpha*op(A)*B or B = alpha*B*op(A)
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
import type { Matrix, MatmulTriangularOptions, MatmulTriangularResult, TransposeOp } from './types.js';

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
 * Compute triangular matrix multiplication B = alpha*op(A)*B or B = alpha*B*op(A).
 *
 * Uses DTRMM from BLAS.
 *
 * @param A - Triangular matrix (n × n)
 * @param B - General matrix (m × n or n × m depending on side)
 * @param options - Multiplication options
 * @returns Result matrix B
 *
 * @example
 * ```typescript
 * // Upper triangular A
 * const A = [[2, 3], [0, 4]];
 * const B = [[1, 2], [3, 4]];
 *
 * // B = A*B (left side, default)
 * const { B: result1 } = matmulTriangular(A, B, { upper: true });
 *
 * // B = B*A (right side)
 * const { B: result2 } = matmulTriangular(A, B, { side: 'R', upper: true });
 *
 * // B = 2*A^T*B
 * const { B: result3 } = matmulTriangular(A, B, { transA: 'T', alpha: 2 });
 * ```
 */
export function matmulTriangular(
  A: Matrix,
  B: Matrix,
  options: MatmulTriangularOptions = {}
): MatmulTriangularResult {
  const Module = getLAPACKModule();

  const {
    side = 'L',
    upper = true,
    transA = 'N',
    unitDiagonal = false,
    alpha = 1.0,
  } = options;

  // Get dimensions
  const [nA] = getMatrixDimensions(A);
  const [m, n] = getMatrixDimensions(B);

  // Check that A is square
  const [rowsA, colsA] = getMatrixDimensions(A);
  if (rowsA !== colsA) {
    throw new Error(`Matrix A must be square, got ${rowsA}×${colsA}`);
  }

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
    // Call DTRMM
    Module._dtrmm_(
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

    // Read result (B is overwritten)
    const result = readDoubles(Module, bPtr, m * n);

    return {
      B: result,
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
