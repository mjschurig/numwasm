/**
 * General Matrix Multiplication (DGEMM)
 *
 * Compute C = alpha*op(A)*op(B) + beta*C
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
import type { Matrix, MatmulOptions, MatmulResult, TransposeOp } from './types.js';

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
 * Compute general matrix multiplication C = alpha*op(A)*op(B) + beta*C.
 *
 * Uses DGEMM from BLAS.
 *
 * @param A - First input matrix
 * @param B - Second input matrix
 * @param options - Multiplication options
 * @returns Result matrix C
 *
 * @example
 * ```typescript
 * const A = [[1, 2], [3, 4]];
 * const B = [[5, 6], [7, 8]];
 *
 * // Simple multiplication C = A*B
 * const { C } = matmul(A, B);
 *
 * // With transpose: C = A^T * B
 * const { C: C2 } = matmul(A, B, { transA: 'T' });
 *
 * // With scaling: C = 2*A*B + 3*C0
 * const C0 = [[1, 0], [0, 1]];
 * const { C: C3 } = matmul(A, B, { alpha: 2, beta: 3, C: C0 });
 * ```
 */
export function matmul(A: Matrix, B: Matrix, options: MatmulOptions = {}): MatmulResult {
  const Module = getLAPACKModule();

  const {
    transA = 'N',
    transB = 'N',
    alpha = 1.0,
    beta = 0.0,
    C,
  } = options;

  // Get dimensions of A and B
  const [rowsA, colsA] = getMatrixDimensions(A);
  const [rowsB, colsB] = getMatrixDimensions(B);

  // Determine effective dimensions after transpose
  const m = transA === 'N' ? rowsA : colsA; // Rows of op(A) and C
  const k = transA === 'N' ? colsA : rowsA; // Cols of op(A) = Rows of op(B)
  const n = transB === 'N' ? colsB : rowsB; // Cols of op(B) and C

  // Check dimension compatibility
  const kB = transB === 'N' ? rowsB : colsB;
  if (k !== kB) {
    throw new Error(
      `Incompatible dimensions: op(A) is ${m}×${k}, op(B) is ${kB}×${n}`
    );
  }

  // Prepare matrices
  const aData = prepareMatrix(A);
  const bData = prepareMatrix(B);

  // Prepare or create C
  let cData: Float64Array;
  if (C !== undefined) {
    const [rowsC, colsC] = getMatrixDimensions(C);
    if (rowsC !== m || colsC !== n) {
      throw new Error(
        `C dimensions ${rowsC}×${colsC} incompatible with result ${m}×${n}`
      );
    }
    cData = prepareMatrix(C);
  } else {
    cData = new Float64Array(m * n);
  }

  // Leading dimensions
  const lda = rowsA;
  const ldb = rowsB;
  const ldc = m;

  // Allocate WASM memory
  const aPtr = allocateDoubles(Module, aData, rowsA * colsA);
  const bPtr = allocateDoubles(Module, bData, rowsB * colsB);
  const cPtr = allocateDoubles(Module, cData, m * n);

  // Parameter pointers
  const transAPtr = allocateInts(Module, [getTransChar(transA)]);
  const transBPtr = allocateInts(Module, [getTransChar(transB)]);
  const mPtr = allocateInts(Module, [m]);
  const nPtr = allocateInts(Module, [n]);
  const kPtr = allocateInts(Module, [k]);
  const alphaPtr = allocateDoubles(Module, [alpha], 1);
  const betaPtr = allocateDoubles(Module, [beta], 1);
  const ldaPtr = allocateInts(Module, [lda]);
  const ldbPtr = allocateInts(Module, [ldb]);
  const ldcPtr = allocateInts(Module, [ldc]);

  try {
    // Call DGEMM
    Module._dgemm_(
      transAPtr,
      transBPtr,
      mPtr,
      nPtr,
      kPtr,
      alphaPtr,
      aPtr,
      ldaPtr,
      bPtr,
      ldbPtr,
      betaPtr,
      cPtr,
      ldcPtr
    );

    // Read result
    const result = readDoubles(Module, cPtr, m * n);

    return {
      C: result,
      m,
      n,
      success: true,
      message: 'Success',
    };
  } finally {
    freeAll(Module, [
      aPtr,
      bPtr,
      cPtr,
      transAPtr,
      transBPtr,
      mPtr,
      nPtr,
      kPtr,
      alphaPtr,
      betaPtr,
      ldaPtr,
      ldbPtr,
      ldcPtr,
    ]);
  }
}
