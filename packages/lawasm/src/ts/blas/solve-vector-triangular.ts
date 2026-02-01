/**
 * Triangular Vector Solve (DTRSV)
 *
 * Solve op(A)*x = b where A is triangular
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
import type { Matrix, Vector, SolveVectorTriangularOptions, SolveVectorTriangularResult, TransposeOp } from './types.js';

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
 * Solve triangular system op(A)*x = b.
 *
 * Uses DTRSV from BLAS. The input vector b is overwritten with the solution x.
 *
 * @param A - Triangular matrix (n × n)
 * @param b - Right-hand side vector (length n)
 * @param options - Solve options
 * @returns Solution vector x
 *
 * @example
 * ```typescript
 * // Upper triangular A
 * const A = [[2, 3], [0, 4]];
 * const b = [1, 2];
 *
 * // Solve A*x = b
 * const { x } = solveVectorTriangular(A, b, { upper: true });
 *
 * // Solve A^T*x = b
 * const { x: x2 } = solveVectorTriangular(A, b, { upper: true, trans: 'T' });
 * ```
 */
export function solveVectorTriangular(
  A: Matrix,
  b: Vector,
  options: SolveVectorTriangularOptions = {}
): SolveVectorTriangularResult {
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

  const bData = prepareVector(b);
  if (bData.length !== n) {
    throw new Error(`Vector b length ${bData.length} must equal matrix dimension ${n}`);
  }

  // Prepare matrix
  const aData = prepareMatrix(A);

  // Leading dimension and increment
  const lda = n;
  const incx = 1;

  // Allocate WASM memory
  const aPtr = allocateDoubles(Module, aData, n * n);
  const xPtr = allocateDoubles(Module, bData, n); // b is overwritten with x

  // Parameter pointers
  const uploPtr = allocateInts(Module, [upper ? CHAR.U : CHAR.L]);
  const transPtr = allocateInts(Module, [getTransChar(trans)]);
  const diagPtr = allocateInts(Module, [unitDiagonal ? CHAR.UNIT : CHAR.N]);
  const nPtr = allocateInts(Module, [n]);
  const ldaPtr = allocateInts(Module, [lda]);
  const incxPtr = allocateInts(Module, [incx]);

  try {
    // Call DTRSV
    Module._dtrsv_(
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
