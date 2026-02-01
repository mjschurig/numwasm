/**
 * Least Squares via QR Factorization
 *
 * Solve least squares problems using QR (overdetermined) or LQ (underdetermined) factorization.
 */

import { getLAPACKModule } from '../loader.js';
import {
  allocateDoubles,
  allocateInts,
  readDoubles,
  readInt,
  freeAll,
  prepareMatrix,
  prepareVector,
  getMatrixDimensions,
  getLapackErrorMessage,
  CHAR,
} from '../helpers.js';
import type { Matrix, Vector, LstSqOptions, LstSqResult } from './types.js';

/**
 * Solve a least squares problem using QR factorization.
 *
 * For overdetermined systems (m >= n): minimize ||Ax - b||_2
 * For underdetermined systems (m < n): find minimum norm solution to Ax = b
 *
 * This function assumes A has full rank. For rank-deficient matrices,
 * use lstsqSVD or lstsqGelsy instead.
 *
 * @param A - Coefficient matrix (m × n)
 * @param b - Right-hand side vector/matrix (m × nrhs)
 * @param options - Computation options
 * @returns Least squares solution
 *
 * @example
 * ```typescript
 * // Overdetermined system (3 equations, 2 unknowns)
 * const A = [[1, 1], [1, 2], [1, 3]];
 * const b = [1, 2, 2];
 * const { x, residuals } = lstsq(A, b);
 * ```
 *
 * @example
 * ```typescript
 * // Underdetermined system (2 equations, 3 unknowns)
 * const A = [[1, 1, 1], [1, 2, 3]];
 * const b = [3, 6];
 * const { x } = lstsq(A, b);
 * // x is the minimum norm solution
 * ```
 */
export function lstsq(
  A: Matrix,
  b: Vector | Matrix,
  options: LstSqOptions = {}
): LstSqResult {
  const Module = getLAPACKModule();

  const { overwriteA = false, overwriteB = false } = options;

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);

  // Determine if b is a vector or matrix
  let bData: Float64Array;
  let nrhs: number;
  let bRows: number;

  if (Array.isArray(b) && Array.isArray(b[0])) {
    // b is a 2D matrix
    const [bm, bn] = getMatrixDimensions(b as Matrix);
    if (bm !== m) {
      throw new Error(`Dimension mismatch: A has ${m} rows but b has ${bm} rows`);
    }
    bData = prepareMatrix(b as Matrix);
    nrhs = bn;
    bRows = bm;
  } else {
    // b is a vector
    const bVec = prepareVector(b as Vector);
    if (bVec.length !== m) {
      throw new Error(`Dimension mismatch: A has ${m} rows but b has ${bVec.length} elements`);
    }
    bData = bVec;
    nrhs = 1;
    bRows = bVec.length;
  }

  // Prepare matrix A
  const aData = prepareMatrix(A);

  // DGELS requires b to have at least max(m, n) rows for the solution
  const ldb = Math.max(m, n);
  const bExtended = new Float64Array(ldb * nrhs);

  // Copy b data (column-major)
  for (let j = 0; j < nrhs; j++) {
    for (let i = 0; i < bRows; i++) {
      bExtended[j * ldb + i] = bData[j * bRows + i];
    }
  }

  // Allocate WASM memory
  const aPtr = allocateDoubles(Module, overwriteA ? null : aData, m * n);
  const bPtr = allocateDoubles(Module, overwriteB ? null : bExtended, ldb * nrhs);
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
  if (!overwriteB) {
    const baseIdx = bPtr >> 3;
    for (let i = 0; i < bExtended.length; i++) {
      Module.HEAPF64[baseIdx + i] = bExtended[i];
    }
  }

  // Parameter pointers
  const transPtr = allocateInts(Module, [CHAR.N]);
  const mPtr = allocateInts(Module, [m]);
  const nPtr = allocateInts(Module, [n]);
  const nrhsPtr = allocateInts(Module, [nrhs]);
  const ldaPtr = allocateInts(Module, [m]);
  const ldbPtr = allocateInts(Module, [ldb]);

  try {
    // Workspace query
    Module._dgels_(
      transPtr, mPtr, nPtr, nrhsPtr, aPtr, ldaPtr,
      bPtr, ldbPtr, workQueryPtr, lworkPtr, infoPtr
    );

    let info = readInt(Module, infoPtr);
    if (info !== 0) {
      return createErrorResult(m, n, nrhs, info, 'DGELS');
    }

    const lwork = Math.ceil(Module.HEAPF64[workQueryPtr >> 3]);
    const workPtr = allocateDoubles(Module, null, lwork);
    Module.HEAP32[lworkPtr >> 2] = lwork;

    // Compute least squares solution
    Module._dgels_(
      transPtr, mPtr, nPtr, nrhsPtr, aPtr, ldaPtr,
      bPtr, ldbPtr, workPtr, lworkPtr, infoPtr
    );

    info = readInt(Module, infoPtr);
    freeAll(Module, [workPtr]);

    if (info !== 0) {
      return createErrorResult(m, n, nrhs, info, 'DGELS');
    }

    // Read solution (first n rows of b for each column)
    const bResult = readDoubles(Module, bPtr, ldb * nrhs);
    const x = new Float64Array(n * nrhs);
    for (let j = 0; j < nrhs; j++) {
      for (let i = 0; i < n; i++) {
        x[j * n + i] = bResult[j * ldb + i];
      }
    }

    const result: LstSqResult = {
      x,
      m,
      n,
      nrhs,
      info: 0,
      success: true,
      message: 'Success',
    };

    // Compute residuals for overdetermined systems
    if (m > n) {
      const residuals = new Float64Array(nrhs);
      for (let j = 0; j < nrhs; j++) {
        let sumSq = 0;
        for (let i = n; i < m; i++) {
          const r = bResult[j * ldb + i];
          sumSq += r * r;
        }
        residuals[j] = sumSq;
      }
      result.residuals = residuals;
    }

    return result;
  } finally {
    freeAll(Module, [
      aPtr, bPtr, infoPtr, lworkPtr, workQueryPtr,
      transPtr, mPtr, nPtr, nrhsPtr, ldaPtr, ldbPtr
    ]);
  }
}

function createErrorResult(m: number, n: number, nrhs: number, info: number, routine: string): LstSqResult {
  return {
    x: new Float64Array(0),
    m,
    n,
    nrhs,
    info,
    success: false,
    message: getLapackErrorMessage(routine, info),
  };
}
