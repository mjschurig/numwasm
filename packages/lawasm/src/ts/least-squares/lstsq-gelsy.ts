/**
 * Least Squares via QR with Column Pivoting
 *
 * Solve least squares problems using complete orthogonal factorization.
 * Generally faster than SVD for rank-deficient problems.
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
} from '../helpers.js';
import type { Matrix, Vector, LstSqGelsyOptions, LstSqGelsyResult } from './types.js';

/**
 * Solve a least squares problem using QR with column pivoting.
 *
 * Computes the minimum-norm solution to: minimize ||Ax - b||_2
 * using a complete orthogonal factorization of A.
 *
 * This method is generally faster than SVD-based methods for rank-deficient problems.
 *
 * @param A - Coefficient matrix (m × n)
 * @param b - Right-hand side vector/matrix (m × nrhs)
 * @param options - Computation options
 * @returns Least squares solution with rank and permutation info
 *
 * @example
 * ```typescript
 * const A = [[1, 2], [2, 4], [3, 6]];
 * const b = [1, 2, 3];
 * const { x, rank } = lstsqGelsy(A, b);
 * ```
 */
export function lstsqGelsy(
  A: Matrix,
  b: Vector | Matrix,
  options: LstSqGelsyOptions = {}
): LstSqGelsyResult {
  const Module = getLAPACKModule();

  const { rcond = -1, overwriteA = false, overwriteB = false } = options;

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);

  // Determine if b is a vector or matrix
  let bData: Float64Array;
  let nrhs: number;
  let bRows: number;

  if (Array.isArray(b) && Array.isArray(b[0])) {
    const [bm, bn] = getMatrixDimensions(b as Matrix);
    if (bm !== m) {
      throw new Error(`Dimension mismatch: A has ${m} rows but b has ${bm} rows`);
    }
    bData = prepareMatrix(b as Matrix);
    nrhs = bn;
    bRows = bm;
  } else {
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

  // DGELSY requires b to have at least max(m, n) rows
  const ldb = Math.max(m, n);
  const bExtended = new Float64Array(ldb * nrhs);
  for (let j = 0; j < nrhs; j++) {
    for (let i = 0; i < bRows; i++) {
      bExtended[j * ldb + i] = bData[j * bRows + i];
    }
  }

  // Allocate WASM memory
  const aPtr = allocateDoubles(Module, overwriteA ? null : aData, m * n);
  const bPtr = allocateDoubles(Module, overwriteB ? null : bExtended, ldb * nrhs);
  const jpvtPtr = allocateInts(Module, null, n); // Initialize to zeros (no initial pivoting preference)
  const rcondPtr = allocateDoubles(Module, [rcond]);
  const rankPtr = allocateInts(Module, [0]);
  const infoPtr = allocateInts(Module, [0]);
  const lworkPtr = allocateInts(Module, [-1]);
  const workQueryPtr = allocateDoubles(Module, null, 1);

  // Write data
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

  // Initialize jpvt to zeros
  const jpvtBaseIdx = jpvtPtr >> 2;
  for (let i = 0; i < n; i++) {
    Module.HEAP32[jpvtBaseIdx + i] = 0;
  }

  // Parameter pointers
  const mPtr = allocateInts(Module, [m]);
  const nPtr = allocateInts(Module, [n]);
  const nrhsPtr = allocateInts(Module, [nrhs]);
  const ldaPtr = allocateInts(Module, [m]);
  const ldbPtr = allocateInts(Module, [ldb]);

  try {
    // Workspace query
    Module._dgelsy_(
      mPtr, nPtr, nrhsPtr, aPtr, ldaPtr, bPtr, ldbPtr,
      jpvtPtr, rcondPtr, rankPtr, workQueryPtr, lworkPtr, infoPtr
    );

    let info = readInt(Module, infoPtr);
    if (info !== 0) {
      return createErrorResult(m, n, nrhs, info, 'DGELSY');
    }

    const lwork = Math.ceil(Module.HEAPF64[workQueryPtr >> 3]);
    const workPtr = allocateDoubles(Module, null, lwork);
    Module.HEAP32[lworkPtr >> 2] = lwork;

    // Compute least squares
    Module._dgelsy_(
      mPtr, nPtr, nrhsPtr, aPtr, ldaPtr, bPtr, ldbPtr,
      jpvtPtr, rcondPtr, rankPtr, workPtr, lworkPtr, infoPtr
    );

    info = readInt(Module, infoPtr);
    freeAll(Module, [workPtr]);

    if (info !== 0) {
      return createErrorResult(m, n, nrhs, info, 'DGELSY');
    }

    // Read results
    const bResult = readDoubles(Module, bPtr, ldb * nrhs);
    const x = new Float64Array(n * nrhs);
    for (let j = 0; j < nrhs; j++) {
      for (let i = 0; i < n; i++) {
        x[j * n + i] = bResult[j * ldb + i];
      }
    }

    const rank = readInt(Module, rankPtr);

    // Read pivot array
    const jpvt = new Int32Array(n);
    for (let i = 0; i < n; i++) {
      jpvt[i] = Module.HEAP32[jpvtBaseIdx + i];
    }

    return {
      x,
      rank,
      jpvt,
      m,
      n,
      nrhs,
      info: 0,
      success: true,
      message: 'Success',
    };
  } finally {
    freeAll(Module, [
      aPtr, bPtr, jpvtPtr, rcondPtr, rankPtr, infoPtr, lworkPtr, workQueryPtr,
      mPtr, nPtr, nrhsPtr, ldaPtr, ldbPtr
    ]);
  }
}

function createErrorResult(m: number, n: number, nrhs: number, info: number, routine: string): LstSqGelsyResult {
  return {
    x: new Float64Array(0),
    rank: 0,
    jpvt: new Int32Array(0),
    m,
    n,
    nrhs,
    info,
    success: false,
    message: getLapackErrorMessage(routine, info),
  };
}
