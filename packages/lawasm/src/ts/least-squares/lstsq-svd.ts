/**
 * Least Squares via SVD
 *
 * Solve least squares problems using singular value decomposition.
 * Safe for rank-deficient matrices.
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
import type { Matrix, Vector, LstSqSVDOptions, LstSqSVDResult } from './types.js';

/**
 * Solve a least squares problem using SVD.
 *
 * Computes the minimum-norm solution to: minimize ||Ax - b||_2
 *
 * This method is safe for rank-deficient matrices. Singular values
 * below rcond * max(s) are treated as zero.
 *
 * @param A - Coefficient matrix (m × n)
 * @param b - Right-hand side vector/matrix (m × nrhs)
 * @param options - Computation options
 * @returns Least squares solution with singular values and rank
 *
 * @example
 * ```typescript
 * // Rank-deficient system
 * const A = [[1, 2], [2, 4], [3, 6]];
 * const b = [1, 2, 3];
 * const { x, rank, s } = lstsqSVD(A, b);
 * // rank = 1 (columns are linearly dependent)
 * ```
 */
export function lstsqSVD(
  A: Matrix,
  b: Vector | Matrix,
  options: LstSqSVDOptions = {}
): LstSqSVDResult {
  const Module = getLAPACKModule();

  const {
    rcond = -1,
    algorithm = 'gelsd',
    overwriteA = false,
    overwriteB = false,
  } = options;

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);
  const minmn = Math.min(m, n);

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

  // DGELSS/DGELSD require b to have at least max(m, n) rows
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
  const sPtr = allocateDoubles(Module, null, minmn);
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

  // Parameter pointers
  const mPtr = allocateInts(Module, [m]);
  const nPtr = allocateInts(Module, [n]);
  const nrhsPtr = allocateInts(Module, [nrhs]);
  const ldaPtr = allocateInts(Module, [m]);
  const ldbPtr = allocateInts(Module, [ldb]);

  try {
    let info: number;

    if (algorithm === 'gelsd') {
      // DGELSD - divide and conquer (faster)
      const liworkPtr = allocateInts(Module, [-1]);
      const iworkQueryPtr = allocateInts(Module, [0]);

      // Workspace query
      Module._dgelsd_(
        mPtr, nPtr, nrhsPtr, aPtr, ldaPtr, bPtr, ldbPtr,
        sPtr, rcondPtr, rankPtr, workQueryPtr, lworkPtr,
        iworkQueryPtr, infoPtr
      );

      info = readInt(Module, infoPtr);
      if (info !== 0) {
        freeAll(Module, [liworkPtr, iworkQueryPtr]);
        return createErrorResult(m, n, nrhs, info, 'DGELSD');
      }

      const lwork = Math.ceil(Module.HEAPF64[workQueryPtr >> 3]);
      const liwork = Module.HEAP32[iworkQueryPtr >> 2];

      const workPtr = allocateDoubles(Module, null, lwork);
      const iworkPtr = allocateInts(Module, null, liwork);

      Module.HEAP32[lworkPtr >> 2] = lwork;

      // Compute least squares
      Module._dgelsd_(
        mPtr, nPtr, nrhsPtr, aPtr, ldaPtr, bPtr, ldbPtr,
        sPtr, rcondPtr, rankPtr, workPtr, lworkPtr,
        iworkPtr, infoPtr
      );

      info = readInt(Module, infoPtr);
      freeAll(Module, [workPtr, iworkPtr, liworkPtr, iworkQueryPtr]);

    } else {
      // DGELSS - standard SVD
      // Workspace query
      Module._dgelss_(
        mPtr, nPtr, nrhsPtr, aPtr, ldaPtr, bPtr, ldbPtr,
        sPtr, rcondPtr, rankPtr, workQueryPtr, lworkPtr, infoPtr
      );

      info = readInt(Module, infoPtr);
      if (info !== 0) {
        return createErrorResult(m, n, nrhs, info, 'DGELSS');
      }

      const lwork = Math.ceil(Module.HEAPF64[workQueryPtr >> 3]);
      const workPtr = allocateDoubles(Module, null, lwork);
      Module.HEAP32[lworkPtr >> 2] = lwork;

      // Compute least squares
      Module._dgelss_(
        mPtr, nPtr, nrhsPtr, aPtr, ldaPtr, bPtr, ldbPtr,
        sPtr, rcondPtr, rankPtr, workPtr, lworkPtr, infoPtr
      );

      info = readInt(Module, infoPtr);
      freeAll(Module, [workPtr]);
    }

    if (info !== 0) {
      return createErrorResult(m, n, nrhs, info, algorithm === 'gelsd' ? 'DGELSD' : 'DGELSS');
    }

    // Read results
    const bResult = readDoubles(Module, bPtr, ldb * nrhs);
    const x = new Float64Array(n * nrhs);
    for (let j = 0; j < nrhs; j++) {
      for (let i = 0; i < n; i++) {
        x[j * n + i] = bResult[j * ldb + i];
      }
    }

    const s = readDoubles(Module, sPtr, minmn);
    const rank = readInt(Module, rankPtr);

    const result: LstSqSVDResult = {
      x,
      s,
      rank,
      m,
      n,
      nrhs,
      info: 0,
      success: true,
      message: 'Success',
    };

    // Compute residuals for overdetermined full-rank systems
    if (m > n && rank === n) {
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
      aPtr, bPtr, sPtr, rcondPtr, rankPtr, infoPtr, lworkPtr, workQueryPtr,
      mPtr, nPtr, nrhsPtr, ldaPtr, ldbPtr
    ]);
  }
}

function createErrorResult(m: number, n: number, nrhs: number, info: number, routine: string): LstSqSVDResult {
  return {
    x: new Float64Array(0),
    s: new Float64Array(0),
    rank: 0,
    m,
    n,
    nrhs,
    info,
    success: false,
    message: getLapackErrorMessage(routine, info),
  };
}
