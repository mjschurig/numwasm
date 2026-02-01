/**
 * Singular Value Decomposition
 *
 * Compute the SVD of a general matrix: A = U * S * V^T
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
import type { Matrix, SVDOptions, SVDResult } from './types.js';

/**
 * Compute the Singular Value Decomposition of a general m×n matrix A.
 *
 * A = U * S * V^T
 *
 * where U is m×m (or m×k) orthogonal, S is a diagonal matrix with singular values,
 * and V is n×n (or k×n) orthogonal, with k = min(m, n).
 *
 * @param A - Input matrix (m × n). Can be 2D array (row-major) or 1D array (column-major).
 * @param options - SVD options
 * @returns SVD result with U, s (singular values), and Vt
 *
 * @example
 * ```typescript
 * const A = [[1, 2], [3, 4], [5, 6]];
 * const { U, s, Vt, success } = svd(A);
 * // U is 3×2, s has 2 elements, Vt is 2×2 (reduced SVD)
 * ```
 *
 * @example
 * ```typescript
 * // Full SVD
 * const { U, s, Vt } = svd(A, { mode: 'full' });
 * // U is 3×3, Vt is 2×2
 * ```
 *
 * @example
 * ```typescript
 * // Only singular values (faster)
 * const { s } = svd(A, { mode: 'none' });
 * ```
 */
export function svd(A: Matrix, options: SVDOptions = {}): SVDResult {
  const Module = getLAPACKModule();

  const {
    mode = 'reduced',
    algorithm = 'gesdd',
    overwriteA = false,
  } = options;

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);
  const k = Math.min(m, n);

  // Prepare matrix
  const aData = prepareMatrix(A);

  // Determine job parameters
  let jobu: number;
  let jobvt: number;
  let jobz: number;
  let uRows: number;
  let uCols: number;
  let vtRows: number;
  let vtCols: number;

  if (mode === 'full') {
    jobu = CHAR.A;   // All columns of U
    jobvt = CHAR.A;  // All rows of Vt
    jobz = CHAR.A;   // All for gesdd
    uRows = m;
    uCols = m;
    vtRows = n;
    vtCols = n;
  } else if (mode === 'reduced') {
    jobu = CHAR.S;   // First k columns of U
    jobvt = CHAR.S;  // First k rows of Vt
    jobz = CHAR.S;   // Some for gesdd
    uRows = m;
    uCols = k;
    vtRows = k;
    vtCols = n;
  } else {
    // mode === 'none'
    jobu = CHAR.N;   // No U
    jobvt = CHAR.N;  // No Vt
    jobz = CHAR.N;   // None for gesdd
    uRows = 1;
    uCols = 1;
    vtRows = 1;
    vtCols = 1;
  }

  // Allocate WASM memory
  const aPtr = allocateDoubles(Module, overwriteA ? null : aData, m * n);
  const sPtr = allocateDoubles(Module, null, k);
  const uPtr = allocateDoubles(Module, null, mode !== 'none' ? uRows * uCols : 1);
  const vtPtr = allocateDoubles(Module, null, mode !== 'none' ? vtRows * vtCols : 1);
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

  // Additional allocations for parameters
  const mPtr = allocateInts(Module, [m]);
  const nPtr = allocateInts(Module, [n]);
  const ldaPtr = allocateInts(Module, [m]);
  const lduPtr = allocateInts(Module, [mode !== 'none' ? uRows : 1]);
  const ldvtPtr = allocateInts(Module, [mode !== 'none' ? vtRows : 1]);

  let iworkPtr = 0;
  if (algorithm === 'gesdd') {
    iworkPtr = allocateInts(Module, null, 8 * k);
  }

  try {
    let info: number;

    if (algorithm === 'gesdd') {
      // Use divide and conquer (faster)
      const jobzPtr = allocateInts(Module, [jobz]);

      // Workspace query
      Module._dgesdd_(
        jobzPtr, mPtr, nPtr, aPtr, ldaPtr,
        sPtr, uPtr, lduPtr, vtPtr, ldvtPtr,
        workQueryPtr, lworkPtr, iworkPtr, infoPtr
      );

      info = readInt(Module, infoPtr);
      if (info !== 0) {
        freeAll(Module, [jobzPtr]);
        return createErrorResult(m, n, k, info, 'DGESDD');
      }

      const lwork = Math.ceil(Module.HEAPF64[workQueryPtr >> 3]);
      const workPtr = allocateDoubles(Module, null, lwork);
      Module.HEAP32[lworkPtr >> 2] = lwork;

      // Compute SVD
      Module._dgesdd_(
        jobzPtr, mPtr, nPtr, aPtr, ldaPtr,
        sPtr, uPtr, lduPtr, vtPtr, ldvtPtr,
        workPtr, lworkPtr, iworkPtr, infoPtr
      );

      info = readInt(Module, infoPtr);
      freeAll(Module, [workPtr, jobzPtr]);
    } else {
      // Use standard algorithm (more accurate)
      const jobuPtr = allocateInts(Module, [jobu]);
      const jobvtPtr = allocateInts(Module, [jobvt]);

      // Workspace query
      Module._dgesvd_(
        jobuPtr, jobvtPtr, mPtr, nPtr, aPtr, ldaPtr,
        sPtr, uPtr, lduPtr, vtPtr, ldvtPtr,
        workQueryPtr, lworkPtr, infoPtr
      );

      info = readInt(Module, infoPtr);
      if (info !== 0) {
        freeAll(Module, [jobuPtr, jobvtPtr]);
        return createErrorResult(m, n, k, info, 'DGESVD');
      }

      const lwork = Math.ceil(Module.HEAPF64[workQueryPtr >> 3]);
      const workPtr = allocateDoubles(Module, null, lwork);
      Module.HEAP32[lworkPtr >> 2] = lwork;

      // Compute SVD
      Module._dgesvd_(
        jobuPtr, jobvtPtr, mPtr, nPtr, aPtr, ldaPtr,
        sPtr, uPtr, lduPtr, vtPtr, ldvtPtr,
        workPtr, lworkPtr, infoPtr
      );

      info = readInt(Module, infoPtr);
      freeAll(Module, [workPtr, jobuPtr, jobvtPtr]);
    }

    if (info !== 0) {
      return createErrorResult(m, n, k, info, algorithm === 'gesdd' ? 'DGESDD' : 'DGESVD');
    }

    // Read results
    const s = readDoubles(Module, sPtr, k);

    const result: SVDResult = {
      s,
      m,
      n,
      k,
      info: 0,
      success: true,
      message: 'Success',
    };

    if (mode !== 'none') {
      result.U = readDoubles(Module, uPtr, uRows * uCols);
      result.Vt = readDoubles(Module, vtPtr, vtRows * vtCols);
    }

    return result;
  } finally {
    const ptrsToFree = [
      aPtr, sPtr, uPtr, vtPtr, infoPtr, lworkPtr, workQueryPtr,
      mPtr, nPtr, ldaPtr, lduPtr, ldvtPtr
    ];
    if (iworkPtr !== 0) {
      ptrsToFree.push(iworkPtr);
    }
    freeAll(Module, ptrsToFree);
  }
}

function createErrorResult(m: number, n: number, k: number, info: number, routine: string): SVDResult {
  return {
    s: new Float64Array(0),
    m,
    n,
    k,
    info,
    success: false,
    message: getLapackErrorMessage(routine, info),
  };
}
