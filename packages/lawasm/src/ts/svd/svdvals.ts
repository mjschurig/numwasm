/**
 * Singular Values Only
 *
 * Compute only the singular values (faster than full SVD).
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
import type { Matrix, SVDValsOptions, SVDValsResult } from './types.js';

/**
 * Compute only the singular values of a general m×n matrix A.
 *
 * This is faster than computing the full SVD when only singular values are needed.
 *
 * @param A - Input matrix (m × n). Can be 2D array (row-major) or 1D array (column-major).
 * @param options - Options
 * @returns Singular values in descending order
 *
 * @example
 * ```typescript
 * const A = [[1, 2], [3, 4], [5, 6]];
 * const { s, success } = svdvals(A);
 * // s contains the 2 singular values in descending order
 * ```
 */
export function svdvals(A: Matrix, options: SVDValsOptions = {}): SVDValsResult {
  const Module = getLAPACKModule();

  const { algorithm = 'gesdd', overwriteA = false } = options;

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);
  const k = Math.min(m, n);

  // Prepare matrix
  const aData = prepareMatrix(A);

  // Allocate WASM memory
  const aPtr = allocateDoubles(Module, overwriteA ? null : aData, m * n);
  const sPtr = allocateDoubles(Module, null, k);
  const uPtr = allocateDoubles(Module, null, 1); // Dummy
  const vtPtr = allocateDoubles(Module, null, 1); // Dummy
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

  const mPtr = allocateInts(Module, [m]);
  const nPtr = allocateInts(Module, [n]);
  const ldaPtr = allocateInts(Module, [m]);
  const lduPtr = allocateInts(Module, [1]);
  const ldvtPtr = allocateInts(Module, [1]);

  let iworkPtr = 0;
  if (algorithm === 'gesdd') {
    iworkPtr = allocateInts(Module, null, 8 * k);
  }

  try {
    let info: number;

    if (algorithm === 'gesdd') {
      const jobzPtr = allocateInts(Module, [CHAR.N]); // No vectors

      // Workspace query
      Module._dgesdd_(
        jobzPtr, mPtr, nPtr, aPtr, ldaPtr,
        sPtr, uPtr, lduPtr, vtPtr, ldvtPtr,
        workQueryPtr, lworkPtr, iworkPtr, infoPtr
      );

      info = readInt(Module, infoPtr);
      if (info !== 0) {
        freeAll(Module, [jobzPtr]);
        return {
          s: new Float64Array(0),
          m,
          n,
          info,
          success: false,
          message: getLapackErrorMessage('DGESDD', info),
        };
      }

      const lwork = Math.ceil(Module.HEAPF64[workQueryPtr >> 3]);
      const workPtr = allocateDoubles(Module, null, lwork);
      Module.HEAP32[lworkPtr >> 2] = lwork;

      // Compute singular values
      Module._dgesdd_(
        jobzPtr, mPtr, nPtr, aPtr, ldaPtr,
        sPtr, uPtr, lduPtr, vtPtr, ldvtPtr,
        workPtr, lworkPtr, iworkPtr, infoPtr
      );

      info = readInt(Module, infoPtr);
      freeAll(Module, [workPtr, jobzPtr]);
    } else {
      const jobuPtr = allocateInts(Module, [CHAR.N]);
      const jobvtPtr = allocateInts(Module, [CHAR.N]);

      // Workspace query
      Module._dgesvd_(
        jobuPtr, jobvtPtr, mPtr, nPtr, aPtr, ldaPtr,
        sPtr, uPtr, lduPtr, vtPtr, ldvtPtr,
        workQueryPtr, lworkPtr, infoPtr
      );

      info = readInt(Module, infoPtr);
      if (info !== 0) {
        freeAll(Module, [jobuPtr, jobvtPtr]);
        return {
          s: new Float64Array(0),
          m,
          n,
          info,
          success: false,
          message: getLapackErrorMessage('DGESVD', info),
        };
      }

      const lwork = Math.ceil(Module.HEAPF64[workQueryPtr >> 3]);
      const workPtr = allocateDoubles(Module, null, lwork);
      Module.HEAP32[lworkPtr >> 2] = lwork;

      // Compute singular values
      Module._dgesvd_(
        jobuPtr, jobvtPtr, mPtr, nPtr, aPtr, ldaPtr,
        sPtr, uPtr, lduPtr, vtPtr, ldvtPtr,
        workPtr, lworkPtr, infoPtr
      );

      info = readInt(Module, infoPtr);
      freeAll(Module, [workPtr, jobuPtr, jobvtPtr]);
    }

    if (info !== 0) {
      return {
        s: new Float64Array(0),
        m,
        n,
        info,
        success: false,
        message: getLapackErrorMessage(algorithm === 'gesdd' ? 'DGESDD' : 'DGESVD', info),
      };
    }

    // Read singular values
    const s = readDoubles(Module, sPtr, k);

    return {
      s,
      m,
      n,
      info: 0,
      success: true,
      message: 'Success',
    };
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
