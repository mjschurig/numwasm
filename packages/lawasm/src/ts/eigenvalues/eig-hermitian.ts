/**
 * Hermitian Eigenvalue Problem
 *
 * Compute eigenvalues and eigenvectors of a complex Hermitian matrix.
 * All eigenvalues of a Hermitian matrix are real.
 */

import { getLAPACKModule } from '../loader.js';
import {
  allocateDoubles,
  allocateInts,
  readDoubles,
  readInt,
  freeAll,
  getLapackErrorMessage,
  CHAR,
} from '../helpers.js';
import type { Matrix, EigHermitianOptions, EigHermitianResult } from './types.js';

/**
 * Compute eigenvalues and optionally eigenvectors of a complex Hermitian n×n matrix A.
 *
 * For a Hermitian matrix (A = A^H), all eigenvalues are real and the eigenvectors
 * form a unitary basis.
 *
 * Input format: The matrix A should be provided as interleaved real/imaginary pairs.
 * For an n×n matrix, provide 2*n*n values where pairs (A[2k], A[2k+1]) represent
 * the real and imaginary parts of complex element k.
 *
 * @param A - Input Hermitian matrix as interleaved complex array.
 *            For 2D array input, each element should be [re, im].
 *            For 1D input, provide interleaved real/imaginary values.
 * @param n - Matrix dimension (must be provided for complex matrices)
 * @param options - Computation options
 * @returns Eigenvalues (real, in ascending order) and optionally complex eigenvectors
 *
 * @example
 * ```typescript
 * // 2x2 Hermitian matrix: [[2, 1+i], [1-i, 3]]
 * // Interleaved format: [2, 0, 1, 1, 1, -1, 3, 0] (column-major)
 * const A = [2, 0, 1, -1, 1, 1, 3, 0];
 * const { values, vectors, success } = eigHermitian(A, 2);
 * ```
 */
export function eigHermitian(
  A: Matrix,
  n: number,
  options: EigHermitianOptions = {}
): EigHermitianResult {
  const Module = getLAPACKModule();

  const {
    computeVectors = true,
    algorithm = 'heevd',
    uplo = 'lower',
    overwriteA = false,
  } = options;

  // For complex matrices, we expect interleaved format
  // Total length should be 2*n*n
  let aData: Float64Array;

  if (Array.isArray(A) && Array.isArray(A[0]) && Array.isArray((A as unknown as number[][][])[0][0])) {
    // 2D array of [re, im] pairs - convert to interleaved
    const matrix = A as unknown as number[][][];
    if (matrix.length !== n || matrix[0].length !== n) {
      throw new Error(`Matrix dimensions don't match n=${n}`);
    }
    aData = new Float64Array(2 * n * n);
    for (let j = 0; j < n; j++) {
      for (let i = 0; i < n; i++) {
        const idx = 2 * (j * n + i);
        aData[idx] = matrix[i][j][0];     // Real part
        aData[idx + 1] = matrix[i][j][1]; // Imaginary part
      }
    }
  } else {
    // Assume already interleaved column-major
    const arr = A as number[] | Float64Array;
    if (arr.length !== 2 * n * n) {
      throw new Error(`Expected ${2 * n * n} elements for ${n}×${n} complex matrix, got ${arr.length}`);
    }
    aData = arr instanceof Float64Array ? arr : new Float64Array(arr);
  }

  // Allocate WASM memory
  const aPtr = allocateDoubles(Module, overwriteA ? null : aData, 2 * n * n);
  const wPtr = allocateDoubles(Module, null, n); // Eigenvalues are real
  const infoPtr = allocateInts(Module, [0]);

  // Write data if not overwriting
  if (!overwriteA) {
    const baseIdx = aPtr >> 3;
    for (let i = 0; i < aData.length; i++) {
      Module.HEAPF64[baseIdx + i] = aData[i];
    }
  }

  // Parameter pointers
  const jobzPtr = allocateInts(Module, [computeVectors ? CHAR.V : CHAR.N]);
  const uploPtr = allocateInts(Module, [uplo === 'upper' ? CHAR.U : CHAR.L]);
  const nPtr = allocateInts(Module, [n]);
  const ldaPtr = allocateInts(Module, [n]);

  try {
    let info: number;

    if (algorithm === 'heev') {
      // Standard algorithm
      const lworkPtr = allocateInts(Module, [-1]);
      const workQueryPtr = allocateDoubles(Module, null, 2); // Complex query value
      const rworkPtr = allocateDoubles(Module, null, Math.max(1, 3 * n - 2));

      // Workspace query
      Module._zheev_(
        jobzPtr, uploPtr, nPtr, aPtr, ldaPtr, wPtr,
        workQueryPtr, lworkPtr, rworkPtr, infoPtr
      );

      info = readInt(Module, infoPtr);
      if (info !== 0) {
        freeAll(Module, [lworkPtr, workQueryPtr, rworkPtr]);
        return createErrorResult(n, info, 'ZHEEV');
      }

      // Complex workspace size is at first element (real part)
      const lwork = Math.ceil(Module.HEAPF64[workQueryPtr >> 3]);
      const workPtr = allocateDoubles(Module, null, 2 * lwork); // Complex workspace
      Module.HEAP32[lworkPtr >> 2] = lwork;

      // Compute eigenvalues and eigenvectors
      Module._zheev_(
        jobzPtr, uploPtr, nPtr, aPtr, ldaPtr, wPtr,
        workPtr, lworkPtr, rworkPtr, infoPtr
      );

      info = readInt(Module, infoPtr);
      freeAll(Module, [workPtr, lworkPtr, workQueryPtr, rworkPtr]);

    } else if (algorithm === 'heevd') {
      // Divide and conquer
      const lworkPtr = allocateInts(Module, [-1]);
      const lrworkPtr = allocateInts(Module, [-1]);
      const liworkPtr = allocateInts(Module, [-1]);
      const workQueryPtr = allocateDoubles(Module, null, 2);
      const rworkQueryPtr = allocateDoubles(Module, null, 1);
      const iworkQueryPtr = allocateInts(Module, [0]);

      // Workspace query
      Module._zheevd_(
        jobzPtr, uploPtr, nPtr, aPtr, ldaPtr, wPtr,
        workQueryPtr, lworkPtr, rworkQueryPtr, lrworkPtr,
        iworkQueryPtr, liworkPtr, infoPtr
      );

      info = readInt(Module, infoPtr);
      if (info !== 0) {
        freeAll(Module, [lworkPtr, lrworkPtr, liworkPtr, workQueryPtr, rworkQueryPtr, iworkQueryPtr]);
        return createErrorResult(n, info, 'ZHEEVD');
      }

      const lwork = Math.ceil(Module.HEAPF64[workQueryPtr >> 3]);
      const lrwork = Math.ceil(Module.HEAPF64[rworkQueryPtr >> 3]);
      const liwork = Module.HEAP32[iworkQueryPtr >> 2];

      const workPtr = allocateDoubles(Module, null, 2 * lwork);
      const rworkPtr = allocateDoubles(Module, null, lrwork);
      const iworkPtr = allocateInts(Module, null, liwork);

      Module.HEAP32[lworkPtr >> 2] = lwork;
      Module.HEAP32[lrworkPtr >> 2] = lrwork;
      Module.HEAP32[liworkPtr >> 2] = liwork;

      // Compute eigenvalues and eigenvectors
      Module._zheevd_(
        jobzPtr, uploPtr, nPtr, aPtr, ldaPtr, wPtr,
        workPtr, lworkPtr, rworkPtr, lrworkPtr,
        iworkPtr, liworkPtr, infoPtr
      );

      info = readInt(Module, infoPtr);
      freeAll(Module, [
        workPtr, rworkPtr, iworkPtr,
        lworkPtr, lrworkPtr, liworkPtr,
        workQueryPtr, rworkQueryPtr, iworkQueryPtr
      ]);

    } else {
      // RRR algorithm (heevr)
      const rangePtr = allocateInts(Module, [CHAR.A]);
      const vlPtr = allocateDoubles(Module, [0]);
      const vuPtr = allocateDoubles(Module, [0]);
      const ilPtr = allocateInts(Module, [0]);
      const iuPtr = allocateInts(Module, [0]);
      const abstolPtr = allocateDoubles(Module, [0]);
      const mPtr = allocateInts(Module, [0]);
      const zPtr = allocateDoubles(Module, null, computeVectors ? 2 * n * n : 2);
      const ldzPtr = allocateInts(Module, [computeVectors ? n : 1]);
      const isuppzPtr = allocateInts(Module, null, 2 * n);

      const lworkPtr = allocateInts(Module, [-1]);
      const lrworkPtr = allocateInts(Module, [-1]);
      const liworkPtr = allocateInts(Module, [-1]);
      const workQueryPtr = allocateDoubles(Module, null, 2);
      const rworkQueryPtr = allocateDoubles(Module, null, 1);
      const iworkQueryPtr = allocateInts(Module, [0]);

      // Workspace query
      Module._zheevr_(
        jobzPtr, rangePtr, uploPtr, nPtr, aPtr, ldaPtr,
        vlPtr, vuPtr, ilPtr, iuPtr, abstolPtr, mPtr,
        wPtr, zPtr, ldzPtr, isuppzPtr,
        workQueryPtr, lworkPtr, rworkQueryPtr, lrworkPtr,
        iworkQueryPtr, liworkPtr, infoPtr
      );

      info = readInt(Module, infoPtr);
      if (info !== 0) {
        freeAll(Module, [
          rangePtr, vlPtr, vuPtr, ilPtr, iuPtr, abstolPtr, mPtr,
          zPtr, ldzPtr, isuppzPtr,
          lworkPtr, lrworkPtr, liworkPtr, workQueryPtr, rworkQueryPtr, iworkQueryPtr
        ]);
        return createErrorResult(n, info, 'ZHEEVR');
      }

      const lwork = Math.ceil(Module.HEAPF64[workQueryPtr >> 3]);
      const lrwork = Math.ceil(Module.HEAPF64[rworkQueryPtr >> 3]);
      const liwork = Module.HEAP32[iworkQueryPtr >> 2];

      const workPtr = allocateDoubles(Module, null, 2 * lwork);
      const rworkPtr = allocateDoubles(Module, null, lrwork);
      const iworkPtr = allocateInts(Module, null, liwork);

      Module.HEAP32[lworkPtr >> 2] = lwork;
      Module.HEAP32[lrworkPtr >> 2] = lrwork;
      Module.HEAP32[liworkPtr >> 2] = liwork;

      // Compute eigenvalues and eigenvectors
      Module._zheevr_(
        jobzPtr, rangePtr, uploPtr, nPtr, aPtr, ldaPtr,
        vlPtr, vuPtr, ilPtr, iuPtr, abstolPtr, mPtr,
        wPtr, zPtr, ldzPtr, isuppzPtr,
        workPtr, lworkPtr, rworkPtr, lrworkPtr,
        iworkPtr, liworkPtr, infoPtr
      );

      info = readInt(Module, infoPtr);

      const result: EigHermitianResult = {
        values: readDoubles(Module, wPtr, n),
        n,
        info: 0,
        success: true,
        message: 'Success',
      };

      if (info !== 0) {
        freeAll(Module, [
          rangePtr, vlPtr, vuPtr, ilPtr, iuPtr, abstolPtr, mPtr,
          zPtr, ldzPtr, isuppzPtr,
          workPtr, rworkPtr, iworkPtr,
          lworkPtr, lrworkPtr, liworkPtr, workQueryPtr, rworkQueryPtr, iworkQueryPtr
        ]);
        return createErrorResult(n, info, 'ZHEEVR');
      }

      if (computeVectors) {
        result.vectors = readDoubles(Module, zPtr, 2 * n * n);
      }

      freeAll(Module, [
        rangePtr, vlPtr, vuPtr, ilPtr, iuPtr, abstolPtr, mPtr,
        zPtr, ldzPtr, isuppzPtr,
        workPtr, rworkPtr, iworkPtr,
        lworkPtr, lrworkPtr, liworkPtr, workQueryPtr, rworkQueryPtr, iworkQueryPtr
      ]);

      return result;
    }

    if (info !== 0) {
      return createErrorResult(n, info, algorithm === 'heev' ? 'ZHEEV' : 'ZHEEVD');
    }

    // Read results
    const result: EigHermitianResult = {
      values: readDoubles(Module, wPtr, n),
      n,
      info: 0,
      success: true,
      message: 'Success',
    };

    if (computeVectors) {
      result.vectors = readDoubles(Module, aPtr, 2 * n * n);
    }

    return result;
  } finally {
    freeAll(Module, [aPtr, wPtr, infoPtr, jobzPtr, uploPtr, nPtr, ldaPtr]);
  }
}

function createErrorResult(n: number, info: number, routine: string): EigHermitianResult {
  return {
    values: new Float64Array(0),
    n,
    info,
    success: false,
    message: getLapackErrorMessage(routine, info),
  };
}
