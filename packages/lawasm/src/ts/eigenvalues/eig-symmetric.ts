/**
 * Symmetric Eigenvalue Problem
 *
 * Compute eigenvalues and eigenvectors of a real symmetric matrix.
 * All eigenvalues of a symmetric matrix are real.
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
import type { Matrix, EigSymmetricOptions, EigSymmetricResult } from './types.js';

/**
 * Compute eigenvalues and optionally eigenvectors of a real symmetric n×n matrix A.
 *
 * For a symmetric matrix, all eigenvalues are real and the eigenvectors
 * form an orthonormal basis.
 *
 * @param A - Input symmetric matrix (n × n). Only the specified triangle is accessed.
 * @param options - Computation options
 * @returns Eigenvalues (in ascending order) and optionally eigenvectors
 *
 * @example
 * ```typescript
 * const A = [[2, 1], [1, 2]];
 * const { values, vectors, success } = eigSymmetric(A);
 * // values = [1, 3] (real eigenvalues in ascending order)
 * // vectors contains orthonormal eigenvectors
 * ```
 *
 * @example
 * ```typescript
 * // Using divide-and-conquer for large matrices
 * const result = eigSymmetric(A, { algorithm: 'syevd' });
 * ```
 */
export function eigSymmetric(A: Matrix, options: EigSymmetricOptions = {}): EigSymmetricResult {
  const Module = getLAPACKModule();

  const {
    computeVectors = true,
    algorithm = 'syevd',
    uplo = 'lower',
    overwriteA = false,
  } = options;

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);
  if (m !== n) {
    throw new Error(`Matrix A must be square, got ${m}×${n}`);
  }

  // Prepare matrix
  const aData = prepareMatrix(A);

  // Allocate WASM memory
  const aPtr = allocateDoubles(Module, overwriteA ? null : aData, n * n);
  const wPtr = allocateDoubles(Module, null, n);
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

    if (algorithm === 'syev') {
      // Standard algorithm
      const lworkPtr = allocateInts(Module, [-1]);
      const workQueryPtr = allocateDoubles(Module, null, 1);

      // Workspace query
      Module._dsyev_(
        jobzPtr, uploPtr, nPtr, aPtr, ldaPtr, wPtr,
        workQueryPtr, lworkPtr, infoPtr
      );

      info = readInt(Module, infoPtr);
      if (info !== 0) {
        freeAll(Module, [lworkPtr, workQueryPtr]);
        return createErrorResult(n, info, 'DSYEV');
      }

      const lwork = Math.ceil(Module.HEAPF64[workQueryPtr >> 3]);
      const workPtr = allocateDoubles(Module, null, lwork);
      Module.HEAP32[lworkPtr >> 2] = lwork;

      // Compute eigenvalues and eigenvectors
      Module._dsyev_(
        jobzPtr, uploPtr, nPtr, aPtr, ldaPtr, wPtr,
        workPtr, lworkPtr, infoPtr
      );

      info = readInt(Module, infoPtr);
      freeAll(Module, [workPtr, lworkPtr, workQueryPtr]);

    } else if (algorithm === 'syevd') {
      // Divide and conquer - faster for large matrices
      const lworkPtr = allocateInts(Module, [-1]);
      const liworkPtr = allocateInts(Module, [-1]);
      const workQueryPtr = allocateDoubles(Module, null, 1);
      const iworkQueryPtr = allocateInts(Module, [0]);

      // Workspace query
      Module._dsyevd_(
        jobzPtr, uploPtr, nPtr, aPtr, ldaPtr, wPtr,
        workQueryPtr, lworkPtr, iworkQueryPtr, liworkPtr, infoPtr
      );

      info = readInt(Module, infoPtr);
      if (info !== 0) {
        freeAll(Module, [lworkPtr, liworkPtr, workQueryPtr, iworkQueryPtr]);
        return createErrorResult(n, info, 'DSYEVD');
      }

      const lwork = Math.ceil(Module.HEAPF64[workQueryPtr >> 3]);
      const liwork = Module.HEAP32[iworkQueryPtr >> 2];
      const workPtr = allocateDoubles(Module, null, lwork);
      const iworkPtr = allocateInts(Module, null, liwork);
      Module.HEAP32[lworkPtr >> 2] = lwork;
      Module.HEAP32[liworkPtr >> 2] = liwork;

      // Compute eigenvalues and eigenvectors
      Module._dsyevd_(
        jobzPtr, uploPtr, nPtr, aPtr, ldaPtr, wPtr,
        workPtr, lworkPtr, iworkPtr, liworkPtr, infoPtr
      );

      info = readInt(Module, infoPtr);
      freeAll(Module, [workPtr, iworkPtr, lworkPtr, liworkPtr, workQueryPtr, iworkQueryPtr]);

    } else {
      // RRR algorithm (syevr) - fastest
      const rangePtr = allocateInts(Module, [CHAR.A]); // All eigenvalues
      const vlPtr = allocateDoubles(Module, [0]);
      const vuPtr = allocateDoubles(Module, [0]);
      const ilPtr = allocateInts(Module, [0]);
      const iuPtr = allocateInts(Module, [0]);
      const abstolPtr = allocateDoubles(Module, [0]); // Use default tolerance
      const mPtr = allocateInts(Module, [0]);
      const zPtr = allocateDoubles(Module, null, computeVectors ? n * n : 1);
      const ldzPtr = allocateInts(Module, [computeVectors ? n : 1]);
      const isuppzPtr = allocateInts(Module, null, 2 * n);
      const lworkPtr = allocateInts(Module, [-1]);
      const liworkPtr = allocateInts(Module, [-1]);
      const workQueryPtr = allocateDoubles(Module, null, 1);
      const iworkQueryPtr = allocateInts(Module, [0]);

      // Workspace query
      Module._dsyevr_(
        jobzPtr, rangePtr, uploPtr, nPtr, aPtr, ldaPtr,
        vlPtr, vuPtr, ilPtr, iuPtr, abstolPtr, mPtr,
        wPtr, zPtr, ldzPtr, isuppzPtr,
        workQueryPtr, lworkPtr, iworkQueryPtr, liworkPtr, infoPtr
      );

      info = readInt(Module, infoPtr);
      if (info !== 0) {
        freeAll(Module, [
          rangePtr, vlPtr, vuPtr, ilPtr, iuPtr, abstolPtr, mPtr,
          zPtr, ldzPtr, isuppzPtr, lworkPtr, liworkPtr, workQueryPtr, iworkQueryPtr
        ]);
        return createErrorResult(n, info, 'DSYEVR');
      }

      const lwork = Math.ceil(Module.HEAPF64[workQueryPtr >> 3]);
      const liwork = Module.HEAP32[iworkQueryPtr >> 2];
      const workPtr = allocateDoubles(Module, null, lwork);
      const iworkPtr = allocateInts(Module, null, liwork);
      Module.HEAP32[lworkPtr >> 2] = lwork;
      Module.HEAP32[liworkPtr >> 2] = liwork;

      // Compute eigenvalues and eigenvectors
      Module._dsyevr_(
        jobzPtr, rangePtr, uploPtr, nPtr, aPtr, ldaPtr,
        vlPtr, vuPtr, ilPtr, iuPtr, abstolPtr, mPtr,
        wPtr, zPtr, ldzPtr, isuppzPtr,
        workPtr, lworkPtr, iworkPtr, liworkPtr, infoPtr
      );

      info = readInt(Module, infoPtr);

      // Read eigenvectors from Z (not A)
      const result: EigSymmetricResult = {
        values: readDoubles(Module, wPtr, n),
        n,
        info: 0,
        success: true,
        message: 'Success',
      };

      if (info !== 0) {
        freeAll(Module, [
          rangePtr, vlPtr, vuPtr, ilPtr, iuPtr, abstolPtr, mPtr,
          zPtr, ldzPtr, isuppzPtr, workPtr, iworkPtr, lworkPtr, liworkPtr, workQueryPtr, iworkQueryPtr
        ]);
        return createErrorResult(n, info, 'DSYEVR');
      }

      if (computeVectors) {
        result.vectors = readDoubles(Module, zPtr, n * n);
      }

      freeAll(Module, [
        rangePtr, vlPtr, vuPtr, ilPtr, iuPtr, abstolPtr, mPtr,
        zPtr, ldzPtr, isuppzPtr, workPtr, iworkPtr, lworkPtr, liworkPtr, workQueryPtr, iworkQueryPtr
      ]);

      return result;
    }

    if (info !== 0) {
      return createErrorResult(n, info, algorithm === 'syev' ? 'DSYEV' : 'DSYEVD');
    }

    // Read results
    const result: EigSymmetricResult = {
      values: readDoubles(Module, wPtr, n),
      n,
      info: 0,
      success: true,
      message: 'Success',
    };

    if (computeVectors) {
      result.vectors = readDoubles(Module, aPtr, n * n);
    }

    return result;
  } finally {
    freeAll(Module, [aPtr, wPtr, infoPtr, jobzPtr, uploPtr, nPtr, ldaPtr]);
  }
}

function createErrorResult(n: number, info: number, routine: string): EigSymmetricResult {
  return {
    values: new Float64Array(0),
    n,
    info,
    success: false,
    message: getLapackErrorMessage(routine, info),
  };
}
