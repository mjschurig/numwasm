/**
 * General Eigenvalue Problem
 *
 * Compute eigenvalues and eigenvectors of a general (non-symmetric) matrix.
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
import type { Matrix, EigOptions, EigResult, Complex } from './types.js';

/**
 * Compute eigenvalues and optionally eigenvectors of a general n×n matrix A.
 *
 * For a general real matrix, eigenvalues may be complex. Complex eigenvalues
 * appear in conjugate pairs.
 *
 * @param A - Input square matrix (n × n). Can be 2D array (row-major) or 1D array (column-major).
 * @param options - Computation options
 * @returns Eigenvalues and optionally eigenvectors
 *
 * @example
 * ```typescript
 * const A = [[1, 2], [3, 4]];
 * const { values, vectors, success } = eig(A);
 * // values is an array of complex numbers
 * // vectors contains the eigenvectors (if computed)
 * ```
 *
 * @example
 * ```typescript
 * // Eigenvalues only (faster)
 * const { values } = eig(A, { computeVectors: false });
 * ```
 */
export function eig(A: Matrix, options: EigOptions = {}): EigResult {
  const Module = getLAPACKModule();

  const {
    computeVectors = true,
    computeLeftVectors = false,
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
  const wrPtr = allocateDoubles(Module, null, n);
  const wiPtr = allocateDoubles(Module, null, n);
  const vlPtr = allocateDoubles(Module, null, computeLeftVectors ? n * n : 1);
  const vrPtr = allocateDoubles(Module, null, computeVectors ? n * n : 1);
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

  // Additional parameter pointers
  const jobvlPtr = allocateInts(Module, [computeLeftVectors ? CHAR.V : CHAR.N]);
  const jobvrPtr = allocateInts(Module, [computeVectors ? CHAR.V : CHAR.N]);
  const nPtr = allocateInts(Module, [n]);
  const ldaPtr = allocateInts(Module, [n]);
  const ldvlPtr = allocateInts(Module, [computeLeftVectors ? n : 1]);
  const ldvrPtr = allocateInts(Module, [computeVectors ? n : 1]);

  try {
    // Workspace query
    Module._dgeev_(
      jobvlPtr, jobvrPtr, nPtr, aPtr, ldaPtr,
      wrPtr, wiPtr, vlPtr, ldvlPtr, vrPtr, ldvrPtr,
      workQueryPtr, lworkPtr, infoPtr
    );

    let info = readInt(Module, infoPtr);
    if (info !== 0) {
      return createErrorResult(n, info, 'DGEEV');
    }

    const lwork = Math.ceil(Module.HEAPF64[workQueryPtr >> 3]);
    const workPtr = allocateDoubles(Module, null, lwork);
    Module.HEAP32[lworkPtr >> 2] = lwork;

    // Compute eigenvalues and eigenvectors
    Module._dgeev_(
      jobvlPtr, jobvrPtr, nPtr, aPtr, ldaPtr,
      wrPtr, wiPtr, vlPtr, ldvlPtr, vrPtr, ldvrPtr,
      workPtr, lworkPtr, infoPtr
    );

    info = readInt(Module, infoPtr);
    freeAll(Module, [workPtr]);

    if (info !== 0) {
      return createErrorResult(n, info, 'DGEEV');
    }

    // Read results
    const wr = readDoubles(Module, wrPtr, n);
    const wi = readDoubles(Module, wiPtr, n);

    // Convert to complex array
    const values: Complex[] = [];
    for (let i = 0; i < n; i++) {
      values.push({ re: wr[i], im: wi[i] });
    }

    const result: EigResult = {
      values,
      wr,
      wi,
      n,
      info: 0,
      success: true,
      message: 'Success',
    };

    if (computeVectors) {
      result.vectors = readDoubles(Module, vrPtr, n * n);
    }

    if (computeLeftVectors) {
      result.leftVectors = readDoubles(Module, vlPtr, n * n);
    }

    return result;
  } finally {
    freeAll(Module, [
      aPtr, wrPtr, wiPtr, vlPtr, vrPtr,
      infoPtr, lworkPtr, workQueryPtr,
      jobvlPtr, jobvrPtr, nPtr, ldaPtr, ldvlPtr, ldvrPtr
    ]);
  }
}

function createErrorResult(n: number, info: number, routine: string): EigResult {
  return {
    values: [],
    wr: new Float64Array(0),
    wi: new Float64Array(0),
    n,
    info,
    success: false,
    message: getLapackErrorMessage(routine, info),
  };
}
