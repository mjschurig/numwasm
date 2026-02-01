/**
 * Eigenvalues Only
 *
 * Compute only the eigenvalues of a general matrix (faster than full eigendecomposition).
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
import type { Matrix, EigvalsOptions, EigvalsResult, Complex } from './types.js';

/**
 * Compute only the eigenvalues of a general n×n matrix A.
 *
 * This is faster than computing the full eigendecomposition when only
 * eigenvalues are needed.
 *
 * @param A - Input square matrix (n × n). Can be 2D array (row-major) or 1D array (column-major).
 * @param options - Computation options
 * @returns Eigenvalues (possibly complex)
 *
 * @example
 * ```typescript
 * const A = [[1, 2], [3, 4]];
 * const { values, wr, wi, success } = eigvals(A);
 * // values is an array of complex numbers
 * // wr contains real parts, wi contains imaginary parts
 * ```
 */
export function eigvals(A: Matrix, options: EigvalsOptions = {}): EigvalsResult {
  const Module = getLAPACKModule();

  const { overwriteA = false } = options;

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
  const vlPtr = allocateDoubles(Module, null, 1); // Dummy
  const vrPtr = allocateDoubles(Module, null, 1); // Dummy
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
  const jobvlPtr = allocateInts(Module, [CHAR.N]); // No left eigenvectors
  const jobvrPtr = allocateInts(Module, [CHAR.N]); // No right eigenvectors
  const nPtr = allocateInts(Module, [n]);
  const ldaPtr = allocateInts(Module, [n]);
  const ldvlPtr = allocateInts(Module, [1]);
  const ldvrPtr = allocateInts(Module, [1]);

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

    // Compute eigenvalues
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

    return {
      values,
      wr,
      wi,
      n,
      info: 0,
      success: true,
      message: 'Success',
    };
  } finally {
    freeAll(Module, [
      aPtr, wrPtr, wiPtr, vlPtr, vrPtr,
      infoPtr, lworkPtr, workQueryPtr,
      jobvlPtr, jobvrPtr, nPtr, ldaPtr, ldvlPtr, ldvrPtr
    ]);
  }
}

function createErrorResult(n: number, info: number, routine: string): EigvalsResult {
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
