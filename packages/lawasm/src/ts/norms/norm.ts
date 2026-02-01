/**
 * Matrix Norm
 *
 * Compute various matrix norms.
 */

import { getLAPACKModule } from '../loader.js';
import {
  allocateDoubles,
  allocateInts,
  freeAll,
  prepareMatrix,
  getMatrixDimensions,
  CHAR,
} from '../helpers.js';
import { svdvals } from '../svd/svdvals.js';
import type { Matrix, NormOptions, NormResult, NormType } from './types.js';

/**
 * Compute a matrix norm.
 *
 * @param A - Input matrix (m × n)
 * @param options - Norm options
 * @returns The computed norm
 *
 * @example
 * ```typescript
 * const A = [[1, 2], [3, 4]];
 *
 * // Frobenius norm (default)
 * const { norm: fro } = norm(A);
 *
 * // 1-norm (max column sum)
 * const { norm: one } = norm(A, { ord: 1 });
 *
 * // 2-norm (spectral norm)
 * const { norm: two } = norm(A, { ord: 2 });
 *
 * // Infinity norm (max row sum)
 * const { norm: inf } = norm(A, { ord: 'inf' });
 * ```
 */
export function norm(A: Matrix, options: NormOptions = {}): NormResult {
  const { ord = 'fro' } = options;

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);

  // For 2-norm, use SVD
  if (ord === 2 || ord === '2') {
    const svdResult = svdvals(A);
    if (!svdResult.success) {
      return {
        norm: NaN,
        ord: 2,
        m,
        n,
        success: false,
        message: svdResult.message,
      };
    }
    return {
      norm: svdResult.s[0], // Largest singular value
      ord: 2,
      m,
      n,
      success: true,
      message: 'Success',
    };
  }

  // For other norms, use DLANGE
  const Module = getLAPACKModule();

  // Prepare matrix
  const aData = prepareMatrix(A);

  // Determine norm character
  let normChar: number;
  let normTypeResult: NormType;

  if (ord === 1 || ord === '1') {
    normChar = CHAR.ONE;
    normTypeResult = 1;
  } else if (ord === 'inf') {
    normChar = CHAR.INF;
    normTypeResult = 'inf';
  } else if (ord === 'fro') {
    normChar = CHAR.FRO;
    normTypeResult = 'fro';
  } else if (ord === 'max') {
    normChar = CHAR.MAX;
    normTypeResult = 'max';
  } else {
    normChar = CHAR.FRO;
    normTypeResult = 'fro';
  }

  // Allocate WASM memory
  const aPtr = allocateDoubles(Module, aData, m * n);
  const workPtr = allocateDoubles(Module, null, m); // Needed for infinity norm

  // Parameter pointers
  const normPtr = allocateInts(Module, [normChar]);
  const mPtr = allocateInts(Module, [m]);
  const nPtr = allocateInts(Module, [n]);
  const ldaPtr = allocateInts(Module, [m]);

  try {
    // Compute norm using DLANGE
    const normValue = Module._dlange_(normPtr, mPtr, nPtr, aPtr, ldaPtr, workPtr);

    return {
      norm: normValue,
      ord: normTypeResult,
      m,
      n,
      success: true,
      message: 'Success',
    };
  } finally {
    freeAll(Module, [aPtr, workPtr, normPtr, mPtr, nPtr, ldaPtr]);
  }
}
