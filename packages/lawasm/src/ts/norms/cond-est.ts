/**
 * Fast Condition Number Estimate
 *
 * Compute a fast estimate of the condition number using LU factorization.
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
import type { Matrix, CondEstOptions, CondEstResult } from './types.js';

/**
 * Compute a fast estimate of the condition number.
 *
 * Uses LU factorization (DGETRF) followed by condition estimation (DGECON).
 * This is faster than computing the exact condition number via SVD.
 *
 * @param A - Input square matrix (n × n)
 * @param options - Computation options
 * @returns Estimated condition number
 *
 * @example
 * ```typescript
 * const A = [[1, 2], [3, 4]];
 * const { cond, rcond } = condEst(A);
 * // cond is approximately 14.9 (exact 2-norm is ~14.93)
 * // rcond is approximately 0.067
 * ```
 */
export function condEst(A: Matrix, options: CondEstOptions = {}): CondEstResult {
  const Module = getLAPACKModule();

  const { norm: normType = '1' } = options;

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);
  if (m !== n) {
    throw new Error(`Matrix A must be square, got ${m}×${n}`);
  }

  // Prepare matrix
  const aData = prepareMatrix(A);

  // Determine norm character
  const normChar = normType === 'inf' ? CHAR.INF : CHAR.ONE;

  // First compute the norm of A using DLANGE
  const aPtr = allocateDoubles(Module, aData, n * n);
  const workNormPtr = allocateDoubles(Module, null, n);
  const normPtr = allocateInts(Module, [normChar]);
  const nPtr = allocateInts(Module, [n]);
  const ldaPtr = allocateInts(Module, [n]);

  let anorm: number;
  try {
    anorm = Module._dlange_(normPtr, nPtr, nPtr, aPtr, ldaPtr, workNormPtr);
  } finally {
    freeAll(Module, [workNormPtr, normPtr]);
  }

  // LU factorization using DGETRF
  const ipivPtr = allocateInts(Module, null, n);
  const infoPtr = allocateInts(Module, [0]);
  const mPtr = allocateInts(Module, [n]);

  try {
    Module._dgetrf_(mPtr, nPtr, aPtr, ldaPtr, ipivPtr, infoPtr);

    let info = readInt(Module, infoPtr);
    if (info !== 0) {
      if (info > 0) {
        // Matrix is singular
        return {
          cond: Infinity,
          rcond: 0,
          norm: normType,
          n,
          info,
          success: true,
          message: 'Matrix is singular',
        };
      }
      return {
        cond: NaN,
        rcond: NaN,
        norm: normType,
        n,
        info,
        success: false,
        message: getLapackErrorMessage('DGETRF', info),
      };
    }

    // Estimate reciprocal condition number using DGECON
    const rcondPtr = allocateDoubles(Module, [0], 1);
    const workPtr = allocateDoubles(Module, null, 4 * n);
    const iworkPtr = allocateInts(Module, null, n);
    const normCharPtr = allocateInts(Module, [normChar]);
    const anormPtr = allocateDoubles(Module, [anorm], 1);

    try {
      Module._dgecon_(
        normCharPtr,
        nPtr,
        aPtr,
        ldaPtr,
        anormPtr,
        rcondPtr,
        workPtr,
        iworkPtr,
        infoPtr
      );

      info = readInt(Module, infoPtr);
      if (info !== 0) {
        return {
          cond: NaN,
          rcond: NaN,
          norm: normType,
          n,
          info,
          success: false,
          message: getLapackErrorMessage('DGECON', info),
        };
      }

      const rcond = readDoubles(Module, rcondPtr, 1)[0];
      const condValue = rcond === 0 ? Infinity : 1 / rcond;

      return {
        cond: condValue,
        rcond,
        norm: normType,
        n,
        info: 0,
        success: true,
        message: 'Success',
      };
    } finally {
      freeAll(Module, [rcondPtr, workPtr, iworkPtr, normCharPtr, anormPtr]);
    }
  } finally {
    freeAll(Module, [aPtr, ipivPtr, infoPtr, mPtr, nPtr, ldaPtr]);
  }
}
