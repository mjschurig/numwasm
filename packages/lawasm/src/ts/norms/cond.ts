/**
 * Matrix Condition Number
 *
 * Compute the condition number of a matrix.
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
import { svdvals } from '../svd/svdvals.js';
import type { Matrix, CondOptions, CondResult } from './types.js';

/**
 * Compute the condition number of a matrix.
 *
 * The condition number κ(A) = ||A|| * ||A^(-1)||
 *
 * For ord=2, this is the ratio of largest to smallest singular value.
 * For ord=1 or 'inf', uses LU factorization and DGECON.
 *
 * @param A - Input square matrix (n × n)
 * @param options - Computation options
 * @returns The condition number
 *
 * @example
 * ```typescript
 * const A = [[1, 2], [3, 4]];
 *
 * // 2-norm condition number (default)
 * const { cond: kappa2 } = cond(A);
 *
 * // 1-norm condition number
 * const { cond: kappa1 } = cond(A, { ord: 1 });
 *
 * // Infinity-norm condition number
 * const { cond: kappaInf } = cond(A, { ord: 'inf' });
 * ```
 */
export function cond(A: Matrix, options: CondOptions = {}): CondResult {
  const { ord = 2 } = options;

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);
  if (m !== n) {
    throw new Error(`Matrix A must be square, got ${m}×${n}`);
  }

  // For 2-norm, use SVD
  if (ord === 2 || ord === '2') {
    const svdResult = svdvals(A);
    if (!svdResult.success) {
      return {
        cond: NaN,
        ord: 2,
        n,
        info: svdResult.info,
        success: false,
        message: svdResult.message,
      };
    }

    const s = svdResult.s;
    const sMax = s[0];
    const sMin = s[s.length - 1];

    // If smallest singular value is zero, matrix is singular
    if (sMin === 0) {
      return {
        cond: Infinity,
        ord: 2,
        n,
        info: 0,
        success: true,
        message: 'Matrix is singular',
      };
    }

    return {
      cond: sMax / sMin,
      ord: 2,
      n,
      info: 0,
      success: true,
      message: 'Success',
    };
  }

  // For 1-norm or infinity-norm, use LU factorization and DGECON
  const Module = getLAPACKModule();

  // Prepare matrix
  const aData = prepareMatrix(A);

  // Determine norm type
  const normOrd: 1 | 'inf' = ord === 'inf' ? 'inf' : 1;
  const normChar = normOrd === 'inf' ? CHAR.INF : CHAR.ONE;

  // First compute the norm of A using DLANGE
  const aPtr = allocateDoubles(Module, aData, n * n);
  const workNormPtr = allocateDoubles(Module, null, n); // Needed for infinity norm
  const normPtr = allocateInts(Module, [normChar]);
  const nPtr = allocateInts(Module, [n]);
  const ldaPtr = allocateInts(Module, [n]);

  let anorm: number;
  try {
    anorm = Module._dlange_(normPtr, nPtr, nPtr, aPtr, ldaPtr, workNormPtr);
  } finally {
    freeAll(Module, [workNormPtr, normPtr]);
  }

  // Now do LU factorization using DGETRF
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
          ord: normOrd,
          n,
          info,
          success: true,
          message: 'Matrix is singular',
        };
      }
      return {
        cond: NaN,
        ord: normOrd,
        n,
        info,
        success: false,
        message: getLapackErrorMessage('DGETRF', info),
      };
    }

    // Now estimate reciprocal condition number using DGECON
    const rcondPtr = allocateDoubles(Module, [0], 1);
    const workPtr = allocateDoubles(Module, null, 4 * n);
    const iworkPtr = allocateInts(Module, null, n);
    const normCharPtr = allocateInts(Module, [normChar]);

    try {
      Module._dgecon_(
        normCharPtr,
        nPtr,
        aPtr,
        ldaPtr,
        allocateDoubles(Module, [anorm], 1), // anorm as pointer
        rcondPtr,
        workPtr,
        iworkPtr,
        infoPtr
      );

      info = readInt(Module, infoPtr);
      if (info !== 0) {
        return {
          cond: NaN,
          ord: normOrd,
          n,
          info,
          success: false,
          message: getLapackErrorMessage('DGECON', info),
        };
      }

      const rcond = readDoubles(Module, rcondPtr, 1)[0];

      // Condition number is 1/rcond
      const condValue = rcond === 0 ? Infinity : 1 / rcond;

      return {
        cond: condValue,
        ord: normOrd,
        n,
        info: 0,
        success: true,
        message: 'Success',
      };
    } finally {
      freeAll(Module, [rcondPtr, workPtr, iworkPtr, normCharPtr]);
    }
  } finally {
    freeAll(Module, [aPtr, ipivPtr, infoPtr, mPtr, nPtr, ldaPtr]);
  }
}
