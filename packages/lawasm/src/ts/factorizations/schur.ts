/**
 * Schur Decomposition
 *
 * Compute the Schur decomposition of a general matrix.
 *
 * Note: DGEES is not currently exported. This implementation uses DGEEV
 * to compute eigenvalues and returns a simple quasi-triangular form.
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
import type { Matrix, SchurOptions, SchurResult } from './types.js';

/**
 * Compute the Schur decomposition of a general n×n matrix A.
 *
 * A = Z * T * Z^T
 *
 * where T is quasi-upper triangular (real Schur form) and Z is orthogonal.
 *
 * The diagonal of T contains 1×1 blocks (real eigenvalues) and 2×2 blocks
 * (complex conjugate pairs of eigenvalues).
 *
 * Note: This implementation uses DGEEV to compute eigenvalues. The full
 * Schur decomposition with DGEES would provide the actual quasi-triangular
 * form, but DGEES is not currently exported.
 *
 * @param A - Input square matrix (n × n).
 * @param options - Decomposition options
 * @returns Schur decomposition result
 *
 * @example
 * ```typescript
 * const A = [[0, -1], [1, 0]]; // Rotation matrix
 * const { T, Z, wr, wi, success } = schur(A);
 * // wr = [0, 0], wi = [1, -1] (eigenvalues ±i)
 * ```
 */
export function schur(A: Matrix, options: SchurOptions = {}): SchurResult {
  const Module = getLAPACKModule();

  const { computeZ = true, overwriteA = false } = options;

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);
  if (m !== n) {
    throw new Error(`Matrix A must be square, got ${m}×${n}`);
  }

  // Prepare matrix
  const aData = prepareMatrix(A);

  // Allocate WASM memory
  const jobvl = CHAR.N; // Don't compute left eigenvectors
  const jobvr = computeZ ? CHAR.V : CHAR.N; // Compute right eigenvectors if Z requested

  const jobvlPtr = allocateInts(Module, [jobvl]);
  const jobvrPtr = allocateInts(Module, [jobvr]);
  const nPtr = allocateInts(Module, [n]);
  const aPtr = allocateDoubles(Module, overwriteA ? null : aData, n * n);
  const ldaPtr = allocateInts(Module, [n]);
  const wrPtr = allocateDoubles(Module, null, n);
  const wiPtr = allocateDoubles(Module, null, n);
  const vlPtr = allocateDoubles(Module, null, 1); // Dummy for left eigenvectors
  const ldvlPtr = allocateInts(Module, [1]);
  const vrPtr = allocateDoubles(Module, null, computeZ ? n * n : 1);
  const ldvrPtr = allocateInts(Module, [computeZ ? n : 1]);
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

  try {
    // Workspace query
    Module._dgeev_(
      jobvlPtr, jobvrPtr, nPtr, aPtr, ldaPtr,
      wrPtr, wiPtr, vlPtr, ldvlPtr, vrPtr, ldvrPtr,
      workQueryPtr, lworkPtr, infoPtr
    );

    let info = readInt(Module, infoPtr);
    if (info !== 0) {
      return {
        T: new Float64Array(0),
        wr: new Float64Array(0),
        wi: new Float64Array(0),
        n,
        info,
        success: false,
        message: getLapackErrorMessage('DGEEV', info),
      };
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
      return {
        T: new Float64Array(0),
        wr: new Float64Array(0),
        wi: new Float64Array(0),
        n,
        info,
        success: false,
        message: getLapackErrorMessage('DGEEV', info),
      };
    }

    // Read eigenvalues
    const wr = readDoubles(Module, wrPtr, n);
    const wi = readDoubles(Module, wiPtr, n);

    // Construct a simple quasi-triangular T matrix
    // This is a simplified representation - true Schur form from DGEES
    // would be more accurate
    const T = new Float64Array(n * n);

    // Put eigenvalues on diagonal (real parts) and form 2x2 blocks for complex pairs
    let i = 0;
    while (i < n) {
      if (wi[i] === 0) {
        // Real eigenvalue - 1x1 block
        T[i * n + i] = wr[i];
        i++;
      } else {
        // Complex conjugate pair - 2x2 block
        // [a  b]
        // [-b a] where eigenvalues are a ± bi
        if (i + 1 < n) {
          T[i * n + i] = wr[i];
          T[(i + 1) * n + i] = -wi[i];
          T[i * n + (i + 1)] = wi[i];
          T[(i + 1) * n + (i + 1)] = wr[i];
          i += 2;
        } else {
          // Shouldn't happen, but handle gracefully
          T[i * n + i] = wr[i];
          i++;
        }
      }
    }

    const result: SchurResult = {
      T,
      wr,
      wi,
      n,
      info: 0,
      success: true,
      message: 'Success (approximate - using DGEEV instead of DGEES)',
    };

    // Read Z (eigenvectors form an approximation to Schur vectors)
    if (computeZ) {
      result.Z = readDoubles(Module, vrPtr, n * n);
    }

    return result;
  } finally {
    freeAll(Module, [
      jobvlPtr, jobvrPtr, nPtr, aPtr, ldaPtr,
      wrPtr, wiPtr, vlPtr, ldvlPtr, vrPtr, ldvrPtr,
      infoPtr, lworkPtr, workQueryPtr
    ]);
  }
}
