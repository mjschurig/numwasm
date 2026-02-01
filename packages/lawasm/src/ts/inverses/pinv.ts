/**
 * Moore-Penrose Pseudoinverse
 *
 * Compute the pseudoinverse of a matrix using SVD.
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
import type { Matrix, PinvOptions, PinvResult } from './types.js';

/**
 * Compute the Moore-Penrose pseudoinverse of an m×n matrix A.
 *
 * The pseudoinverse A^+ satisfies:
 * - A * A^+ * A = A
 * - A^+ * A * A^+ = A^+
 * - (A * A^+)^T = A * A^+
 * - (A^+ * A)^T = A^+ * A
 *
 * For full-rank matrices:
 * - If m >= n (overdetermined): A^+ = (A^T*A)^(-1) * A^T
 * - If m < n (underdetermined): A^+ = A^T * (A*A^T)^(-1)
 *
 * This implementation uses SVD: A = U * S * V^T => A^+ = V * S^+ * U^T
 * where S^+ is the pseudoinverse of the diagonal singular value matrix.
 *
 * @param A - Input matrix (m × n)
 * @param options - Computation options
 * @returns Pseudoinverse (n × m)
 *
 * @example
 * ```typescript
 * const A = [[1, 2], [3, 4], [5, 6]];
 * const { pinv, rank } = pinv(A);
 * // pinv is the 2×3 pseudoinverse
 * ```
 */
export function pinv(A: Matrix, options: PinvOptions = {}): PinvResult {
  const Module = getLAPACKModule();

  const { rcond = -1, overwriteA = false } = options;

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);
  const k = Math.min(m, n);

  // Prepare matrix
  const aData = prepareMatrix(A);

  // Allocate WASM memory for SVD
  const aPtr = allocateDoubles(Module, overwriteA ? null : aData, m * n);
  const sPtr = allocateDoubles(Module, null, k);
  const uPtr = allocateDoubles(Module, null, m * m); // Full U
  const vtPtr = allocateDoubles(Module, null, n * n); // Full Vt
  const infoPtr = allocateInts(Module, [0]);
  const lworkPtr = allocateInts(Module, [-1]);
  const workQueryPtr = allocateDoubles(Module, null, 1);
  const iworkPtr = allocateInts(Module, null, 8 * k);

  // Write data if not overwriting
  if (!overwriteA) {
    const baseIdx = aPtr >> 3;
    for (let i = 0; i < aData.length; i++) {
      Module.HEAPF64[baseIdx + i] = aData[i];
    }
  }

  // Parameter pointers for DGESDD
  const jobzPtr = allocateInts(Module, [CHAR.A]); // All of U and Vt
  const mPtr = allocateInts(Module, [m]);
  const nPtr = allocateInts(Module, [n]);
  const ldaPtr = allocateInts(Module, [m]);
  const lduPtr = allocateInts(Module, [m]);
  const ldvtPtr = allocateInts(Module, [n]);

  try {
    // Workspace query for DGESDD
    Module._dgesdd_(
      jobzPtr, mPtr, nPtr, aPtr, ldaPtr,
      sPtr, uPtr, lduPtr, vtPtr, ldvtPtr,
      workQueryPtr, lworkPtr, iworkPtr, infoPtr
    );

    let info = readInt(Module, infoPtr);
    if (info !== 0) {
      return createErrorResult(m, n, info, 'DGESDD');
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
    freeAll(Module, [workPtr]);

    if (info !== 0) {
      return createErrorResult(m, n, info, 'DGESDD');
    }

    // Read SVD results
    const s = readDoubles(Module, sPtr, k);
    const U = readDoubles(Module, uPtr, m * m);
    const Vt = readDoubles(Module, vtPtr, n * n);

    // Determine tolerance for rank
    const eps = Number.EPSILON;
    const tol = rcond < 0 ? Math.max(m, n) * eps * s[0] : rcond * s[0];

    // Compute effective rank
    let rank = 0;
    for (let i = 0; i < k; i++) {
      if (s[i] > tol) {
        rank++;
      }
    }

    // Compute pseudoinverse: A^+ = V * S^+ * U^T
    // S^+ is k×k diagonal with 1/s[i] for s[i] > tol, 0 otherwise
    // Result is n×m

    // First compute V * S^+
    // V is n×n (columns of Vt transposed)
    // S^+ is k×k diagonal
    // VS_plus is n×k
    const VS_plus = new Float64Array(n * k);
    for (let j = 0; j < k; j++) {
      const sinv = s[j] > tol ? 1.0 / s[j] : 0.0;
      for (let i = 0; i < n; i++) {
        // V[i,j] = Vt[j,i] = Vt[j*n + i]
        // Actually Vt is stored column-major, so Vt[i,j] is at Vt[j*n + i]
        // So V[i,j] = Vt[j,i] which is at index i*n + j in Vt (row-major)
        // But Vt is column-major: Vt[col][row] at col*n + row
        // Vt[i][j] at j*n + i
        // V = Vt^T, so V[i][j] = Vt[j][i] at i*n + j
        VS_plus[j * n + i] = Vt[i * n + j] * sinv;
      }
    }

    // Now compute (V * S^+) * U^T = pinv (n×m)
    // U is m×m, U^T is m×m
    // (V * S^+) is n×k, U^T restricted to first k rows is k×m
    // Result is n×m
    const pinvMatrix = new Float64Array(n * m);
    for (let j = 0; j < m; j++) {
      for (let i = 0; i < n; i++) {
        let sum = 0;
        for (let l = 0; l < k; l++) {
          // VS_plus[i,l] at l*n + i
          // U^T[l,j] = U[j,l] at l*m + j
          sum += VS_plus[l * n + i] * U[l * m + j];
        }
        // pinv[i,j] at j*n + i (column-major)
        pinvMatrix[j * n + i] = sum;
      }
    }

    return {
      pinv: pinvMatrix,
      s,
      rank,
      m,
      n,
      info: 0,
      success: true,
      message: 'Success',
    };
  } finally {
    freeAll(Module, [
      aPtr, sPtr, uPtr, vtPtr, infoPtr, lworkPtr, workQueryPtr, iworkPtr,
      jobzPtr, mPtr, nPtr, ldaPtr, lduPtr, ldvtPtr
    ]);
  }
}

function createErrorResult(m: number, n: number, info: number, routine: string): PinvResult {
  return {
    pinv: new Float64Array(0),
    s: new Float64Array(0),
    rank: 0,
    m,
    n,
    info,
    success: false,
    message: getLapackErrorMessage(routine, info),
  };
}
