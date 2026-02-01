/**
 * Rank-Revealing QR Decomposition
 *
 * QR decomposition with column pivoting that reveals numerical rank.
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
} from '../helpers.js';
import type { Matrix, RRQROptions, RRQRResult } from './types.js';

/**
 * Compute rank-revealing QR decomposition with column pivoting.
 *
 * Computes A*P = Q*R where P is a permutation matrix chosen to
 * reveal the numerical rank of A.
 *
 * Uses DGEQP3 from LAPACK.
 *
 * @param A - Input matrix (m × n)
 * @param options - RRQR options
 * @returns Rank-revealing QR decomposition
 *
 * @example
 * ```typescript
 * // Rank-deficient matrix
 * const A = [[1, 2, 3], [2, 4, 6], [1, 1, 1]];
 * const { Q, R, P, rank } = rrqr(A);
 * // rank = 2 (rows 1 and 2 are linearly dependent)
 * // A[:, P] = Q * R
 * ```
 */
export function rrqr(A: Matrix, options: RRQROptions = {}): RRQRResult {
  const Module = getLAPACKModule();

  const [m, n] = getMatrixDimensions(A);

  if (m === 0 || n === 0) {
    return {
      Q: new Float64Array(0),
      R: new Float64Array(0),
      P: new Int32Array(0),
      rank: 0,
      m,
      n,
      info: 0,
      success: true,
      message: 'Success',
    };
  }

  const k = Math.min(m, n);

  // Prepare matrix
  const aData = prepareMatrix(A);

  // Allocate WASM memory
  const aPtr = allocateDoubles(Module, aData, m * n);
  const tauPtr = allocateDoubles(Module, null, k);
  const jpvtPtr = allocateInts(Module, new Array(n).fill(0), n); // Initialize to 0 for free pivoting
  const infoPtr = allocateInts(Module, [0]);

  // Parameter pointers
  const mPtr = allocateInts(Module, [m]);
  const nPtr = allocateInts(Module, [n]);
  const ldaPtr = allocateInts(Module, [m]);

  // Workspace query
  const workQuery = allocateDoubles(Module, [0], 1);
  const lworkQuery = allocateInts(Module, [-1]);

  try {
    // Query optimal workspace size
    Module._dgeqp3_(
      mPtr,
      nPtr,
      aPtr,
      ldaPtr,
      jpvtPtr,
      tauPtr,
      workQuery,
      lworkQuery,
      infoPtr
    );

    const optimalLwork = Math.max(1, Math.ceil(readDoubles(Module, workQuery, 1)[0]));
    freeAll(Module, [workQuery, lworkQuery]);

    // Allocate workspace
    const workPtr = allocateDoubles(Module, null, optimalLwork);
    const lworkPtr = allocateInts(Module, [optimalLwork]);

    // Compute QR with column pivoting
    Module._dgeqp3_(
      mPtr,
      nPtr,
      aPtr,
      ldaPtr,
      jpvtPtr,
      tauPtr,
      workPtr,
      lworkPtr,
      infoPtr
    );

    const info = readInt(Module, infoPtr);
    freeAll(Module, [workPtr, lworkPtr]);

    if (info !== 0) {
      return {
        Q: new Float64Array(0),
        R: new Float64Array(0),
        P: new Int32Array(0),
        rank: 0,
        m,
        n,
        info,
        success: false,
        message: getLapackErrorMessage('DGEQP3', info),
      };
    }

    // Read results
    const qrData = readDoubles(Module, aPtr, m * n);
    // tau contains Householder reflector scalars (used for reconstructing Q)
    const _tau = readDoubles(Module, tauPtr, k);
    void _tau; // Available for future Q reconstruction

    // Read pivot indices (1-based from LAPACK)
    const jpvt = new Int32Array(n);
    const baseIdx = jpvtPtr >> 2;
    for (let i = 0; i < n; i++) {
      jpvt[i] = Module.HEAP32[baseIdx + i] - 1; // Convert to 0-based
    }

    // Determine numerical rank from R diagonal
    const eps = 2.220446049250313e-16;

    // Compute ||A|| using Frobenius norm
    let normA = 0;
    for (let i = 0; i < m * n; i++) {
      normA += aData[i] * aData[i];
    }
    normA = Math.sqrt(normA);

    const defaultTol = Math.max(m, n) * eps * normA;
    const tol = options.tol !== undefined ? options.tol : defaultTol;

    let rank = 0;
    for (let i = 0; i < k; i++) {
      // R[i,i] is at qrData[i + i*m] in column-major
      if (Math.abs(qrData[i + i * m]) > tol) {
        rank++;
      } else {
        break; // R is sorted, so remaining diagonals are smaller
      }
    }

    // Extract R (upper triangular part, k × n)
    const R = new Float64Array(k * n);
    for (let i = 0; i < k; i++) {
      for (let j = i; j < n; j++) {
        R[i + j * k] = qrData[i + j * m];
      }
    }

    // Build Q from Householder reflectors using DORGQR
    // Q is m × k
    const qPtr = allocateDoubles(Module, qrData, m * n);
    const kPtr = allocateInts(Module, [k]);

    // Query workspace for DORGQR
    const workQuery2 = allocateDoubles(Module, [0], 1);
    const lworkQuery2 = allocateInts(Module, [-1]);

    Module._dorgqr_(
      mPtr,
      kPtr,
      kPtr,
      qPtr,
      ldaPtr,
      tauPtr,
      workQuery2,
      lworkQuery2,
      infoPtr
    );

    const optimalLwork2 = Math.max(1, Math.ceil(readDoubles(Module, workQuery2, 1)[0]));
    freeAll(Module, [workQuery2, lworkQuery2]);

    const workPtr2 = allocateDoubles(Module, null, optimalLwork2);
    const lworkPtr2 = allocateInts(Module, [optimalLwork2]);

    Module._dorgqr_(
      mPtr,
      kPtr,
      kPtr,
      qPtr,
      ldaPtr,
      tauPtr,
      workPtr2,
      lworkPtr2,
      infoPtr
    );

    freeAll(Module, [workPtr2, lworkPtr2]);

    // Read Q (m × k)
    const qData = readDoubles(Module, qPtr, m * n);
    const Q = new Float64Array(m * k);
    for (let i = 0; i < m; i++) {
      for (let j = 0; j < k; j++) {
        Q[i + j * m] = qData[i + j * m];
      }
    }

    freeAll(Module, [qPtr, kPtr]);

    return {
      Q,
      R,
      P: jpvt,
      rank,
      m,
      n,
      info: 0,
      success: true,
      message: 'Success',
    };
  } finally {
    freeAll(Module, [aPtr, tauPtr, jpvtPtr, infoPtr, mPtr, nPtr, ldaPtr]);
  }
}
