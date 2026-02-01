/**
 * QR Factorization with Column Pivoting
 *
 * Compute the QR factorization with column pivoting of a general matrix.
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
import type { Matrix, QRPivotedOptions, QRPivotedResult } from './types.js';

/**
 * Compute the QR factorization with column pivoting of a general m×n matrix A.
 *
 * A * P = Q * R
 *
 * where P is a permutation matrix, Q is orthogonal, and R is upper triangular
 * with diagonal elements of decreasing magnitude.
 *
 * This is useful for rank-revealing factorization: the numerical rank of A
 * can be estimated from the diagonal of R.
 *
 * @param A - Input matrix (m × n). Can be 2D array (row-major) or 1D array (column-major).
 * @param options - Factorization options
 * @returns QR factorization result with Q, R, P, and estimated rank
 *
 * @example
 * ```typescript
 * const A = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]; // Rank 2 matrix
 * const { Q, R, P, rank } = qrPivoted(A);
 * // rank ≈ 2 (estimated from R diagonal)
 * ```
 */
export function qrPivoted(A: Matrix, options: QRPivotedOptions = {}): QRPivotedResult {
  const Module = getLAPACKModule();

  const { mode = 'reduced', overwriteA = false } = options;

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);
  const k = Math.min(m, n);

  // Prepare matrix
  const aData = prepareMatrix(A);

  // Allocate WASM memory
  const aPtr = allocateDoubles(Module, overwriteA ? null : aData, m * n);
  const jpvtPtr = allocateInts(Module, null, n); // Initialize to 0 for free pivoting
  const tauPtr = allocateDoubles(Module, null, k);
  const mPtr = allocateInts(Module, [m]);
  const nPtr = allocateInts(Module, [n]);
  const ldaPtr = allocateInts(Module, [m]);
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

  // Initialize jpvt to 0 (free pivoting)
  const jpvtBase = jpvtPtr >> 2;
  for (let i = 0; i < n; i++) {
    Module.HEAP32[jpvtBase + i] = 0;
  }

  try {
    // Workspace query for DGEQP3
    Module._dgeqp3_(mPtr, nPtr, aPtr, ldaPtr, jpvtPtr, tauPtr, workQueryPtr, lworkPtr, infoPtr);

    let info = readInt(Module, infoPtr);
    if (info !== 0) {
      return {
        R: new Float64Array(0),
        P: new Int32Array(0),
        m,
        n,
        rank: 0,
        info,
        success: false,
        message: getLapackErrorMessage('DGEQP3', info),
      };
    }

    const lwork = Math.ceil(Module.HEAPF64[workQueryPtr >> 3]);
    const workPtr = allocateDoubles(Module, null, lwork);

    // Write lwork
    Module.HEAP32[lworkPtr >> 2] = lwork;

    // Compute QR factorization with pivoting
    Module._dgeqp3_(mPtr, nPtr, aPtr, ldaPtr, jpvtPtr, tauPtr, workPtr, lworkPtr, infoPtr);

    info = readInt(Module, infoPtr);
    if (info !== 0) {
      freeAll(Module, [workPtr]);
      return {
        R: new Float64Array(0),
        P: new Int32Array(0),
        m,
        n,
        rank: 0,
        info,
        success: false,
        message: getLapackErrorMessage('DGEQP3', info),
      };
    }

    // Read the result
    const qrData = readDoubles(Module, aPtr, m * n);

    // Read pivot indices (convert from 1-based to 0-based)
    const P = new Int32Array(n);
    for (let i = 0; i < n; i++) {
      P[i] = Module.HEAP32[jpvtBase + i] - 1; // Convert to 0-based
    }

    // Estimate rank from R diagonal
    // Count diagonal elements above a tolerance
    const tol = Math.max(m, n) * Number.EPSILON * Math.abs(qrData[0]); // Use |R[0,0]| as reference
    let rank = 0;
    for (let i = 0; i < k; i++) {
      if (Math.abs(qrData[i * m + i]) > tol) {
        rank++;
      } else {
        break; // Diagonal elements are in decreasing order
      }
    }

    // Extract R
    let rRows: number;
    if (mode === 'complete') {
      rRows = m;
    } else {
      rRows = k;
    }

    const R = new Float64Array(rRows * n);
    for (let j = 0; j < n; j++) {
      for (let i = 0; i < rRows; i++) {
        if (i <= j) {
          R[j * rRows + i] = qrData[j * m + i];
        }
      }
    }

    // If mode is 'r', we're done
    if (mode === 'r') {
      freeAll(Module, [workPtr]);
      return {
        R,
        P,
        m,
        n,
        rank,
        info: 0,
        success: true,
        message: 'Success',
      };
    }

    // Compute Q using DORGQR
    let qCols: number;
    if (mode === 'complete') {
      qCols = m;
    } else {
      qCols = k;
    }

    // For DORGQR, we need to expand A if mode is 'complete' and m > n
    let qPtr: number;
    if (mode === 'complete' && m > n) {
      qPtr = allocateDoubles(Module, null, m * m);
      const qBase = qPtr >> 3;
      for (let j = 0; j < n; j++) {
        for (let i = 0; i < m; i++) {
          Module.HEAPF64[qBase + j * m + i] = qrData[j * m + i];
        }
      }
      for (let j = n; j < m; j++) {
        for (let i = 0; i < m; i++) {
          Module.HEAPF64[qBase + j * m + i] = 0;
        }
      }
    } else {
      qPtr = aPtr;
    }

    // Workspace query for DORGQR
    const kPtr = allocateInts(Module, [k]);
    const ldqPtr = allocateInts(Module, [m]);
    const qColsPtr = allocateInts(Module, [qCols]);

    Module.HEAP32[lworkPtr >> 2] = -1;
    Module._dorgqr_(mPtr, qColsPtr, kPtr, qPtr, ldqPtr, tauPtr, workQueryPtr, lworkPtr, infoPtr);

    info = readInt(Module, infoPtr);
    if (info !== 0) {
      if (mode === 'complete' && m > n) {
        Module._free(qPtr);
      }
      freeAll(Module, [workPtr, kPtr, ldqPtr, qColsPtr]);
      return {
        R,
        P,
        m,
        n,
        rank,
        info,
        success: false,
        message: getLapackErrorMessage('DORGQR', info),
      };
    }

    const lwork2 = Math.ceil(Module.HEAPF64[workQueryPtr >> 3]);
    if (lwork2 > lwork) {
      freeAll(Module, [workPtr]);
      const workPtr2 = allocateDoubles(Module, null, lwork2);
      Module.HEAP32[lworkPtr >> 2] = lwork2;
      Module._dorgqr_(mPtr, qColsPtr, kPtr, qPtr, ldqPtr, tauPtr, workPtr2, lworkPtr, infoPtr);
      freeAll(Module, [workPtr2]);
    } else {
      Module.HEAP32[lworkPtr >> 2] = lwork;
      Module._dorgqr_(mPtr, qColsPtr, kPtr, qPtr, ldqPtr, tauPtr, workPtr, lworkPtr, infoPtr);
      freeAll(Module, [workPtr]);
    }

    info = readInt(Module, infoPtr);
    if (info !== 0) {
      if (mode === 'complete' && m > n) {
        Module._free(qPtr);
      }
      freeAll(Module, [kPtr, ldqPtr, qColsPtr]);
      return {
        R,
        P,
        m,
        n,
        rank,
        info,
        success: false,
        message: getLapackErrorMessage('DORGQR', info),
      };
    }

    // Read Q
    const Q = readDoubles(Module, qPtr, m * qCols);

    if (mode === 'complete' && m > n) {
      Module._free(qPtr);
    }
    freeAll(Module, [kPtr, ldqPtr, qColsPtr]);

    return {
      Q,
      R,
      P,
      m,
      n,
      rank,
      info: 0,
      success: true,
      message: 'Success',
    };
  } finally {
    freeAll(Module, [aPtr, jpvtPtr, tauPtr, mPtr, nPtr, ldaPtr, infoPtr, lworkPtr, workQueryPtr]);
  }
}
