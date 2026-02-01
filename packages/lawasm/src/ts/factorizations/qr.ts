/**
 * QR Factorization
 *
 * Compute the QR factorization of a general matrix.
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
import type { Matrix, QROptions, QRResult } from './types.js';

/**
 * Compute the QR factorization of a general m×n matrix A.
 *
 * A = Q * R
 *
 * where Q is orthogonal (m×m or m×k) and R is upper triangular (m×n or k×n),
 * with k = min(m, n).
 *
 * @param A - Input matrix (m × n). Can be 2D array (row-major) or 1D array (column-major).
 * @param options - Factorization options
 * @returns QR factorization result with Q and R
 *
 * @example
 * ```typescript
 * const A = [[1, 2], [3, 4], [5, 6]];
 * const { Q, R, success } = qr(A);
 * // Q is 3×2 orthogonal matrix (reduced QR)
 * // R is 2×2 upper triangular
 * ```
 *
 * @example
 * ```typescript
 * // Full QR factorization
 * const { Q, R } = qr(A, { mode: 'complete' });
 * // Q is 3×3 orthogonal matrix
 * // R is 3×2 upper triangular
 * ```
 *
 * @example
 * ```typescript
 * // Only compute R (faster)
 * const { R } = qr(A, { mode: 'r' });
 * ```
 */
export function qr(A: Matrix, options: QROptions = {}): QRResult {
  const Module = getLAPACKModule();

  const { mode = 'reduced', overwriteA = false } = options;

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);
  const k = Math.min(m, n);

  // Prepare matrix
  const aData = prepareMatrix(A);

  // Allocate WASM memory
  const aPtr = allocateDoubles(Module, overwriteA ? null : aData, m * n);
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

  try {
    // Workspace query for DGEQRF
    Module._dgeqrf_(mPtr, nPtr, aPtr, ldaPtr, tauPtr, workQueryPtr, lworkPtr, infoPtr);

    let info = readInt(Module, infoPtr);
    if (info !== 0) {
      return {
        R: new Float64Array(0),
        m,
        n,
        k,
        info,
        success: false,
        message: getLapackErrorMessage('DGEQRF', info),
      };
    }

    const lwork = Math.ceil(Module.HEAPF64[workQueryPtr >> 3]);
    const workPtr = allocateDoubles(Module, null, lwork);

    // Write lwork
    Module.HEAP32[lworkPtr >> 2] = lwork;

    // Compute QR factorization
    Module._dgeqrf_(mPtr, nPtr, aPtr, ldaPtr, tauPtr, workPtr, lworkPtr, infoPtr);

    info = readInt(Module, infoPtr);
    if (info !== 0) {
      freeAll(Module, [workPtr]);
      return {
        R: new Float64Array(0),
        m,
        n,
        k,
        info,
        success: false,
        message: getLapackErrorMessage('DGEQRF', info),
      };
    }

    // Read the result (A now contains R in upper triangle and reflectors below)
    const qrData = readDoubles(Module, aPtr, m * n);
    const tau = readDoubles(Module, tauPtr, k);

    // Extract R
    let R: Float64Array;
    let rRows: number;
    if (mode === 'complete') {
      rRows = m;
    } else {
      rRows = k;
    }

    R = new Float64Array(rRows * n);
    for (let j = 0; j < n; j++) {
      for (let i = 0; i < rRows; i++) {
        if (i <= j) {
          R[j * rRows + i] = qrData[j * m + i];
        }
        // Below diagonal is 0 (already initialized)
      }
    }

    // If mode is 'r' or 'raw', we're done
    if (mode === 'r') {
      freeAll(Module, [workPtr]);
      return {
        R,
        m,
        n,
        k,
        info: 0,
        success: true,
        message: 'Success',
      };
    }

    if (mode === 'raw') {
      freeAll(Module, [workPtr]);
      return {
        R,
        tau,
        m,
        n,
        k,
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
      // Need to allocate m×m and copy the reflectors
      qPtr = allocateDoubles(Module, null, m * m);
      const qBase = qPtr >> 3;
      // Copy the first n columns (containing reflectors)
      for (let j = 0; j < n; j++) {
        for (let i = 0; i < m; i++) {
          Module.HEAPF64[qBase + j * m + i] = qrData[j * m + i];
        }
      }
      // Zero the remaining columns
      for (let j = n; j < m; j++) {
        for (let i = 0; i < m; i++) {
          Module.HEAPF64[qBase + j * m + i] = 0;
        }
      }
    } else {
      qPtr = aPtr; // Reuse A which already contains reflectors
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
        m,
        n,
        k,
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
        m,
        n,
        k,
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
      m,
      n,
      k,
      info: 0,
      success: true,
      message: 'Success',
    };
  } finally {
    freeAll(Module, [aPtr, tauPtr, mPtr, nPtr, ldaPtr, infoPtr, lworkPtr, workQueryPtr]);
  }
}
