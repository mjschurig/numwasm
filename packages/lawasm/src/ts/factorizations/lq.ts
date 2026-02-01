/**
 * LQ Factorization
 *
 * Compute the LQ factorization of a general matrix.
 *
 * Note: DGELQF is not currently exported. This implementation uses QR of A^T.
 */

import type { Matrix, LQOptions, LQResult } from './types.js';
import { qr } from './qr.js';
import { getMatrixDimensions, prepareMatrix } from '../helpers.js';

/**
 * Compute the LQ factorization of a general m×n matrix A.
 *
 * A = L * Q
 *
 * where L is lower triangular (m×k) and Q is orthogonal (k×n),
 * with k = min(m, n).
 *
 * Implementation note: This uses QR of A^T since DGELQF is not exported.
 * A = L * Q implies A^T = Q^T * L^T, so we compute QR of A^T.
 *
 * @param A - Input matrix (m × n). Can be 2D array (row-major) or 1D array (column-major).
 * @param options - Factorization options
 * @returns LQ factorization result with L and Q
 *
 * @example
 * ```typescript
 * const A = [[1, 2, 3], [4, 5, 6]];
 * const { L, Q, success } = lq(A);
 * // L is 2×2 lower triangular
 * // Q is 2×3 orthogonal
 * ```
 */
export function lq(A: Matrix, options: LQOptions = {}): LQResult {
  const { mode = 'reduced' } = options;

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);
  const k = Math.min(m, n);

  // Prepare matrix and compute transpose
  const aData = prepareMatrix(A);

  // A is m×n in column-major, A^T is n×m in column-major
  const atData = new Float64Array(n * m);
  for (let j = 0; j < n; j++) {
    for (let i = 0; i < m; i++) {
      // A[i,j] in column-major = aData[j*m + i]
      // A^T[j,i] in column-major = atData[i*n + j]
      atData[i * n + j] = aData[j * m + i];
    }
  }

  // Compute QR of A^T (n×m matrix)
  // A^T = Qt * Rt where Qt is n×k and Rt is k×m
  const qrResult = qr(atData, {
    mode: mode === 'complete' ? 'complete' : 'reduced',
    overwriteA: true,
  });

  if (!qrResult.success) {
    return {
      L: new Float64Array(0),
      Q: new Float64Array(0),
      m,
      n,
      info: qrResult.info,
      success: false,
      message: qrResult.message,
    };
  }

  // A = L * Q means A^T = Q^T * L^T
  // So Qt = Q^T and Rt = L^T
  // Therefore Q = Qt^T and L = Rt^T

  // Qt is n×qCols, so Q = Qt^T is qCols×n
  const qCols = mode === 'complete' ? n : k;
  const Q = new Float64Array(qCols * n);
  if (qrResult.Q) {
    // Transpose Qt (n×qCols) to get Q (qCols×n)
    for (let j = 0; j < n; j++) {
      for (let i = 0; i < qCols; i++) {
        // Qt[j,i] in column-major = qrResult.Q[i*n + j]
        // Q[i,j] in column-major = Q[j*qCols + i]
        Q[j * qCols + i] = qrResult.Q[i * n + j];
      }
    }
  }

  // Rt is k×m (or n×m for complete), so L = Rt^T is m×k (or m×n)
  const lCols = mode === 'complete' ? n : k;
  const L = new Float64Array(m * lCols);
  // R from QR is stored as rRows×m where rRows = k (or n for complete)
  const rRows = qrResult.R.length / m;
  // Transpose Rt (rRows×m) to get L (m×rRows)
  for (let j = 0; j < m; j++) {
    for (let i = 0; i < Math.min(rRows, lCols); i++) {
      // Rt[i,j] in column-major = qrResult.R[j*rRows + i]
      // L[j,i] in column-major = L[i*m + j]
      L[i * m + j] = qrResult.R[j * rRows + i];
    }
  }

  return {
    L,
    Q,
    m,
    n,
    info: 0,
    success: true,
    message: 'Success',
  };
}
