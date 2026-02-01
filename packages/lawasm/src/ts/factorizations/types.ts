/**
 * Matrix Factorization Types
 *
 * Type definitions for matrix factorization functions.
 */

import type { RealArray } from '../helpers.js';

// ============================================================
// INPUT TYPES
// ============================================================

/**
 * Matrix input type - supports both 2D arrays (row-major) and 1D typed arrays (column-major).
 */
export type Matrix = number[][] | RealArray;

// ============================================================
// LU FACTORIZATION
// ============================================================

/**
 * Options for LU factorization.
 */
export interface LUOptions {
  /**
   * Whether to overwrite the input matrix A.
   * @default false
   */
  overwriteA?: boolean;
}

/**
 * Result from LU factorization.
 * A = P * L * U where P is a permutation matrix, L is lower triangular with unit diagonal,
 * and U is upper triangular.
 */
export interface LUResult {
  /**
   * Lower triangular factor L with unit diagonal.
   * Stored in column-major order.
   */
  L: Float64Array;

  /**
   * Upper triangular factor U.
   * Stored in column-major order.
   */
  U: Float64Array;

  /**
   * Permutation vector P.
   * Row i of A was interchanged with row P[i].
   */
  P: Int32Array;

  /**
   * Combined LU factors in packed format (as returned by DGETRF).
   * L and U are stored in the same array: L is below diagonal, U is on and above diagonal.
   */
  LU: Float64Array;

  /**
   * Number of rows.
   */
  m: number;

  /**
   * Number of columns.
   */
  n: number;

  /**
   * LAPACK info code.
   */
  info: number;

  /**
   * Whether the computation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

// ============================================================
// CHOLESKY FACTORIZATION
// ============================================================

/**
 * Options for Cholesky factorization.
 */
export interface CholeskyOptions {
  /**
   * Whether to compute upper (true) or lower (false) triangular factor.
   * If upper=true: A = U^T * U
   * If upper=false: A = L * L^T
   * @default false (lower triangular)
   */
  upper?: boolean;

  /**
   * Whether to overwrite the input matrix A.
   * @default false
   */
  overwriteA?: boolean;
}

/**
 * Result from Cholesky factorization.
 */
export interface CholeskyResult {
  /**
   * Cholesky factor (L or U depending on options).
   * If upper=false: L such that A = L * L^T
   * If upper=true: U such that A = U^T * U
   */
  factor: Float64Array;

  /**
   * Whether the upper or lower factor was computed.
   */
  upper: boolean;

  /**
   * Matrix dimension.
   */
  n: number;

  /**
   * LAPACK info code.
   */
  info: number;

  /**
   * Whether the computation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

// ============================================================
// QR FACTORIZATION
// ============================================================

/**
 * Options for QR factorization.
 */
export interface QROptions {
  /**
   * What to compute:
   * - 'reduced': Q is m×k and R is k×n where k = min(m,n) (default)
   * - 'complete': Q is m×m and R is m×n
   * - 'r': Only compute R (faster)
   * - 'raw': Return the raw LAPACK output (tau and reflectors)
   * @default 'reduced'
   */
  mode?: 'reduced' | 'complete' | 'r' | 'raw';

  /**
   * Whether to overwrite the input matrix A.
   * @default false
   */
  overwriteA?: boolean;
}

/**
 * Result from QR factorization.
 * A = Q * R where Q is orthogonal and R is upper triangular.
 */
export interface QRResult {
  /**
   * Orthogonal matrix Q.
   * For 'reduced' mode: m×k where k = min(m,n)
   * For 'complete' mode: m×m
   * Not present if mode='r' or mode='raw'
   */
  Q?: Float64Array;

  /**
   * Upper triangular matrix R.
   * For 'reduced' mode: k×n where k = min(m,n)
   * For 'complete' mode: m×n
   */
  R: Float64Array;

  /**
   * Scalar factors of elementary reflectors (for mode='raw').
   */
  tau?: Float64Array;

  /**
   * Number of rows of input matrix.
   */
  m: number;

  /**
   * Number of columns of input matrix.
   */
  n: number;

  /**
   * Rank (min(m, n)).
   */
  k: number;

  /**
   * LAPACK info code.
   */
  info: number;

  /**
   * Whether the computation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

// ============================================================
// QR WITH COLUMN PIVOTING
// ============================================================

/**
 * Options for QR factorization with column pivoting.
 */
export interface QRPivotedOptions {
  /**
   * What to compute:
   * - 'reduced': Q is m×k and R is k×n where k = min(m,n) (default)
   * - 'complete': Q is m×m and R is m×n
   * - 'r': Only compute R (faster)
   * @default 'reduced'
   */
  mode?: 'reduced' | 'complete' | 'r';

  /**
   * Whether to overwrite the input matrix A.
   * @default false
   */
  overwriteA?: boolean;
}

/**
 * Result from QR factorization with column pivoting.
 * A * P = Q * R where P is a permutation matrix.
 */
export interface QRPivotedResult {
  /**
   * Orthogonal matrix Q (if computed).
   */
  Q?: Float64Array;

  /**
   * Upper triangular matrix R.
   */
  R: Float64Array;

  /**
   * Permutation vector.
   * Column j of A*P is column P[j] of A.
   */
  P: Int32Array;

  /**
   * Number of rows.
   */
  m: number;

  /**
   * Number of columns.
   */
  n: number;

  /**
   * Effective rank (number of independent columns).
   */
  rank: number;

  /**
   * LAPACK info code.
   */
  info: number;

  /**
   * Whether the computation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

// ============================================================
// LQ FACTORIZATION
// ============================================================

/**
 * Options for LQ factorization.
 */
export interface LQOptions {
  /**
   * What to compute:
   * - 'reduced': L is m×k and Q is k×n where k = min(m,n) (default)
   * - 'complete': L is m×n and Q is n×n
   * @default 'reduced'
   */
  mode?: 'reduced' | 'complete';

  /**
   * Whether to overwrite the input matrix A.
   * @default false
   */
  overwriteA?: boolean;
}

/**
 * Result from LQ factorization.
 * A = L * Q where L is lower triangular and Q is orthogonal.
 */
export interface LQResult {
  /**
   * Lower triangular matrix L.
   */
  L: Float64Array;

  /**
   * Orthogonal matrix Q.
   */
  Q: Float64Array;

  /**
   * Number of rows.
   */
  m: number;

  /**
   * Number of columns.
   */
  n: number;

  /**
   * LAPACK info code.
   */
  info: number;

  /**
   * Whether the computation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

// ============================================================
// LDL FACTORIZATION
// ============================================================

/**
 * Options for LDL factorization.
 */
export interface LDLOptions {
  /**
   * Whether to use upper (true) or lower (false) triangle.
   * If upper=true: A = U * D * U^T
   * If upper=false: A = L * D * L^T
   * @default false (lower triangular)
   */
  upper?: boolean;

  /**
   * Whether to overwrite the input matrix A.
   * @default false
   */
  overwriteA?: boolean;
}

/**
 * Result from LDL^T factorization.
 * A = L * D * L^T (or U * D * U^T) where L (or U) is unit triangular
 * and D is block diagonal with 1x1 and 2x2 blocks.
 */
export interface LDLResult {
  /**
   * Unit triangular factor (L or U).
   */
  factor: Float64Array;

  /**
   * Block diagonal matrix D.
   * Stored as a vector of diagonal elements.
   * For 2x2 blocks, consecutive elements represent the block.
   */
  D: Float64Array;

  /**
   * Pivot indices indicating the structure of D.
   * If ipiv[k] > 0, D[k,k] is a 1x1 block.
   * If ipiv[k] < 0, D[k:k+1, k:k+1] is a 2x2 block.
   */
  ipiv: Int32Array;

  /**
   * Whether upper or lower factor was computed.
   */
  upper: boolean;

  /**
   * Matrix dimension.
   */
  n: number;

  /**
   * LAPACK info code.
   */
  info: number;

  /**
   * Whether the computation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

// ============================================================
// SCHUR DECOMPOSITION
// ============================================================

/**
 * Options for Schur decomposition.
 */
export interface SchurOptions {
  /**
   * Whether to compute the Schur vectors Z.
   * @default true
   */
  computeZ?: boolean;

  /**
   * Whether to sort eigenvalues.
   * @default false
   */
  sort?: boolean;

  /**
   * Selection function for sorting (if sort=true).
   * Takes real and imaginary parts of eigenvalue, returns true to select.
   */
  select?: (wr: number, wi: number) => boolean;

  /**
   * Whether to overwrite the input matrix A.
   * @default false
   */
  overwriteA?: boolean;
}

/**
 * Result from Schur decomposition.
 * A = Z * T * Z^T where T is quasi-upper triangular (real Schur form)
 * and Z is orthogonal.
 */
export interface SchurResult {
  /**
   * Schur form T (quasi-upper triangular).
   * Diagonal blocks are 1x1 (real eigenvalues) or 2x2 (complex conjugate pairs).
   */
  T: Float64Array;

  /**
   * Orthogonal matrix Z of Schur vectors (if computed).
   */
  Z?: Float64Array;

  /**
   * Real parts of eigenvalues.
   */
  wr: Float64Array;

  /**
   * Imaginary parts of eigenvalues.
   */
  wi: Float64Array;

  /**
   * Matrix dimension.
   */
  n: number;

  /**
   * Number of eigenvalues for which the condition was satisfied (if sorted).
   */
  sdim?: number;

  /**
   * LAPACK info code.
   */
  info: number;

  /**
   * Whether the computation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

// ============================================================
// HESSENBERG REDUCTION
// ============================================================

/**
 * Options for Hessenberg reduction.
 */
export interface HessenbergOptions {
  /**
   * Whether to compute the orthogonal matrix Q.
   * @default true
   */
  computeQ?: boolean;

  /**
   * Whether to overwrite the input matrix A.
   * @default false
   */
  overwriteA?: boolean;
}

/**
 * Result from Hessenberg reduction.
 * A = Q * H * Q^T where H is upper Hessenberg and Q is orthogonal.
 */
export interface HessenbergResult {
  /**
   * Upper Hessenberg matrix H.
   */
  H: Float64Array;

  /**
   * Orthogonal matrix Q (if computed).
   */
  Q?: Float64Array;

  /**
   * Matrix dimension.
   */
  n: number;

  /**
   * LAPACK info code.
   */
  info: number;

  /**
   * Whether the computation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}
