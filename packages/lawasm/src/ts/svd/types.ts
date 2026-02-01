/**
 * Singular Value Decomposition Types
 *
 * Type definitions for SVD functions.
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
// SVD OPTIONS
// ============================================================

/**
 * Options for full SVD computation.
 */
export interface SVDOptions {
  /**
   * What to compute:
   * - 'full': Full SVD with U (m×m) and Vt (n×n)
   * - 'reduced': Reduced SVD with U (m×k) and Vt (k×n) where k = min(m,n)
   * - 'none': Only singular values, no vectors
   * @default 'reduced'
   */
  mode?: 'full' | 'reduced' | 'none';

  /**
   * Which algorithm to use:
   * - 'gesvd': Standard algorithm (DGESVD) - more accurate
   * - 'gesdd': Divide and conquer (DGESDD) - faster for large matrices
   * @default 'gesdd'
   */
  algorithm?: 'gesvd' | 'gesdd';

  /**
   * Whether to overwrite the input matrix A.
   * @default false
   */
  overwriteA?: boolean;
}

/**
 * Options for computing only singular values.
 */
export interface SVDValsOptions {
  /**
   * Which algorithm to use.
   * @default 'gesdd'
   */
  algorithm?: 'gesvd' | 'gesdd';

  /**
   * Whether to overwrite the input matrix A.
   * @default false
   */
  overwriteA?: boolean;
}

/**
 * Options for SVD-based rank computation.
 */
export interface SVDRankOptions {
  /**
   * Tolerance for determining rank.
   * Singular values below tol * max(s) are considered zero.
   * @default max(m,n) * eps * max(s)
   */
  tol?: number;

  /**
   * Which algorithm to use.
   * @default 'gesdd'
   */
  algorithm?: 'gesvd' | 'gesdd';
}

// ============================================================
// SVD RESULTS
// ============================================================

/**
 * Result from SVD computation.
 */
export interface SVDResult {
  /**
   * Left singular vectors U.
   * - For 'full' mode: m×m orthogonal matrix
   * - For 'reduced' mode: m×k matrix where k = min(m,n)
   * - Not present for 'none' mode
   *
   * Stored in column-major order.
   */
  U?: Float64Array;

  /**
   * Singular values in descending order.
   * Length is min(m, n).
   */
  s: Float64Array;

  /**
   * Right singular vectors V^T (transposed).
   * - For 'full' mode: n×n orthogonal matrix
   * - For 'reduced' mode: k×n matrix where k = min(m,n)
   * - Not present for 'none' mode
   *
   * Stored in column-major order.
   */
  Vt?: Float64Array;

  /**
   * Number of rows of input matrix.
   */
  m: number;

  /**
   * Number of columns of input matrix.
   */
  n: number;

  /**
   * Number of singular values (min(m, n)).
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

/**
 * Result from singular values only computation.
 */
export interface SVDValsResult {
  /**
   * Singular values in descending order.
   */
  s: Float64Array;

  /**
   * Number of rows of input matrix.
   */
  m: number;

  /**
   * Number of columns of input matrix.
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

/**
 * Result from SVD-based rank computation.
 */
export interface SVDRankResult {
  /**
   * Estimated numerical rank.
   */
  rank: number;

  /**
   * Singular values in descending order.
   */
  s: Float64Array;

  /**
   * Tolerance used for rank determination.
   */
  tol: number;

  /**
   * Number of rows of input matrix.
   */
  m: number;

  /**
   * Number of columns of input matrix.
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
