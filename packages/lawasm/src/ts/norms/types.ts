/**
 * Matrix Norms, Condition Numbers, Determinants & Rank Types
 *
 * Type definitions for matrix analysis functions.
 */

import type { RealArray } from '../helpers.js';

// ============================================================
// INPUT TYPES
// ============================================================

/**
 * Matrix input type - supports both 2D arrays (row-major) and 1D typed arrays (column-major).
 */
export type Matrix = number[][] | RealArray;

/**
 * Norm type for matrices.
 * - 1 or '1': 1-norm (max column sum)
 * - 2 or '2': 2-norm (spectral norm, largest singular value)
 * - Infinity or 'inf': infinity-norm (max row sum)
 * - 'fro': Frobenius norm
 * - 'max': max absolute element
 */
export type NormType = 1 | 2 | '1' | '2' | 'inf' | 'fro' | 'max';

// ============================================================
// NORM
// ============================================================

/**
 * Options for matrix norm computation.
 */
export interface NormOptions {
  /**
   * Which norm to compute.
   * @default 'fro'
   */
  ord?: NormType;
}

/**
 * Result from matrix norm computation.
 */
export interface NormResult {
  /**
   * The computed norm value.
   */
  norm: number;

  /**
   * Which norm was computed.
   */
  ord: NormType;

  /**
   * Number of rows.
   */
  m: number;

  /**
   * Number of columns.
   */
  n: number;

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
// CONDITION NUMBER
// ============================================================

/**
 * Options for condition number computation.
 */
export interface CondOptions {
  /**
   * Which norm to use for condition number.
   * Only 1, 2, and 'inf' are supported.
   * @default 2
   */
  ord?: 1 | 2 | '1' | '2' | 'inf';
}

/**
 * Result from condition number computation.
 */
export interface CondResult {
  /**
   * The condition number κ(A) = ||A|| * ||A^(-1)||.
   * A large condition number indicates an ill-conditioned matrix.
   */
  cond: number;

  /**
   * Which norm was used.
   */
  ord: 1 | 2 | 'inf';

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
// CONDITION ESTIMATE
// ============================================================

/**
 * Options for fast condition number estimate.
 */
export interface CondEstOptions {
  /**
   * Which norm to use ('1' or 'inf').
   * @default '1'
   */
  norm?: '1' | 'inf';
}

/**
 * Result from condition number estimate.
 */
export interface CondEstResult {
  /**
   * Estimated condition number.
   */
  cond: number;

  /**
   * Reciprocal condition number (1/cond).
   */
  rcond: number;

  /**
   * Which norm was used.
   */
  norm: '1' | 'inf';

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
// RECIPROCAL CONDITION NUMBER
// ============================================================

/**
 * Options for reciprocal condition number.
 */
export interface RcondOptions {
  /**
   * Which norm to use ('1' or 'inf').
   * @default '1'
   */
  norm?: '1' | 'inf';
}

/**
 * Result from reciprocal condition number computation.
 */
export interface RcondResult {
  /**
   * Reciprocal of the condition number.
   * If rcond is small (e.g., < machine_epsilon * n), the matrix is nearly singular.
   */
  rcond: number;

  /**
   * Which norm was used.
   */
  norm: '1' | 'inf';

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
// DETERMINANT
// ============================================================

/**
 * Options for determinant computation.
 */
export interface DetOptions {
  /**
   * Whether to overwrite the input matrix A.
   * @default false
   */
  overwriteA?: boolean;
}

/**
 * Result from determinant computation.
 */
export interface DetResult {
  /**
   * The determinant of A.
   */
  det: number;

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
// LOG DETERMINANT
// ============================================================

/**
 * Options for log-determinant computation.
 */
export interface LogDetOptions {
  /**
   * Whether to overwrite the input matrix A.
   * @default false
   */
  overwriteA?: boolean;
}

/**
 * Result from log-determinant computation.
 */
export interface LogDetResult {
  /**
   * The natural logarithm of the absolute value of the determinant.
   * log(|det(A)|)
   */
  logdet: number;

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
// SIGN AND LOG DETERMINANT
// ============================================================

/**
 * Options for sign and log-determinant computation.
 */
export interface SlogDetOptions {
  /**
   * Whether to overwrite the input matrix A.
   * @default false
   */
  overwriteA?: boolean;
}

/**
 * Result from sign and log-determinant computation.
 */
export interface SlogDetResult {
  /**
   * The sign of the determinant: -1, 0, or 1.
   * det(A) = sign * exp(logabsdet)
   */
  sign: number;

  /**
   * The natural logarithm of the absolute value of the determinant.
   */
  logabsdet: number;

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
// MATRIX RANK
// ============================================================

/**
 * Options for matrix rank computation.
 */
export interface RankOptions {
  /**
   * Tolerance for determining rank.
   * Singular values below tol * max(s) are considered zero.
   * @default max(m,n) * eps * max(s)
   */
  tol?: number;
}

/**
 * Result from matrix rank computation.
 */
export interface RankResult {
  /**
   * The numerical rank of A.
   */
  rank: number;

  /**
   * Singular values of A (for reference).
   */
  s: Float64Array;

  /**
   * Tolerance used for rank determination.
   */
  tol: number;

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
