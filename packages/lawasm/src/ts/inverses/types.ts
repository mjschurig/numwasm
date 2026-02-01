/**
 * Matrix Inverse Types
 *
 * Type definitions for matrix inversion functions.
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
// GENERAL INVERSE
// ============================================================

/**
 * Options for general matrix inverse.
 */
export interface InvOptions {
  /**
   * Whether to overwrite the input matrix A.
   * @default false
   */
  overwriteA?: boolean;
}

/**
 * Result from general matrix inverse.
 */
export interface InvResult {
  /**
   * Inverse matrix A^(-1).
   * Stored in column-major order.
   */
  inv: Float64Array;

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
// TRIANGULAR INVERSE
// ============================================================

/**
 * Options for triangular matrix inverse.
 */
export interface InvTriangularOptions {
  /**
   * Whether A is upper triangular (true) or lower triangular (false).
   * @default true
   */
  upper?: boolean;

  /**
   * Whether A has unit diagonal (diagonal elements are all 1).
   * @default false
   */
  unitDiagonal?: boolean;

  /**
   * Whether to overwrite the input matrix A.
   * @default false
   */
  overwriteA?: boolean;
}

/**
 * Result from triangular matrix inverse.
 */
export interface InvTriangularResult {
  /**
   * Inverse matrix A^(-1).
   * Stored in column-major order.
   */
  inv: Float64Array;

  /**
   * Matrix dimension.
   */
  n: number;

  /**
   * Whether upper or lower triangular.
   */
  upper: boolean;

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
// SYMMETRIC POSITIVE DEFINITE INVERSE
// ============================================================

/**
 * Options for symmetric positive definite matrix inverse.
 */
export interface InvSymmetricOptions {
  /**
   * Which triangle of A contains the data.
   * @default 'lower'
   */
  uplo?: 'upper' | 'lower';

  /**
   * Whether to overwrite the input matrix A.
   * @default false
   */
  overwriteA?: boolean;
}

/**
 * Result from symmetric positive definite matrix inverse.
 */
export interface InvSymmetricResult {
  /**
   * Inverse matrix A^(-1).
   * Only the triangle specified by uplo is meaningful.
   * Stored in column-major order.
   */
  inv: Float64Array;

  /**
   * Matrix dimension.
   */
  n: number;

  /**
   * Which triangle contains the result.
   */
  uplo: 'upper' | 'lower';

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
// PSEUDOINVERSE
// ============================================================

/**
 * Options for Moore-Penrose pseudoinverse.
 */
export interface PinvOptions {
  /**
   * Reciprocal condition number threshold for determining rank.
   * Singular values s[i] <= rcond * s[0] are treated as zero.
   * If rcond < 0, machine precision is used.
   * @default -1 (machine precision)
   */
  rcond?: number;

  /**
   * Whether to overwrite the input matrix A.
   * @default false
   */
  overwriteA?: boolean;
}

/**
 * Result from pseudoinverse computation.
 */
export interface PinvResult {
  /**
   * Pseudoinverse matrix A^+.
   * For m×n input, the result is n×m.
   * Stored in column-major order.
   */
  pinv: Float64Array;

  /**
   * Singular values of A in descending order.
   */
  s: Float64Array;

  /**
   * Effective rank of A.
   */
  rank: number;

  /**
   * Number of rows of input matrix A.
   */
  m: number;

  /**
   * Number of columns of input matrix A.
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
