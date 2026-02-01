/**
 * Least Squares Types
 *
 * Type definitions for least squares and minimum norm functions.
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
 * Vector input type.
 */
export type Vector = RealArray;

// ============================================================
// LSTSQ (QR-based)
// ============================================================

/**
 * Options for QR-based least squares.
 */
export interface LstSqOptions {
  /**
   * Whether to overwrite the input matrix A.
   * @default false
   */
  overwriteA?: boolean;

  /**
   * Whether to overwrite the input vector/matrix b.
   * @default false
   */
  overwriteB?: boolean;
}

/**
 * Result from QR-based least squares.
 */
export interface LstSqResult {
  /**
   * Solution vector/matrix x that minimizes ||Ax - b||.
   * For overdetermined systems (m > n), this is the least squares solution.
   * For underdetermined systems (m < n), this is the minimum norm solution.
   */
  x: Float64Array;

  /**
   * Residual sum of squares for each column of b.
   * Only computed for overdetermined systems (m > n).
   */
  residuals?: Float64Array;

  /**
   * Number of rows of input matrix A.
   */
  m: number;

  /**
   * Number of columns of input matrix A.
   */
  n: number;

  /**
   * Number of right-hand sides.
   */
  nrhs: number;

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
// LSTSQ SVD (rank-deficient safe)
// ============================================================

/**
 * Options for SVD-based least squares.
 */
export interface LstSqSVDOptions {
  /**
   * Reciprocal condition number threshold for determining rank.
   * Singular values s[i] <= rcond * s[0] are treated as zero.
   * If rcond < 0, machine precision is used.
   * @default -1 (machine precision)
   */
  rcond?: number;

  /**
   * Which SVD algorithm to use:
   * - 'gelss': Standard SVD (DGELSS)
   * - 'gelsd': Divide and conquer SVD (DGELSD) - faster for large matrices
   * @default 'gelsd'
   */
  algorithm?: 'gelss' | 'gelsd';

  /**
   * Whether to overwrite the input matrix A.
   * @default false
   */
  overwriteA?: boolean;

  /**
   * Whether to overwrite the input vector/matrix b.
   * @default false
   */
  overwriteB?: boolean;
}

/**
 * Result from SVD-based least squares.
 */
export interface LstSqSVDResult {
  /**
   * Solution vector/matrix x.
   */
  x: Float64Array;

  /**
   * Singular values of A in descending order.
   */
  s: Float64Array;

  /**
   * Effective rank of A.
   * Number of singular values > rcond * max(s).
   */
  rank: number;

  /**
   * Residual sum of squares for each column of b.
   * Only computed for overdetermined systems (m > n) with rank = n.
   */
  residuals?: Float64Array;

  /**
   * Number of rows of input matrix A.
   */
  m: number;

  /**
   * Number of columns of input matrix A.
   */
  n: number;

  /**
   * Number of right-hand sides.
   */
  nrhs: number;

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
// LSTSQ GELSY (QR with pivoting)
// ============================================================

/**
 * Options for QR with pivoting least squares.
 */
export interface LstSqGelsyOptions {
  /**
   * Reciprocal condition number threshold for determining rank.
   * @default machine precision
   */
  rcond?: number;

  /**
   * Whether to overwrite the input matrix A.
   * @default false
   */
  overwriteA?: boolean;

  /**
   * Whether to overwrite the input vector/matrix b.
   * @default false
   */
  overwriteB?: boolean;
}

/**
 * Result from QR with pivoting least squares.
 */
export interface LstSqGelsyResult {
  /**
   * Solution vector/matrix x.
   */
  x: Float64Array;

  /**
   * Effective rank of A.
   */
  rank: number;

  /**
   * Column permutation from pivoting.
   * jpvt[i] = k means column k of A was moved to position i.
   */
  jpvt: Int32Array;

  /**
   * Number of rows of input matrix A.
   */
  m: number;

  /**
   * Number of columns of input matrix A.
   */
  n: number;

  /**
   * Number of right-hand sides.
   */
  nrhs: number;

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
// CONSTRAINED LEAST SQUARES
// ============================================================

/**
 * Options for equality-constrained least squares.
 */
export interface ConstrainedLstSqOptions {
  /**
   * Whether to overwrite input matrices.
   * @default false
   */
  overwriteInputs?: boolean;
}

/**
 * Result from equality-constrained least squares.
 * Solves: minimize ||Ax - b||_2 subject to Cx = d
 */
export interface ConstrainedLstSqResult {
  /**
   * Solution vector x.
   */
  x: Float64Array;

  /**
   * Residual norm ||Ax - b||.
   */
  residualNorm: number;

  /**
   * Number of variables.
   */
  n: number;

  /**
   * Number of equations in A.
   */
  m: number;

  /**
   * Number of equality constraints.
   */
  p: number;

  /**
   * LAPACK info code (or custom error code).
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
// GENERALIZED LEAST SQUARES
// ============================================================

/**
 * Options for generalized least squares.
 */
export interface GeneralizedLstSqOptions {
  /**
   * Whether to overwrite input matrices.
   * @default false
   */
  overwriteInputs?: boolean;
}

/**
 * Result from generalized least squares.
 * Solves: minimize ||y||_2 subject to Ax + By = d
 */
export interface GeneralizedLstSqResult {
  /**
   * Solution vector x.
   */
  x: Float64Array;

  /**
   * Residual vector y such that Ax + By = d.
   */
  y: Float64Array;

  /**
   * Number of columns of A (dimension of x).
   */
  n: number;

  /**
   * Number of columns of B (dimension of y).
   */
  m: number;

  /**
   * Number of rows (equations).
   */
  p: number;

  /**
   * LAPACK info code (or custom error code).
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
