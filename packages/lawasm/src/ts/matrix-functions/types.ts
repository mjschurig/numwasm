/**
 * Matrix Functions Types
 *
 * Type definitions for matrix functions (expm, logm, sqrtm, etc.)
 * and specialized decompositions.
 */

import type { RealArray } from '../helpers.js';

// ============================================================
// INPUT TYPES
// ============================================================

/**
 * Matrix input type.
 */
export type Matrix = number[][] | RealArray;

// ============================================================
// MATRIX FUNCTIONS
// ============================================================

/**
 * Result from matrix exponential.
 */
export interface ExpmResult {
  /**
   * Matrix exponential e^A.
   */
  E: Float64Array;

  /**
   * Matrix dimension.
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

/**
 * Result from matrix logarithm.
 */
export interface LogmResult {
  /**
   * Matrix logarithm log(A).
   */
  L: Float64Array;

  /**
   * Matrix dimension.
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

/**
 * Result from matrix square root.
 */
export interface SqrtmResult {
  /**
   * Matrix square root sqrt(A) such that S*S = A.
   */
  S: Float64Array;

  /**
   * Matrix dimension.
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

/**
 * Result from matrix power.
 */
export interface PowmResult {
  /**
   * Matrix power A^p.
   */
  P: Float64Array;

  /**
   * Matrix dimension.
   */
  n: number;

  /**
   * The power used.
   */
  p: number;

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
 * Result from general matrix function.
 */
export interface FunmResult {
  /**
   * Result f(A).
   */
  F: Float64Array;

  /**
   * Matrix dimension.
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
// SPECIALIZED DECOMPOSITIONS
// ============================================================

/**
 * Result from polar decomposition.
 */
export interface PolarDecompositionResult {
  /**
   * Unitary (orthogonal for real) matrix U.
   */
  U: Float64Array;

  /**
   * Positive semidefinite matrix P such that A = U*P.
   */
  P: Float64Array;

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

/**
 * Options for rank-revealing QR.
 */
export interface RRQROptions {
  /**
   * Tolerance for determining rank.
   * @default machine epsilon * max(m,n) * ||A||
   */
  tol?: number;
}

/**
 * Result from rank-revealing QR.
 */
export interface RRQRResult {
  /**
   * Orthogonal matrix Q (m × m or m × k depending on mode).
   */
  Q: Float64Array;

  /**
   * Upper triangular matrix R.
   */
  R: Float64Array;

  /**
   * Column permutation (0-based indices).
   */
  P: Int32Array;

  /**
   * Numerical rank.
   */
  rank: number;

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

/**
 * Result from Cosine-Sine decomposition.
 */
export interface CSDResult {
  /**
   * Left unitary matrix U1 (p × p).
   */
  U1: Float64Array;

  /**
   * Left unitary matrix U2 ((m-p) × (m-p)).
   */
  U2: Float64Array;

  /**
   * Right unitary matrix V1 (q × q).
   */
  V1: Float64Array;

  /**
   * Right unitary matrix V2 ((m-q) × (m-q)).
   */
  V2: Float64Array;

  /**
   * Cosine values.
   */
  theta: Float64Array;

  /**
   * Partition size p.
   */
  p: number;

  /**
   * Partition size q.
   */
  q: number;

  /**
   * Matrix dimension m.
   */
  m: number;

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
// COMPLEX-SPECIFIC FUNCTIONS
// ============================================================

/**
 * Complex number type.
 */
export interface Complex {
  re: number;
  im: number;
}

/**
 * Result from Hermitian eigenvalue problem.
 */
export interface EigHermitianResult {
  /**
   * Eigenvalues (real, in ascending order).
   */
  w: Float64Array;

  /**
   * Eigenvectors (column-major, complex interleaved).
   */
  V?: Float64Array;

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

/**
 * Result from Hermitian Cholesky.
 */
export interface CholeskyHermitianResult {
  /**
   * Cholesky factor (complex, interleaved format).
   */
  L: Float64Array;

  /**
   * Matrix dimension.
   */
  n: number;

  /**
   * Which triangle was computed.
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

/**
 * Result from Hermitian solve.
 */
export interface SolveHermitianResult {
  /**
   * Solution vector/matrix (complex, interleaved format).
   */
  x: Float64Array;

  /**
   * Matrix dimension.
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
