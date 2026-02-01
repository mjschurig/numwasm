/**
 * Matrix Utilities & Properties Types
 *
 * Type definitions for matrix utility functions and property tests.
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
 * Complex matrix input type - array of [real, imag] pairs or interleaved real/imag.
 */
export type ComplexMatrix = [number, number][][] | RealArray;

// ============================================================
// MATRIX UTILITIES
// ============================================================

/**
 * Result from matrix transpose.
 */
export interface TransposeResult {
  /**
   * Transposed matrix (n × m).
   */
  T: Float64Array;

  /**
   * Number of rows of result.
   */
  m: number;

  /**
   * Number of columns of result.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Result from complex conjugate.
 */
export interface ConjugateResult {
  /**
   * Conjugated matrix.
   */
  conj: Float64Array;

  /**
   * Number of rows.
   */
  m: number;

  /**
   * Number of columns.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Result from Hermitian (conjugate transpose).
 */
export interface HermitianResult {
  /**
   * Conjugate transpose (n × m).
   */
  H: Float64Array;

  /**
   * Number of rows of result.
   */
  m: number;

  /**
   * Number of columns of result.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Options for triangular extraction.
 */
export interface TriuOptions {
  /**
   * Diagonal offset. k=0 is the main diagonal, k>0 is above, k<0 is below.
   * @default 0
   */
  k?: number;
}

/**
 * Result from upper triangular extraction.
 */
export interface TriuResult {
  /**
   * Upper triangular part of the matrix.
   */
  U: Float64Array;

  /**
   * Number of rows.
   */
  m: number;

  /**
   * Number of columns.
   */
  n: number;

  /**
   * Diagonal offset used.
   */
  k: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Options for lower triangular extraction.
 */
export interface TrilOptions {
  /**
   * Diagonal offset. k=0 is the main diagonal, k>0 is above, k<0 is below.
   * @default 0
   */
  k?: number;
}

/**
 * Result from lower triangular extraction.
 */
export interface TrilResult {
  /**
   * Lower triangular part of the matrix.
   */
  L: Float64Array;

  /**
   * Number of rows.
   */
  m: number;

  /**
   * Number of columns.
   */
  n: number;

  /**
   * Diagonal offset used.
   */
  k: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Result from diagonal extraction/creation.
 */
export interface DiagResult {
  /**
   * Diagonal elements (if extracting) or diagonal matrix (if creating).
   */
  diag: Float64Array;

  /**
   * Number of elements or matrix dimension.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Result from matrix trace.
 */
export interface TraceResult {
  /**
   * The trace (sum of diagonal elements).
   */
  trace: number;

  /**
   * Matrix dimension.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Result from matrix balancing.
 */
export interface BalanceResult {
  /**
   * Balanced matrix.
   */
  B: Float64Array;

  /**
   * Diagonal scaling factors.
   */
  scale: Float64Array;

  /**
   * Index of first balanced row/column.
   */
  ilo: number;

  /**
   * Index of last balanced row/column.
   */
  ihi: number;

  /**
   * Matrix dimension.
   */
  n: number;

  /**
   * LAPACK info code.
   */
  info: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

// ============================================================
// MATRIX PROPERTIES & TESTS
// ============================================================

/**
 * Options for symmetry test.
 */
export interface IsSymmetricOptions {
  /**
   * Tolerance for comparison.
   * @default 1e-10
   */
  tol?: number;
}

/**
 * Result from symmetry test.
 */
export interface IsSymmetricResult {
  /**
   * Whether the matrix is symmetric.
   */
  isSymmetric: boolean;

  /**
   * Maximum deviation from symmetry (max |A[i,j] - A[j,i]|).
   */
  maxDeviation: number;

  /**
   * Matrix dimension.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Options for Hermitian test.
 */
export interface IsHermitianOptions {
  /**
   * Tolerance for comparison.
   * @default 1e-10
   */
  tol?: number;
}

/**
 * Result from Hermitian test.
 */
export interface IsHermitianResult {
  /**
   * Whether the matrix is Hermitian.
   */
  isHermitian: boolean;

  /**
   * Maximum deviation from Hermitian property.
   */
  maxDeviation: number;

  /**
   * Matrix dimension.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Result from positive definiteness test.
 */
export interface IsPositiveDefiniteResult {
  /**
   * Whether the matrix is positive definite.
   */
  isPositiveDefinite: boolean;

  /**
   * Smallest eigenvalue (if computed).
   */
  minEigenvalue?: number;

  /**
   * Matrix dimension.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Options for orthogonality test.
 */
export interface IsOrthogonalOptions {
  /**
   * Tolerance for comparison.
   * @default 1e-10
   */
  tol?: number;
}

/**
 * Result from orthogonality test.
 */
export interface IsOrthogonalResult {
  /**
   * Whether the matrix is orthogonal (Q^T * Q = I).
   */
  isOrthogonal: boolean;

  /**
   * Maximum deviation from identity in Q^T * Q.
   */
  maxDeviation: number;

  /**
   * Matrix dimension.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Options for unitarity test.
 */
export interface IsUnitaryOptions {
  /**
   * Tolerance for comparison.
   * @default 1e-10
   */
  tol?: number;
}

/**
 * Result from unitarity test.
 */
export interface IsUnitaryResult {
  /**
   * Whether the matrix is unitary (Q^H * Q = I).
   */
  isUnitary: boolean;

  /**
   * Maximum deviation from identity in Q^H * Q.
   */
  maxDeviation: number;

  /**
   * Matrix dimension.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Options for singularity test.
 */
export interface IsSingularOptions {
  /**
   * Tolerance for considering a value as zero.
   * @default machine epsilon * max dimension
   */
  tol?: number;
}

/**
 * Result from singularity test.
 */
export interface IsSingularResult {
  /**
   * Whether the matrix is singular (non-invertible).
   */
  isSingular: boolean;

  /**
   * Smallest singular value.
   */
  minSingularValue: number;

  /**
   * Reciprocal condition number (small = ill-conditioned).
   */
  rcond: number;

  /**
   * Matrix dimension.
   */
  n: number;

  /**
   * Whether the operation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}
