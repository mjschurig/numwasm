/**
 * Linear System Solver Types
 *
 * Type definitions for high-level linear system solving functions.
 */

import type { RealArray, ComplexArray } from '../helpers.js';

// ============================================================
// MATRIX INPUT TYPES
// ============================================================

/**
 * Matrix input type - supports both 2D arrays (row-major) and 1D typed arrays (column-major).
 *
 * For 2D arrays: A[i][j] = element at row i, column j
 * For 1D arrays: Assumed to be in column-major (Fortran) order
 */
export type Matrix = number[][] | RealArray;

/**
 * Vector input type.
 */
export type Vector = RealArray;

/**
 * Complex matrix input type - interleaved real/imaginary pairs.
 * For a complex value z = a + bi, stored as [a, b, ...].
 * Total length = 2 * (number of elements).
 */
export type ComplexMatrix = ComplexArray;

/**
 * Complex vector input type - interleaved real/imaginary pairs.
 */
export type ComplexVector = ComplexArray;

// ============================================================
// SOLVE OPTIONS
// ============================================================

/**
 * Options for the general solve function.
 */
export interface SolveOptions {
  /**
   * Whether to overwrite the input matrix A with the LU factors.
   * If false (default), a copy is made.
   * @default false
   */
  overwriteA?: boolean;

  /**
   * Whether to overwrite the input matrix B with the solution.
   * If false (default), a copy is made.
   * @default false
   */
  overwriteB?: boolean;

  /**
   * Whether to also return the LU factorization.
   * @default false
   */
  returnFactorization?: boolean;
}

/**
 * Options for triangular solve.
 */
export interface SolveTriangularOptions {
  /**
   * Whether A is upper triangular ('U') or lower triangular ('L').
   * @default 'U'
   */
  upper?: boolean;

  /**
   * Whether to transpose A before solving.
   * 'N' = no transpose, 'T' = transpose, 'C' = conjugate transpose (complex)
   * @default 'N'
   */
  trans?: 'N' | 'T' | 'C';

  /**
   * Whether A has unit diagonal (all diagonal elements are 1).
   * If true, diagonal elements are not referenced.
   * @default false
   */
  unitDiagonal?: boolean;

  /**
   * Whether to overwrite B with the solution.
   * @default false
   */
  overwriteB?: boolean;
}

/**
 * Options for symmetric solve.
 */
export interface SolveSymmetricOptions {
  /**
   * Whether the upper ('U') or lower ('L') triangle of A is stored.
   * @default 'U'
   */
  upper?: boolean;

  /**
   * Whether to overwrite A with the Cholesky factor.
   * @default false
   */
  overwriteA?: boolean;

  /**
   * Whether to overwrite B with the solution.
   * @default false
   */
  overwriteB?: boolean;

  /**
   * Whether to also return the Cholesky factorization.
   * @default false
   */
  returnFactorization?: boolean;
}

/**
 * Options for Hermitian solve (complex symmetric positive definite).
 */
export interface SolveHermitianOptions {
  /**
   * Whether the upper ('U') or lower ('L') triangle of A is stored.
   * @default 'U'
   */
  upper?: boolean;

  /**
   * Whether to overwrite A with the Cholesky factor.
   * @default false
   */
  overwriteA?: boolean;

  /**
   * Whether to overwrite B with the solution.
   * @default false
   */
  overwriteB?: boolean;
}

/**
 * Options for banded solve.
 */
export interface SolveBandedOptions {
  /**
   * Whether to overwrite the banded matrix with factors.
   * @default false
   */
  overwriteAB?: boolean;

  /**
   * Whether to overwrite B with the solution.
   * @default false
   */
  overwriteB?: boolean;
}

// ============================================================
// RESULT TYPES
// ============================================================

/**
 * Result from linear system solve.
 */
export interface SolveResult {
  /**
   * The solution matrix/vector X where AX = B.
   * Has the same shape as B.
   */
  x: Float64Array;

  /**
   * Number of right-hand sides (columns in X).
   */
  nrhs: number;

  /**
   * LAPACK info return code.
   * 0 = success, < 0 = illegal argument, > 0 = singular matrix
   */
  info: number;

  /**
   * Whether the computation was successful.
   */
  success: boolean;

  /**
   * Human-readable status message.
   */
  message: string;

  /**
   * LU factorization (if returnFactorization was true).
   * Contains L and U factors packed in the same array.
   */
  lu?: Float64Array;

  /**
   * Pivot indices from LU factorization (if returnFactorization was true).
   */
  ipiv?: Int32Array;
}

/**
 * Result from symmetric/Hermitian solve.
 */
export interface SolveSymmetricResult {
  /**
   * The solution matrix/vector X.
   */
  x: Float64Array;

  /**
   * Number of right-hand sides.
   */
  nrhs: number;

  /**
   * LAPACK info code.
   */
  info: number;

  /**
   * Success flag.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;

  /**
   * Cholesky factor (if returnFactorization was true).
   */
  cholesky?: Float64Array;
}

/**
 * Result from complex solve.
 */
export interface SolveComplexResult {
  /**
   * The solution - interleaved real/imaginary pairs.
   */
  x: Float64Array;

  /**
   * Number of right-hand sides.
   */
  nrhs: number;

  /**
   * LAPACK info code.
   */
  info: number;

  /**
   * Success flag.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Result from triangular solve.
 */
export interface SolveTriangularResult {
  /**
   * The solution matrix/vector X.
   */
  x: Float64Array;

  /**
   * Number of right-hand sides.
   */
  nrhs: number;

  /**
   * LAPACK info code.
   */
  info: number;

  /**
   * Success flag.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Result from tridiagonal solve.
 */
export interface SolveTridiagonalResult {
  /**
   * The solution matrix/vector X.
   */
  x: Float64Array;

  /**
   * Number of right-hand sides.
   */
  nrhs: number;

  /**
   * LAPACK info code.
   */
  info: number;

  /**
   * Success flag.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}
