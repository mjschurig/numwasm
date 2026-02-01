/**
 * High-Level Type Definitions for SuperLU WASM
 *
 * These types provide a user-friendly interface for working with sparse
 * matrices and linear solvers, abstracting away the low-level WASM details.
 *
 * @module high-level-types
 */

import type {
  SuperLUColPerm,
  SuperLURowPerm,
  SuperLUILUDrop,
  SuperLUMiluT,
  SuperLUNorm,
} from './types.js';

// ============================================================
// Numeric Array Types
// ============================================================

/**
 * Numeric array types accepted for real-valued data.
 * Supports both JavaScript arrays and typed arrays for flexibility.
 */
export type RealArray = number[] | Float64Array | Float32Array;

/**
 * Complex array type - interleaved real/imaginary pairs.
 * For a complex vector of length n, the array has 2n elements:
 * [re_0, im_0, re_1, im_1, ..., re_{n-1}, im_{n-1}]
 */
export type ComplexArray = number[] | Float64Array | Float32Array;

/**
 * Integer array types for indices and pointers.
 */
export type IntArray = number[] | Int32Array;

// ============================================================
// Data Type Enumeration
// ============================================================

/**
 * Supported numeric data types for sparse matrices.
 */
export type DataType = 'float32' | 'float64' | 'complex64' | 'complex128';

// ============================================================
// Sparse Matrix Types
// ============================================================

/**
 * Sparse matrix in Compressed Sparse Column (CSC) format.
 *
 * CSC is the native format for SuperLU. For a matrix with m rows, n columns,
 * and nnz nonzeros:
 * - `colPointers[j]` to `colPointers[j+1]-1` gives the indices into
 *   `rowIndices` and `values` for column j
 * - `rowIndices[k]` is the row index of the k-th nonzero
 * - `values[k]` is the value of the k-th nonzero
 *
 * @example
 * ```typescript
 * // 3x3 matrix:
 * // [ 1  0  2 ]
 * // [ 0  3  0 ]
 * // [ 4  0  5 ]
 *
 * const A: SparseMatrixCSC = {
 *   m: 3,
 *   n: 3,
 *   values: new Float64Array([1, 4, 3, 2, 5]),
 *   rowIndices: new Int32Array([0, 2, 1, 0, 2]),
 *   colPointers: new Int32Array([0, 2, 3, 5]),
 *   dtype: 'float64'
 * };
 * ```
 */
export interface SparseMatrixCSC {
  /** Number of rows */
  m: number;
  /** Number of columns */
  n: number;
  /** Nonzero values (length = nnz) */
  values: RealArray | ComplexArray;
  /** Row indices for each nonzero (length = nnz) */
  rowIndices: IntArray;
  /** Column pointers (length = n + 1) */
  colPointers: IntArray;
  /** Data type of values */
  dtype?: DataType;
}

/**
 * Sparse matrix in Compressed Sparse Row (CSR) format.
 *
 * CSR is an alternative format that is row-oriented. For a matrix with
 * m rows, n columns, and nnz nonzeros:
 * - `rowPointers[i]` to `rowPointers[i+1]-1` gives the indices into
 *   `colIndices` and `values` for row i
 * - `colIndices[k]` is the column index of the k-th nonzero
 * - `values[k]` is the value of the k-th nonzero
 *
 * @example
 * ```typescript
 * // Same 3x3 matrix in CSR:
 * const A: SparseMatrixCSR = {
 *   m: 3,
 *   n: 3,
 *   values: new Float64Array([1, 2, 3, 4, 5]),
 *   colIndices: new Int32Array([0, 2, 1, 0, 2]),
 *   rowPointers: new Int32Array([0, 2, 3, 5]),
 *   dtype: 'float64'
 * };
 * ```
 */
export interface SparseMatrixCSR {
  /** Number of rows */
  m: number;
  /** Number of columns */
  n: number;
  /** Nonzero values (length = nnz) */
  values: RealArray | ComplexArray;
  /** Column indices for each nonzero (length = nnz) */
  colIndices: IntArray;
  /** Row pointers (length = m + 1) */
  rowPointers: IntArray;
  /** Data type of values */
  dtype?: DataType;
}

/**
 * Union type for any supported sparse matrix format.
 */
export type SparseMatrix = SparseMatrixCSC | SparseMatrixCSR;

// ============================================================
// Solve Options
// ============================================================

/**
 * Options for sparse linear solvers.
 *
 * These options control the behavior of the solve operation including
 * permutation strategies, equilibration, and iterative refinement.
 */
export interface SolveOptions {
  /**
   * Column permutation strategy for fill-in reduction.
   * @default SuperLUColPerm.COLAMD
   */
  columnPermutation?: SuperLUColPerm;

  /**
   * Row permutation strategy for numerical stability.
   * @default SuperLURowPerm.LargeDiag_MC64
   */
  rowPermutation?: SuperLURowPerm;

  /**
   * Whether to equilibrate (scale) the matrix for numerical stability.
   * Equilibration scales rows and columns to have similar magnitudes.
   * @default true
   */
  equilibrate?: boolean;

  /**
   * Whether to perform iterative refinement to improve accuracy.
   * Iterative refinement uses the computed solution to iteratively
   * improve accuracy using the original matrix.
   * @default false
   */
  iterativeRefinement?: boolean;

  /**
   * Maximum iterations for iterative refinement.
   * @default 5
   */
  maxRefinementIterations?: number;

  /**
   * Whether to compute and return the condition number estimate.
   * Computing the condition number adds some overhead.
   * @default false
   */
  computeConditionNumber?: boolean;

  /**
   * Whether to compute forward and backward error bounds.
   * @default false
   */
  computeErrorBounds?: boolean;

  /**
   * Norm type to use for condition number estimation.
   * @default SuperLUNorm.ONE_NORM
   */
  normType?: SuperLUNorm;

  /**
   * Whether to print statistics (timing, memory usage).
   * @default false
   */
  printStatistics?: boolean;
}

/**
 * Extended options for the expert solver interface.
 * Provides full control over all SuperLU parameters.
 */
export interface ExpertSolveOptions extends SolveOptions {
  /**
   * User-supplied column permutation array.
   * Only used when columnPermutation is MY_PERMC.
   */
  userColumnPerm?: IntArray;

  /**
   * User-supplied row permutation array.
   * Only used when rowPermutation is MY_PERMR.
   */
  userRowPerm?: IntArray;

  /**
   * Pre-computed L and U factors for reuse.
   * When provided, skips factorization and uses these factors.
   */
  reuseFactorization?: LUFactorization;

  /**
   * Workspace size in bytes. If 0, SuperLU allocates automatically.
   * @default 0
   */
  workspaceSize?: number;
}

// ============================================================
// LU Factorization Options
// ============================================================

/**
 * Options for LU factorization.
 */
export interface LUFactorizationOptions {
  /**
   * Column permutation strategy for fill-in reduction.
   * @default SuperLUColPerm.COLAMD
   */
  columnPermutation?: SuperLUColPerm;

  /**
   * Row permutation strategy for numerical stability.
   * @default SuperLURowPerm.LargeDiag_MC64
   */
  rowPermutation?: SuperLURowPerm;

  /**
   * Whether to equilibrate the matrix before factorization.
   * @default true
   */
  equilibrate?: boolean;

  /**
   * Panel size for supernodal factorization.
   * @default (automatic)
   */
  panelSize?: number;

  /**
   * Relaxation parameter for supernode amalgamation.
   * @default (automatic)
   */
  relaxation?: number;

  /**
   * Whether to print statistics.
   * @default false
   */
  printStatistics?: boolean;
}

/**
 * Options for Incomplete LU (ILU) factorization.
 */
export interface ILUFactorizationOptions extends LUFactorizationOptions {
  /**
   * Drop tolerance for ILU. Elements smaller than this
   * (relative to the pivot) are dropped.
   * @default 1e-4
   */
  dropTolerance?: number;

  /**
   * Drop rule controlling which elements are dropped.
   * @default SuperLUILUDrop.DROP_BASIC
   */
  dropRule?: SuperLUILUDrop;

  /**
   * Modified ILU type.
   * @default SuperLUMiluT.SILU
   */
  miluType?: SuperLUMiluT;

  /**
   * Maximum fill-in factor. The ILU factors will have at most
   * fillFactor * nnz(A) nonzeros.
   * @default 10.0
   */
  fillFactor?: number;

  /**
   * Gamma parameter for ILU (used with some drop rules).
   * @default 1.0
   */
  gamma?: number;
}

// ============================================================
// Result Types
// ============================================================

/**
 * Statistics from a solve or factorization operation.
 */
export interface SolveStatistics {
  /** Time spent in factorization (milliseconds) */
  factorizationTime?: number;
  /** Time spent in triangular solve (milliseconds) */
  solveTime?: number;
  /** Time spent in iterative refinement (milliseconds) */
  refinementTime?: number;
  /** Total operation time (milliseconds) */
  totalTime: number;
  /** Number of nonzeros in L factor */
  nnzL?: number;
  /** Number of nonzeros in U factor */
  nnzU?: number;
  /** Memory used by L and U factors (bytes) */
  memoryUsed?: number;
  /** Number of floating point operations */
  flops?: number;
  /** Number of refinement iterations performed */
  refinementIterations?: number;
}

/**
 * Result of solving a sparse linear system.
 */
export interface SolveResult {
  /** Solution vector x */
  x: Float64Array | Float32Array;

  /** Operation statistics */
  statistics: SolveStatistics;

  /**
   * Reciprocal condition number estimate (1/κ(A)).
   * Values close to 0 indicate ill-conditioning.
   * Only present if computeConditionNumber was true.
   */
  rcond?: number;

  /**
   * Forward error bound for each right-hand side.
   * ferr[j] is an upper bound on |x_true - x_computed| / |x_true|
   * for the j-th RHS.
   * Only present if computeErrorBounds was true.
   */
  forwardError?: Float64Array;

  /**
   * Backward error bound for each right-hand side.
   * berr[j] is the relative residual |A*x - b| / (|A|*|x| + |b|)
   * for the j-th RHS.
   * Only present if computeErrorBounds was true.
   */
  backwardError?: Float64Array;

  /**
   * Reciprocal pivot growth factor.
   * Small values may indicate numerical instability.
   */
  pivotGrowth?: number;

  /**
   * The computed row permutation.
   */
  rowPermutation?: Int32Array;

  /**
   * The computed column permutation.
   */
  columnPermutation?: Int32Array;

  /**
   * Equilibration type applied:
   * - 'N': No equilibration
   * - 'R': Row equilibration
   * - 'C': Column equilibration
   * - 'B': Both row and column equilibration
   */
  equilibration?: 'N' | 'R' | 'C' | 'B';
}

/**
 * LU factorization result.
 *
 * Contains the L and U factors along with permutation information
 * that can be reused for solving multiple systems with the same matrix.
 */
export interface LUFactorization {
  /**
   * Lower triangular factor L in SuperLU's internal format.
   * This is an opaque handle - use with solveReusedFactorization.
   */
  L: LUFactorHandle;

  /**
   * Upper triangular factor U in SuperLU's internal format.
   * This is an opaque handle - use with solveReusedFactorization.
   */
  U: LUFactorHandle;

  /** Row permutation vector */
  rowPermutation: Int32Array;

  /** Column permutation vector */
  columnPermutation: Int32Array;

  /** Elimination tree */
  eliminationTree: Int32Array;

  /** Row scaling factors (if equilibration was applied) */
  rowScale?: Float64Array;

  /** Column scaling factors (if equilibration was applied) */
  columnScale?: Float64Array;

  /** Equilibration type applied */
  equilibration: 'N' | 'R' | 'C' | 'B';

  /** Matrix dimensions */
  m: number;
  n: number;

  /** Data type of the factorization */
  dtype: DataType;

  /** Statistics from factorization */
  statistics: SolveStatistics;

  /**
   * Free the memory used by this factorization.
   * Call this when the factorization is no longer needed.
   */
  dispose: () => void;
}

/**
 * Opaque handle to an LU factor stored in WASM memory.
 * Do not access the internal pointer directly.
 */
export interface LUFactorHandle {
  /** @internal WASM pointer - do not use directly */
  readonly _ptr: number;
  /** @internal Whether this handle has been disposed */
  readonly _disposed: boolean;
}

/**
 * ILU factorization result.
 * Similar to LUFactorization but for incomplete factorization.
 */
export interface ILUFactorization extends Omit<LUFactorization, 'L' | 'U'> {
  /** Incomplete lower triangular factor */
  L: LUFactorHandle;

  /** Incomplete upper triangular factor */
  U: LUFactorHandle;

  /** Actual fill-in achieved (nnz(L+U) / nnz(A)) */
  fillIn: number;
}

// ============================================================
// Type Guards
// ============================================================

/**
 * Check if a sparse matrix is in CSC format.
 */
export function isSparseMatrixCSC(matrix: SparseMatrix): matrix is SparseMatrixCSC {
  return 'colPointers' in matrix && 'rowIndices' in matrix;
}

/**
 * Check if a sparse matrix is in CSR format.
 */
export function isSparseMatrixCSR(matrix: SparseMatrix): matrix is SparseMatrixCSR {
  return 'rowPointers' in matrix && 'colIndices' in matrix;
}

/**
 * Get the number of nonzeros in a sparse matrix.
 */
export function getNnz(matrix: SparseMatrix): number {
  return matrix.values.length;
}

/**
 * Infer the data type from an array.
 */
export function inferDataType(arr: RealArray | ComplexArray): DataType {
  if (arr instanceof Float32Array) return 'float32';
  if (arr instanceof Float64Array) return 'float64';
  // For number[], default to float64
  return 'float64';
}
