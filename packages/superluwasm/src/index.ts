/**
 * # SuperLUWasm
 *
 * WebAssembly bindings for SuperLU, a library for the direct solution of
 * large, sparse, nonsymmetric systems of linear equations.
 *
 * ## Overview
 *
 * SuperLU uses LU factorization with partial pivoting and triangular
 * solves to compute the solution to sparse linear systems Ax = b.
 * Key features include:
 *
 * - **Sparse LU Factorization**: Efficient factorization of sparse matrices
 * - **Multiple Precisions**: Single (float32), double (float64), and complex variants
 * - **Fill-Reducing Orderings**: Column permutations to minimize fill-in
 * - **Numerical Stability**: Row permutations and equilibration
 * - **Iterative Refinement**: Improve solution accuracy
 * - **ILU Preconditioning**: Incomplete LU for iterative solvers
 *
 * ## Quick Start
 *
 * ```typescript
 * import { loadSuperLUModule, solveSparseCSC } from 'superluwasm';
 *
 * // 1. Load the WASM module (do this once at startup)
 * await loadSuperLUModule();
 *
 * // 2. Create a sparse matrix in CSC format
 * const A = {
 *   m: 3, n: 3,
 *   values: new Float64Array([4, 1, 1, 3, 2, 5]),
 *   rowIndices: new Int32Array([0, 1, 1, 2, 0, 2]),
 *   colPointers: new Int32Array([0, 2, 4, 6]),
 *   dtype: 'float64' as const
 * };
 *
 * // 3. Solve Ax = b
 * const b = new Float64Array([1, 2, 3]);
 * const { x, statistics } = solveSparseCSC(A, b);
 *
 * console.log('Solution:', x);
 * console.log('Time:', statistics.totalTime, 'ms');
 * ```
 *
 * ## Matrix Formats
 *
 * SuperLU uses compressed sparse formats:
 *
 * - **CSC (Compressed Sparse Column)**: Native format, most efficient
 * - **CSR (Compressed Sparse Row)**: Supported via automatic conversion
 *
 * ## High-Level API
 *
 * ### Solvers
 * - {@link solveSparseCSC} - Solve Ax=b with CSC matrix
 * - {@link solveSparseCSR} - Solve Ax=b with CSR matrix
 * - {@link solveSparseTranspose} - Solve A^T x = b
 * - {@link solveSparseConjugateTranspose} - Solve A^H x = b (complex)
 * - {@link solveSparseExpert} - Full control over all options
 *
 * ### Factorization
 * - {@link sparseLU} - Compute LU factorization
 * - {@link sparseILU} - Compute incomplete LU for preconditioning
 *
 * ### Module Loading
 * - {@link loadSuperLUModule} - Initialize the WASM module
 * - {@link getSuperLUModule} - Get the loaded module
 * - {@link isSuperLULoaded} - Check if loaded
 * - {@link configureSuperLU} - Configure WASM location
 *
 * @packageDocumentation
 * @module superluwasm
 */

// ============================================================
// Module Loading API
// ============================================================

export {
  loadSuperLUModule,
  getSuperLUModule,
  isSuperLULoaded,
  isSuperLULoading,
  resetSuperLUModule,
  configureSuperLU,
  getSuperLUConfig,
  type SuperLULoadConfig,
} from './ts/core/loader.js';

// ============================================================
// High-Level Solvers
// ============================================================

export {
  solveSparseCSC,
  solveSparseCSR,
  solveSparseTranspose,
  solveSparseConjugateTranspose,
} from './ts/solvers/direct.js';

export { solveSparseExpert } from './ts/solvers/expert.js';

// ============================================================
// Factorization
// ============================================================

export { sparseLU } from './ts/factorization/lu.js';
export { sparseILU } from './ts/factorization/ilu.js';

// ============================================================
// High-Level Types
// ============================================================

export type {
  // Array types
  RealArray,
  ComplexArray,
  IntArray,
  DataType,

  // Sparse matrix types
  SparseMatrixCSC,
  SparseMatrixCSR,
  SparseMatrix,

  // Options
  SolveOptions,
  ExpertSolveOptions,
  LUFactorizationOptions,
  ILUFactorizationOptions,

  // Results
  SolveResult,
  SolveStatistics,
  LUFactorization,
  ILUFactorization,
  LUFactorHandle,
} from './ts/high-level-types.js';

export {
  isSparseMatrixCSC,
  isSparseMatrixCSR,
  getNnz,
  inferDataType,
} from './ts/high-level-types.js';

// ============================================================
// Low-Level Module Types
// ============================================================

export type {
  SuperLUModule,
  SuperLUModuleFactory,
  SuperLUModuleOptions,
} from './ts/types.js';

// ============================================================
// Enumerations
// ============================================================

export {
  SuperLUFactOption,
  SuperLURowPerm,
  SuperLUColPerm,
  SuperLUTrans,
  SuperLUILUDrop,
  SuperLUMiluT,
  SuperLUNorm,
  SuperLUStype,
  SuperLUDtype,
  SuperLUMtype,
} from './ts/types.js';

// ============================================================
// Submodule Namespaces
// ============================================================

export * as core from './ts/core/index.js';
export * as solvers from './ts/solvers/index.js';
export * as factorization from './ts/factorization/index.js';
