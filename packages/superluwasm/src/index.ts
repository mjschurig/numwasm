/**
 * SuperLUWasm - SuperLU WebAssembly Module
 *
 * This package provides a WebAssembly build of SuperLU, a library for
 * solving sparse linear systems using LU factorization with partial pivoting.
 *
 * SuperLU provides routines for:
 * - Sparse LU factorization
 * - Solving sparse linear systems Ax = b
 * - Column and row permutations for numerical stability
 * - Memory-efficient storage of sparse matrices
 *
 * @example
 * ```typescript
 * import { loadSuperLUModule, getSuperLUModule } from 'superluwasm';
 *
 * // Load the WASM module
 * await loadSuperLUModule();
 *
 * // Get the module instance
 * const superlu = getSuperLUModule();
 *
 * // Use SuperLU routines via the module interface
 * ```
 *
 * @packageDocumentation
 */

// Re-export loader functions
export {
  loadSuperLUModule,
  getSuperLUModule,
  isSuperLULoaded,
  resetSuperLUModule,
  configureSuperLU,
  type SuperLULoadConfig,
} from './ts/loader.js';

// Re-export types
export type {
  SuperLUModule,
  SuperLUModuleFactory,
} from './ts/types.js';

// Re-export enums for configuration
export {
  SuperLUFactOption,
  SuperLURowPerm,
  SuperLUColPerm,
  SuperLUTrans,
  SuperLUILUDrop,
  SuperLUMiluT,
  SuperLUStype,
  SuperLUDtype,
  SuperLUMtype,
} from './ts/types.js';
