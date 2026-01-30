/**
 * LAWasm - LAPACK + BLAS WebAssembly Module
 *
 * This package provides a WebAssembly build of LAPACK (Linear Algebra PACKage)
 * with BLAS (Basic Linear Algebra Subprograms) bundled internally.
 *
 * LAPACK provides routines for:
 * - Solving systems of linear equations
 * - Least squares solutions
 * - Eigenvalue problems
 * - Singular value decomposition
 * - Matrix factorizations (LU, Cholesky, QR, etc.)
 *
 * @example
 * ```typescript
 * import { loadLAPACKModule, getLAPACKModule } from 'lawasm';
 *
 * // Load the WASM module
 * await loadLAPACKModule();
 *
 * // Get the module instance
 * const lapack = getLAPACKModule();
 *
 * // Use LAPACK routines via the module interface
 * // Note: All operations use column-major order (Fortran convention)
 * ```
 *
 * @packageDocumentation
 */

// Re-export loader functions
export {
  loadLAPACKModule,
  getLAPACKModule,
  isLAPACKLoaded,
  resetLAPACKModule,
  configureLAPACK,
  type LAPACKLoadConfig,
} from './ts/loader.js';

// Re-export types
export type {
  LAPACKModule,
  LAPACKModuleFactory,
} from './ts/types.js';
