/**
 * ARWasm - ARPACK WebAssembly Module
 *
 * This package provides a WebAssembly build of ARPACK (ARnoldi PACKage),
 * a collection of Fortran77 subroutines designed to solve large scale
 * eigenvalue problems.
 *
 * ARPACK provides routines for:
 * - Real symmetric eigenvalue problems (dsaupd/dseupd, ssaupd/sseupd)
 * - Real nonsymmetric eigenvalue problems (dnaupd/dneupd, snaupd/sneupd)
 * - Complex eigenvalue problems (cnaupd/cneupd, znaupd/zneupd)
 *
 * ## High-Level API
 *
 * For most use cases, use the high-level functions that handle all the
 * complexity of the reverse communication interface:
 *
 * @example
 * ```typescript
 * import { eigs } from 'arwasm';
 *
 * // Find the 6 smallest eigenvalues of a symmetric matrix
 * const n = 100;
 * const matvec = (x: Float64Array): Float64Array => {
 *   const y = new Float64Array(n);
 *   for (let i = 0; i < n; i++) {
 *     y[i] = 2 * x[i];
 *     if (i > 0) y[i] -= x[i - 1];
 *     if (i < n - 1) y[i] -= x[i + 1];
 *   }
 *   return y;
 * };
 *
 * const result = await eigs(matvec, n, 6, { which: 'SM' });
 * console.log('Smallest eigenvalues:', result.eigenvalues);
 * ```
 *
 * ## Low-Level API
 *
 * For advanced use cases, you can use the raw WASM module interface:
 *
 * @example
 * ```typescript
 * import { loadARPACKModule, getARPACKModule } from 'arwasm';
 *
 * await loadARPACKModule();
 * const arpack = getARPACKModule();
 * // Use reverse communication interface directly
 * ```
 *
 * @packageDocumentation
 */

// ============================================================
// HIGH-LEVEL API (recommended for most users)
// ============================================================

// High-level eigenvalue solvers
export { eigs } from './ts/eigs.js';
export { eign } from './ts/eign.js';
export { eigsh, isEigsResult, isEignResult } from './ts/eigsh.js';

// High-level type definitions
export type {
  RealArray,
  MatVecFunction,
  OperatorFunction,
  BMatVecFunction,
  EigOptionsBase,
  EigsOptions,
  EignOptions,
  EigsResult,
  EignResult,
  EigshOptions,
  ProblemType,
} from './ts/high-level-types.js';

// Helper utilities
export {
  dsaupdWorklSize,
  dnaupdWorklSize,
  dneupdWorkevSize,
  defaultNcv,
  getDsaupdMessage,
  getDnaupdMessage,
  getDseupdMessage,
  getDneupdMessage,
} from './ts/helpers.js';

// ============================================================
// LOW-LEVEL API (for advanced usage)
// ============================================================

// Module loader functions
export {
  loadARPACKModule,
  getARPACKModule,
  isARPACKLoaded,
  resetARPACKModule,
  configureARPACK,
  type ARPACKLoadConfig,
} from './ts/loader.js';

// Low-level module types
export type {
  ARPACKModule,
  ARPACKModuleFactory,
  ARPACKModuleOptions,
  WhichSymmetric,
  WhichNonSymmetric,
} from './ts/types.js';

// Error code mappings and constants
export {
  DSAUPD_ERRORS,
  DSEUPD_ERRORS,
  DNAUPD_ERRORS,
  DNEUPD_ERRORS,
  WHICH_SYMMETRIC,
  WHICH_NONSYMMETRIC,
} from './ts/types.js';
