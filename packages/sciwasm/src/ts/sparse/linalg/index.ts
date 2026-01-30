/**
 * Sparse linear algebra module
 *
 * Provides iterative and direct solvers, eigensolvers, and utilities
 * for sparse linear systems.
 */

// Types
export type {
  IterativeSolverResult,
  IterativeSolverOptions,
  GMRESOptions,
  LinearOperatorLike,
  LinearOperator as ILinearOperator,
  // Eigenvalue types
  EigshOptions,
  EigsOptions,
  EigsResult,
  SvdsOptions,
  SvdsResult,
  // Matrix exponential types
  ExpmOptions,
  ExpmMultiplyOptions,
  ExpmMultiplyResult,
} from './types.js';

// LinearOperator
export {
  LinearOperator,
  IdentityOperator,
  aslinearoperator,
  isLinearOperator,
} from './interface.js';

// Norm
export { norm } from './norm.js';

// Iterative solvers
export { cg, bicgstab, gmres } from './iterative.js';

// Direct solvers (SuperLU-based)
export {
  spsolve,
  splu,
  spilu,
  inv,
} from './direct.js';

export type {
  SpsolveOptions,
  SpLUOptions,
  SpLUResult,
  SpILUOptions,
  SpILUResult,
  InvOptions,
  ColPermSpec,
  TransSpec,
} from './direct.js';

// Eigenvalue solvers (ARPACK-based)
export { eigsh, eigs, svds } from './eigen.js';

// Matrix exponential
export { expm, expm_multiply } from './expm.js';
