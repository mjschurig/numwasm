/**
 * Sparse linear algebra module
 *
 * Provides iterative solvers and utilities for sparse linear systems.
 */

// Types
export type {
  IterativeSolverResult,
  IterativeSolverOptions,
  GMRESOptions,
  LinearOperatorLike,
  LinearOperator as ILinearOperator,
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
