/**
 * NumJS Linear Algebra Module
 *
 * Re-exports from submodules for convenience.
 */

// Types and utilities
export { LinAlgError, getLinalgDtype, prepareMatrix, assertSquare } from "./types.js";
export type {
  EigResult,
  SVDResult,
  QRResult,
  LstsqResult,
  SlogdetResult,
} from "./types.js";

// Matrix products
export {
  matmul,
  dot,
  vdot,
  inner,
  outer,
  tensordot,
  multi_dot,
  kron,
  cross,
} from "./products.js";

// Decompositions
export {
  cholesky,
  qr,
  svd,
  svdvals,
  eig,
  eigvals,
  eigh,
  eigvalsh,
} from "./decomposition.js";

// Solving and inverting
export {
  solve,
  inv,
  pinv,
  lstsq,
  matrix_power,
  tensorsolve,
  tensorinv,
} from "./solve.js";

// Norms and numbers
export {
  norm,
  det,
  slogdet,
  matrix_rank,
  trace,
  cond,
  matrix_norm,
  vector_norm,
} from "./norms.js";

// Re-export linalg namespace object
import { LinAlgError } from "./types.js";
import { matmul, dot, vdot, inner, outer, tensordot, multi_dot, kron, cross } from "./products.js";
import { cholesky, qr, svd, svdvals, eig, eigvals, eigh, eigvalsh } from "./decomposition.js";
import { solve, inv, pinv, lstsq, matrix_power, tensorsolve, tensorinv } from "./solve.js";
import { norm, det, slogdet, matrix_rank, trace, cond, matrix_norm, vector_norm } from "./norms.js";

export const linalg = {
  // Error class
  LinAlgError,

  // Matrix products
  dot,
  vdot,
  inner,
  outer,
  matmul,

  // Decompositions
  cholesky,
  qr,
  svd,
  svdvals,

  // Eigenvalues
  eig,
  eigh,
  eigvals,
  eigvalsh,

  // Norms & Numbers
  norm,
  det,
  slogdet,
  matrix_rank,
  trace,
  cond,

  // Solving & Inverting
  solve,
  lstsq,
  inv,
  pinv,

  // Matrix Operations
  matrix_power,

  // Advanced Linear Algebra
  tensordot,
  multi_dot,
  kron,
  cross,
  tensorsolve,
  tensorinv,
  matrix_norm,
  vector_norm,
};
