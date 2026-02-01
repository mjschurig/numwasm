/**
 * NumJS Linear Algebra Module
 *
 * Re-exports from linalg/ for backwards compatibility.
 */

// Types and utilities
export { LinAlgError } from "./linalg/types.js";
export type {
  EigResult,
  SVDResult,
  QRResult,
  LstsqResult,
  SlogdetResult,
} from "./linalg/types.js";

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
  // NumPy 2.0
  matrix_transpose,
  vecdot,
  matvec,
  vecmat,
} from "./linalg/products.js";

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
} from "./linalg/decomposition.js";

// Solving and inverting
export {
  solve,
  inv,
  pinv,
  lstsq,
  matrix_power,
  tensorsolve,
  tensorinv,
} from "./linalg/solve.js";

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
} from "./linalg/norms.js";

// linalg namespace object
export { linalg } from "./linalg/index.js";
