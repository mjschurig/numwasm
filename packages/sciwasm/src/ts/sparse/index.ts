// Import sparse formats first to trigger factory registration
export { CSRMatrix, csr_matrix } from './csr.js';
export { CSCMatrix, csc_matrix } from './csc.js';
export { COOMatrix, coo_matrix } from './coo.js';
export { DOKMatrix, dok_matrix } from './dok.js';
export { LILMatrix, lil_matrix } from './lil.js';
export { DIAMatrix, dia_matrix } from './dia.js';
export { BSRMatrix, bsr_matrix } from './bsr.js';

export { SparseMatrix } from './base.js';
export { CompressedSparseMatrix } from './compressed.js';
export {
  eye,
  diags,
  issparse,
  tril,
  triu,
  hstack,
  vstack,
  block_diag,
  kron,
  kronsum,
  random,
} from './construct.js';
export type {
  SparseFormat,
  SparseConstructorArrays,
  COOConstructorArrays,
  LILRow,
  LILConstructorArrays,
  DIAConstructorArrays,
  BSRConstructorArrays,
  SparseMatrixOptions,
} from './types.js';

// Linear algebra submodule
export * as linalg from './linalg/index.js';
export {
  LinearOperator,
  IdentityOperator,
  aslinearoperator,
  isLinearOperator,
  norm as linalg_norm,
  cg,
  bicgstab,
  gmres,
} from './linalg/index.js';
export type {
  IterativeSolverResult,
  IterativeSolverOptions,
  GMRESOptions,
  LinearOperatorLike,
  ILinearOperator,
} from './linalg/index.js';
