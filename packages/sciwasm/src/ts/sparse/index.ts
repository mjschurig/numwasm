// Import sparse formats first to trigger factory registration
export { CSRMatrix, csr_matrix } from './csr.js';
export { CSCMatrix, csc_matrix } from './csc.js';
export { COOMatrix, coo_matrix } from './coo.js';

export { SparseMatrix } from './base.js';
export { CompressedSparseMatrix } from './compressed.js';
export { eye, diags } from './construct.js';
export type { SparseFormat, SparseConstructorArrays, COOConstructorArrays, SparseMatrixOptions } from './types.js';
