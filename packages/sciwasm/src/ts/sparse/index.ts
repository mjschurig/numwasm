/**
 * Sparse matrix data structures and operations.
 * @module sparse
 */

import { NotImplementedError } from '../errors.js';

/** Sparse matrix interface. */
export interface SparseMatrix {
  readonly shape: [number, number];
  readonly nnz: number;
  toarray(): number[][];
  todense(): number[][];
  dot(other: SparseMatrix | number[][]): SparseMatrix;
}

/**
 * Compressed Sparse Row matrix.
 * Mirrors scipy.sparse.csr_matrix.
 */
export function csr_matrix(
  _data: number[][] | { data: number[]; indices: number[]; indptr: number[]; shape: [number, number] },
): SparseMatrix {
  throw new NotImplementedError('sciwasm.sparse.csr_matrix');
}

/**
 * Compressed Sparse Column matrix.
 * Mirrors scipy.sparse.csc_matrix.
 */
export function csc_matrix(
  _data: number[][] | { data: number[]; indices: number[]; indptr: number[]; shape: [number, number] },
): SparseMatrix {
  throw new NotImplementedError('sciwasm.sparse.csc_matrix');
}

/**
 * Sparse identity matrix.
 * Mirrors scipy.sparse.eye.
 */
export function eye(_m: number, _n?: number, _k?: number): SparseMatrix {
  throw new NotImplementedError('sciwasm.sparse.eye');
}

/**
 * Construct a sparse matrix from diagonals.
 * Mirrors scipy.sparse.diags.
 */
export function diags(_diagonals: number[][], _offsets?: number[], _shape?: [number, number]): SparseMatrix {
  throw new NotImplementedError('sciwasm.sparse.diags');
}
