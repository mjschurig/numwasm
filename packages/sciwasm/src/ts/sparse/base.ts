import type { SparseFormat } from './types.js';

export abstract class SparseMatrix {
  abstract get shape(): [number, number];
  abstract get format(): SparseFormat;
  abstract get nnz(): number;

  get ndim(): number {
    return 2;
  }

  get dtype(): string {
    return 'float64';
  }

  abstract toArray(): number[][];
  abstract tocsr(): SparseMatrix;
  abstract tocsc(): SparseMatrix;
  abstract diagonal(k?: number): Float64Array;
  abstract copy(): SparseMatrix;
  abstract transpose(): SparseMatrix;
}
