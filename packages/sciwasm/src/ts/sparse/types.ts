export type SparseFormat = 'csr' | 'csc';

export interface SparseConstructorArrays {
  data: Float64Array;
  indices: Int32Array;
  indptr: Int32Array;
}

export interface SparseMatrixOptions {
  shape?: [number, number];
  copy?: boolean;
}
