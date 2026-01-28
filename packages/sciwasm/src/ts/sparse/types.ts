export type SparseFormat = 'csr' | 'csc' | 'coo';

export interface SparseConstructorArrays {
  data: Float64Array;
  indices: Int32Array;
  indptr: Int32Array;
}

export interface COOConstructorArrays {
  row: Int32Array;
  col: Int32Array;
  data: Float64Array;
}

export interface SparseMatrixOptions {
  shape?: [number, number];
  copy?: boolean;
}
