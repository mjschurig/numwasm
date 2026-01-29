export type SparseFormat = 'csr' | 'csc' | 'coo' | 'lil' | 'dok' | 'dia' | 'bsr';

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

// LIL (List of Lists) format types
export interface LILRow {
  indices: number[];
  data: number[];
}

export interface LILConstructorArrays {
  rows: LILRow[];
}

// DIA (Diagonal) format types
export interface DIAConstructorArrays {
  data: Float64Array[] | number[][];
  offsets: Int32Array | number[];
}

// BSR (Block Sparse Row) format types
export interface BSRConstructorArrays {
  data: Float64Array;
  indices: Int32Array;
  indptr: Int32Array;
  blocksize: [number, number];
}
