import type { CompressedSparseMatrix } from './compressed.js';
import { SparseMatrix } from './base.js';
import { createCSR, createCOO } from './_factory.js';
import type { SparseFormat } from './types.js';
import { getWasmModule } from '../wasm-loader.js';
import type { SciWasmModule } from '../wasm-types.js';

/**
 * Helper: copy a typed array into WASM heap and return the pointer.
 */
function toWasm(wasm: SciWasmModule, arr: Int32Array | Float64Array): number {
  const ptr = wasm._malloc(arr.byteLength);
  if (arr instanceof Float64Array) {
    wasm.HEAPF64.set(arr, ptr >> 3);
  } else {
    wasm.HEAP32.set(arr, ptr >> 2);
  }
  return ptr;
}

function fromWasmF64(wasm: SciWasmModule, ptr: number, len: number): Float64Array {
  const result = new Float64Array(len);
  result.set(wasm.HEAPF64.subarray(ptr >> 3, (ptr >> 3) + len));
  return result;
}

function fromWasmI32(wasm: SciWasmModule, ptr: number, len: number): Int32Array {
  const result = new Int32Array(len);
  result.set(wasm.HEAP32.subarray(ptr >> 2, (ptr >> 2) + len));
  return result;
}

/**
 * Convert input to COO format
 */
function toCOO(input: SparseMatrix | number[][]): { row: Int32Array; col: Int32Array; data: Float64Array; shape: [number, number] } {
  if (input instanceof SparseMatrix) {
    // Check if it's already a COO matrix
    if (input.format === 'coo') {
      return {
        row: (input as any)._row,
        col: (input as any)._col,
        data: (input as any)._data,
        shape: input.shape,
      };
    }

    // Convert to CSR first, then expand to COO arrays
    const csr = input.tocsr();
    const indptr = (csr as any)._indptr as Int32Array;
    const indices = (csr as any)._indices as Int32Array;
    const data = (csr as any)._data as Float64Array;
    const [nrow] = csr.shape;
    const nnz = csr.nnz;

    // Expand CSR to COO (row indices from indptr)
    const row = new Int32Array(nnz);
    let idx = 0;
    for (let i = 0; i < nrow; i++) {
      const rowStart = indptr[i];
      const rowEnd = indptr[i + 1];
      for (let j = rowStart; j < rowEnd; j++) {
        row[idx++] = i;
      }
    }

    return {
      row,
      col: new Int32Array(indices),
      data: new Float64Array(data),
      shape: csr.shape,
    };
  }
  // Dense 2D array
  const dense = input;
  const nrow = dense.length;
  const ncol = nrow > 0 ? dense[0].length : 0;
  const rows: number[] = [];
  const cols: number[] = [];
  const vals: number[] = [];
  for (let i = 0; i < nrow; i++) {
    for (let j = 0; j < ncol; j++) {
      if (dense[i][j] !== 0) {
        rows.push(i);
        cols.push(j);
        vals.push(dense[i][j]);
      }
    }
  }
  return {
    row: new Int32Array(rows),
    col: new Int32Array(cols),
    data: new Float64Array(vals),
    shape: [nrow, ncol] as [number, number],
  };
}

/**
 * Convert COO arrays to requested format
 */
function cooToFormat(
  row: Int32Array,
  col: Int32Array,
  data: Float64Array,
  shape: [number, number],
  format: SparseFormat
): SparseMatrix {
  const coo = createCOO({ row, col, data }, shape);
  switch (format) {
    case 'coo':
      return coo;
    case 'csr':
      return coo.tocsr();
    case 'csc':
      return coo.tocsc();
    default:
      return coo.tocsr();
  }
}

/**
 * Sparse identity matrix.
 *
 * @param m - Number of rows
 * @param n - Number of columns (defaults to m)
 * @param k - Diagonal offset (default 0 = main diagonal, positive = above, negative = below)
 * @param format - Output format: 'csr' or 'csc' (default 'csr')
 * @returns Sparse matrix with ones on the k-th diagonal
 */
export function eye(
  m: number,
  n?: number,
  k: number = 0,
  format: 'csr' | 'csc' = 'csr'
): CompressedSparseMatrix {
  if (n === undefined) n = m;

  const firstRow = k >= 0 ? 0 : -k;
  const firstCol = k >= 0 ? k : 0;
  const diagLen = Math.max(0, Math.min(m - firstRow, n - firstCol));

  // Build CSR arrays directly
  const indptr = new Int32Array(m + 1);
  const indices = new Int32Array(diagLen);
  const data = new Float64Array(diagLen);

  let nnzIdx = 0;
  for (let i = 0; i < m; i++) {
    indptr[i] = nnzIdx;
    const col = i + k;
    if (col >= 0 && col < n) {
      indices[nnzIdx] = col;
      data[nnzIdx] = 1.0;
      nnzIdx++;
    }
  }
  indptr[m] = nnzIdx;

  const result = createCSR(
    { data, indices, indptr },
    [m, n]
  );

  if (format === 'csc') {
    return result.tocsc() as CompressedSparseMatrix;
  }
  return result;
}

/**
 * Construct a sparse matrix from diagonals.
 *
 * @param diagonals - Array of diagonal values. Each element can be:
 *   - A number (scalar broadcast to fill the diagonal)
 *   - A number[] or Float64Array (explicit diagonal values)
 * @param offsets - Diagonal offsets (default [0]). Can be a single number or array.
 * @param shape - Matrix shape [m, n]. If not given, the matrix is square with
 *   size = len(diagonals[0]) + abs(offsets[0])
 * @param format - Output format: 'csr' or 'csc' (default 'csr')
 */
export function diags(
  diagonals: (number | number[] | Float64Array)[],
  offsets?: number | number[],
  shape?: [number, number],
  format: 'csr' | 'csc' = 'csr'
): CompressedSparseMatrix {
  // Normalize offsets
  let offsetArr: number[];
  if (offsets === undefined) {
    offsetArr = [0];
  } else if (typeof offsets === 'number') {
    offsetArr = [offsets];
  } else {
    offsetArr = offsets;
  }

  if (diagonals.length !== offsetArr.length) {
    throw new Error(
      `Number of diagonals (${diagonals.length}) does not match number of offsets (${offsetArr.length})`
    );
  }

  // Determine shape
  let m: number, n: number;
  if (shape) {
    [m, n] = shape;
  } else {
    // Infer from first diagonal
    const d0 = diagonals[0];
    const len0 = typeof d0 === 'number' ? 1 : d0.length;
    const absK = Math.abs(offsetArr[0]);
    m = len0 + absK;
    n = m;
  }

  // Collect COO entries
  const rows: number[] = [];
  const cols: number[] = [];
  const vals: number[] = [];

  for (let d = 0; d < diagonals.length; d++) {
    const k = offsetArr[d];
    const firstRow = k >= 0 ? 0 : -k;
    const firstCol = k >= 0 ? k : 0;
    const diagLen = Math.max(0, Math.min(m - firstRow, n - firstCol));

    const diag = diagonals[d];

    for (let i = 0; i < diagLen; i++) {
      let val: number;
      if (typeof diag === 'number') {
        val = diag;
      } else {
        val = diag[i];
      }
      if (val !== 0) {
        rows.push(firstRow + i);
        cols.push(firstCol + i);
        vals.push(val);
      }
    }
  }

  // Build CSR directly
  const nnz = vals.length;
  const indptr = new Int32Array(m + 1);
  const indices = new Int32Array(nnz);
  const data = new Float64Array(nnz);

  // Count per row
  for (let i = 0; i < nnz; i++) {
    indptr[rows[i] + 1]++;
  }
  for (let i = 1; i <= m; i++) {
    indptr[i] += indptr[i - 1];
  }

  // Fill (using temp offsets)
  const offset = new Int32Array(m);
  for (let i = 0; i < nnz; i++) {
    const row = rows[i];
    const dest = indptr[row] + offset[row];
    indices[dest] = cols[i];
    data[dest] = vals[i];
    offset[row]++;
  }

  const result = createCSR(
    { data, indices, indptr },
    [m, n]
  );

  if (format === 'csc') {
    return result.tocsc() as CompressedSparseMatrix;
  }
  return result;
}

/**
 * Check if x is a sparse matrix.
 *
 * @param x - Object to check
 * @returns True if x is a SparseMatrix instance
 */
export function issparse(x: unknown): x is SparseMatrix {
  return x instanceof SparseMatrix;
}

/**
 * Return the lower triangular portion of a sparse matrix.
 *
 * @param A - Sparse matrix or dense 2D array
 * @param k - Diagonal offset (default 0). k=0 is main diagonal, k>0 above, k<0 below.
 * @param format - Output format: 'coo', 'csr', or 'csc' (default 'coo')
 * @returns Sparse matrix with elements on and below the k-th diagonal
 */
export function tril(
  A: SparseMatrix | number[][],
  k: number = 0,
  format: SparseFormat = 'coo'
): SparseMatrix {
  const wasm = getWasmModule();
  const coo = toCOO(A);
  const { row, col, data, shape } = coo;
  const nnz = row.length;

  if (nnz === 0) {
    return cooToFormat(row, col, data, shape, format);
  }

  const rowPtr = toWasm(wasm, row);
  const colPtr = toWasm(wasm, col);
  const dataPtr = toWasm(wasm, data);
  const outRowPtr = wasm._malloc(nnz * 4);
  const outColPtr = wasm._malloc(nnz * 4);
  const outDataPtr = wasm._malloc(nnz * 8);

  try {
    const newNnz = wasm._sp_coo_tril_f64_entry(
      nnz, k,
      rowPtr, colPtr, dataPtr,
      outRowPtr, outColPtr, outDataPtr
    );

    const outRow = fromWasmI32(wasm, outRowPtr, newNnz);
    const outCol = fromWasmI32(wasm, outColPtr, newNnz);
    const outData = fromWasmF64(wasm, outDataPtr, newNnz);

    return cooToFormat(outRow, outCol, outData, shape, format);
  } finally {
    wasm._free(rowPtr);
    wasm._free(colPtr);
    wasm._free(dataPtr);
    wasm._free(outRowPtr);
    wasm._free(outColPtr);
    wasm._free(outDataPtr);
  }
}

/**
 * Return the upper triangular portion of a sparse matrix.
 *
 * @param A - Sparse matrix or dense 2D array
 * @param k - Diagonal offset (default 0). k=0 is main diagonal, k>0 above, k<0 below.
 * @param format - Output format: 'coo', 'csr', or 'csc' (default 'coo')
 * @returns Sparse matrix with elements on and above the k-th diagonal
 */
export function triu(
  A: SparseMatrix | number[][],
  k: number = 0,
  format: SparseFormat = 'coo'
): SparseMatrix {
  const wasm = getWasmModule();
  const coo = toCOO(A);
  const { row, col, data, shape } = coo;
  const nnz = row.length;

  if (nnz === 0) {
    return cooToFormat(row, col, data, shape, format);
  }

  const rowPtr = toWasm(wasm, row);
  const colPtr = toWasm(wasm, col);
  const dataPtr = toWasm(wasm, data);
  const outRowPtr = wasm._malloc(nnz * 4);
  const outColPtr = wasm._malloc(nnz * 4);
  const outDataPtr = wasm._malloc(nnz * 8);

  try {
    const newNnz = wasm._sp_coo_triu_f64_entry(
      nnz, k,
      rowPtr, colPtr, dataPtr,
      outRowPtr, outColPtr, outDataPtr
    );

    const outRow = fromWasmI32(wasm, outRowPtr, newNnz);
    const outCol = fromWasmI32(wasm, outColPtr, newNnz);
    const outData = fromWasmF64(wasm, outDataPtr, newNnz);

    return cooToFormat(outRow, outCol, outData, shape, format);
  } finally {
    wasm._free(rowPtr);
    wasm._free(colPtr);
    wasm._free(dataPtr);
    wasm._free(outRowPtr);
    wasm._free(outColPtr);
    wasm._free(outDataPtr);
  }
}

/**
 * Stack sparse matrices horizontally (column-wise).
 *
 * @param blocks - Array of sparse matrices to stack
 * @param format - Output format: 'coo', 'csr', or 'csc' (default 'csr')
 * @returns Horizontally stacked sparse matrix
 */
export function hstack(
  blocks: (SparseMatrix | number[][])[],
  format: SparseFormat = 'csr'
): SparseMatrix {
  if (blocks.length === 0) {
    throw new Error('Need at least one block to stack');
  }

  if (blocks.length === 1) {
    const coo = toCOO(blocks[0]);
    return cooToFormat(coo.row, coo.col, coo.data, coo.shape, format);
  }

  // Convert all to COO and validate row counts
  const coos = blocks.map(b => toCOO(b));
  const nRow = coos[0].shape[0];

  for (let i = 1; i < coos.length; i++) {
    if (coos[i].shape[0] !== nRow) {
      throw new Error(
        `Incompatible row dimensions: block 0 has ${nRow} rows, block ${i} has ${coos[i].shape[0]} rows`
      );
    }
  }

  const wasm = getWasmModule();
  const nBlocks = coos.length;

  // Prepare arrays
  const nColArr = new Int32Array(nBlocks);
  const nnzPerBlock = new Int32Array(nBlocks);
  let totalNnz = 0;
  let totalCols = 0;

  for (let i = 0; i < nBlocks; i++) {
    nColArr[i] = coos[i].shape[1];
    nnzPerBlock[i] = coos[i].row.length;
    totalNnz += nnzPerBlock[i];
    totalCols += nColArr[i];
  }

  // Concatenate arrays
  const rowCat = new Int32Array(totalNnz);
  const colCat = new Int32Array(totalNnz);
  const dataCat = new Float64Array(totalNnz);
  let offset = 0;

  for (let i = 0; i < nBlocks; i++) {
    rowCat.set(coos[i].row, offset);
    colCat.set(coos[i].col, offset);
    dataCat.set(coos[i].data, offset);
    offset += nnzPerBlock[i];
  }

  // Allocate WASM memory
  const nColPtr = toWasm(wasm, nColArr);
  const nnzPtr = toWasm(wasm, nnzPerBlock);
  const rowCatPtr = toWasm(wasm, rowCat);
  const colCatPtr = toWasm(wasm, colCat);
  const dataCatPtr = toWasm(wasm, dataCat);
  const outRowPtr = wasm._malloc(totalNnz * 4);
  const outColPtr = wasm._malloc(totalNnz * 4);
  const outDataPtr = wasm._malloc(totalNnz * 8);

  try {
    wasm._sp_coo_hstack_f64_entry(
      nBlocks, nColPtr, nnzPtr,
      rowCatPtr, colCatPtr, dataCatPtr,
      outRowPtr, outColPtr, outDataPtr
    );

    const outRow = fromWasmI32(wasm, outRowPtr, totalNnz);
    const outCol = fromWasmI32(wasm, outColPtr, totalNnz);
    const outData = fromWasmF64(wasm, outDataPtr, totalNnz);

    return cooToFormat(outRow, outCol, outData, [nRow, totalCols], format);
  } finally {
    wasm._free(nColPtr);
    wasm._free(nnzPtr);
    wasm._free(rowCatPtr);
    wasm._free(colCatPtr);
    wasm._free(dataCatPtr);
    wasm._free(outRowPtr);
    wasm._free(outColPtr);
    wasm._free(outDataPtr);
  }
}

/**
 * Stack sparse matrices vertically (row-wise).
 *
 * @param blocks - Array of sparse matrices to stack
 * @param format - Output format: 'coo', 'csr', or 'csc' (default 'csr')
 * @returns Vertically stacked sparse matrix
 */
export function vstack(
  blocks: (SparseMatrix | number[][])[],
  format: SparseFormat = 'csr'
): SparseMatrix {
  if (blocks.length === 0) {
    throw new Error('Need at least one block to stack');
  }

  if (blocks.length === 1) {
    const coo = toCOO(blocks[0]);
    return cooToFormat(coo.row, coo.col, coo.data, coo.shape, format);
  }

  // Convert all to COO and validate column counts
  const coos = blocks.map(b => toCOO(b));
  const nCol = coos[0].shape[1];

  for (let i = 1; i < coos.length; i++) {
    if (coos[i].shape[1] !== nCol) {
      throw new Error(
        `Incompatible column dimensions: block 0 has ${nCol} columns, block ${i} has ${coos[i].shape[1]} columns`
      );
    }
  }

  const wasm = getWasmModule();
  const nBlocks = coos.length;

  // Prepare arrays
  const nRowArr = new Int32Array(nBlocks);
  const nnzPerBlock = new Int32Array(nBlocks);
  let totalNnz = 0;
  let totalRows = 0;

  for (let i = 0; i < nBlocks; i++) {
    nRowArr[i] = coos[i].shape[0];
    nnzPerBlock[i] = coos[i].row.length;
    totalNnz += nnzPerBlock[i];
    totalRows += nRowArr[i];
  }

  // Concatenate arrays
  const rowCat = new Int32Array(totalNnz);
  const colCat = new Int32Array(totalNnz);
  const dataCat = new Float64Array(totalNnz);
  let offset = 0;

  for (let i = 0; i < nBlocks; i++) {
    rowCat.set(coos[i].row, offset);
    colCat.set(coos[i].col, offset);
    dataCat.set(coos[i].data, offset);
    offset += nnzPerBlock[i];
  }

  // Allocate WASM memory
  const nRowPtr = toWasm(wasm, nRowArr);
  const nnzPtr = toWasm(wasm, nnzPerBlock);
  const rowCatPtr = toWasm(wasm, rowCat);
  const colCatPtr = toWasm(wasm, colCat);
  const dataCatPtr = toWasm(wasm, dataCat);
  const outRowPtr = wasm._malloc(totalNnz * 4);
  const outColPtr = wasm._malloc(totalNnz * 4);
  const outDataPtr = wasm._malloc(totalNnz * 8);

  try {
    wasm._sp_coo_vstack_f64_entry(
      nBlocks, nRowPtr, nnzPtr,
      rowCatPtr, colCatPtr, dataCatPtr,
      outRowPtr, outColPtr, outDataPtr
    );

    const outRow = fromWasmI32(wasm, outRowPtr, totalNnz);
    const outCol = fromWasmI32(wasm, outColPtr, totalNnz);
    const outData = fromWasmF64(wasm, outDataPtr, totalNnz);

    return cooToFormat(outRow, outCol, outData, [totalRows, nCol], format);
  } finally {
    wasm._free(nRowPtr);
    wasm._free(nnzPtr);
    wasm._free(rowCatPtr);
    wasm._free(colCatPtr);
    wasm._free(dataCatPtr);
    wasm._free(outRowPtr);
    wasm._free(outColPtr);
    wasm._free(outDataPtr);
  }
}

/**
 * Build a block diagonal sparse matrix from provided matrices.
 *
 * @param mats - Array of sparse matrices or dense 2D arrays
 * @param format - Output format: 'coo', 'csr', or 'csc' (default 'csr')
 * @returns Block diagonal sparse matrix
 */
export function block_diag(
  mats: (SparseMatrix | number[][])[],
  format: SparseFormat = 'csr'
): SparseMatrix {
  if (mats.length === 0) {
    return cooToFormat(new Int32Array(0), new Int32Array(0), new Float64Array(0), [0, 0], format);
  }

  if (mats.length === 1) {
    const coo = toCOO(mats[0]);
    return cooToFormat(coo.row, coo.col, coo.data, coo.shape, format);
  }

  const wasm = getWasmModule();
  const coos = mats.map(m => toCOO(m));
  const nBlocks = coos.length;

  // Prepare arrays
  const nRowArr = new Int32Array(nBlocks);
  const nColArr = new Int32Array(nBlocks);
  const nnzPerBlock = new Int32Array(nBlocks);
  let totalNnz = 0;
  let totalRows = 0;
  let totalCols = 0;

  for (let i = 0; i < nBlocks; i++) {
    nRowArr[i] = coos[i].shape[0];
    nColArr[i] = coos[i].shape[1];
    nnzPerBlock[i] = coos[i].row.length;
    totalNnz += nnzPerBlock[i];
    totalRows += nRowArr[i];
    totalCols += nColArr[i];
  }

  // Concatenate arrays
  const rowCat = new Int32Array(totalNnz);
  const colCat = new Int32Array(totalNnz);
  const dataCat = new Float64Array(totalNnz);
  let offset = 0;

  for (let i = 0; i < nBlocks; i++) {
    rowCat.set(coos[i].row, offset);
    colCat.set(coos[i].col, offset);
    dataCat.set(coos[i].data, offset);
    offset += nnzPerBlock[i];
  }

  // Allocate WASM memory
  const nRowPtr = toWasm(wasm, nRowArr);
  const nColPtr = toWasm(wasm, nColArr);
  const nnzPtr = toWasm(wasm, nnzPerBlock);
  const rowCatPtr = toWasm(wasm, rowCat);
  const colCatPtr = toWasm(wasm, colCat);
  const dataCatPtr = toWasm(wasm, dataCat);
  const outRowPtr = wasm._malloc(totalNnz * 4);
  const outColPtr = wasm._malloc(totalNnz * 4);
  const outDataPtr = wasm._malloc(totalNnz * 8);

  try {
    wasm._sp_coo_block_diag_f64_entry(
      nBlocks, nRowPtr, nColPtr, nnzPtr,
      rowCatPtr, colCatPtr, dataCatPtr,
      outRowPtr, outColPtr, outDataPtr
    );

    const outRow = fromWasmI32(wasm, outRowPtr, totalNnz);
    const outCol = fromWasmI32(wasm, outColPtr, totalNnz);
    const outData = fromWasmF64(wasm, outDataPtr, totalNnz);

    return cooToFormat(outRow, outCol, outData, [totalRows, totalCols], format);
  } finally {
    wasm._free(nRowPtr);
    wasm._free(nColPtr);
    wasm._free(nnzPtr);
    wasm._free(rowCatPtr);
    wasm._free(colCatPtr);
    wasm._free(dataCatPtr);
    wasm._free(outRowPtr);
    wasm._free(outColPtr);
    wasm._free(outDataPtr);
  }
}

/**
 * Kronecker product of sparse matrices.
 *
 * @param A - First sparse matrix or dense 2D array
 * @param B - Second sparse matrix or dense 2D array
 * @param format - Output format: 'coo', 'csr', or 'csc' (default 'csr')
 * @returns Kronecker product A ⊗ B
 */
export function kron(
  A: SparseMatrix | number[][],
  B: SparseMatrix | number[][],
  format: SparseFormat = 'csr'
): SparseMatrix {
  const wasm = getWasmModule();
  const cooA = toCOO(A);
  const cooB = toCOO(B);

  const nnzA = cooA.row.length;
  const nnzB = cooB.row.length;
  const outNnz = nnzA * nnzB;

  if (outNnz === 0) {
    const outShape: [number, number] = [
      cooA.shape[0] * cooB.shape[0],
      cooA.shape[1] * cooB.shape[1]
    ];
    return cooToFormat(new Int32Array(0), new Int32Array(0), new Float64Array(0), outShape, format);
  }

  // Allocate WASM memory
  const aRowPtr = toWasm(wasm, cooA.row);
  const aColPtr = toWasm(wasm, cooA.col);
  const aDataPtr = toWasm(wasm, cooA.data);
  const bRowPtr = toWasm(wasm, cooB.row);
  const bColPtr = toWasm(wasm, cooB.col);
  const bDataPtr = toWasm(wasm, cooB.data);
  const outRowPtr = wasm._malloc(outNnz * 4);
  const outColPtr = wasm._malloc(outNnz * 4);
  const outDataPtr = wasm._malloc(outNnz * 8);

  try {
    wasm._sp_coo_kron_f64_entry(
      nnzA, nnzB,
      cooB.shape[0], cooB.shape[1],
      aRowPtr, aColPtr, aDataPtr,
      bRowPtr, bColPtr, bDataPtr,
      outRowPtr, outColPtr, outDataPtr
    );

    const outRow = fromWasmI32(wasm, outRowPtr, outNnz);
    const outCol = fromWasmI32(wasm, outColPtr, outNnz);
    const outData = fromWasmF64(wasm, outDataPtr, outNnz);

    const outShape: [number, number] = [
      cooA.shape[0] * cooB.shape[0],
      cooA.shape[1] * cooB.shape[1]
    ];

    return cooToFormat(outRow, outCol, outData, outShape, format);
  } finally {
    wasm._free(aRowPtr);
    wasm._free(aColPtr);
    wasm._free(aDataPtr);
    wasm._free(bRowPtr);
    wasm._free(bColPtr);
    wasm._free(bDataPtr);
    wasm._free(outRowPtr);
    wasm._free(outColPtr);
    wasm._free(outDataPtr);
  }
}

/**
 * Kronecker sum of sparse matrices.
 *
 * kronsum(A, B) = kron(A, I_b) + kron(I_a, B)
 *
 * where I_a and I_b are identity matrices of shape A and B respectively.
 * Both A and B must be square matrices.
 *
 * @param A - First square sparse matrix or dense 2D array
 * @param B - Second square sparse matrix or dense 2D array
 * @param format - Output format: 'coo', 'csr', or 'csc' (default 'csr')
 * @returns Kronecker sum
 */
export function kronsum(
  A: SparseMatrix | number[][],
  B: SparseMatrix | number[][],
  format: SparseFormat = 'csr'
): SparseMatrix {
  const cooA = toCOO(A);
  const cooB = toCOO(B);

  if (cooA.shape[0] !== cooA.shape[1]) {
    throw new Error(`A must be square, got shape [${cooA.shape[0]}, ${cooA.shape[1]}]`);
  }
  if (cooB.shape[0] !== cooB.shape[1]) {
    throw new Error(`B must be square, got shape [${cooB.shape[0]}, ${cooB.shape[1]}]`);
  }

  const na = cooA.shape[0];
  const nb = cooB.shape[0];

  // kron(A, I_b) + kron(I_a, B)
  const Ib = eye(nb, nb);
  const Ia = eye(na, na);

  const term1 = kron(A, Ib, 'csr') as CompressedSparseMatrix;
  const term2 = kron(Ia, B, 'csr') as CompressedSparseMatrix;

  const result = term1.add(term2);

  if (format === 'coo') {
    // Convert CSR result to COO using our helper, then wrap in cooToFormat
    const coo = toCOO(result);
    return cooToFormat(coo.row, coo.col, coo.data, coo.shape, 'coo');
  } else if (format === 'csc') {
    return result.tocsc();
  }
  return result;
}

/**
 * Generate a sparse matrix with uniformly distributed random values.
 *
 * @param m - Number of rows
 * @param n - Number of columns
 * @param density - Density of non-zero values (0 to 1, default 0.01)
 * @param format - Output format: 'coo', 'csr', or 'csc' (default 'csr')
 * @param rng - Optional random number generator function (default Math.random)
 * @returns Random sparse matrix
 */
export function random(
  m: number,
  n: number,
  density: number = 0.01,
  format: SparseFormat = 'csr',
  rng: () => number = Math.random
): SparseMatrix {
  if (density < 0 || density > 1) {
    throw new Error(`density must be between 0 and 1, got ${density}`);
  }

  const totalSize = m * n;
  const targetNnz = Math.round(totalSize * density);

  if (targetNnz === 0) {
    return cooToFormat(new Int32Array(0), new Int32Array(0), new Float64Array(0), [m, n], format);
  }

  // Generate unique random flat indices using Fisher-Yates shuffle variant
  const nnz = Math.min(targetNnz, totalSize);

  // For small density, use reservoir sampling
  // For large density, use shuffle
  let flatIndices: Int32Array;

  if (density <= 0.5) {
    // Use set-based sampling for low density
    const indexSet = new Set<number>();
    while (indexSet.size < nnz) {
      const idx = Math.floor(rng() * totalSize);
      indexSet.add(idx);
    }
    flatIndices = new Int32Array([...indexSet]);
  } else {
    // Use shuffle for high density
    const allIndices = new Int32Array(totalSize);
    for (let i = 0; i < totalSize; i++) {
      allIndices[i] = i;
    }
    // Partial Fisher-Yates
    for (let i = 0; i < nnz; i++) {
      const j = i + Math.floor(rng() * (totalSize - i));
      const temp = allIndices[i];
      allIndices[i] = allIndices[j];
      allIndices[j] = temp;
    }
    flatIndices = allIndices.slice(0, nnz);
  }

  // Generate random values
  const randomValues = new Float64Array(nnz);
  for (let i = 0; i < nnz; i++) {
    randomValues[i] = rng();
  }

  const wasm = getWasmModule();

  const flatPtr = toWasm(wasm, flatIndices);
  const valPtr = toWasm(wasm, randomValues);
  const outRowPtr = wasm._malloc(nnz * 4);
  const outColPtr = wasm._malloc(nnz * 4);
  const outDataPtr = wasm._malloc(nnz * 8);

  try {
    wasm._sp_coo_random_f64_entry(
      nnz, m, n,
      flatPtr, valPtr,
      outRowPtr, outColPtr, outDataPtr
    );

    const outRow = fromWasmI32(wasm, outRowPtr, nnz);
    const outCol = fromWasmI32(wasm, outColPtr, nnz);
    const outData = fromWasmF64(wasm, outDataPtr, nnz);

    return cooToFormat(outRow, outCol, outData, [m, n], format);
  } finally {
    wasm._free(flatPtr);
    wasm._free(valPtr);
    wasm._free(outRowPtr);
    wasm._free(outColPtr);
    wasm._free(outDataPtr);
  }
}

