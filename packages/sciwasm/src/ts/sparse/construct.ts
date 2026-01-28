import type { CompressedSparseMatrix } from './compressed.js';
import { createCSR } from './_factory.js';

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

