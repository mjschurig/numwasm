/**
 * Create sparse matrix from COO format
 * @param rows - Number of rows
 * @param cols - Number of columns
 * @param data - COO data {row: number[], col: number[], val: number[]}
 * @returns SparseMatrix instance
 */
export function fromCOO(
  rows: number,
  cols: number,
  data: { row: number[]; col: number[]; val: number[] }
): unknown {
  throw new Error('Not implemented');
}
