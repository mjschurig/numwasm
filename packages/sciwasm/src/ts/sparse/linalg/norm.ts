/**
 * Sparse matrix norms
 *
 * Based on scipy.sparse.linalg.norm
 */

import type { SparseMatrix } from '../base.js';
import type { CompressedSparseMatrix } from '../compressed.js';

/**
 * Sparse matrix or vector norm
 *
 * @param x - Sparse matrix
 * @param ord - Order of the norm:
 *   - 'fro' or undefined: Frobenius norm (default)
 *   - 1: max column sum
 *   - -1: min column sum
 *   - Infinity: max row sum
 *   - -Infinity: min row sum
 *   - 2: largest singular value (not yet implemented)
 *   - -2: smallest singular value (not yet implemented)
 * @returns The norm value
 *
 * @example
 * ```typescript
 * const A = csr_matrix([[1, 2], [3, 4]]);
 * norm(A);        // Frobenius norm: sqrt(1 + 4 + 9 + 16) = sqrt(30)
 * norm(A, 1);     // max column sum: max(|1|+|3|, |2|+|4|) = 6
 * norm(A, Infinity); // max row sum: max(|1|+|2|, |3|+|4|) = 7
 * ```
 */
export function norm(
  x: SparseMatrix,
  ord?: 'fro' | 1 | -1 | 2 | -2 | number
): number {

  // Default to Frobenius norm
  if (ord === undefined || ord === 'fro') {
    return frobeniusNorm(x);
  }

  if (ord === 1) {
    // Max column sum of absolute values
    return maxColumnSum(x);
  }

  if (ord === -1) {
    // Min column sum of absolute values
    return minColumnSum(x);
  }

  if (ord === Infinity) {
    // Max row sum of absolute values
    return maxRowSum(x);
  }

  if (ord === -Infinity) {
    // Min row sum of absolute values
    return minRowSum(x);
  }

  if (ord === 2 || ord === -2) {
    // Requires svds - defer to Phase 5
    throw new Error(`norm with ord=${ord} requires svds, which is not yet implemented`);
  }

  throw new Error(`Invalid norm order: ${ord}`);
}

/**
 * Frobenius norm: sqrt(sum(|a_ij|^2))
 */
function frobeniusNorm(x: SparseMatrix): number {
  // Get data array from any sparse format
  const csr = x.tocsr() as CompressedSparseMatrix;
  const data = csr.data;

  let sumSq = 0;
  for (let i = 0; i < data.length; i++) {
    sumSq += data[i] * data[i];
  }
  return Math.sqrt(sumSq);
}

/**
 * Max column sum: max_j(sum_i |a_ij|)
 * Convert to CSC and sum each column
 */
function maxColumnSum(x: SparseMatrix): number {
  const csc = x.tocsc() as CompressedSparseMatrix;
  const [, ncol] = x.shape;
  const data = csc.data;
  const indptr = csc.indptr;

  if (ncol === 0) return 0;

  let maxSum = 0;
  for (let j = 0; j < ncol; j++) {
    let colSum = 0;
    for (let k = indptr[j]; k < indptr[j + 1]; k++) {
      colSum += Math.abs(data[k]);
    }
    if (colSum > maxSum) {
      maxSum = colSum;
    }
  }
  return maxSum;
}

/**
 * Min column sum: min_j(sum_i |a_ij|)
 */
function minColumnSum(x: SparseMatrix): number {
  const csc = x.tocsc() as CompressedSparseMatrix;
  const [, ncol] = x.shape;
  const data = csc.data;
  const indptr = csc.indptr;

  if (ncol === 0) return 0;

  let minSum = Infinity;
  for (let j = 0; j < ncol; j++) {
    let colSum = 0;
    for (let k = indptr[j]; k < indptr[j + 1]; k++) {
      colSum += Math.abs(data[k]);
    }
    if (colSum < minSum) {
      minSum = colSum;
    }
  }
  return minSum === Infinity ? 0 : minSum;
}

/**
 * Max row sum: max_i(sum_j |a_ij|)
 */
function maxRowSum(x: SparseMatrix): number {
  const csr = x.tocsr() as CompressedSparseMatrix;
  const [nrow] = x.shape;
  const data = csr.data;
  const indptr = csr.indptr;

  if (nrow === 0) return 0;

  let maxSum = 0;
  for (let i = 0; i < nrow; i++) {
    let rowSum = 0;
    for (let k = indptr[i]; k < indptr[i + 1]; k++) {
      rowSum += Math.abs(data[k]);
    }
    if (rowSum > maxSum) {
      maxSum = rowSum;
    }
  }
  return maxSum;
}

/**
 * Min row sum: min_i(sum_j |a_ij|)
 */
function minRowSum(x: SparseMatrix): number {
  const csr = x.tocsr() as CompressedSparseMatrix;
  const [nrow] = x.shape;
  const data = csr.data;
  const indptr = csr.indptr;

  if (nrow === 0) return 0;

  let minSum = Infinity;
  for (let i = 0; i < nrow; i++) {
    let rowSum = 0;
    for (let k = indptr[i]; k < indptr[i + 1]; k++) {
      rowSum += Math.abs(data[k]);
    }
    if (rowSum < minSum) {
      minSum = rowSum;
    }
  }
  return minSum === Infinity ? 0 : minSum;
}
