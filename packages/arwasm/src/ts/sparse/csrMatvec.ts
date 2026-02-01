/**
 * CSR (Compressed Sparse Row) Matrix-Vector Product
 *
 * Creates a matvec function from CSR sparse matrix format.
 * CSR is efficient for row slicing and matrix-vector products.
 */

import type { MatVecFunction } from '../high-level-types.js';

/**
 * Create a matrix-vector product function from CSR sparse matrix format.
 *
 * CSR format stores a sparse matrix using three arrays:
 * - indptr: Row pointers (length m+1). indptr[i] to indptr[i+1] gives the
 *   range of indices in data/indices for row i.
 * - indices: Column indices for each non-zero entry
 * - data: Non-zero values
 *
 * @param indptr - Row pointer array (length m+1)
 * @param indices - Column indices array
 * @param data - Non-zero values array
 * @param shape - Matrix dimensions [m, n] (rows, cols)
 * @returns Function computing y = A*x
 *
 * @example
 * ```ts
 * // Matrix: [[1, 0, 2],
 * //          [0, 3, 0],
 * //          [4, 0, 5]]
 * const indptr = new Int32Array([0, 2, 3, 5]);
 * const indices = new Int32Array([0, 2, 1, 0, 2]);
 * const data = new Float64Array([1, 2, 3, 4, 5]);
 *
 * const matvec = csrMatvec(indptr, indices, data, [3, 3]);
 * const x = new Float64Array([1, 1, 1]);
 * const y = matvec(x); // [3, 3, 9]
 * ```
 */
export function csrMatvec(
  indptr: Int32Array,
  indices: Int32Array,
  data: Float64Array,
  shape: [number, number]
): MatVecFunction {
  const [m, n] = shape;

  // Validate inputs
  if (indptr.length !== m + 1) {
    throw new Error(`indptr length must be m+1 (${m + 1}), got ${indptr.length}`);
  }
  if (indices.length !== data.length) {
    throw new Error(`indices and data must have same length`);
  }
  if (indptr[m] !== data.length) {
    throw new Error(`indptr[m] must equal nnz (${data.length}), got ${indptr[m]}`);
  }

  return (x: Float64Array): Float64Array => {
    if (x.length !== n) {
      throw new Error(`Input vector length must be ${n}, got ${x.length}`);
    }

    const y = new Float64Array(m);

    for (let i = 0; i < m; i++) {
      let sum = 0;
      const rowStart = indptr[i];
      const rowEnd = indptr[i + 1];

      for (let k = rowStart; k < rowEnd; k++) {
        sum += data[k] * x[indices[k]];
      }

      y[i] = sum;
    }

    return y;
  };
}

/**
 * Create a transpose matrix-vector product function from CSR format.
 *
 * Computes y = A^T * x efficiently using CSR format.
 *
 * @param indptr - Row pointer array
 * @param indices - Column indices array
 * @param data - Non-zero values array
 * @param shape - Matrix dimensions [m, n]
 * @returns Function computing y = A^T * x
 */
export function csrMatvecT(
  indptr: Int32Array,
  indices: Int32Array,
  data: Float64Array,
  shape: [number, number]
): MatVecFunction {
  const [m, n] = shape;

  return (x: Float64Array): Float64Array => {
    if (x.length !== m) {
      throw new Error(`Input vector length must be ${m}, got ${x.length}`);
    }

    const y = new Float64Array(n);

    for (let i = 0; i < m; i++) {
      const xi = x[i];
      const rowStart = indptr[i];
      const rowEnd = indptr[i + 1];

      for (let k = rowStart; k < rowEnd; k++) {
        y[indices[k]] += data[k] * xi;
      }
    }

    return y;
  };
}
