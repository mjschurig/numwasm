/**
 * CSC (Compressed Sparse Column) Matrix-Vector Product
 *
 * Creates a matvec function from CSC sparse matrix format.
 * CSC is efficient for column slicing and A^T*x products.
 */

import type { MatVecFunction } from '../high-level-types.js';

/**
 * Create a matrix-vector product function from CSC sparse matrix format.
 *
 * CSC format stores a sparse matrix using three arrays:
 * - indptr: Column pointers (length n+1). indptr[j] to indptr[j+1] gives the
 *   range of indices in data/indices for column j.
 * - indices: Row indices for each non-zero entry
 * - data: Non-zero values
 *
 * @param indptr - Column pointer array (length n+1)
 * @param indices - Row indices array
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
 * const data = new Float64Array([1, 4, 3, 2, 5]);
 *
 * const matvec = cscMatvec(indptr, indices, data, [3, 3]);
 * const x = new Float64Array([1, 1, 1]);
 * const y = matvec(x); // [3, 3, 9]
 * ```
 */
export function cscMatvec(
  indptr: Int32Array,
  indices: Int32Array,
  data: Float64Array,
  shape: [number, number]
): MatVecFunction {
  const [m, n] = shape;

  // Validate inputs
  if (indptr.length !== n + 1) {
    throw new Error(`indptr length must be n+1 (${n + 1}), got ${indptr.length}`);
  }
  if (indices.length !== data.length) {
    throw new Error(`indices and data must have same length`);
  }
  if (indptr[n] !== data.length) {
    throw new Error(`indptr[n] must equal nnz (${data.length}), got ${indptr[n]}`);
  }

  return (x: Float64Array): Float64Array => {
    if (x.length !== n) {
      throw new Error(`Input vector length must be ${n}, got ${x.length}`);
    }

    const y = new Float64Array(m);

    for (let j = 0; j < n; j++) {
      const xj = x[j];
      if (xj === 0) continue; // Skip zero entries for efficiency

      const colStart = indptr[j];
      const colEnd = indptr[j + 1];

      for (let k = colStart; k < colEnd; k++) {
        y[indices[k]] += data[k] * xj;
      }
    }

    return y;
  };
}

/**
 * Create a transpose matrix-vector product function from CSC format.
 *
 * Computes y = A^T * x efficiently using CSC format.
 * This is equivalent to CSR matvec on the original matrix.
 *
 * @param indptr - Column pointer array
 * @param indices - Row indices array
 * @param data - Non-zero values array
 * @param shape - Matrix dimensions [m, n]
 * @returns Function computing y = A^T * x
 */
export function cscMatvecT(
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

    for (let j = 0; j < n; j++) {
      let sum = 0;
      const colStart = indptr[j];
      const colEnd = indptr[j + 1];

      for (let k = colStart; k < colEnd; k++) {
        sum += data[k] * x[indices[k]];
      }

      y[j] = sum;
    }

    return y;
  };
}
