/**
 * COO (Coordinate) Matrix-Vector Product
 *
 * Creates a matvec function from COO sparse matrix format.
 * COO is simple to construct but less efficient for repeated matvec operations.
 */

import type { MatVecFunction } from '../high-level-types.js';

/**
 * Create a matrix-vector product function from COO sparse matrix format.
 *
 * COO format stores a sparse matrix using three arrays:
 * - rows: Row indices for each non-zero entry
 * - cols: Column indices for each non-zero entry
 * - values: Non-zero values
 *
 * @param rows - Row indices array
 * @param cols - Column indices array
 * @param values - Non-zero values array
 * @param shape - Matrix dimensions [m, n] (rows, cols)
 * @returns Function computing y = A*x
 *
 * @example
 * ```ts
 * // Matrix: [[1, 0, 2],
 * //          [0, 3, 0],
 * //          [4, 0, 5]]
 * const rows = new Int32Array([0, 0, 1, 2, 2]);
 * const cols = new Int32Array([0, 2, 1, 0, 2]);
 * const values = new Float64Array([1, 2, 3, 4, 5]);
 *
 * const matvec = cooMatvec(rows, cols, values, [3, 3]);
 * const x = new Float64Array([1, 1, 1]);
 * const y = matvec(x); // [3, 3, 9]
 * ```
 */
export function cooMatvec(
  rows: Int32Array,
  cols: Int32Array,
  values: Float64Array,
  shape: [number, number]
): MatVecFunction {
  const [m, n] = shape;
  const nnz = values.length;

  // Validate inputs
  if (rows.length !== nnz || cols.length !== nnz) {
    throw new Error(`rows, cols, and values must have same length`);
  }

  return (x: Float64Array): Float64Array => {
    if (x.length !== n) {
      throw new Error(`Input vector length must be ${n}, got ${x.length}`);
    }

    const y = new Float64Array(m);

    for (let k = 0; k < nnz; k++) {
      y[rows[k]] += values[k] * x[cols[k]];
    }

    return y;
  };
}

/**
 * Create a transpose matrix-vector product function from COO format.
 *
 * Computes y = A^T * x using COO format.
 *
 * @param rows - Row indices array
 * @param cols - Column indices array
 * @param values - Non-zero values array
 * @param shape - Matrix dimensions [m, n]
 * @returns Function computing y = A^T * x
 */
export function cooMatvecT(
  rows: Int32Array,
  cols: Int32Array,
  values: Float64Array,
  shape: [number, number]
): MatVecFunction {
  const [m, n] = shape;
  const nnz = values.length;

  return (x: Float64Array): Float64Array => {
    if (x.length !== m) {
      throw new Error(`Input vector length must be ${m}, got ${x.length}`);
    }

    const y = new Float64Array(n);

    for (let k = 0; k < nnz; k++) {
      y[cols[k]] += values[k] * x[rows[k]];
    }

    return y;
  };
}
