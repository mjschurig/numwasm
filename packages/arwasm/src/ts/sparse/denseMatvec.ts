/**
 * Dense Matrix-Vector Product
 *
 * Creates a matvec function from a dense matrix stored as a flat array.
 * Useful for small problems or as a reference implementation.
 */

import type { MatVecFunction } from '../high-level-types.js';

/**
 * Create a matrix-vector product function from a dense matrix.
 *
 * The matrix is stored as a flat Float64Array in either row-major
 * (C-style) or column-major (Fortran-style) order.
 *
 * @param matrix - Flat array containing matrix elements (length m*n)
 * @param m - Number of rows
 * @param n - Number of columns
 * @param rowMajor - If true, matrix is row-major (default). If false, column-major.
 * @returns Function computing y = A*x
 *
 * @example
 * ```ts
 * // Matrix: [[1, 2, 3],
 * //          [4, 5, 6]]
 * // Row-major storage: [1, 2, 3, 4, 5, 6]
 * const matrix = new Float64Array([1, 2, 3, 4, 5, 6]);
 * const matvec = denseMatvec(matrix, 2, 3, true);
 *
 * const x = new Float64Array([1, 1, 1]);
 * const y = matvec(x); // [6, 15]
 * ```
 */
export function denseMatvec(
  matrix: Float64Array,
  m: number,
  n: number,
  rowMajor: boolean = true
): MatVecFunction {
  // Validate inputs
  if (matrix.length !== m * n) {
    throw new Error(`Matrix length must be m*n (${m * n}), got ${matrix.length}`);
  }

  if (rowMajor) {
    // Row-major: A[i,j] = matrix[i*n + j]
    return (x: Float64Array): Float64Array => {
      if (x.length !== n) {
        throw new Error(`Input vector length must be ${n}, got ${x.length}`);
      }

      const y = new Float64Array(m);

      for (let i = 0; i < m; i++) {
        let sum = 0;
        const rowOffset = i * n;
        for (let j = 0; j < n; j++) {
          sum += matrix[rowOffset + j] * x[j];
        }
        y[i] = sum;
      }

      return y;
    };
  } else {
    // Column-major: A[i,j] = matrix[i + j*m]
    return (x: Float64Array): Float64Array => {
      if (x.length !== n) {
        throw new Error(`Input vector length must be ${n}, got ${x.length}`);
      }

      const y = new Float64Array(m);

      for (let j = 0; j < n; j++) {
        const xj = x[j];
        if (xj === 0) continue;

        const colOffset = j * m;
        for (let i = 0; i < m; i++) {
          y[i] += matrix[colOffset + i] * xj;
        }
      }

      return y;
    };
  }
}

/**
 * Create a transpose matrix-vector product function from a dense matrix.
 *
 * @param matrix - Flat array containing matrix elements
 * @param m - Number of rows
 * @param n - Number of columns
 * @param rowMajor - If true, matrix is row-major. If false, column-major.
 * @returns Function computing y = A^T * x
 */
export function denseMatvecT(
  matrix: Float64Array,
  m: number,
  n: number,
  rowMajor: boolean = true
): MatVecFunction {
  if (rowMajor) {
    // Row-major: A[i,j] = matrix[i*n + j], so A^T[j,i] = matrix[i*n + j]
    return (x: Float64Array): Float64Array => {
      if (x.length !== m) {
        throw new Error(`Input vector length must be ${m}, got ${x.length}`);
      }

      const y = new Float64Array(n);

      for (let i = 0; i < m; i++) {
        const xi = x[i];
        if (xi === 0) continue;

        const rowOffset = i * n;
        for (let j = 0; j < n; j++) {
          y[j] += matrix[rowOffset + j] * xi;
        }
      }

      return y;
    };
  } else {
    // Column-major: A[i,j] = matrix[i + j*m], so A^T[j,i] = matrix[i + j*m]
    return (x: Float64Array): Float64Array => {
      if (x.length !== m) {
        throw new Error(`Input vector length must be ${m}, got ${x.length}`);
      }

      const y = new Float64Array(n);

      for (let j = 0; j < n; j++) {
        let sum = 0;
        const colOffset = j * m;
        for (let i = 0; i < m; i++) {
          sum += matrix[colOffset + i] * x[i];
        }
        y[j] = sum;
      }

      return y;
    };
  }
}
