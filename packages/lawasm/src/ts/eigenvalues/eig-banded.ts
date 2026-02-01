/**
 * Banded Matrix Eigenvalue Problem
 *
 * Compute eigenvalues and eigenvectors of a symmetric banded matrix.
 *
 * Note: DSBEV is not currently exported from LAPACK WASM.
 * This provides a pure TypeScript implementation by converting to dense form.
 */

import type { Matrix, EigBandedOptions, EigBandedResult } from './types.js';
import { eigSymmetric } from './eig-symmetric.js';

/**
 * Compute eigenvalues and eigenvectors of a symmetric banded n×n matrix.
 *
 * The matrix is stored in banded format where:
 * - A[i][j] represents element (i, j-i+kd) of the original matrix
 * - kd is the number of super-diagonals (or sub-diagonals)
 *
 * For a symmetric matrix, only the upper or lower triangle needs to be stored.
 *
 * @param A - Banded matrix in LAPACK band storage format
 * @param n - Matrix dimension
 * @param kd - Number of super-diagonals (bandwidth - 1)
 * @param options - Computation options
 * @returns Eigenvalues (in ascending order) and optionally eigenvectors
 *
 * @example
 * ```typescript
 * // Tridiagonal matrix:
 * // [2 1 0]
 * // [1 2 1]
 * // [0 1 2]
 * // In upper banded format (kd=1):
 * // [* 1 1]  <- super-diagonal
 * // [2 2 2]  <- diagonal
 * const AB = [[0, 1, 1], [2, 2, 2]]; // (kd+1) x n
 * const { values, vectors } = eigBanded(AB, 3, 1);
 * ```
 */
export function eigBanded(
  A: Matrix,
  n: number,
  kd: number,
  options: EigBandedOptions = {}
): EigBandedResult {
  const { computeVectors = true, uplo = 'lower' } = options;

  try {
    // Convert banded format to dense symmetric matrix
    const dense = new Float64Array(n * n);

    // Get data from input
    let bandData: number[];
    if (Array.isArray(A) && Array.isArray(A[0])) {
      // 2D array - flatten
      const rows = A as number[][];
      if (rows.length !== kd + 1) {
        throw new Error(`Expected ${kd + 1} rows in banded format, got ${rows.length}`);
      }
      bandData = [];
      for (let i = 0; i < rows.length; i++) {
        for (let j = 0; j < rows[i].length; j++) {
          bandData.push(rows[i][j]);
        }
      }
    } else {
      // 1D array
      const arr = A as number[] | Float64Array;
      bandData = Array.from(arr);
    }

    // Convert banded to dense format
    // LAPACK band storage: AB(kd+1-i+j, i) = A(i, j) for max(1,j-kd) <= i <= j (upper)
    // or AB(1+i-j, j) = A(i, j) for j <= i <= min(n, j+kd) (lower)

    if (uplo === 'upper') {
      // Upper triangular banded storage
      for (let j = 0; j < n; j++) {
        for (let i = Math.max(0, j - kd); i <= j; i++) {
          // AB(kd + i - j, j) = A(i, j)
          const bandIdx = (kd + i - j) * n + j;
          if (bandIdx < bandData.length) {
            const val = bandData[bandIdx];
            // Column-major: A[i,j] at index j*n + i
            dense[j * n + i] = val;
            dense[i * n + j] = val; // Symmetric
          }
        }
      }
    } else {
      // Lower triangular banded storage
      for (let j = 0; j < n; j++) {
        for (let i = j; i <= Math.min(n - 1, j + kd); i++) {
          // AB(i - j, j) = A(i, j)
          const bandIdx = (i - j) * n + j;
          if (bandIdx < bandData.length) {
            const val = bandData[bandIdx];
            dense[j * n + i] = val;
            dense[i * n + j] = val; // Symmetric
          }
        }
      }
    }

    // Use eigSymmetric for the dense matrix
    const result = eigSymmetric(dense, {
      computeVectors,
      uplo: 'lower', // Dense matrix is fully populated
    });

    return {
      values: result.values,
      vectors: result.vectors,
      n,
      info: result.info,
      success: result.success,
      message: result.message,
    };
  } catch (e) {
    return {
      values: new Float64Array(0),
      n,
      info: -1,
      success: false,
      message: `Error: ${e instanceof Error ? e.message : String(e)}`,
    };
  }
}
