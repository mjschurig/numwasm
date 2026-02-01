/**
 * Banded System Solver
 *
 * Solve Ax = b where A is a general banded matrix.
 */

import type { RealArray } from '../helpers.js';
import type { Matrix, SolveResult } from './types.js';
import { solve } from './solve.js';

/**
 * Solve a banded system Ax = b where A is a general banded matrix.
 *
 * A banded matrix has non-zero elements only within kl diagonals below
 * and ku diagonals above the main diagonal.
 *
 * The banded matrix is stored in band storage format:
 * - Row i of the band array contains diagonal (i - kl)
 * - The band array has (2*kl + ku + 1) rows and n columns
 *
 * Note: This function converts the banded format to dense and uses
 * the general solver. For large sparse banded systems, consider using
 * specialized routines.
 *
 * @param AB - Banded matrix in LAPACK band storage format
 * @param b - Right-hand side vector or matrix
 * @param n - Matrix dimension
 * @param kl - Number of sub-diagonals
 * @param ku - Number of super-diagonals
 * @returns Solution x and diagnostic information
 *
 * @example
 * ```typescript
 * // Tridiagonal matrix (kl=1, ku=1):
 * // [2  1  0  0]
 * // [1  2  1  0]
 * // [0  1  2  1]
 * // [0  0  1  2]
 * //
 * // Band storage (rows are diagonals, bottom to top):
 * // Row 0 (sub-diagonal): [*, 1, 1, 1]
 * // Row 1 (main diagonal): [2, 2, 2, 2]
 * // Row 2 (super-diagonal): [1, 1, 1, *]
 *
 * const AB = [
 *   0, 1, 1, 1,  // sub-diagonal (padded at start)
 *   2, 2, 2, 2,  // main diagonal
 *   1, 1, 1, 0,  // super-diagonal (padded at end)
 * ];
 * const b = [1, 2, 2, 1];
 * const { x } = solveBanded(AB, b, 4, 1, 1);
 * ```
 */
export function solveBanded(
  AB: RealArray,
  b: RealArray | Matrix,
  n: number,
  kl: number,
  ku: number
): SolveResult {
  // Convert banded to dense format and use general solver
  const ldab = 2 * kl + ku + 1;
  const abData = Array.isArray(AB) ? AB : Array.from(AB);

  if (abData.length !== ldab * n) {
    throw new Error(
      `Banded matrix AB should have ${ldab * n} elements, got ${abData.length}`
    );
  }

  // Convert band storage to dense matrix (column-major)
  const A = new Float64Array(n * n);

  for (let j = 0; j < n; j++) {
    for (let i = Math.max(0, j - ku); i <= Math.min(n - 1, j + kl); i++) {
      // In band storage: AB[kl + ku + i - j, j] contains A[i, j]
      const bandRow = kl + ku + i - j;
      const bandIdx = bandRow + j * ldab;
      const denseIdx = i + j * n; // column-major
      A[denseIdx] = abData[bandIdx];
    }
  }

  // Use general solver
  return solve(A, b as Matrix, { overwriteA: true });
}
