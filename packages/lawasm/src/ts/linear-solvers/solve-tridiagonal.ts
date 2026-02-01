/**
 * Tridiagonal System Solver
 *
 * Solve Ax = b where A is a tridiagonal matrix using the Thomas algorithm.
 */

import type { RealArray } from '../helpers.js';
import { is2DArray } from '../helpers.js';
import type { Matrix, SolveTridiagonalResult } from './types.js';

/**
 * Solve a tridiagonal system Ax = b where A is tridiagonal.
 *
 * A tridiagonal matrix has non-zero elements only on the main diagonal
 * and the diagonals immediately above and below it.
 *
 * Uses the Thomas algorithm (tridiagonal matrix algorithm) for efficient O(n) solution.
 *
 * @param dl - Sub-diagonal elements (length n-1)
 * @param d - Main diagonal elements (length n)
 * @param du - Super-diagonal elements (length n-1)
 * @param b - Right-hand side vector or matrix
 * @returns Solution x and diagnostic information
 *
 * @example
 * ```typescript
 * // Solve tridiagonal system:
 * // [2  1  0] [x0]   [1]
 * // [1  2  1] [x1] = [2]
 * // [0  1  2] [x2]   [1]
 *
 * const dl = [1, 1];     // sub-diagonal
 * const d = [2, 2, 2];   // main diagonal
 * const du = [1, 1];     // super-diagonal
 * const b = [1, 2, 1];
 *
 * const { x } = solveTridiagonal(dl, d, du, b);
 * ```
 */
export function solveTridiagonal(
  dl: RealArray,
  d: RealArray,
  du: RealArray,
  b: RealArray | Matrix
): SolveTridiagonalResult {
  // Convert to arrays
  const dlArr = Array.isArray(dl) ? dl : Array.from(dl);
  const dArr = Array.isArray(d) ? d : Array.from(d);
  const duArr = Array.isArray(du) ? du : Array.from(du);

  const n = dArr.length;

  // Validate dimensions
  if (dlArr.length !== n - 1) {
    throw new Error(
      `Sub-diagonal dl should have ${n - 1} elements, got ${dlArr.length}`
    );
  }
  if (duArr.length !== n - 1) {
    throw new Error(
      `Super-diagonal du should have ${n - 1} elements, got ${duArr.length}`
    );
  }

  // Prepare b
  let bData: Float64Array;
  let nrhs: number;

  if (is2DArray(b)) {
    // Multiple right-hand sides
    nrhs = (b as number[][])[0].length;
    bData = new Float64Array(n * nrhs);
    // Convert to column-major
    for (let j = 0; j < nrhs; j++) {
      for (let i = 0; i < n; i++) {
        bData[j * n + i] = (b as number[][])[i][j];
      }
    }
  } else {
    const bArr = Array.isArray(b) ? b : Array.from(b);
    if (bArr.length === n) {
      nrhs = 1;
      bData = new Float64Array(bArr);
    } else {
      nrhs = Math.floor(bArr.length / n);
      bData = new Float64Array(bArr);
    }
  }

  // Thomas algorithm (tridiagonal matrix algorithm)
  // This is a direct O(n) solver for tridiagonal systems

  // Make copies to avoid modifying input
  const c = new Float64Array(n - 1); // modified super-diagonal
  const dMod = new Float64Array(n); // modified diagonal
  const x = new Float64Array(n * nrhs);

  try {
    for (let rhs = 0; rhs < nrhs; rhs++) {
      const bOffset = rhs * n;
      const xOffset = rhs * n;

      // Copy b for this RHS
      const bCol = new Float64Array(n);
      for (let i = 0; i < n; i++) {
        bCol[i] = bData[bOffset + i];
      }

      // Forward elimination
      dMod[0] = dArr[0];
      if (Math.abs(dMod[0]) < 1e-15) {
        return {
          x: new Float64Array(0),
          nrhs,
          info: 1,
          success: false,
          message: 'Zero pivot encountered at row 1',
        };
      }

      c[0] = duArr[0] / dMod[0];
      bCol[0] = bCol[0] / dMod[0];

      for (let i = 1; i < n; i++) {
        const m = dlArr[i - 1];
        dMod[i] = dArr[i] - m * c[i - 1];

        if (Math.abs(dMod[i]) < 1e-15) {
          return {
            x: new Float64Array(0),
            nrhs,
            info: i + 1,
            success: false,
            message: `Zero pivot encountered at row ${i + 1}`,
          };
        }

        if (i < n - 1) {
          c[i] = duArr[i] / dMod[i];
        }
        bCol[i] = (bCol[i] - m * bCol[i - 1]) / dMod[i];
      }

      // Back substitution
      x[xOffset + n - 1] = bCol[n - 1];
      for (let i = n - 2; i >= 0; i--) {
        x[xOffset + i] = bCol[i] - c[i] * x[xOffset + i + 1];
      }
    }

    return {
      x,
      nrhs,
      info: 0,
      success: true,
      message: 'Success',
    };
  } catch (e) {
    return {
      x: new Float64Array(0),
      nrhs,
      info: -1,
      success: false,
      message: `Error during solve: ${e instanceof Error ? e.message : String(e)}`,
    };
  }
}
