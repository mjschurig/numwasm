/**
 * Matrix Square Root
 *
 * Compute the principal square root of a matrix.
 */

import {
  prepareMatrix,
  getMatrixDimensions,
} from '../helpers.js';
import type { Matrix, SqrtmResult } from './types.js';

/**
 * Compute the principal matrix square root.
 *
 * Finds S such that S*S = A, where S has eigenvalues with positive real parts.
 * Uses the Denman-Beavers iteration.
 *
 * @param A - Input square matrix (n × n)
 * @returns Matrix square root S such that S*S = A
 *
 * @example
 * ```typescript
 * const A = [[4, 0], [0, 9]];
 * const { S } = sqrtm(A);
 * // S ≈ [[2, 0], [0, 3]]
 *
 * // For identity matrix
 * const I = [[1, 0], [0, 1]];
 * const { S: SI } = sqrtm(I);
 * // SI = [[1, 0], [0, 1]]
 * ```
 */
export function sqrtm(A: Matrix): SqrtmResult {
  // Get dimensions
  const [m, n] = getMatrixDimensions(A);

  if (m !== n) {
    throw new Error(`Matrix must be square, got ${m}×${n}`);
  }

  if (n === 0) {
    return {
      S: new Float64Array(0),
      n: 0,
      success: true,
      message: 'Success',
    };
  }

  // Prepare matrix (column-major)
  const aData = prepareMatrix(A);

  // Denman-Beavers iteration for matrix square root:
  // Y_0 = A, Z_0 = I
  // Y_{k+1} = (Y_k + Z_k^{-1}) / 2
  // Z_{k+1} = (Z_k + Y_k^{-1}) / 2
  // Y_k -> sqrt(A), Z_k -> sqrt(A)^{-1}

  let Y = new Float64Array(n * n);
  let Z = new Float64Array(n * n);
  const Ynew = new Float64Array(n * n);
  const Znew = new Float64Array(n * n);

  // Initialize: Y = A, Z = I
  for (let i = 0; i < n * n; i++) {
    Y[i] = aData[i];
    Z[i] = 0;
  }
  for (let i = 0; i < n; i++) {
    Z[i + i * n] = 1;
  }

  const maxIter = 100;
  const eps = 1e-14;

  // Helper: compute matrix inverse using Gauss-Jordan elimination
  const matrixInverse = (M: Float64Array): Float64Array | null => {
    const aug = new Float64Array(n * 2 * n);

    // Create augmented matrix [M | I]
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        aug[i + j * n] = M[i + j * n];
        aug[i + (j + n) * n] = i === j ? 1 : 0;
      }
    }

    // Gauss-Jordan elimination
    for (let col = 0; col < n; col++) {
      // Find pivot
      let maxVal = Math.abs(aug[col + col * n]);
      let maxRow = col;
      for (let row = col + 1; row < n; row++) {
        if (Math.abs(aug[row + col * n]) > maxVal) {
          maxVal = Math.abs(aug[row + col * n]);
          maxRow = row;
        }
      }

      if (maxVal < 1e-15) {
        return null; // Singular matrix
      }

      // Swap rows
      if (maxRow !== col) {
        for (let j = 0; j < 2 * n; j++) {
          const temp = aug[col + j * n];
          aug[col + j * n] = aug[maxRow + j * n];
          aug[maxRow + j * n] = temp;
        }
      }

      // Scale pivot row
      const pivot = aug[col + col * n];
      for (let j = 0; j < 2 * n; j++) {
        aug[col + j * n] /= pivot;
      }

      // Eliminate column
      for (let row = 0; row < n; row++) {
        if (row !== col) {
          const factor = aug[row + col * n];
          for (let j = 0; j < 2 * n; j++) {
            aug[row + j * n] -= factor * aug[col + j * n];
          }
        }
      }
    }

    // Extract inverse
    const inv = new Float64Array(n * n);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        inv[i + j * n] = aug[i + (j + n) * n];
      }
    }

    return inv;
  };

  for (let iter = 0; iter < maxIter; iter++) {
    // Compute Z^{-1} and Y^{-1}
    const Zinv = matrixInverse(Z);
    const Yinv = matrixInverse(Y);

    if (!Zinv || !Yinv) {
      return {
        S: Y, // Return best approximation
        n,
        success: false,
        message: 'Matrix became singular during iteration',
      };
    }

    // Y_{k+1} = (Y_k + Z_k^{-1}) / 2
    // Z_{k+1} = (Z_k + Y_k^{-1}) / 2
    let maxChange = 0;
    for (let i = 0; i < n * n; i++) {
      Ynew[i] = (Y[i] + Zinv[i]) / 2;
      Znew[i] = (Z[i] + Yinv[i]) / 2;
      const change = Math.abs(Ynew[i] - Y[i]);
      if (change > maxChange) maxChange = change;
    }

    // Update
    for (let i = 0; i < n * n; i++) {
      Y[i] = Ynew[i];
      Z[i] = Znew[i];
    }

    // Check convergence
    if (maxChange < eps) {
      break;
    }
  }

  return {
    S: Y,
    n,
    success: true,
    message: 'Success',
  };
}
