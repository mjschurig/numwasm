/**
 * Hessenberg Reduction
 *
 * Reduce a general matrix to upper Hessenberg form.
 *
 * Note: DGEHRD is not currently exported. This provides a pure TypeScript
 * implementation using Householder reflections.
 */

import type { Matrix, HessenbergOptions, HessenbergResult } from './types.js';
import { getMatrixDimensions, prepareMatrix } from '../helpers.js';

/**
 * Reduce a general n×n matrix A to upper Hessenberg form.
 *
 * A = Q * H * Q^T
 *
 * where H is upper Hessenberg (zero below the first subdiagonal)
 * and Q is orthogonal.
 *
 * Upper Hessenberg form is useful as a preprocessing step for
 * eigenvalue computations.
 *
 * @param A - Input square matrix (n × n).
 * @param options - Reduction options
 * @returns Hessenberg reduction result
 *
 * @example
 * ```typescript
 * const A = [[1, 2, 3], [4, 5, 6], [7, 8, 9]];
 * const { H, Q, success } = hessenberg(A);
 * // H is upper Hessenberg (H[i,j] = 0 for i > j+1)
 * ```
 */
export function hessenberg(A: Matrix, options: HessenbergOptions = {}): HessenbergResult {
  const { computeQ = true } = options;

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);
  if (m !== n) {
    throw new Error(`Matrix A must be square, got ${m}×${n}`);
  }

  // Prepare matrix - work with a copy
  const aData = prepareMatrix(A);
  const H = new Float64Array(n * n);
  for (let i = 0; i < n * n; i++) {
    H[i] = aData[i];
  }

  // Initialize Q as identity if needed
  let Q: Float64Array | undefined;
  if (computeQ) {
    Q = new Float64Array(n * n);
    for (let i = 0; i < n; i++) {
      Q[i * n + i] = 1.0;
    }
  }

  try {
    // Householder reduction to Hessenberg form
    for (let k = 0; k < n - 2; k++) {
      // Compute Householder vector for column k, rows k+1:n
      const colStart = k + 1;
      const colLen = n - colStart;

      // Extract the column segment
      const x = new Float64Array(colLen);
      for (let i = 0; i < colLen; i++) {
        x[i] = H[k * n + colStart + i];
      }

      // Compute norm
      let norm = 0;
      for (let i = 0; i < colLen; i++) {
        norm += x[i] * x[i];
      }
      norm = Math.sqrt(norm);

      if (norm < 1e-15) {
        continue; // Skip if column is already zero
      }

      // Compute Householder vector v such that (I - 2*v*v^T) * x = ||x|| * e1
      const sign = x[0] >= 0 ? 1 : -1;
      const alpha = -sign * norm;

      // v = x - alpha * e1, normalized
      x[0] -= alpha;
      let vnorm = 0;
      for (let i = 0; i < colLen; i++) {
        vnorm += x[i] * x[i];
      }
      vnorm = Math.sqrt(vnorm);

      if (vnorm < 1e-15) {
        continue;
      }

      for (let i = 0; i < colLen; i++) {
        x[i] /= vnorm;
      }

      // Apply Householder transformation from the left: H = (I - 2*v*v^T) * H
      // Only affects rows k+1:n
      for (let j = k; j < n; j++) {
        // Compute v^T * H[:,j] for rows k+1:n
        let dot = 0;
        for (let i = 0; i < colLen; i++) {
          dot += x[i] * H[j * n + colStart + i];
        }
        // H[:,j] -= 2 * dot * v
        for (let i = 0; i < colLen; i++) {
          H[j * n + colStart + i] -= 2 * dot * x[i];
        }
      }

      // Apply Householder transformation from the right: H = H * (I - 2*v*v^T)
      // Affects all rows, columns k+1:n
      for (let i = 0; i < n; i++) {
        // Compute H[i,:] * v for columns k+1:n
        let dot = 0;
        for (let j = 0; j < colLen; j++) {
          dot += H[(colStart + j) * n + i] * x[j];
        }
        // H[i,:] -= 2 * dot * v^T
        for (let j = 0; j < colLen; j++) {
          H[(colStart + j) * n + i] -= 2 * dot * x[j];
        }
      }

      // Accumulate Q if requested
      if (Q) {
        // Q = Q * (I - 2*v*v^T)
        for (let i = 0; i < n; i++) {
          let dot = 0;
          for (let j = 0; j < colLen; j++) {
            dot += Q[(colStart + j) * n + i] * x[j];
          }
          for (let j = 0; j < colLen; j++) {
            Q[(colStart + j) * n + i] -= 2 * dot * x[j];
          }
        }
      }
    }

    // Clean up numerical noise below subdiagonal
    for (let j = 0; j < n - 2; j++) {
      for (let i = j + 2; i < n; i++) {
        H[j * n + i] = 0;
      }
    }

    return {
      H,
      Q,
      n,
      info: 0,
      success: true,
      message: 'Success',
    };
  } catch (e) {
    return {
      H: new Float64Array(0),
      n,
      info: -1,
      success: false,
      message: `Error during reduction: ${e instanceof Error ? e.message : String(e)}`,
    };
  }
}
