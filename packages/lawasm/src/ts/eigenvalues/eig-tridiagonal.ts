/**
 * Tridiagonal Eigenvalue Problem
 *
 * Compute eigenvalues and eigenvectors of a symmetric tridiagonal matrix.
 *
 * Note: DSTEV is not currently exported from LAPACK WASM.
 * This provides a pure TypeScript implementation using the symmetric QR algorithm.
 */

import type { EigTridiagonalOptions, EigTridiagonalResult } from './types.js';
import type { RealArray } from '../helpers.js';

/**
 * Compute eigenvalues and eigenvectors of a symmetric tridiagonal matrix.
 *
 * The matrix is defined by:
 * - d: diagonal elements [d0, d1, ..., d(n-1)]
 * - e: off-diagonal elements [e0, e1, ..., e(n-2)]
 *
 * The matrix T is:
 * [d0  e0  0   0  ...]
 * [e0  d1  e1  0  ...]
 * [0   e1  d2  e2 ...]
 * [.................]
 *
 * @param d - Diagonal elements (length n)
 * @param e - Off-diagonal elements (length n-1)
 * @param options - Computation options
 * @returns Eigenvalues (in ascending order) and optionally eigenvectors
 *
 * @example
 * ```typescript
 * // Tridiagonal matrix:
 * // [2 1 0]
 * // [1 2 1]
 * // [0 1 2]
 * const d = [2, 2, 2];
 * const e = [1, 1];
 * const { values, vectors } = eigTridiagonal(d, e);
 * ```
 */
export function eigTridiagonal(
  d: RealArray,
  e: RealArray,
  options: EigTridiagonalOptions = {}
): EigTridiagonalResult {
  const { computeVectors = true } = options;

  const n = d.length;
  if (e.length !== n - 1 && n > 1) {
    throw new Error(`Expected ${n - 1} off-diagonal elements, got ${e.length}`);
  }

  if (n === 0) {
    return {
      values: new Float64Array(0),
      n: 0,
      info: 0,
      success: true,
      message: 'Success',
    };
  }

  if (n === 1) {
    const result: EigTridiagonalResult = {
      values: new Float64Array([d[0]]),
      n: 1,
      info: 0,
      success: true,
      message: 'Success',
    };
    if (computeVectors) {
      result.vectors = new Float64Array([1]);
    }
    return result;
  }

  try {
    // Work with copies
    const diag = new Float64Array(d);
    const offdiag = new Float64Array(n - 1);
    for (let i = 0; i < n - 1; i++) {
      offdiag[i] = e[i];
    }

    // Initialize eigenvectors as identity matrix (column-major)
    let Z: Float64Array | undefined;
    if (computeVectors) {
      Z = new Float64Array(n * n);
      for (let i = 0; i < n; i++) {
        Z[i * n + i] = 1;
      }
    }

    // Symmetric QR algorithm with implicit shifts
    const maxIterations = 30 * n;
    const eps = Number.EPSILON;
    let m = n - 1;
    let iterations = 0;

    while (m > 0 && iterations < maxIterations) {
      iterations++;

      // Find the largest off-diagonal element in the unreduced block
      let l = m;
      while (l > 0) {
        const test = Math.abs(diag[l - 1]) + Math.abs(diag[l]);
        if (Math.abs(offdiag[l - 1]) <= eps * test) {
          offdiag[l - 1] = 0;
          break;
        }
        l--;
      }

      if (l === m) {
        // Eigenvalue found
        m--;
        continue;
      }

      // Wilkinson shift
      const d0 = diag[m - 1];
      const d1 = diag[m];
      const e0 = offdiag[m - 1];
      const delta = (d0 - d1) / 2;
      const sign = delta >= 0 ? 1 : -1;
      const shift = d1 - (e0 * e0) / (delta + sign * Math.sqrt(delta * delta + e0 * e0));

      // Apply shift
      let x = diag[l] - shift;
      let z = offdiag[l];

      // Chase the bulge
      for (let k = l; k < m; k++) {
        // Compute Givens rotation
        let c: number, s: number;
        if (Math.abs(z) > Math.abs(x)) {
          const t = -x / z;
          s = 1 / Math.sqrt(1 + t * t);
          c = s * t;
        } else if (Math.abs(x) > 0) {
          const t = -z / x;
          c = 1 / Math.sqrt(1 + t * t);
          s = c * t;
        } else {
          c = 1;
          s = 0;
        }

        // Apply rotation to tridiagonal matrix
        if (k > l) {
          offdiag[k - 1] = c * offdiag[k - 1] - s * z;
        }

        const dk = diag[k];
        const dk1 = diag[k + 1];
        const ek = offdiag[k];

        const p = c * dk - s * ek;
        const q = c * ek - s * dk1;
        const r = s * dk + c * ek;

        diag[k] = c * p - s * q;
        diag[k + 1] = s * r + c * dk1;
        offdiag[k] = s * p + c * q;

        if (k < m - 1) {
          x = offdiag[k];
          z = -s * offdiag[k + 1];
          offdiag[k + 1] = c * offdiag[k + 1];
        }

        // Apply rotation to eigenvector matrix
        if (Z) {
          for (let i = 0; i < n; i++) {
            const zik = Z[k * n + i];
            const zik1 = Z[(k + 1) * n + i];
            Z[k * n + i] = c * zik - s * zik1;
            Z[(k + 1) * n + i] = s * zik + c * zik1;
          }
        }
      }
    }

    if (iterations >= maxIterations) {
      return {
        values: new Float64Array(0),
        n,
        info: 1,
        success: false,
        message: 'Failed to converge',
      };
    }

    // Sort eigenvalues (and eigenvectors) in ascending order
    const indices = Array.from({ length: n }, (_, i) => i);
    indices.sort((a, b) => diag[a] - diag[b]);

    const sortedValues = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      sortedValues[i] = diag[indices[i]];
    }

    let sortedVectors: Float64Array | undefined;
    if (Z) {
      sortedVectors = new Float64Array(n * n);
      for (let j = 0; j < n; j++) {
        const srcCol = indices[j];
        for (let i = 0; i < n; i++) {
          sortedVectors[j * n + i] = Z[srcCol * n + i];
        }
      }
    }

    return {
      values: sortedValues,
      vectors: sortedVectors,
      n,
      info: 0,
      success: true,
      message: 'Success',
    };
  } catch (err) {
    return {
      values: new Float64Array(0),
      n,
      info: -1,
      success: false,
      message: `Error: ${err instanceof Error ? err.message : String(err)}`,
    };
  }
}
