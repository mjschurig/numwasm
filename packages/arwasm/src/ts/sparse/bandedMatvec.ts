/**
 * Banded Matrix-Vector Product
 *
 * Creates a matvec function from a general banded matrix.
 * Handles matrices with multiple diagonals at various offsets.
 */

import type { MatVecFunction } from '../high-level-types.js';

/**
 * Create a matrix-vector product function from a banded matrix.
 *
 * A banded matrix has non-zero entries only on certain diagonals.
 * Each diagonal is specified by its offset from the main diagonal:
 * - offset = 0: main diagonal
 * - offset > 0: super-diagonals (above main)
 * - offset < 0: sub-diagonals (below main)
 *
 * @param bands - Array of diagonal entries, one Float64Array per band
 * @param offsets - Offset for each band (0 = main diagonal, positive = super, negative = sub)
 * @param n - Matrix dimension (square matrix)
 * @returns Function computing y = A*x
 *
 * @example
 * ```ts
 * // Pentadiagonal matrix (5-point stencil for 2D Laplacian flattened):
 * // Offsets: -n, -1, 0, 1, n where n is the grid size
 * const gridSize = 10;
 * const n = gridSize * gridSize;
 *
 * // Main diagonal
 * const diag = new Float64Array(n).fill(4);
 *
 * // Adjacent diagonals (offset +-1)
 * const adj1 = new Float64Array(n - 1).fill(-1);
 *
 * // Grid row diagonals (offset +-gridSize)
 * const adjN = new Float64Array(n - gridSize).fill(-1);
 *
 * const bands = [adjN, adj1, diag, adj1, adjN];
 * const offsets = [-gridSize, -1, 0, 1, gridSize];
 *
 * const matvec = bandedMatvec(bands, offsets, n);
 * ```
 */
export function bandedMatvec(
  bands: Float64Array[],
  offsets: number[],
  n: number
): MatVecFunction {
  // Validate inputs
  if (bands.length !== offsets.length) {
    throw new Error(`bands and offsets must have same length`);
  }

  // Validate band lengths
  for (let k = 0; k < bands.length; k++) {
    const offset = offsets[k];
    const expectedLength = n - Math.abs(offset);
    if (bands[k].length !== expectedLength) {
      throw new Error(
        `Band at offset ${offset} should have length ${expectedLength}, got ${bands[k].length}`
      );
    }
  }

  return (x: Float64Array): Float64Array => {
    if (x.length !== n) {
      throw new Error(`Input vector length must be ${n}, got ${x.length}`);
    }

    const y = new Float64Array(n);

    for (let k = 0; k < bands.length; k++) {
      const band = bands[k];
      const offset = offsets[k];

      if (offset >= 0) {
        // Super-diagonal or main diagonal
        // A[i, i+offset] = band[i] for i = 0, 1, ..., n-1-offset
        for (let i = 0; i < n - offset; i++) {
          y[i] += band[i] * x[i + offset];
        }
      } else {
        // Sub-diagonal
        // A[i, i+offset] = band[i+offset] for i = -offset, ..., n-1
        // Or equivalently: A[i-offset, i] = band[i] for i = 0, ..., n-1+offset
        const absOffset = -offset;
        for (let i = 0; i < n - absOffset; i++) {
          y[i + absOffset] += band[i] * x[i];
        }
      }
    }

    return y;
  };
}

/**
 * Create a symmetric banded matrix-vector product function.
 *
 * For symmetric banded matrices, only the upper (or lower) bands need
 * to be specified. This function mirrors them automatically.
 *
 * @param bands - Upper bands (including main diagonal), indexed by offset (0, 1, 2, ...)
 * @param n - Matrix dimension
 * @returns Function computing y = A*x
 *
 * @example
 * ```ts
 * // Symmetric pentadiagonal: offsets 0, 1, 2
 * const diag = new Float64Array(n).fill(6);
 * const off1 = new Float64Array(n - 1).fill(-4);
 * const off2 = new Float64Array(n - 2).fill(1);
 *
 * // bands[0] = main diagonal, bands[1] = first off-diagonal, etc.
 * const bands = [diag, off1, off2];
 *
 * const matvec = symBandedMatvec(bands, n);
 * ```
 */
export function symBandedMatvec(
  bands: Float64Array[],
  n: number
): MatVecFunction {
  // Validate inputs
  for (let k = 0; k < bands.length; k++) {
    const expectedLength = n - k;
    if (bands[k].length !== expectedLength) {
      throw new Error(
        `Band at offset ${k} should have length ${expectedLength}, got ${bands[k].length}`
      );
    }
  }

  return (x: Float64Array): Float64Array => {
    if (x.length !== n) {
      throw new Error(`Input vector length must be ${n}, got ${x.length}`);
    }

    const y = new Float64Array(n);

    // Main diagonal
    const diag = bands[0];
    for (let i = 0; i < n; i++) {
      y[i] = diag[i] * x[i];
    }

    // Off-diagonals (both upper and lower due to symmetry)
    for (let k = 1; k < bands.length; k++) {
      const band = bands[k];

      for (let i = 0; i < n - k; i++) {
        // Upper: A[i, i+k] = band[i]
        y[i] += band[i] * x[i + k];
        // Lower: A[i+k, i] = band[i] (symmetry)
        y[i + k] += band[i] * x[i];
      }
    }

    return y;
  };
}

/**
 * Create a Toeplitz banded matrix-vector product function.
 *
 * For banded matrices where each diagonal has a constant value.
 *
 * @param bandValues - Value for each band
 * @param offsets - Offset for each band
 * @param n - Matrix dimension
 * @returns Function computing y = A*x
 *
 * @example
 * ```ts
 * // 1D Laplacian with periodic boundary (circulant approximation):
 * const values = [-1, 2, -1];
 * const offsets = [-1, 0, 1];
 * const matvec = toeplitzBandedMatvec(values, offsets, 100);
 * ```
 */
export function toeplitzBandedMatvec(
  bandValues: number[],
  offsets: number[],
  n: number
): MatVecFunction {
  if (bandValues.length !== offsets.length) {
    throw new Error(`bandValues and offsets must have same length`);
  }

  return (x: Float64Array): Float64Array => {
    if (x.length !== n) {
      throw new Error(`Input vector length must be ${n}, got ${x.length}`);
    }

    const y = new Float64Array(n);

    for (let k = 0; k < bandValues.length; k++) {
      const val = bandValues[k];
      if (val === 0) continue;

      const offset = offsets[k];

      if (offset >= 0) {
        // Super-diagonal or main diagonal
        for (let i = 0; i < n - offset; i++) {
          y[i] += val * x[i + offset];
        }
      } else {
        // Sub-diagonal
        const absOffset = -offset;
        for (let i = 0; i < n - absOffset; i++) {
          y[i + absOffset] += val * x[i];
        }
      }
    }

    return y;
  };
}
