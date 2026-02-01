/**
 * Diagonal Matrix-Vector Product
 *
 * Creates a matvec function from a diagonal matrix.
 * Very efficient O(n) implementation.
 */

import type { MatVecFunction } from '../high-level-types.js';

/**
 * Create a matrix-vector product function from a diagonal matrix.
 *
 * For a diagonal matrix D with diagonal entries d[i], computes:
 * y[i] = d[i] * x[i]
 *
 * @param diagonal - Diagonal entries (length n)
 * @returns Function computing y = D*x
 *
 * @example
 * ```ts
 * // Diagonal matrix: diag([2, 3, 4])
 * const diag = new Float64Array([2, 3, 4]);
 * const matvec = diagMatvec(diag);
 *
 * const x = new Float64Array([1, 2, 3]);
 * const y = matvec(x); // [2, 6, 12]
 * ```
 */
export function diagMatvec(diagonal: Float64Array): MatVecFunction {
  const n = diagonal.length;

  return (x: Float64Array): Float64Array => {
    if (x.length !== n) {
      throw new Error(`Input vector length must be ${n}, got ${x.length}`);
    }

    const y = new Float64Array(n);

    for (let i = 0; i < n; i++) {
      y[i] = diagonal[i] * x[i];
    }

    return y;
  };
}

/**
 * Create an inverse diagonal matrix-vector product function.
 *
 * Computes y = D^{-1} * x where D is diagonal.
 * Useful for preconditioning.
 *
 * @param diagonal - Diagonal entries (must be non-zero)
 * @param tol - Tolerance for "zero" entries (default 1e-14)
 * @returns Function computing y = D^{-1}*x
 * @throws Error if any diagonal entry is smaller than tol
 */
export function diagMatvecInv(
  diagonal: Float64Array,
  tol: number = 1e-14
): MatVecFunction {
  const n = diagonal.length;

  // Pre-compute inverse diagonal
  const invDiag = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    if (Math.abs(diagonal[i]) < tol) {
      throw new Error(`Diagonal entry ${i} is near-zero: ${diagonal[i]}`);
    }
    invDiag[i] = 1 / diagonal[i];
  }

  return (x: Float64Array): Float64Array => {
    if (x.length !== n) {
      throw new Error(`Input vector length must be ${n}, got ${x.length}`);
    }

    const y = new Float64Array(n);

    for (let i = 0; i < n; i++) {
      y[i] = invDiag[i] * x[i];
    }

    return y;
  };
}

/**
 * Create a sqrt-diagonal matrix-vector product function.
 *
 * Computes y = D^{1/2} * x where D is diagonal with non-negative entries.
 * Useful for symmetric normalization.
 *
 * @param diagonal - Diagonal entries (must be non-negative)
 * @returns Function computing y = D^{1/2}*x
 */
export function diagMatvecSqrt(diagonal: Float64Array): MatVecFunction {
  const n = diagonal.length;

  // Pre-compute sqrt diagonal
  const sqrtDiag = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    if (diagonal[i] < 0) {
      throw new Error(`Diagonal entry ${i} is negative: ${diagonal[i]}`);
    }
    sqrtDiag[i] = Math.sqrt(diagonal[i]);
  }

  return (x: Float64Array): Float64Array => {
    if (x.length !== n) {
      throw new Error(`Input vector length must be ${n}, got ${x.length}`);
    }

    const y = new Float64Array(n);

    for (let i = 0; i < n; i++) {
      y[i] = sqrtDiag[i] * x[i];
    }

    return y;
  };
}

/**
 * Create an inverse-sqrt-diagonal matrix-vector product function.
 *
 * Computes y = D^{-1/2} * x where D is diagonal with positive entries.
 * Useful for symmetric normalization of Laplacians.
 *
 * @param diagonal - Diagonal entries (must be positive)
 * @param tol - Tolerance for "zero" entries (default 1e-14)
 * @returns Function computing y = D^{-1/2}*x
 */
export function diagMatvecInvSqrt(
  diagonal: Float64Array,
  tol: number = 1e-14
): MatVecFunction {
  const n = diagonal.length;

  // Pre-compute inverse sqrt diagonal
  const invSqrtDiag = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    if (diagonal[i] < tol) {
      throw new Error(`Diagonal entry ${i} is non-positive or near-zero: ${diagonal[i]}`);
    }
    invSqrtDiag[i] = 1 / Math.sqrt(diagonal[i]);
  }

  return (x: Float64Array): Float64Array => {
    if (x.length !== n) {
      throw new Error(`Input vector length must be ${n}, got ${x.length}`);
    }

    const y = new Float64Array(n);

    for (let i = 0; i < n; i++) {
      y[i] = invSqrtDiag[i] * x[i];
    }

    return y;
  };
}
