/**
 * EIGSDEFLATED - Deflated Eigenvalue Solver
 *
 * Computes eigenvalues while excluding (deflating) known eigenvectors.
 * This is useful for:
 * 1. Computing additional eigenvalues after an initial solve
 * 2. Finding eigenvalues orthogonal to a known subspace
 * 3. Avoiding repeated computation of known eigenpairs
 *
 * The deflation is performed by modifying the operator:
 *   A_deflated = A - Σ λ_i * v_i * v_i^T
 *
 * or equivalently, by orthogonalizing against known eigenvectors in the
 * Arnoldi/Lanczos iteration.
 */

import { eigs } from '../core/eigs.js';
import { eign } from '../core/eign.js';
import type {
  MatVecFunction,
  EigsResult,
  EignResult,
} from '../high-level-types.js';
import type { WhichSymmetric, WhichNonSymmetric } from '../types.js';

/**
 * Options for deflated eigenvalue computation.
 */
export interface EigsDeflatedOptions {
  /**
   * Which eigenvalues to compute.
   * @default 'LM'
   */
  which?: WhichSymmetric;

  /**
   * Convergence tolerance.
   * @default 0 (machine precision)
   */
  tol?: number;

  /**
   * Number of Lanczos vectors.
   */
  ncv?: number;

  /**
   * Maximum iterations.
   */
  maxiter?: number;

  /**
   * Whether to compute eigenvectors.
   * @default true
   */
  return_eigenvectors?: boolean;
}

/**
 * Compute the 2-norm of a vector.
 */
function norm2(v: Float64Array): number {
  let sum = 0;
  for (let i = 0; i < v.length; i++) {
    sum += v[i] * v[i];
  }
  return Math.sqrt(sum);
}

/**
 * Compute dot product of two vectors.
 */
function dot(a: Float64Array, b: Float64Array): number {
  let sum = 0;
  const n = a.length;
  for (let i = 0; i < n; i++) {
    sum += a[i] * b[i];
  }
  return sum;
}

/**
 * Orthogonalize vector v against a set of orthonormal vectors.
 * Returns the modified vector (in place).
 */
function orthogonalize(v: Float64Array, Q: Float64Array[]): Float64Array {
  for (const q of Q) {
    const coeff = dot(v, q);
    for (let i = 0; i < v.length; i++) {
      v[i] -= coeff * q[i];
    }
  }
  return v;
}

/**
 * Compute eigenvalues with deflation against known eigenvectors.
 *
 * This modifies the matrix operator to project out the known eigenvectors,
 * ensuring that the new eigenvalues are orthogonal to the known subspace.
 *
 * @param matvec - Symmetric matrix operator A*x
 * @param n - Matrix dimension
 * @param nev - Number of new eigenvalues to compute
 * @param knownEigenvectors - Known eigenvectors to deflate
 * @param knownEigenvalues - Eigenvalues corresponding to known eigenvectors (optional)
 * @param options - Solver options
 * @returns New eigenvalues orthogonal to the known eigenvectors
 *
 * @example
 * ```ts
 * // First, compute 6 eigenvalues
 * const result1 = await eigs(matvec, n, 6);
 *
 * // Then compute 6 more, excluding the first 6
 * const result2 = await eigsDeflated(
 *   matvec, n, 6,
 *   result1.eigenvectors!,
 *   result1.eigenvalues
 * );
 * ```
 */
export async function eigsDeflated(
  matvec: MatVecFunction,
  n: number,
  nev: number,
  knownEigenvectors: Float64Array[],
  knownEigenvalues?: Float64Array,
  options?: EigsDeflatedOptions
): Promise<EigsResult> {
  const {
    which = 'LM',
    tol = 0,
    ncv,
    maxiter,
    return_eigenvectors = true,
  } = options ?? {};

  // Ensure eigenvectors are normalized
  const Q: Float64Array[] = [];
  for (const v of knownEigenvectors) {
    const normV = norm2(v);
    if (normV > 1e-14) {
      const normalized = new Float64Array(v.length);
      for (let i = 0; i < v.length; i++) {
        normalized[i] = v[i] / normV;
      }
      Q.push(normalized);
    }
  }

  // Create deflated operator
  // If we have eigenvalues, use explicit deflation: A - Σ λ_i * v_i * v_i^T
  // Otherwise, use projection: (I - Q*Q^T) * A * (I - Q*Q^T)
  let deflatedMatvec: MatVecFunction;

  if (knownEigenvalues && knownEigenvalues.length === Q.length) {
    // Explicit deflation
    deflatedMatvec = (x: Float64Array): Float64Array => {
      const Ax = matvec(x);
      const AxArr = Ax instanceof Float64Array ? Ax : new Float64Array(Ax);
      const result = new Float64Array(AxArr);

      for (let i = 0; i < Q.length; i++) {
        const lambda = knownEigenvalues[i];
        const viTx = dot(Q[i], x);
        for (let j = 0; j < n; j++) {
          result[j] -= lambda * viTx * Q[i][j];
        }
      }

      return result;
    };
  } else {
    // Projection-based deflation
    deflatedMatvec = (x: Float64Array): Float64Array => {
      // Project x: x' = x - Q*Q^T*x
      const xProj = new Float64Array(x);
      orthogonalize(xProj, Q);

      // Apply A
      const AxProj = matvec(xProj);
      const AxProjArr = AxProj instanceof Float64Array ? AxProj : new Float64Array(AxProj);

      // Project result: y = (I - Q*Q^T) * A * x'
      orthogonalize(AxProjArr, Q);

      return AxProjArr;
    };
  }

  // Compute eigenvalues of deflated operator
  const result = await eigs(deflatedMatvec, n, nev, {
    which,
    tol,
    ncv,
    maxiter,
    return_eigenvectors,
  });

  // Ensure computed eigenvectors are orthogonal to known ones
  if (result.eigenvectors) {
    for (const v of result.eigenvectors) {
      orthogonalize(v, Q);
      // Re-normalize
      const normV = norm2(v);
      if (normV > 1e-14) {
        for (let i = 0; i < v.length; i++) {
          v[i] /= normV;
        }
      }
    }
  }

  return result;
}

/**
 * Continue eigenvalue computation from a previous result.
 *
 * Computes additional eigenvalues while excluding those already found.
 *
 * @param matvec - Matrix operator A*x
 * @param n - Matrix dimension
 * @param additionalNev - Number of additional eigenvalues to compute
 * @param previousResult - Result from previous eigs call
 * @param options - Solver options
 * @returns Combined result with all eigenvalues
 */
export async function eigsContinue(
  matvec: MatVecFunction,
  n: number,
  additionalNev: number,
  previousResult: EigsResult,
  options?: EigsDeflatedOptions
): Promise<EigsResult> {
  if (!previousResult.eigenvectors) {
    throw new Error('Previous result must include eigenvectors for continuation');
  }

  const newResult = await eigsDeflated(
    matvec,
    n,
    additionalNev,
    previousResult.eigenvectors,
    previousResult.eigenvalues,
    options
  );

  if (!newResult.success) {
    return newResult;
  }

  // Combine results
  const totalConv = previousResult.nconv + newResult.nconv;
  const combinedEigenvalues = new Float64Array(totalConv);
  const combinedEigenvectors: Float64Array[] = [];

  // Copy previous results
  for (let i = 0; i < previousResult.nconv; i++) {
    combinedEigenvalues[i] = previousResult.eigenvalues[i];
    if (previousResult.eigenvectors) {
      combinedEigenvectors.push(previousResult.eigenvectors[i]);
    }
  }

  // Copy new results
  for (let i = 0; i < newResult.nconv; i++) {
    combinedEigenvalues[previousResult.nconv + i] = newResult.eigenvalues[i];
    if (newResult.eigenvectors) {
      combinedEigenvectors.push(newResult.eigenvectors[i]);
    }
  }

  return {
    eigenvalues: combinedEigenvalues,
    eigenvectors: combinedEigenvectors.length > 0 ? combinedEigenvectors : undefined,
    niter: previousResult.niter + newResult.niter,
    nops: previousResult.nops + newResult.nops,
    nconv: totalConv,
    info: 0,
    success: true,
    message: `Computed ${totalConv} eigenvalues total (${previousResult.nconv} + ${newResult.nconv})`,
  };
}

/**
 * Options for non-symmetric deflated eigenvalue computation.
 */
export interface EignDeflatedOptions {
  /**
   * Which eigenvalues to compute.
   * @default 'LM'
   */
  which?: WhichNonSymmetric;

  /**
   * Convergence tolerance.
   * @default 0 (machine precision)
   */
  tol?: number;

  /**
   * Number of Arnoldi vectors.
   */
  ncv?: number;

  /**
   * Maximum iterations.
   */
  maxiter?: number;

  /**
   * Whether to compute eigenvectors.
   * @default true
   */
  return_eigenvectors?: boolean;
}

/**
 * Compute non-symmetric eigenvalues with deflation.
 *
 * @param matvec - Non-symmetric matrix operator A*x
 * @param n - Matrix dimension
 * @param nev - Number of eigenvalues to compute
 * @param knownEigenvectors - Known eigenvectors to deflate
 * @param options - Solver options
 * @returns New eigenvalues
 */
export async function eignDeflated(
  matvec: MatVecFunction,
  n: number,
  nev: number,
  knownEigenvectors: Float64Array[],
  options?: EignDeflatedOptions
): Promise<EignResult> {
  const {
    which = 'LM',
    tol = 0,
    ncv,
    maxiter,
    return_eigenvectors = true,
  } = options ?? {};

  // Normalize eigenvectors
  const Q: Float64Array[] = [];
  for (const v of knownEigenvectors) {
    const normV = norm2(v);
    if (normV > 1e-14) {
      const normalized = new Float64Array(v.length);
      for (let i = 0; i < v.length; i++) {
        normalized[i] = v[i] / normV;
      }
      Q.push(normalized);
    }
  }

  // Projection-based deflation for non-symmetric case
  const deflatedMatvec: MatVecFunction = (x: Float64Array): Float64Array => {
    const xProj = new Float64Array(x);
    orthogonalize(xProj, Q);

    const AxProj = matvec(xProj);
    const AxProjArr = AxProj instanceof Float64Array ? AxProj : new Float64Array(AxProj);

    orthogonalize(AxProjArr, Q);

    return AxProjArr;
  };

  return eign(deflatedMatvec, n, nev, {
    which,
    tol,
    ncv,
    maxiter,
    return_eigenvectors,
  });
}
