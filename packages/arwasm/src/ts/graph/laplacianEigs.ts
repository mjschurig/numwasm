/**
 * LAPLACIANEIGS - Graph Laplacian Eigenvalues
 *
 * Computes eigenvalues and eigenvectors of the graph Laplacian matrix
 * for spectral clustering and graph analysis applications.
 *
 * The graph Laplacian L = D - A where:
 * - A is the adjacency matrix
 * - D is the degree matrix (diagonal with row sums of A)
 *
 * For normalized Laplacian: L_norm = D^(-1/2) * L * D^(-1/2) = I - D^(-1/2) * A * D^(-1/2)
 */

import { eigs } from '../core/eigs.js';
import type { MatVecFunction, EigsResult } from '../high-level-types.js';

/**
 * Options for graph Laplacian eigenvalue computation.
 */
export interface LaplacianEigsOptions {
  /**
   * Use normalized Laplacian (recommended for spectral clustering).
   * @default true
   */
  normalized?: boolean;

  /**
   * Which eigenvalues to compute.
   * 'SM' is typically used for spectral clustering (smallest eigenvalues).
   * @default 'SM'
   */
  which?: 'SM' | 'LM';

  /**
   * Drop the first trivial eigenvalue (0 for connected graphs).
   * @default true
   */
  dropFirst?: boolean;

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
}

/**
 * Result from Laplacian eigenvalue computation.
 */
export interface LaplacianEigsResult extends EigsResult {
  /**
   * Fiedler value (second smallest eigenvalue).
   * Indicates algebraic connectivity of the graph.
   */
  fiedlerValue?: number;

  /**
   * Fiedler vector (eigenvector for second smallest eigenvalue).
   * Used for graph bipartitioning.
   */
  fiedlerVector?: Float64Array;
}

/**
 * Compute eigenvalues of the graph Laplacian for spectral clustering.
 *
 * This function computes the smallest eigenvalues of the graph Laplacian,
 * which are used for spectral clustering and graph partitioning.
 *
 * @param adjacencyMatvec - Function computing y = A*x where A is the adjacency matrix
 * @param degrees - Degree of each node (row sums of adjacency matrix), or function to compute D*x
 * @param n - Number of nodes in the graph
 * @param nev - Number of eigenvalues to compute
 * @param options - Solver options
 * @returns Laplacian eigenvalues and eigenvectors
 *
 * @example
 * ```ts
 * import { laplacianEigs } from 'arwasm';
 *
 * // For a simple graph with adjacency matrix
 * const n = 100;
 * const adjacency = ...; // Your sparse adjacency matrix
 *
 * const adjacencyMatvec = (x: Float64Array): Float64Array => {
 *   // Compute y = A*x
 *   return y;
 * };
 *
 * // Compute degrees (row sums)
 * const degrees = new Float64Array(n);
 * for (let i = 0; i < n; i++) {
 *   degrees[i] = adjacency.rowSum(i);
 * }
 *
 * // Find 10 smallest Laplacian eigenvalues
 * const result = await laplacianEigs(adjacencyMatvec, degrees, n, 10);
 * console.log('Fiedler value:', result.fiedlerValue);
 * ```
 */
export async function laplacianEigs(
  adjacencyMatvec: MatVecFunction,
  degrees: Float64Array | MatVecFunction,
  n: number,
  nev: number,
  options?: LaplacianEigsOptions
): Promise<LaplacianEigsResult> {
  const {
    normalized = true,
    which = 'SM',
    dropFirst = true,
    tol = 0,
    ncv,
    maxiter,
  } = options ?? {};

  // Convert degrees to array if it's a function
  let degreeArray: Float64Array;
  if (degrees instanceof Float64Array) {
    degreeArray = degrees;
  } else {
    // Compute degrees by applying D to a ones vector
    const ones = new Float64Array(n).fill(1);
    const result = degrees(ones);
    degreeArray = result instanceof Float64Array ? result : new Float64Array(result);
  }

  // Compute inverse sqrt of degrees for normalized Laplacian
  const invSqrtDeg = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    invSqrtDeg[i] = degreeArray[i] > 0 ? 1 / Math.sqrt(degreeArray[i]) : 0;
  }

  // Create Laplacian matvec
  // Unnormalized: L*x = D*x - A*x
  // Normalized: L_norm*x = x - D^(-1/2) * A * D^(-1/2) * x
  const laplacianMatvec: MatVecFunction = normalized
    ? (x: Float64Array) => {
        // Compute D^(-1/2) * x
        const scaled = new Float64Array(n);
        for (let i = 0; i < n; i++) {
          scaled[i] = invSqrtDeg[i] * x[i];
        }

        // Compute A * (D^(-1/2) * x)
        const Ax = adjacencyMatvec(scaled);
        const AxArr = Ax instanceof Float64Array ? Ax : new Float64Array(Ax);

        // Compute D^(-1/2) * A * D^(-1/2) * x
        const result = new Float64Array(n);
        for (let i = 0; i < n; i++) {
          // L_norm * x = x - D^(-1/2) * A * D^(-1/2) * x
          result[i] = x[i] - invSqrtDeg[i] * AxArr[i];
        }
        return result;
      }
    : (x: Float64Array) => {
        // Compute A*x
        const Ax = adjacencyMatvec(x);
        const AxArr = Ax instanceof Float64Array ? Ax : new Float64Array(Ax);

        // L*x = D*x - A*x
        const result = new Float64Array(n);
        for (let i = 0; i < n; i++) {
          result[i] = degreeArray[i] * x[i] - AxArr[i];
        }
        return result;
      };

  // Request one extra eigenvalue if we're dropping the first
  const requestedNev = dropFirst ? nev + 1 : nev;

  const result = await eigs(laplacianMatvec, n, requestedNev, {
    which,
    tol,
    ncv,
    maxiter,
    return_eigenvectors: true,
  });

  if (!result.success) {
    return {
      ...result,
      fiedlerValue: undefined,
      fiedlerVector: undefined,
    };
  }

  // Process results
  let eigenvalues = result.eigenvalues;
  let eigenvectors = result.eigenvectors;

  // Drop first eigenvalue (should be ~0 for connected graphs)
  if (dropFirst && result.nconv > 1) {
    eigenvalues = eigenvalues.slice(1);
    eigenvectors = eigenvectors?.slice(1);
  }

  // Extract Fiedler value and vector
  const fiedlerValue = eigenvalues.length > 0 ? eigenvalues[0] : undefined;
  const fiedlerVector = eigenvectors && eigenvectors.length > 0 ? eigenvectors[0] : undefined;

  return {
    eigenvalues,
    eigenvectors,
    niter: result.niter,
    nops: result.nops,
    nconv: eigenvalues.length,
    info: result.info,
    success: true,
    message: result.message,
    fiedlerValue,
    fiedlerVector,
  };
}
