/**
 * SPECTRALEMBEDDING - Laplacian Eigenmaps for Manifold Learning
 *
 * Projects high-dimensional data to a lower-dimensional space using
 * the eigenvectors of the graph Laplacian. This preserves local
 * neighborhood structure and is useful for non-linear dimensionality
 * reduction.
 *
 * Also known as Laplacian Eigenmaps.
 */

import { eigs } from '../core/eigs.js';
import type { MatVecFunction } from '../high-level-types.js';

/**
 * Options for spectral embedding.
 */
export interface SpectralEmbeddingOptions {
  /**
   * Use normalized Laplacian.
   * @default true
   */
  normalized?: boolean;

  /**
   * Drop the first trivial eigenvector.
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
 * Result from spectral embedding.
 */
export interface SpectralEmbeddingResult {
  /**
   * Embedded coordinates (nComponents x n).
   * embedding[i] is the i-th component for all n points.
   */
  embedding: Float64Array[];

  /**
   * Eigenvalues corresponding to each component.
   */
  eigenvalues: Float64Array;

  /**
   * Number of iterations.
   */
  niter: number;

  /**
   * Whether computation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Compute spectral embedding (Laplacian Eigenmaps) for dimensionality reduction.
 *
 * Projects data points to a lower-dimensional space by computing the
 * smallest eigenvectors of the graph Laplacian of an affinity matrix.
 *
 * @param affinityMatvec - Function computing y = W*x where W is the affinity matrix
 * @param degrees - Degree of each node (row sums of W), or function computing D*x
 * @param n - Number of data points
 * @param nComponents - Number of dimensions in the embedding
 * @param options - Solver options
 * @returns Low-dimensional embedding
 *
 * @example
 * ```ts
 * import { spectralEmbedding } from 'arwasm';
 *
 * // Build k-nearest neighbors affinity matrix
 * const n = 1000;
 * const affinity = buildKNNAffinity(data, k=10);
 *
 * const affinityMatvec = (x: Float64Array): Float64Array => {
 *   // Compute y = W*x
 *   return affinity.matvec(x);
 * };
 *
 * const degrees = affinity.rowSums();
 *
 * // Embed into 2D
 * const result = await spectralEmbedding(affinityMatvec, degrees, n, 2);
 * const [x, y] = result.embedding;
 * // Plot x vs y for visualization
 * ```
 */
export async function spectralEmbedding(
  affinityMatvec: MatVecFunction,
  degrees: Float64Array | MatVecFunction,
  n: number,
  nComponents: number,
  options?: SpectralEmbeddingOptions
): Promise<SpectralEmbeddingResult> {
  const {
    normalized = true,
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
  const laplacianMatvec: MatVecFunction = normalized
    ? (x: Float64Array) => {
        const scaled = new Float64Array(n);
        for (let i = 0; i < n; i++) {
          scaled[i] = invSqrtDeg[i] * x[i];
        }

        const Wx = affinityMatvec(scaled);
        const WxArr = Wx instanceof Float64Array ? Wx : new Float64Array(Wx);

        const result = new Float64Array(n);
        for (let i = 0; i < n; i++) {
          result[i] = x[i] - invSqrtDeg[i] * WxArr[i];
        }
        return result;
      }
    : (x: Float64Array) => {
        const Wx = affinityMatvec(x);
        const WxArr = Wx instanceof Float64Array ? Wx : new Float64Array(Wx);

        const result = new Float64Array(n);
        for (let i = 0; i < n; i++) {
          result[i] = degreeArray[i] * x[i] - WxArr[i];
        }
        return result;
      };

  // Request extra eigenvectors if dropping first
  const requestedNev = dropFirst ? nComponents + 1 : nComponents;

  const result = await eigs(laplacianMatvec, n, requestedNev, {
    which: 'SM', // Smallest eigenvalues for embedding
    tol,
    ncv,
    maxiter,
    return_eigenvectors: true,
  });

  if (!result.success || !result.eigenvectors) {
    return {
      embedding: [],
      eigenvalues: new Float64Array(0),
      niter: result.niter,
      success: false,
      message: result.message,
    };
  }

  // Extract embedding from eigenvectors
  let eigenvectors = result.eigenvectors;
  let eigenvalues = result.eigenvalues;

  // Drop first trivial eigenvector if requested
  if (dropFirst && eigenvectors.length > nComponents) {
    eigenvectors = eigenvectors.slice(1, nComponents + 1);
    eigenvalues = eigenvalues.slice(1, nComponents + 1);
  } else {
    eigenvectors = eigenvectors.slice(0, nComponents);
    eigenvalues = eigenvalues.slice(0, nComponents);
  }

  return {
    embedding: eigenvectors,
    eigenvalues,
    niter: result.niter,
    success: true,
    message: 'Converged successfully',
  };
}
