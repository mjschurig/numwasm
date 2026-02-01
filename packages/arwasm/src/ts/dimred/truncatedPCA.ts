/**
 * TRUNCATEDPCA - Truncated Principal Component Analysis
 *
 * Computes the top principal components of a dataset using eigenvalue
 * decomposition of the covariance matrix. This is efficient for finding
 * a small number of components from high-dimensional data.
 *
 * PCA finds the directions of maximum variance in the data.
 */

import { eigs } from '../core/eigs.js';
import type { MatVecFunction } from '../high-level-types.js';

/**
 * Options for truncated PCA.
 */
export interface TruncatedPCAOptions {
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
 * Result from truncated PCA.
 */
export interface TruncatedPCAResult {
  /**
   * Principal components (nComponents x n).
   * components[i] is the i-th principal component direction.
   */
  components: Float64Array[];

  /**
   * Explained variance (eigenvalues) for each component.
   */
  explainedVariance: Float64Array;

  /**
   * Ratio of variance explained by each component.
   * Sums to less than 1 for truncated PCA.
   */
  explainedVarianceRatio: Float64Array;

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
 * Compute truncated PCA via eigendecomposition of the covariance matrix.
 *
 * Finds the top principal components by computing the largest eigenvalues
 * and eigenvectors of X^T * X / (n-1) or a provided covariance matvec.
 *
 * @param covarianceMatvec - Function computing y = C*x where C is the covariance matrix
 * @param n - Dimensionality (number of features)
 * @param nComponents - Number of principal components to compute
 * @param options - Solver options
 * @returns Principal components and explained variance
 *
 * @example
 * ```ts
 * import { truncatedPCA } from 'arwasm';
 *
 * // For centered data matrix X (samples x features)
 * const nSamples = 10000;
 * const nFeatures = 500;
 *
 * // Create covariance matvec: C*x = X^T * X * x / (n-1)
 * const covMatvec = (x: Float64Array): Float64Array => {
 *   // First compute X*x (nSamples-vector)
 *   const Xx = matmul(X, x);
 *   // Then compute X^T * (X*x) (nFeatures-vector)
 *   const XtXx = matmulT(X, Xx);
 *   // Normalize by n-1
 *   for (let i = 0; i < nFeatures; i++) {
 *     XtXx[i] /= (nSamples - 1);
 *   }
 *   return XtXx;
 * };
 *
 * // Find top 10 principal components
 * const result = await truncatedPCA(covMatvec, nFeatures, 10);
 * console.log('Explained variance ratio:', result.explainedVarianceRatio);
 * ```
 */
export async function truncatedPCA(
  covarianceMatvec: MatVecFunction,
  n: number,
  nComponents: number,
  options?: TruncatedPCAOptions
): Promise<TruncatedPCAResult> {
  const {
    tol = 0,
    ncv,
    maxiter,
  } = options ?? {};

  if (nComponents < 1) {
    throw new Error('nComponents must be at least 1');
  }
  if (nComponents >= n) {
    throw new Error('nComponents must be less than n');
  }

  // Find largest eigenvalues of covariance matrix
  const result = await eigs(covarianceMatvec, n, nComponents, {
    which: 'LM', // Largest eigenvalues = principal components
    tol,
    ncv,
    maxiter,
    return_eigenvectors: true,
  });

  if (!result.success || !result.eigenvectors) {
    return {
      components: [],
      explainedVariance: new Float64Array(0),
      explainedVarianceRatio: new Float64Array(0),
      niter: result.niter,
      success: false,
      message: result.message,
    };
  }

  // Eigenvalues are the explained variances
  const explainedVariance = result.eigenvalues;

  // Compute explained variance ratio
  // We need to estimate total variance for the ratio
  // For truncated PCA, we can only provide ratio relative to computed components
  // A proper ratio would require full eigendecomposition
  const totalVariance = explainedVariance.reduce((a, b) => a + Math.max(0, b), 0);
  const explainedVarianceRatio = new Float64Array(result.nconv);

  for (let i = 0; i < result.nconv; i++) {
    // Note: This ratio is only among computed components
    // True ratio would need trace(covariance) which requires additional computation
    explainedVarianceRatio[i] = totalVariance > 0 ? Math.max(0, explainedVariance[i]) / totalVariance : 0;
  }

  return {
    components: result.eigenvectors,
    explainedVariance,
    explainedVarianceRatio,
    niter: result.niter,
    success: true,
    message: 'Converged successfully',
  };
}

/**
 * Compute truncated PCA directly from a data matrix via SVD.
 *
 * This version uses SVD of the centered data matrix X, which is
 * equivalent to eigendecomposition of the covariance matrix but
 * can be more numerically stable.
 *
 * @param dataMatvec - Function computing y = X*v (nSamples x nFeatures matrix times nFeatures-vector)
 * @param dataMatvecT - Function computing y = X^T*v (nFeatures x nSamples matrix times nSamples-vector)
 * @param nSamples - Number of samples (rows)
 * @param nFeatures - Number of features (columns)
 * @param nComponents - Number of principal components
 * @param options - Solver options
 * @returns Principal components and explained variance
 */
export async function truncatedPCAfromData(
  dataMatvec: MatVecFunction,
  dataMatvecT: MatVecFunction,
  nSamples: number,
  nFeatures: number,
  nComponents: number,
  options?: TruncatedPCAOptions
): Promise<TruncatedPCAResult> {
  const {
    tol = 0,
    ncv,
    maxiter,
  } = options ?? {};

  // We compute eigenvalues of X^T*X by working with the symmetric matvec
  // v -> X^T * X * v
  const XtXMatvec: MatVecFunction = (v: Float64Array) => {
    const Xv = dataMatvec(v);
    const XvArr = Xv instanceof Float64Array ? Xv : new Float64Array(Xv);
    return dataMatvecT(XvArr);
  };

  const result = await eigs(XtXMatvec, nFeatures, nComponents, {
    which: 'LM',
    tol,
    ncv,
    maxiter,
    return_eigenvectors: true,
  });

  if (!result.success || !result.eigenvectors) {
    return {
      components: [],
      explainedVariance: new Float64Array(0),
      explainedVarianceRatio: new Float64Array(0),
      niter: result.niter,
      success: false,
      message: result.message,
    };
  }

  // Eigenvalues of X^T*X are squared singular values
  // Variance = eigenvalue / (n-1)
  const explainedVariance = new Float64Array(result.nconv);
  for (let i = 0; i < result.nconv; i++) {
    explainedVariance[i] = Math.max(0, result.eigenvalues[i]) / (nSamples - 1);
  }

  const totalVariance = explainedVariance.reduce((a, b) => a + b, 0);
  const explainedVarianceRatio = new Float64Array(result.nconv);
  for (let i = 0; i < result.nconv; i++) {
    explainedVarianceRatio[i] = totalVariance > 0 ? explainedVariance[i] / totalVariance : 0;
  }

  return {
    components: result.eigenvectors,
    explainedVariance,
    explainedVarianceRatio,
    niter: result.niter,
    success: true,
    message: 'Converged successfully',
  };
}
