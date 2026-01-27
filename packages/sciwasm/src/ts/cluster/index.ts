/**
 * Clustering algorithms.
 * @module cluster
 */

import { NotImplementedError } from '../errors.js';

/** K-means clustering result. */
export interface KMeansResult {
  centroids: number[][];
  labels: number[];
  inertia: number;
  n_iter: number;
}

/**
 * K-means clustering.
 * Mirrors scipy.cluster.vq.kmeans2.
 */
export function kmeans(
  _data: number[][],
  _k: number,
  _options?: { maxiter?: number; tol?: number }
): KMeansResult {
  throw new NotImplementedError('sciwasm.cluster.kmeans');
}

/** Hierarchical clustering namespace. */
export const hierarchy = {
  /**
   * Perform hierarchical/agglomerative clustering.
   * Mirrors scipy.cluster.hierarchy.linkage.
   */
  linkage(
    _y: number[] | number[][],
    _method?: 'single' | 'complete' | 'average' | 'ward'
  ): number[][] {
    throw new NotImplementedError('sciwasm.cluster.hierarchy.linkage');
  },

  /**
   * Form flat clusters from the hierarchical clustering.
   * Mirrors scipy.cluster.hierarchy.fcluster.
   */
  fcluster(
    _Z: number[][],
    _t: number,
    _criterion?: 'inconsistent' | 'distance' | 'maxclust'
  ): number[] {
    throw new NotImplementedError('sciwasm.cluster.hierarchy.fcluster');
  },
};
