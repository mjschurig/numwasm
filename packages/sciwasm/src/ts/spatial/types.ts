/**
 * Options for KDTree.query() method
 */
export interface QueryOptions {
  /**
   * Approximation factor. Return approximate nearest neighbors;
   * the k-th returned value is guaranteed to be no further than (1+eps) times
   * the distance to the real k-th nearest neighbor.
   * Default: 0 (exact search)
   */
  eps?: number;

  /**
   * Minkowski p-norm to use for distance calculation.
   * - 1: Manhattan distance
   * - 2: Euclidean distance (default)
   * - Infinity: Chebyshev distance
   */
  p?: number;

  /**
   * Return only neighbors within this distance.
   * Neighbors beyond this distance will have index -1 and distance Infinity.
   */
  distance_upper_bound?: number;
}

/**
 * Options for KDTree.query_ball_point() method
 */
export interface BallQueryOptions {
  /**
   * Approximation factor for branch pruning.
   * Default: 0 (exact search)
   */
  eps?: number;

  /**
   * Minkowski p-norm to use for distance calculation.
   * - 1: Manhattan distance
   * - 2: Euclidean distance (default)
   * - Infinity: Chebyshev distance
   */
  p?: number;

  /**
   * Sort the returned neighbors by distance.
   * Default: false
   */
  return_sorted?: boolean;

  /**
   * Return only the count of neighbors instead of the indices.
   * Default: false
   */
  return_length?: boolean;
}

/**
 * Result from KDTree.query() for a single query point
 */
export interface QueryResultSingle {
  /** Distances to the k nearest neighbors */
  distances: number[];
  /** Indices of the k nearest neighbors in the original dataset */
  indices: number[];
}

/**
 * Result from KDTree.query() for multiple query points
 */
export interface QueryResultBatch {
  /** Distances to the k nearest neighbors for each query point */
  distances: number[][];
  /** Indices of the k nearest neighbors for each query point */
  indices: number[][];
}

/**
 * Result from KDTree.query() - can be single or batch
 */
export type QueryResult = QueryResultSingle | QueryResultBatch;
