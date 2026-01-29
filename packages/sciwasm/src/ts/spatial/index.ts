/**
 * Spatial data structures and algorithms
 *
 * This module provides efficient spatial data structures for nearest neighbor
 * searches and range queries in k-dimensional space.
 */

export { KDTree } from './kdtree.js';
export type {
  QueryOptions,
  BallQueryOptions,
  QueryResult,
  QueryResultSingle,
  QueryResultBatch
} from './types.js';
