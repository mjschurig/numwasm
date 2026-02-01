/**
 * Singular Value Decomposition
 *
 * High-level functions for computing SVD.
 */

// Export all SVD functions
export { svd } from './svd.js';
export { svdvals } from './svdvals.js';
export { svdCompact } from './svd-compact.js';
export { svdRank } from './svd-rank.js';

// Export all types
export type {
  Matrix,
  SVDOptions,
  SVDValsOptions,
  SVDRankOptions,
  SVDResult,
  SVDValsResult,
  SVDRankResult,
} from './types.js';
