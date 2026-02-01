/**
 * Least Squares & Minimum Norm
 *
 * High-level functions for solving least squares problems.
 */

// Export all least squares functions
export { lstsq } from './lstsq.js';
export { lstsqSVD } from './lstsq-svd.js';
export { lstsqGelsy } from './lstsq-gelsy.js';
export { constrainedLstSq } from './constrained-lstsq.js';
export { generalizedLstSq } from './generalized-lstsq.js';

// Export all types
export type {
  Matrix,
  Vector,
  LstSqOptions,
  LstSqResult,
  LstSqSVDOptions,
  LstSqSVDResult,
  LstSqGelsyOptions,
  LstSqGelsyResult,
  ConstrainedLstSqOptions,
  ConstrainedLstSqResult,
  GeneralizedLstSqOptions,
  GeneralizedLstSqResult,
} from './types.js';
