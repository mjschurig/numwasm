/**
 * Matrix Norms, Condition Numbers, Determinants & Rank
 *
 * This module provides functions for computing matrix norms, condition numbers,
 * determinants, and numerical rank.
 *
 * @packageDocumentation
 */

// Matrix Norms
export { norm } from './norm.js';

// Condition Numbers
export { cond } from './cond.js';
export { condEst } from './cond-est.js';
export { rcond } from './rcond.js';

// Determinants
export { det } from './det.js';
export { logdet } from './logdet.js';
export { slogdet } from './slogdet.js';

// Matrix Rank
export { rank } from './rank.js';

// Types
export type {
  Matrix,
  NormType,
  NormOptions,
  NormResult,
  CondOptions,
  CondResult,
  CondEstOptions,
  CondEstResult,
  RcondOptions,
  RcondResult,
  DetOptions,
  DetResult,
  LogDetOptions,
  LogDetResult,
  SlogDetOptions,
  SlogDetResult,
  RankOptions,
  RankResult,
} from './types.js';
