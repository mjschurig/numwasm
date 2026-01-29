/**
 * Expression simplification functions.
 * @module simplify
 */

// Re-export implemented functions from core (Phase 2.4)
export {
  expand,
  simplify,
  numer,
  denom,
  trigsimp,
  radsimp,
  powsimp,
  rewrite_as_exp,
  rewrite_as_sin,
  rewrite_as_cos,
  as_real_imag,
  expand_trig,
  expand_complex,
} from '../core/index.js';
