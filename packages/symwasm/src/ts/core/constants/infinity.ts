/**
 * Infinity constant.
 * @module core/constants/infinity
 */

import { Infinity_ } from '../classes/infinity.js';
import { lazyConstantProxy } from './lazy-proxy.js';

let _infinity: Infinity_ | null = null;

/** Positive infinity. */
export const oo: Infinity_ = lazyConstantProxy(() => {
  if (!_infinity) _infinity = Infinity_.positive();
  return _infinity;
});
