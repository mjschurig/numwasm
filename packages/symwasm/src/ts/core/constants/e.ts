/**
 * Euler's number constant.
 * @module core/constants/e
 */

import { Constant } from '../classes/constant.js';
import { lazyConstantProxy } from './lazy-proxy.js';

let _E: Constant | null = null;

/** Euler's number (e). */
export const E: Constant = lazyConstantProxy(() => {
  if (!_E) _E = Constant.E();
  return _E;
});
