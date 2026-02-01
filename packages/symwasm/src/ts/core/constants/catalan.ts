/**
 * Catalan's constant.
 * @module core/constants/catalan
 */

import { Constant } from '../classes/constant.js';
import { lazyConstantProxy } from './lazy-proxy.js';

let _catalan: Constant | null = null;

/** Catalan's constant (≈ 0.9159). */
export const Catalan: Constant = lazyConstantProxy(() => {
  if (!_catalan) _catalan = Constant.Catalan();
  return _catalan;
});
