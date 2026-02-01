/**
 * Golden ratio constant.
 * @module core/constants/golden-ratio
 */

import { Constant } from '../classes/constant.js';
import { lazyConstantProxy } from './lazy-proxy.js';

let _goldenRatio: Constant | null = null;

/** Golden ratio (φ ≈ 1.618). */
export const GoldenRatio: Constant = lazyConstantProxy(() => {
  if (!_goldenRatio) _goldenRatio = Constant.GoldenRatio();
  return _goldenRatio;
});
