/**
 * Euler-Mascheroni constant.
 * @module core/constants/euler-gamma
 */

import { Constant } from '../classes/constant.js';
import { lazyConstantProxy } from './lazy-proxy.js';

let _eulerGamma: Constant | null = null;

/** Euler-Mascheroni constant (γ ≈ 0.5772). */
export const EulerGamma: Constant = lazyConstantProxy(() => {
  if (!_eulerGamma) _eulerGamma = Constant.EulerGamma();
  return _eulerGamma;
});
