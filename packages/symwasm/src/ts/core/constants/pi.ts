/**
 * Pi constant.
 * @module core/constants/pi
 */

import { Constant } from '../classes/constant.js';
import { lazyConstantProxy } from './lazy-proxy.js';

let _pi: Constant | null = null;

/** Pi (π). */
export const pi: Constant = lazyConstantProxy(() => {
  if (!_pi) _pi = Constant.pi();
  return _pi;
});
