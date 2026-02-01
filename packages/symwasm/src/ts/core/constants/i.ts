/**
 * Imaginary unit constant.
 * @module core/constants/i
 */

import { ImaginaryUnit } from '../classes/imaginary-unit.js';
import { lazyConstantProxy } from './lazy-proxy.js';

let _I: ImaginaryUnit | null = null;

/** Imaginary unit (i). */
export const I: ImaginaryUnit = lazyConstantProxy(() => {
  if (!_I) _I = ImaginaryUnit.create();
  return _I;
});
