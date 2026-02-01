/**
 * Shared utilities for XSF function modules
 */

import { getXSFModule, isXSFLoaded } from './loader.js';
import type { XSFModule } from './types.js';

export type ArrayInput = number[] | Float64Array;

export function ensureLoaded(): XSFModule {
  if (!isXSFLoaded()) {
    throw new Error(
      'XSF WASM module not loaded. Call "await loadXSFModule()" before using special functions.'
    );
  }
  return getXSFModule();
}

export function toFloat64Array(x: ArrayInput): Float64Array {
  return x instanceof Float64Array ? x : new Float64Array(x);
}
