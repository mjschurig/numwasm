/**
 * Internal utilities for the stats module.
 */

import { NDArray } from 'numwasm';

/**
 * Convert a 0-d NDArray to a plain number, otherwise return as-is.
 * Disposes the NDArray if it was converted.
 */
export function scalarize(value: NDArray | number | boolean): NDArray | number {
  if (value instanceof NDArray && value.ndim === 0) {
    const v = value.item();
    value.dispose();
    return v;
  }
  return value as NDArray | number;
}
