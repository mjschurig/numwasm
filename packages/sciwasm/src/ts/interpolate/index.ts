/**
 * Interpolation functions and spline classes.
 * @module interpolate
 */

import { NotImplementedError } from '../errors.js';

/** Interpolation function that can be called with new x values. */
export interface Interpolator {
  /** Evaluate the interpolant at the given points. */
  evaluate(x: number | number[]): number | number[];
}

/**
 * Interpolate a 1-D function.
 * Mirrors scipy.interpolate.interp1d.
 */
export function interp1d(
  _x: number[],
  _y: number[],
  _options?: { kind?: 'linear' | 'nearest' | 'cubic' | 'quadratic'; fill_value?: number | 'extrapolate' }
): Interpolator {
  throw new NotImplementedError('sciwasm.interpolate.interp1d');
}

/**
 * Cubic spline interpolation.
 * Mirrors scipy.interpolate.CubicSpline.
 */
export class CubicSpline implements Interpolator {
  constructor(
    _x: number[],
    _y: number[],
    _options?: { bc_type?: 'not-a-knot' | 'clamped' | 'natural' }
  ) {
    throw new NotImplementedError('sciwasm.interpolate.CubicSpline');
  }

  evaluate(_x: number | number[]): number | number[] {
    throw new NotImplementedError('sciwasm.interpolate.CubicSpline.evaluate');
  }
}

/**
 * PCHIP 1-D monotonic cubic interpolation.
 * Mirrors scipy.interpolate.PchipInterpolator.
 */
export class PchipInterpolator implements Interpolator {
  constructor(_x: number[], _y: number[]) {
    throw new NotImplementedError('sciwasm.interpolate.PchipInterpolator');
  }

  evaluate(_x: number | number[]): number | number[] {
    throw new NotImplementedError('sciwasm.interpolate.PchipInterpolator.evaluate');
  }
}

/**
 * Interpolate unstructured N-D data.
 * Mirrors scipy.interpolate.griddata.
 */
export function griddata(
  _points: number[][],
  _values: number[],
  _xi: number[][],
  _method?: 'nearest' | 'linear' | 'cubic'
): number[] {
  throw new NotImplementedError('sciwasm.interpolate.griddata');
}
