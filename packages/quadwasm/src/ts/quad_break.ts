/**
 * QUAD_BREAK - Adaptive Integration with User-Specified Break Points
 *
 * Adaptive integration using QUADPACK's DQAGP routine, which subdivides
 * at user-specified points where the integrand may have difficulties.
 */

import { loadQUADPACKModule } from './loader.js';
import type { IntegrandFunction, QuadBreakOptions, QuadResult } from './high-level-types.js';
import {
  quadBreakWorkSize,
  quadBreakIworkSize,
  getQuadpackMessage,
  allocateDouble,
  allocateDoubles,
  allocateInt,
  allocateInts,
  readDouble,
  readInt,
  registerIntegrand,
  unregisterIntegrand,
  freeAll,
} from './helpers.js';

/**
 * Compute the integral of f(x) over [a, b] with user-specified break points.
 *
 * Uses QUADPACK's DQAGP routine. This is designed for integrands that have
 * known difficulties (discontinuities, peaks, rapid changes) at specific
 * points within the integration interval.
 *
 * The algorithm subdivides at the specified break points, ensuring that
 * no subinterval crosses these difficult locations.
 *
 * @param f - The integrand function f(x)
 * @param a - Lower limit of integration
 * @param b - Upper limit of integration
 * @param options - Integration options including break points
 * @returns Integration result with value, error estimate, and diagnostics
 *
 * @example
 * ```ts
 * import { quad_break } from 'quadwasm';
 *
 * // Integrate a step function with a jump at x=0.5
 * const step = (x: number) => x < 0.5 ? 0 : 1;
 * const result = await quad_break(step, 0, 1, {
 *   points: [0.5]  // Tell the integrator about the discontinuity
 * });
 *
 * // Integrate a function with known peaks
 * const peaky = (x: number) => Math.exp(-100 * (x - 0.3) ** 2) + Math.exp(-100 * (x - 0.7) ** 2);
 * const result2 = await quad_break(peaky, 0, 1, {
 *   points: [0.3, 0.7]  // Peaks at these locations
 * });
 * ```
 */
export async function quad_break(
  f: IntegrandFunction,
  a: number,
  b: number,
  options: QuadBreakOptions
): Promise<QuadResult> {
  const Module = await loadQUADPACKModule();

  const {
    points,
    epsabs = 1e-8,
    epsrel = 1e-8,
    limit = 50,
  } = options;

  // Validate points
  if (!points || points.length === 0) {
    throw new Error('At least one break point is required. Use quad() for general integration.');
  }

  // Convert points to array if Float64Array
  const pointsArray = points instanceof Float64Array ? Array.from(points) : points;

  // Validate break points are within (a, b)
  for (const p of pointsArray) {
    if (p <= a || p >= b) {
      throw new Error(`Break point ${p} must be strictly within the interval (${a}, ${b})`);
    }
  }

  const npts = pointsArray.length;
  const npts2 = npts + 2;  // QUADPACK includes a and b
  const lenw = quadBreakWorkSize(limit, npts);
  const leniw = quadBreakIworkSize(limit, npts);

  // Register the integrand callback
  const fPtr = registerIntegrand(Module, f);

  // Allocate WASM memory
  const aPtr = allocateDouble(Module, a);
  const bPtr = allocateDouble(Module, b);
  const npts2Ptr = allocateInt(Module, npts2);
  const pointsPtr = allocateDoubles(Module, pointsArray);
  const epsabsPtr = allocateDouble(Module, epsabs);
  const epsrelPtr = allocateDouble(Module, epsrel);
  const resultPtr = allocateDouble(Module, 0);
  const abserrPtr = allocateDouble(Module, 0);
  const nevalPtr = allocateInt(Module, 0);
  const ierPtr = allocateInt(Module, 0);
  const leniwPtr = allocateInt(Module, leniw);
  const lenwPtr = allocateInt(Module, lenw);
  const lastPtr = allocateInt(Module, 0);
  const iworkPtr = allocateInts(Module, null, leniw);
  const workPtr = allocateDoubles(Module, null, lenw);

  const allPtrs = [
    aPtr,
    bPtr,
    npts2Ptr,
    pointsPtr,
    epsabsPtr,
    epsrelPtr,
    resultPtr,
    abserrPtr,
    nevalPtr,
    ierPtr,
    leniwPtr,
    lenwPtr,
    lastPtr,
    iworkPtr,
    workPtr,
  ];

  try {
    // Call DQAGP
    Module._dqagp_(
      fPtr,
      aPtr,
      bPtr,
      npts2Ptr,
      pointsPtr,
      epsabsPtr,
      epsrelPtr,
      resultPtr,
      abserrPtr,
      nevalPtr,
      ierPtr,
      leniwPtr,
      lenwPtr,
      lastPtr,
      iworkPtr,
      workPtr
    );

    // Read results
    const result = readDouble(Module, resultPtr);
    const abserr = readDouble(Module, abserrPtr);
    const neval = readInt(Module, nevalPtr);
    const ier = readInt(Module, ierPtr);
    const last = readInt(Module, lastPtr);

    return {
      result,
      abserr,
      neval,
      subdivisions: last,
      ier,
      success: ier === 0,
      message: getQuadpackMessage(ier),
    };
  } finally {
    unregisterIntegrand(Module, fPtr);
    freeAll(Module, allPtrs);
  }
}
