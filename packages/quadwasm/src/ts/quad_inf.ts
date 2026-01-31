/**
 * QUAD_INF - Infinite Interval Integration
 *
 * Computes integrals over infinite or semi-infinite intervals using
 * QUADPACK's DQAGI routine.
 */

import { loadQUADPACKModule } from './loader.js';
import type { QUADPACKModule } from './types.js';
import { QUADPACKInf } from './types.js';
import type { IntegrandFunction, QuadInfOptions, QuadResult } from './high-level-types.js';
import {
  quadWorkSize,
  quadIworkSize,
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
 * Compute the integral of f(x) over an infinite or semi-infinite interval.
 *
 * Uses QUADPACK's DQAGI routine. The interval type is determined
 * automatically from the bounds:
 * - [a, Infinity]: Semi-infinite upper (∫[a,∞) f(x) dx)
 * - [-Infinity, b]: Semi-infinite lower (∫(-∞,b] f(x) dx)
 * - [-Infinity, Infinity]: Doubly infinite (∫(-∞,∞) f(x) dx)
 *
 * @param f - The integrand function f(x)
 * @param a - Lower limit (use -Infinity for lower infinite)
 * @param b - Upper limit (use Infinity for upper infinite)
 * @param options - Integration options
 * @returns Integration result with value, error estimate, and diagnostics
 *
 * @example
 * ```ts
 * import { quad_inf } from 'quadwasm';
 *
 * // Integrate exp(-x^2) from 0 to infinity (result is √π/2 ≈ 0.886)
 * const result = await quad_inf(x => Math.exp(-x * x), 0, Infinity);
 * console.log(result.result);
 *
 * // Integrate 1/(1+x^2) from -infinity to +infinity (result is π)
 * const result2 = await quad_inf(x => 1 / (1 + x * x), -Infinity, Infinity);
 * console.log(result2.result);
 *
 * // Integrate exp(x) from -infinity to 0 (result is 1)
 * const result3 = await quad_inf(Math.exp, -Infinity, 0);
 * console.log(result3.result);
 * ```
 */
export async function quad_inf(
  f: IntegrandFunction,
  a: number,
  b: number,
  options?: QuadInfOptions
): Promise<QuadResult> {
  const Module = await loadQUADPACKModule();

  const { epsabs = 1e-8, epsrel = 1e-8, limit = 50 } = options ?? {};

  // Determine interval type and bound
  let inf: number;
  let bound: number;

  const aIsInf = !Number.isFinite(a);
  const bIsInf = !Number.isFinite(b);

  if (!aIsInf && bIsInf) {
    // [a, +∞)
    inf = QUADPACKInf.UPPER;
    bound = a;
  } else if (aIsInf && !bIsInf) {
    // (-∞, b]
    inf = QUADPACKInf.LOWER;
    bound = b;
  } else if (aIsInf && bIsInf) {
    // (-∞, +∞)
    inf = QUADPACKInf.BOTH;
    bound = 0; // Not used, but must be provided
  } else {
    throw new Error(
      'For finite intervals, use quad instead of quad_inf'
    );
  }

  return quadInfInternal(Module, f, bound, inf, epsabs, epsrel, limit);
}

/**
 * Internal implementation using DQAGI.
 */
async function quadInfInternal(
  Module: QUADPACKModule,
  f: IntegrandFunction,
  bound: number,
  inf: number,
  epsabs: number,
  epsrel: number,
  limit: number
): Promise<QuadResult> {
  const lenw = quadWorkSize(limit);
  const leniw = quadIworkSize(limit);

  // Register the integrand callback
  const fPtr = registerIntegrand(Module, f);

  // Allocate WASM memory
  const boundPtr = allocateDouble(Module, bound);
  const infPtr = allocateInt(Module, inf);
  const epsabsPtr = allocateDouble(Module, epsabs);
  const epsrelPtr = allocateDouble(Module, epsrel);
  const resultPtr = allocateDouble(Module, 0);
  const abserrPtr = allocateDouble(Module, 0);
  const nevalPtr = allocateInt(Module, 0);
  const ierPtr = allocateInt(Module, 0);
  const limitPtr = allocateInt(Module, limit);
  const lenwPtr = allocateInt(Module, lenw);
  const lastPtr = allocateInt(Module, 0);
  const iworkPtr = allocateInts(Module, null, leniw);
  const workPtr = allocateDoubles(Module, null, lenw);

  const allPtrs = [
    boundPtr,
    infPtr,
    epsabsPtr,
    epsrelPtr,
    resultPtr,
    abserrPtr,
    nevalPtr,
    ierPtr,
    limitPtr,
    lenwPtr,
    lastPtr,
    iworkPtr,
    workPtr,
  ];

  try {
    // Call DQAGI
    Module._dqagi_(
      fPtr,
      boundPtr,
      infPtr,
      epsabsPtr,
      epsrelPtr,
      resultPtr,
      abserrPtr,
      nevalPtr,
      ierPtr,
      limitPtr,
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
