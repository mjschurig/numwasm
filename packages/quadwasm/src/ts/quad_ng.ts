/**
 * QUAD_NG - Non-Adaptive Integration
 *
 * Simple non-adaptive Gauss-Kronrod integration using QUADPACK's DQNG routine.
 * Best for smooth, well-behaved functions where adaptive subdivision is not needed.
 */

import { loadQUADPACKModule } from './loader.js';
import type { IntegrandFunction, QuadNgOptions, QuadResult } from './high-level-types.js';
import {
  getQuadpackMessage,
  allocateDouble,
  allocateInt,
  readDouble,
  readInt,
  registerIntegrand,
  unregisterIntegrand,
  freeAll,
} from './helpers.js';

/**
 * Compute the integral of f(x) over [a, b] using non-adaptive integration.
 *
 * Uses QUADPACK's DQNG routine, which applies progressively higher-order
 * Gauss-Kronrod rules (10, 21, 43, 87 points) until the requested accuracy
 * is achieved or all rules are exhausted.
 *
 * This is the fastest integration routine but only suitable for smooth,
 * well-behaved functions without singularities or rapid oscillations.
 *
 * @param f - The integrand function f(x)
 * @param a - Lower limit of integration
 * @param b - Upper limit of integration
 * @param options - Integration options
 * @returns Integration result with value, error estimate, and diagnostics
 *
 * @example
 * ```ts
 * import { quad_ng } from 'quadwasm';
 *
 * // Integrate a smooth polynomial
 * const result = await quad_ng(x => x * x, 0, 1);
 * console.log(result.result);  // ~0.333...
 *
 * // Integrate with tighter tolerance
 * const result2 = await quad_ng(x => Math.sin(x), 0, Math.PI, {
 *   epsabs: 1e-12,
 *   epsrel: 1e-12
 * });
 * console.log(result2.result);  // ~2.0
 * ```
 */
export async function quad_ng(
  f: IntegrandFunction,
  a: number,
  b: number,
  options?: QuadNgOptions
): Promise<QuadResult> {
  const Module = await loadQUADPACKModule();

  const { epsabs = 1e-8, epsrel = 1e-8 } = options ?? {};

  // Register the integrand callback
  const fPtr = registerIntegrand(Module, f);

  // Allocate WASM memory - DQNG is simple, no work arrays needed
  const aPtr = allocateDouble(Module, a);
  const bPtr = allocateDouble(Module, b);
  const epsabsPtr = allocateDouble(Module, epsabs);
  const epsrelPtr = allocateDouble(Module, epsrel);
  const resultPtr = allocateDouble(Module, 0);
  const abserrPtr = allocateDouble(Module, 0);
  const nevalPtr = allocateInt(Module, 0);
  const ierPtr = allocateInt(Module, 0);

  const allPtrs = [
    aPtr,
    bPtr,
    epsabsPtr,
    epsrelPtr,
    resultPtr,
    abserrPtr,
    nevalPtr,
    ierPtr,
  ];

  try {
    // Call DQNG
    Module._dqng_(
      fPtr,
      aPtr,
      bPtr,
      epsabsPtr,
      epsrelPtr,
      resultPtr,
      abserrPtr,
      nevalPtr,
      ierPtr
    );

    // Read results
    const result = readDouble(Module, resultPtr);
    const abserr = readDouble(Module, abserrPtr);
    const neval = readInt(Module, nevalPtr);
    const ier = readInt(Module, ierPtr);

    return {
      result,
      abserr,
      neval,
      subdivisions: 0, // Non-adaptive, no subdivisions
      ier,
      success: ier === 0,
      message: getQuadpackMessage(ier),
    };
  } finally {
    unregisterIntegrand(Module, fPtr);
    freeAll(Module, allPtrs);
  }
}
