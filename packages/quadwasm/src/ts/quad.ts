/**
 * QUAD - General Adaptive Integration
 *
 * Computes definite integrals using QUADPACK's adaptive quadrature
 * with extrapolation (DQAGS).
 */

import { loadQUADPACKModule } from './loader.js';
import type { QUADPACKModule } from './types.js';
import type { IntegrandFunction, QuadOptions, QuadResult } from './high-level-types.js';
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
 * Compute the definite integral of f(x) over [a, b].
 *
 * Uses QUADPACK's DQAGS routine, which is the best general-purpose
 * adaptive integrator. It handles most integrands well, including
 * those with integrable singularities.
 *
 * @param f - The integrand function f(x)
 * @param a - Lower limit of integration
 * @param b - Upper limit of integration
 * @param options - Integration options
 * @returns Integration result with value, error estimate, and diagnostics
 *
 * @example
 * ```ts
 * import { quad } from 'quadwasm';
 *
 * // Integrate x^2 from 0 to 1 (result should be 1/3)
 * const result = await quad(x => x * x, 0, 1);
 * console.log(result.result);  // ~0.333...
 *
 * // Integrate sin(x) from 0 to π with higher precision
 * const result2 = await quad(Math.sin, 0, Math.PI, {
 *   epsabs: 1e-12,
 *   epsrel: 1e-12
 * });
 * console.log(result2.result);  // ~2.0
 *
 * // Handle integrable singularity: 1/sqrt(x) from 0 to 1
 * const result3 = await quad(x => 1 / Math.sqrt(x), 0, 1);
 * console.log(result3.result);  // ~2.0
 * ```
 */
export async function quad(
  f: IntegrandFunction,
  a: number,
  b: number,
  options?: QuadOptions
): Promise<QuadResult> {
  const Module = await loadQUADPACKModule();

  const { epsabs = 1e-8, epsrel = 1e-8, limit = 50 } = options ?? {};

  // Validate inputs
  if (!Number.isFinite(a) || !Number.isFinite(b)) {
    throw new Error(
      'For infinite intervals, use quad_inf instead of quad'
    );
  }

  return quadInternal(Module, f, a, b, epsabs, epsrel, limit);
}

/**
 * Internal implementation of quad using DQAGS.
 */
async function quadInternal(
  Module: QUADPACKModule,
  f: IntegrandFunction,
  a: number,
  b: number,
  epsabs: number,
  epsrel: number,
  limit: number
): Promise<QuadResult> {
  const lenw = quadWorkSize(limit);
  const leniw = quadIworkSize(limit);

  // Register the integrand callback
  const fPtr = registerIntegrand(Module, f);

  // Allocate WASM memory
  const aPtr = allocateDouble(Module, a);
  const bPtr = allocateDouble(Module, b);
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
    aPtr,
    bPtr,
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
    // Call DQAGS
    Module._dqags_(
      fPtr,
      aPtr,
      bPtr,
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
