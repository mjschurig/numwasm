/**
 * QUAD_RULE - Adaptive Integration with Explicit Rule Selection
 *
 * Adaptive integration using QUADPACK's DQAG routine with user-specified
 * Gauss-Kronrod rule.
 */

import { loadQUADPACKModule } from './loader.js';
import type { IntegrandFunction, QuadRuleOptions, QuadResult } from './high-level-types.js';
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
 * Compute the integral of f(x) over [a, b] with explicit Gauss-Kronrod rule.
 *
 * Uses QUADPACK's DQAG routine. This is similar to quad() but allows you
 * to specify which Gauss-Kronrod rule to use for each subinterval.
 *
 * Higher-order rules are more efficient for smooth functions:
 * - Rule 1 (GK15): 15-point rule, fastest per subdivision
 * - Rule 2 (GK21): 21-point rule
 * - Rule 3 (GK31): 31-point rule, good default
 * - Rule 4 (GK41): 41-point rule
 * - Rule 5 (GK51): 51-point rule
 * - Rule 6 (GK61): 61-point rule, most accurate per subdivision
 *
 * @param f - The integrand function f(x)
 * @param a - Lower limit of integration
 * @param b - Upper limit of integration
 * @param options - Integration options including rule selection
 * @returns Integration result with value, error estimate, and diagnostics
 *
 * @example
 * ```ts
 * import { quad_rule } from 'quadwasm';
 *
 * // Use 61-point rule for a smooth function
 * const result = await quad_rule(x => Math.exp(-x * x), -5, 5, {
 *   rule: 6  // GK61
 * });
 *
 * // Use 15-point rule for faster evaluation
 * const fast = await quad_rule(x => x * x, 0, 1, {
 *   rule: 1  // GK15
 * });
 * ```
 */
export async function quad_rule(
  f: IntegrandFunction,
  a: number,
  b: number,
  options?: QuadRuleOptions
): Promise<QuadResult> {
  const Module = await loadQUADPACKModule();

  const {
    rule = 1,
    epsabs = 1e-8,
    epsrel = 1e-8,
    limit = 50,
  } = options ?? {};

  // Validate rule
  if (rule < 1 || rule > 6) {
    throw new Error(`Invalid rule: ${rule}. Must be 1-6 (GK15, GK21, GK31, GK41, GK51, GK61).`);
  }

  const lenw = quadWorkSize(limit);
  const leniw = quadIworkSize(limit);

  // Register the integrand callback
  const fPtr = registerIntegrand(Module, f);

  // Allocate WASM memory
  const aPtr = allocateDouble(Module, a);
  const bPtr = allocateDouble(Module, b);
  const epsabsPtr = allocateDouble(Module, epsabs);
  const epsrelPtr = allocateDouble(Module, epsrel);
  const keyPtr = allocateInt(Module, rule);
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
    keyPtr,
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
    // Call DQAG
    Module._dqag_(
      fPtr,
      aPtr,
      bPtr,
      epsabsPtr,
      epsrelPtr,
      keyPtr,
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
