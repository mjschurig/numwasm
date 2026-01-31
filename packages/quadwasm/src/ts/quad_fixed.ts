/**
 * QUAD_FIXED - Fixed-Order Gauss-Kronrod Rules
 *
 * Provides direct access to QUADPACK's DQK15-DQK61 fixed-order Gauss-Kronrod rules.
 * These perform a single quadrature rule evaluation without adaptation.
 */

import { loadQUADPACKModule } from './loader.js';
import type { QUADPACKModule } from './types.js';
import type { IntegrandFunction, QuadFixedResult } from './high-level-types.js';
import {
  allocateDouble,
  readDouble,
  registerIntegrand,
  unregisterIntegrand,
  freeAll,
} from './helpers.js';

/**
 * Gauss-Kronrod rule order type.
 */
export type GKRule = 15 | 21 | 31 | 41 | 51 | 61;

/**
 * Evaluate the integral of f(x) over [a, b] using a fixed Gauss-Kronrod rule.
 *
 * This function applies a single Gauss-Kronrod quadrature rule without
 * any adaptive subdivision. It's useful for:
 * - Building custom integration strategies
 * - Quick estimates when you know the function is smooth
 * - Performance-critical code where overhead must be minimized
 *
 * Available rules:
 * - 15-point: 7 Gauss + 8 Kronrod points (fastest)
 * - 21-point: 10 Gauss + 11 Kronrod points
 * - 31-point: 15 Gauss + 16 Kronrod points
 * - 41-point: 20 Gauss + 21 Kronrod points
 * - 51-point: 25 Gauss + 26 Kronrod points
 * - 61-point: 30 Gauss + 31 Kronrod points (most accurate)
 *
 * @param f - The integrand function f(x)
 * @param a - Lower limit of integration
 * @param b - Upper limit of integration
 * @param rule - Gauss-Kronrod rule order (15, 21, 31, 41, 51, or 61)
 * @returns Result with integral, error estimate, and diagnostic values
 *
 * @example
 * ```ts
 * import { quad_fixed } from 'quadwasm';
 *
 * // Use 15-point rule for quick estimate
 * const fast = await quad_fixed(x => x * x, 0, 1, 15);
 * console.log(fast.result);  // ~0.333...
 *
 * // Use 61-point rule for higher accuracy
 * const accurate = await quad_fixed(x => x * x, 0, 1, 61);
 * console.log(accurate.result);  // More precise 0.333...
 * console.log(accurate.abserr);  // Much smaller error estimate
 * ```
 */
export async function quad_fixed(
  f: IntegrandFunction,
  a: number,
  b: number,
  rule: GKRule
): Promise<QuadFixedResult> {
  const Module = await loadQUADPACKModule();

  // Register the integrand callback
  const fPtr = registerIntegrand(Module, f);

  // Allocate WASM memory
  const aPtr = allocateDouble(Module, a);
  const bPtr = allocateDouble(Module, b);
  const resultPtr = allocateDouble(Module, 0);
  const abserrPtr = allocateDouble(Module, 0);
  const resabsPtr = allocateDouble(Module, 0);
  const resascPtr = allocateDouble(Module, 0);

  const allPtrs = [aPtr, bPtr, resultPtr, abserrPtr, resabsPtr, resascPtr];

  try {
    // Call the appropriate DQK routine
    callDqk(Module, fPtr, aPtr, bPtr, resultPtr, abserrPtr, resabsPtr, resascPtr, rule);

    // Read results
    const result = readDouble(Module, resultPtr);
    const abserr = readDouble(Module, abserrPtr);
    const resabs = readDouble(Module, resabsPtr);
    const resasc = readDouble(Module, resascPtr);

    return {
      result,
      abserr,
      resabs,
      resasc,
    };
  } finally {
    unregisterIntegrand(Module, fPtr);
    freeAll(Module, allPtrs);
  }
}

/**
 * Call the appropriate DQK routine based on rule order.
 */
function callDqk(
  Module: QUADPACKModule,
  fPtr: number,
  aPtr: number,
  bPtr: number,
  resultPtr: number,
  abserrPtr: number,
  resabsPtr: number,
  resascPtr: number,
  rule: GKRule
): void {
  switch (rule) {
    case 15:
      Module._dqk15_(fPtr, aPtr, bPtr, resultPtr, abserrPtr, resabsPtr, resascPtr);
      break;
    case 21:
      Module._dqk21_(fPtr, aPtr, bPtr, resultPtr, abserrPtr, resabsPtr, resascPtr);
      break;
    case 31:
      Module._dqk31_(fPtr, aPtr, bPtr, resultPtr, abserrPtr, resabsPtr, resascPtr);
      break;
    case 41:
      Module._dqk41_(fPtr, aPtr, bPtr, resultPtr, abserrPtr, resabsPtr, resascPtr);
      break;
    case 51:
      Module._dqk51_(fPtr, aPtr, bPtr, resultPtr, abserrPtr, resabsPtr, resascPtr);
      break;
    case 61:
      Module._dqk61_(fPtr, aPtr, bPtr, resultPtr, abserrPtr, resabsPtr, resascPtr);
      break;
    default:
      throw new Error(`Invalid Gauss-Kronrod rule: ${rule}. Must be 15, 21, 31, 41, 51, or 61.`);
  }
}
