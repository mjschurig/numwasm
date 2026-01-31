/**
 * QUAD_SINGULAR - Integration with Singular Weight Functions
 *
 * Computes integrals with algebraic-logarithmic singularities and
 * Cauchy principal value integrals using QUADPACK's DQAWS and DQAWC.
 */

import { loadQUADPACKModule } from './loader.js';
import type { QUADPACKModule } from './types.js';
import type {
  IntegrandFunction,
  QuadSingularOptions,
  QuadCauchyOptions,
  QuadResult,
} from './high-level-types.js';
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
 * Compute the integral of f(x) * w(x) over [a, b] with singular weight.
 *
 * Uses QUADPACK's DQAWS routine. Computes integrals of the form:
 * ∫[a,b] f(x) * (x-a)^α * (b-x)^β * v(x) dx
 *
 * where v(x) is an optional logarithmic factor:
 * - wgtfunc=1: v(x) = 1
 * - wgtfunc=2: v(x) = log(x-a)
 * - wgtfunc=3: v(x) = log(b-x)
 * - wgtfunc=4: v(x) = log(x-a) * log(b-x)
 *
 * @param f - The integrand function f(x) (without the singular weight)
 * @param a - Lower limit of integration
 * @param b - Upper limit of integration
 * @param options - Integration options (alfa and beta are required)
 * @returns Integration result with value, error estimate, and diagnostics
 *
 * @example
 * ```ts
 * import { quad_singular } from 'quadwasm';
 *
 * // Integrate f(x) * x^(-0.5) from 0 to 1 (square root singularity at 0)
 * // This is ∫ f(x) * (x-0)^(-0.5) * (1-x)^0 dx
 * const result = await quad_singular(x => Math.sin(x), 0, 1, {
 *   alfa: -0.5,  // (x-a)^alfa at left endpoint
 *   beta: 0,     // (b-x)^beta at right endpoint
 * });
 *
 * // Integrate f(x) * log(x) from 0 to 1
 * const result2 = await quad_singular(x => x, 0, 1, {
 *   alfa: 0,
 *   beta: 0,
 *   wgtfunc: 2  // log(x-a) = log(x)
 * });
 *
 * // Integrate f(x) / sqrt(x*(1-x)) from 0 to 1 (singularities at both ends)
 * const result3 = await quad_singular(x => Math.cos(x), 0, 1, {
 *   alfa: -0.5,
 *   beta: -0.5
 * });
 * ```
 */
export async function quad_singular(
  f: IntegrandFunction,
  a: number,
  b: number,
  options: QuadSingularOptions
): Promise<QuadResult> {
  const Module = await loadQUADPACKModule();

  const {
    alfa,
    beta,
    wgtfunc = 1,
    epsabs = 1e-8,
    epsrel = 1e-8,
    limit = 50,
  } = options;

  if (alfa === undefined || beta === undefined) {
    throw new Error('alfa and beta are required for singular integration');
  }

  if (alfa <= -1 || beta <= -1) {
    throw new Error('alfa and beta must be greater than -1');
  }

  return quadSingularInternal(
    Module,
    f,
    a,
    b,
    alfa,
    beta,
    wgtfunc,
    epsabs,
    epsrel,
    limit
  );
}

/**
 * Internal implementation using DQAWS.
 */
async function quadSingularInternal(
  Module: QUADPACKModule,
  f: IntegrandFunction,
  a: number,
  b: number,
  alfa: number,
  beta: number,
  integr: number,
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
  const alfaPtr = allocateDouble(Module, alfa);
  const betaPtr = allocateDouble(Module, beta);
  const integrPtr = allocateInt(Module, integr);
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
    alfaPtr,
    betaPtr,
    integrPtr,
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
    // Call DQAWS
    Module._dqaws_(
      fPtr,
      aPtr,
      bPtr,
      alfaPtr,
      betaPtr,
      integrPtr,
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

/**
 * Compute the Cauchy principal value of f(x)/(x-c) over [a, b].
 *
 * Uses QUADPACK's DQAWC routine. Computes:
 * P.V. ∫[a,b] f(x) / (x-c) dx
 *
 * where c is a point in the open interval (a, b).
 *
 * @param f - The integrand function f(x) (the numerator)
 * @param a - Lower limit of integration
 * @param b - Upper limit of integration
 * @param options - Integration options (c is required)
 * @returns Integration result with value, error estimate, and diagnostics
 *
 * @example
 * ```ts
 * import { quad_cauchy } from 'quadwasm';
 *
 * // Compute P.V. ∫[-1,1] 1/(x-0.5) dx
 * // This tests the singularity at x = 0.5
 * const result = await quad_cauchy(x => 1, -1, 1, { c: 0.5 });
 *
 * // Compute P.V. ∫[0,2] sin(x)/(x-1) dx
 * const result2 = await quad_cauchy(Math.sin, 0, 2, { c: 1 });
 * ```
 */
export async function quad_cauchy(
  f: IntegrandFunction,
  a: number,
  b: number,
  options: QuadCauchyOptions
): Promise<QuadResult> {
  const Module = await loadQUADPACKModule();

  const { c, epsabs = 1e-8, epsrel = 1e-8, limit = 50 } = options;

  if (c === undefined) {
    throw new Error('c (singular point) is required for Cauchy integration');
  }

  if (c <= a || c >= b) {
    throw new Error('Singular point c must satisfy a < c < b');
  }

  return quadCauchyInternal(Module, f, a, b, c, epsabs, epsrel, limit);
}

/**
 * Internal implementation using DQAWC.
 */
async function quadCauchyInternal(
  Module: QUADPACKModule,
  f: IntegrandFunction,
  a: number,
  b: number,
  c: number,
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
  const cPtr = allocateDouble(Module, c);
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
    cPtr,
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
    // Call DQAWC
    Module._dqawc_(
      fPtr,
      aPtr,
      bPtr,
      cPtr,
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
