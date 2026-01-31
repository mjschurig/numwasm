/**
 * QUAD_OSC - Oscillatory Integration
 *
 * Computes integrals of oscillatory functions using QUADPACK's
 * DQAWO and DQAWF routines.
 */

import { loadQUADPACKModule } from './loader.js';
import type { QUADPACKModule } from './types.js';
import { QUADPACKOscillatory } from './types.js';
import type {
  IntegrandFunction,
  QuadOscOptions,
  QuadFourierOptions,
  QuadResult,
} from './high-level-types.js';
import {
  quadOscWorkSize,
  quadOscIworkSize,
  quadFourierWorkSize,
  quadFourierIworkSize,
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
 * Compute the integral of f(x) * w(x) over [a, b] where w(x) is oscillatory.
 *
 * Uses QUADPACK's DQAWO routine. Computes integrals of the form:
 * - ∫[a,b] f(x) * cos(ω*x) dx  (when weight='cos')
 * - ∫[a,b] f(x) * sin(ω*x) dx  (when weight='sin')
 *
 * @param f - The integrand function f(x) (without the oscillatory weight)
 * @param a - Lower limit of integration
 * @param b - Upper limit of integration
 * @param options - Integration options (omega is required)
 * @returns Integration result with value, error estimate, and diagnostics
 *
 * @example
 * ```ts
 * import { quad_osc } from 'quadwasm';
 *
 * // Integrate x * cos(10*x) from 0 to 2π
 * const result = await quad_osc(x => x, 0, 2 * Math.PI, {
 *   omega: 10,
 *   weight: 'cos'
 * });
 *
 * // Integrate exp(-x) * sin(50*x) from 0 to 1
 * const result2 = await quad_osc(x => Math.exp(-x), 0, 1, {
 *   omega: 50,
 *   weight: 'sin'
 * });
 * ```
 */
export async function quad_osc(
  f: IntegrandFunction,
  a: number,
  b: number,
  options: QuadOscOptions
): Promise<QuadResult> {
  const Module = await loadQUADPACKModule();

  const {
    omega,
    weight = 'cos',
    epsabs = 1e-8,
    epsrel = 1e-8,
    limit = 50,
    maxp1 = 21,
  } = options;

  if (omega === undefined) {
    throw new Error('omega is required for oscillatory integration');
  }

  // Validate inputs
  if (!Number.isFinite(a) || !Number.isFinite(b)) {
    throw new Error(
      'For oscillatory integration over infinite intervals, use quad_fourier'
    );
  }

  const integr = weight === 'cos' ? QUADPACKOscillatory.COS : QUADPACKOscillatory.SIN;

  return quadOscInternal(Module, f, a, b, omega, integr, epsabs, epsrel, limit, maxp1);
}

/**
 * Internal implementation using DQAWO.
 */
async function quadOscInternal(
  Module: QUADPACKModule,
  f: IntegrandFunction,
  a: number,
  b: number,
  omega: number,
  integr: number,
  epsabs: number,
  epsrel: number,
  limit: number,
  maxp1: number
): Promise<QuadResult> {
  const lenw = quadOscWorkSize(limit, maxp1);
  const leniw = quadOscIworkSize(limit);

  // Register the integrand callback
  const fPtr = registerIntegrand(Module, f);

  // Allocate WASM memory
  const aPtr = allocateDouble(Module, a);
  const bPtr = allocateDouble(Module, b);
  const omegaPtr = allocateDouble(Module, omega);
  const integrPtr = allocateInt(Module, integr);
  const epsabsPtr = allocateDouble(Module, epsabs);
  const epsrelPtr = allocateDouble(Module, epsrel);
  const resultPtr = allocateDouble(Module, 0);
  const abserrPtr = allocateDouble(Module, 0);
  const nevalPtr = allocateInt(Module, 0);
  const ierPtr = allocateInt(Module, 0);
  const leniwPtr = allocateInt(Module, leniw);
  const maxp1Ptr = allocateInt(Module, maxp1);
  const lenwPtr = allocateInt(Module, lenw);
  const lastPtr = allocateInt(Module, 0);
  const iworkPtr = allocateInts(Module, null, leniw);
  const workPtr = allocateDoubles(Module, null, lenw);

  const allPtrs = [
    aPtr,
    bPtr,
    omegaPtr,
    integrPtr,
    epsabsPtr,
    epsrelPtr,
    resultPtr,
    abserrPtr,
    nevalPtr,
    ierPtr,
    leniwPtr,
    maxp1Ptr,
    lenwPtr,
    lastPtr,
    iworkPtr,
    workPtr,
  ];

  try {
    // Call DQAWO
    Module._dqawo_(
      fPtr,
      aPtr,
      bPtr,
      omegaPtr,
      integrPtr,
      epsabsPtr,
      epsrelPtr,
      resultPtr,
      abserrPtr,
      nevalPtr,
      ierPtr,
      leniwPtr,
      maxp1Ptr,
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
 * Compute the Fourier integral of f(x) * w(x) over [a, ∞).
 *
 * Uses QUADPACK's DQAWF routine. Computes integrals of the form:
 * - ∫[a,∞) f(x) * cos(ω*x) dx  (when weight='cos')
 * - ∫[a,∞) f(x) * sin(ω*x) dx  (when weight='sin')
 *
 * The function f(x) should decay sufficiently fast as x → ∞.
 *
 * @param f - The integrand function f(x) (without the oscillatory weight)
 * @param a - Lower limit of integration
 * @param options - Integration options (omega is required)
 * @returns Integration result with value, error estimate, and diagnostics
 *
 * @example
 * ```ts
 * import { quad_fourier } from 'quadwasm';
 *
 * // Fourier cosine transform of exp(-x) from 0 to infinity
 * // Result: 1/(1+ω²)
 * const result = await quad_fourier(x => Math.exp(-x), 0, {
 *   omega: 1,
 *   weight: 'cos'
 * });
 * console.log(result.result);  // ~0.5
 *
 * // Fourier sine transform of exp(-x) from 0 to infinity
 * // Result: ω/(1+ω²)
 * const result2 = await quad_fourier(x => Math.exp(-x), 0, {
 *   omega: 1,
 *   weight: 'sin'
 * });
 * console.log(result2.result);  // ~0.5
 * ```
 */
export async function quad_fourier(
  f: IntegrandFunction,
  a: number,
  options: QuadFourierOptions
): Promise<QuadResult> {
  const Module = await loadQUADPACKModule();

  const {
    omega,
    weight = 'cos',
    epsabs = 1e-8,
    limit = 50,
    limlst = 50,
    maxp1 = 21,
  } = options;

  if (omega === undefined) {
    throw new Error('omega is required for Fourier integration');
  }

  const integr = weight === 'cos' ? QUADPACKOscillatory.COS : QUADPACKOscillatory.SIN;

  return quadFourierInternal(Module, f, a, omega, integr, epsabs, limlst, limit, maxp1);
}

/**
 * Internal implementation using DQAWF.
 */
async function quadFourierInternal(
  Module: QUADPACKModule,
  f: IntegrandFunction,
  a: number,
  omega: number,
  integr: number,
  epsabs: number,
  limlst: number,
  limit: number,
  maxp1: number
): Promise<QuadResult> {
  const lenw = quadFourierWorkSize(limlst, limit, maxp1);
  const leniw = quadFourierIworkSize(limlst, limit);

  // Register the integrand callback
  const fPtr = registerIntegrand(Module, f);

  // Allocate WASM memory
  const aPtr = allocateDouble(Module, a);
  const omegaPtr = allocateDouble(Module, omega);
  const integrPtr = allocateInt(Module, integr);
  const epsabsPtr = allocateDouble(Module, epsabs);
  const resultPtr = allocateDouble(Module, 0);
  const abserrPtr = allocateDouble(Module, 0);
  const nevalPtr = allocateInt(Module, 0);
  const ierPtr = allocateInt(Module, 0);
  const limlstPtr = allocateInt(Module, limlst);
  const lstPtr = allocateInt(Module, 0);
  const leniwPtr = allocateInt(Module, leniw);
  const maxp1Ptr = allocateInt(Module, maxp1);
  const lenwPtr = allocateInt(Module, lenw);
  const iworkPtr = allocateInts(Module, null, leniw);
  const workPtr = allocateDoubles(Module, null, lenw);

  const allPtrs = [
    aPtr,
    omegaPtr,
    integrPtr,
    epsabsPtr,
    resultPtr,
    abserrPtr,
    nevalPtr,
    ierPtr,
    limlstPtr,
    lstPtr,
    leniwPtr,
    maxp1Ptr,
    lenwPtr,
    iworkPtr,
    workPtr,
  ];

  try {
    // Call DQAWF
    Module._dqawf_(
      fPtr,
      aPtr,
      omegaPtr,
      integrPtr,
      epsabsPtr,
      resultPtr,
      abserrPtr,
      nevalPtr,
      ierPtr,
      limlstPtr,
      lstPtr,
      leniwPtr,
      maxp1Ptr,
      lenwPtr,
      iworkPtr,
      workPtr
    );

    // Read results
    const result = readDouble(Module, resultPtr);
    const abserr = readDouble(Module, abserrPtr);
    const neval = readInt(Module, nevalPtr);
    const ier = readInt(Module, ierPtr);
    const lst = readInt(Module, lstPtr);

    return {
      result,
      abserr,
      neval,
      subdivisions: lst, // For DQAWF, lst is the number of cycles used
      ier,
      success: ier === 0,
      message: getQuadpackMessage(ier),
    };
  } finally {
    unregisterIntegrand(Module, fPtr);
    freeAll(Module, allPtrs);
  }
}
