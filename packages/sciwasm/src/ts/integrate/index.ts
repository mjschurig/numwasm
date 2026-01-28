/**
 * Integration module: numerical quadrature and ODE solvers.
 * Mirrors scipy.integrate.
 * @module integrate
 */

export { quad } from './quad.js';
export type {
  IntegrandFunction,
  QuadOptions,
  QuadResult,
  QuadFullResult,
  QuadFullResultWithMessage,
  QuadInfoDict,
} from './types.js';
