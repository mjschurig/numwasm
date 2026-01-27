/**
 * Special mathematical functions.
 * @module special
 */

import { NotImplementedError } from '../errors.js';

/** Gamma function. Mirrors scipy.special.gamma. */
export function gamma(_z: number): number {
  throw new NotImplementedError('sciwasm.special.gamma');
}

/** Natural log of the absolute value of the gamma function. */
export function gammaln(_x: number): number {
  throw new NotImplementedError('sciwasm.special.gammaln');
}

/** Beta function. Mirrors scipy.special.beta. */
export function beta(_a: number, _b: number): number {
  throw new NotImplementedError('sciwasm.special.beta');
}

/** Error function. Mirrors scipy.special.erf. */
export function erf(_z: number): number {
  throw new NotImplementedError('sciwasm.special.erf');
}

/** Complementary error function. */
export function erfc(_x: number): number {
  throw new NotImplementedError('sciwasm.special.erfc');
}

/** Bessel function of the first kind, order 0. */
export function j0(_x: number): number {
  throw new NotImplementedError('sciwasm.special.j0');
}

/** Bessel function of the first kind, order 1. */
export function j1(_x: number): number {
  throw new NotImplementedError('sciwasm.special.j1');
}

/** Factorial. Mirrors scipy.special.factorial. */
export function factorial(_n: number, _exact?: boolean): number {
  throw new NotImplementedError('sciwasm.special.factorial');
}

/** Binomial coefficient (n choose k). */
export function comb(_N: number, _k: number, _exact?: boolean): number {
  throw new NotImplementedError('sciwasm.special.comb');
}

/** Permutations (n P k). */
export function perm(_N: number, _k: number, _exact?: boolean): number {
  throw new NotImplementedError('sciwasm.special.perm');
}
