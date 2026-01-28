/**
 * Special mathematical functions module
 *
 * Provides special functions similar to scipy.special, including:
 * - Gamma functions and factorials
 * - Combinatorial functions (binomial coefficients, permutations)
 * - (Future: error functions, Bessel functions, etc.)
 */

export { gamma, gammaln, rgamma } from './gamma.js';
export { factorial, type FactorialOptions } from './factorial.js';
export { factorial2, type Factorial2Options } from './factorial2.js';
export { factorialk, type FactorialkOptions } from './factorialk.js';
export { comb, type CombOptions } from './comb.js';
export { perm, type PermOptions } from './perm.js';
