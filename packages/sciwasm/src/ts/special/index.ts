/**
 * Special mathematical functions module
 *
 * Provides special functions similar to scipy.special, including:
 * - Gamma functions and factorials
 * - Beta functions
 * - Error functions (erf, erfc, erfcx, erfi)
 * - Bessel functions (J, Y, I, K variants)
 * - Combinatorial functions (binomial coefficients, permutations)
 */

// Gamma functions
export { gamma, gammaln, rgamma } from './gamma.js';

// Beta functions
export { beta, betaln } from './beta.js';

// Error functions
export { erf, erfc, erfcx, erfi } from './erf.js';

// Bessel functions
export { j0, j1, jv, y0, y1, yv, i0, i1, iv, k0, k1 } from './bessel.js';

// Factorials
export { factorial, type FactorialOptions } from './factorial.js';
export { factorial2, type Factorial2Options } from './factorial2.js';
export { factorialk, type FactorialkOptions } from './factorialk.js';

// Combinatorial
export { comb, type CombOptions } from './comb.js';
export { perm, type PermOptions } from './perm.js';
