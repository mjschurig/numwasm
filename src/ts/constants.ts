/**
 * NumJS Constants
 *
 * Mathematical and special constants matching NumPy's constant exports.
 * Reference: numpy/_core/numeric.py, numpy/_core/umath.py
 */

/* ============ Mathematical Constants ============ */

/**
 * Euler's number, the base of natural logarithms.
 *
 * `e = 2.718281828459045`
 *
 * @example
 * import { e, exp, NDArray } from 'numjs';
 * console.log(e);  // 2.718281828459045
 *
 * @see https://numpy.org/doc/stable/reference/constants.html#numpy.e
 */
export const e: number = 2.718281828459045;

/**
 * Pi, the ratio of a circle's circumference to its diameter.
 *
 * `π = 3.141592653589793`
 *
 * @example
 * import { pi, sin, NDArray } from 'numjs';
 * console.log(pi);  // 3.141592653589793
 *
 * @see https://numpy.org/doc/stable/reference/constants.html#numpy.pi
 */
export const pi: number = 3.141592653589793;

/**
 * Euler-Mascheroni constant (γ).
 *
 * `γ = lim(n→∞) [1 + 1/2 + 1/3 + ... + 1/n - ln(n)]`
 * `γ ≈ 0.5772156649015329`
 *
 * @example
 * import { euler_gamma } from 'numjs';
 * console.log(euler_gamma);  // 0.5772156649015329
 *
 * @see https://numpy.org/doc/stable/reference/constants.html#numpy.euler_gamma
 */
export const euler_gamma: number = 0.5772156649015329;

/* ============ Special Floating-Point Values ============ */

/**
 * IEEE 754 positive infinity.
 *
 * @example
 * import { inf, isinf } from 'numjs';
 * console.log(inf);           // Infinity
 * console.log(1 / 0 === inf); // true
 *
 * @see https://numpy.org/doc/stable/reference/constants.html#numpy.inf
 */
export const inf: number = Infinity;

/**
 * IEEE 754 positive infinity (alias).
 *
 * @see inf
 */
export const PINF: number = Infinity;

/**
 * IEEE 754 negative infinity.
 *
 * @example
 * import { NINF, isneginf } from 'numjs';
 * console.log(NINF);            // -Infinity
 * console.log(-1 / 0 === NINF); // true
 *
 * @see https://numpy.org/doc/stable/reference/constants.html#numpy.NINF
 */
export const NINF: number = -Infinity;

/**
 * IEEE 754 quiet Not-a-Number.
 *
 * NaN is used to represent undefined or unrepresentable results
 * in floating-point calculations (e.g., 0/0, sqrt(-1)).
 *
 * Note: NaN is not equal to anything, including itself.
 * Use `isnan()` to check for NaN values.
 *
 * @example
 * import { nan, isnan } from 'numjs';
 * console.log(nan);              // NaN
 * console.log(nan === nan);      // false (NaN never equals NaN)
 * console.log(Number.isNaN(nan)); // true
 *
 * @see https://numpy.org/doc/stable/reference/constants.html#numpy.nan
 */
export const nan: number = NaN;

/**
 * IEEE 754 quiet Not-a-Number (alias).
 *
 * @see nan
 */
export const NAN: number = NaN;

/**
 * IEEE 754 positive zero.
 *
 * JavaScript/IEEE 754 distinguishes between +0 and -0.
 * They compare equal, but have different behavior in edge cases.
 *
 * @example
 * import { PZERO, NZERO } from 'numjs';
 * console.log(PZERO);          // 0
 * console.log(1 / PZERO);      // Infinity
 * console.log(1 / NZERO);      // -Infinity
 *
 * @see https://numpy.org/doc/stable/reference/constants.html#numpy.PZERO
 */
export const PZERO: number = 0.0;

/**
 * IEEE 754 negative zero.
 *
 * @see PZERO for examples
 * @see https://numpy.org/doc/stable/reference/constants.html#numpy.NZERO
 */
export const NZERO: number = -0.0;

/* ============ Indexing Constants ============ */

/**
 * Alias for newaxis, used to expand array dimensions.
 *
 * When used in an index expression, adds a new axis of length 1
 * at that position.
 *
 * @example
 * import { newaxis, NDArray, slice } from 'numjs';
 *
 * const arr = await NDArray.fromArray([1, 2, 3]); // shape: [3]
 * const row = arr.slice([newaxis, slice()]); // shape: [1, 3]
 * const col = arr.slice([slice(), newaxis]); // shape: [3, 1]
 *
 * @see https://numpy.org/doc/stable/reference/constants.html#numpy.newaxis
 */
export { newaxis } from './slice.js';
