/**
 * Calculus operations for symbolic mathematics.
 *
 * This module provides fundamental calculus operations including
 * differentiation, integration, limits, series expansion, and summation.
 * All functions work on symbolic expressions and return symbolic results.
 *
 * @module calculus
 *
 * @example
 * ```typescript
 * import { Symbol, diff, series, sin, cos, pow } from 'symwasm';
 *
 * const x = new Symbol('x');
 *
 * // Differentiation
 * diff(sin(x), x);           // cos(x)
 * diff(pow(x, 3), x);        // 3*x²
 * diff(sin(x), x, 2);        // -sin(x) (second derivative)
 *
 * // Series expansion
 * series(sin(x), x);         // x - x³/6 + x⁵/120
 * series(cos(x), x, 0, 4);   // 1 - x²/2
 * ```
 */

export { diff } from './diff.js';
export { series } from './series.js';
export { integrate } from './integrate.js';
export { limit } from './limit.js';
export { summation } from './summation.js';
