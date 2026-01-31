/**
 * QUADPACK High-Level TypeScript Interfaces
 *
 * User-friendly API for numerical integration using QUADPACK's
 * adaptive quadrature routines.
 */

// ============================================================
// FUNCTION TYPES
// ============================================================

/**
 * Integrand function type.
 *
 * A function that takes a single real number and returns a real number.
 * This represents the function f(x) to be integrated.
 *
 * @param x - The evaluation point
 * @returns The function value f(x)
 *
 * @example
 * ```ts
 * // Simple polynomial
 * const f: IntegrandFunction = x => x * x;
 *
 * // Trigonometric function
 * const g: IntegrandFunction = x => Math.sin(x);
 *
 * // Function with parameters via closure
 * const makeGaussian = (mu: number, sigma: number): IntegrandFunction =>
 *   x => Math.exp(-((x - mu) ** 2) / (2 * sigma ** 2));
 * ```
 */
export type IntegrandFunction = (x: number) => number;

// ============================================================
// OPTIONS INTERFACES
// ============================================================

/**
 * Base options shared by all integration routines.
 */
export interface QuadOptionsBase {
  /**
   * Absolute error tolerance.
   * The integration stops when the estimated absolute error is below this value.
   * @default 1e-8
   */
  epsabs?: number;

  /**
   * Relative error tolerance.
   * The integration stops when the estimated relative error is below this value.
   * @default 1e-8
   */
  epsrel?: number;

  /**
   * Maximum number of subdivisions allowed.
   * Higher values allow more accurate results but use more memory.
   * @default 50
   */
  limit?: number;
}

/**
 * Options for general adaptive integration (quad).
 *
 * Uses DQAGS internally, which is the best general-purpose routine for
 * most integrands, including those with integrable singularities.
 *
 * @example
 * ```ts
 * const options: QuadOptions = {
 *   epsabs: 1e-10,
 *   epsrel: 1e-10,
 *   limit: 100
 * };
 * ```
 */
export interface QuadOptions extends QuadOptionsBase {
  // No additional options for basic usage
}

/**
 * Options for adaptive integration with explicit rule selection (quad_rule).
 *
 * Uses DQAG internally, which allows selecting the Gauss-Kronrod rule.
 * Higher-order rules are more efficient for smooth functions.
 */
export interface QuadRuleOptions extends QuadOptionsBase {
  /**
   * Gauss-Kronrod rule to use:
   * - 1: 15-point (7 Gauss + 8 Kronrod) - fastest
   * - 2: 21-point
   * - 3: 31-point - good default
   * - 4: 41-point
   * - 5: 51-point
   * - 6: 61-point - most accurate per subdivision
   * @default 1
   */
  rule?: 1 | 2 | 3 | 4 | 5 | 6;
}

/**
 * Options for infinite interval integration (quad_inf).
 *
 * Uses DQAGI internally. The interval type is determined automatically
 * from the bounds passed to the function:
 * - [a, Infinity]: Semi-infinite upper
 * - [-Infinity, b]: Semi-infinite lower
 * - [-Infinity, Infinity]: Doubly infinite
 */
export interface QuadInfOptions extends QuadOptionsBase {
  // Interval type is auto-detected from bounds
}

/**
 * Options for oscillatory integration (quad_osc).
 *
 * Uses DQAWO internally. Computes integrals of the form:
 * ∫ f(x) * cos(ω*x) dx  or  ∫ f(x) * sin(ω*x) dx
 */
export interface QuadOscOptions extends QuadOptionsBase {
  /**
   * Angular frequency of the oscillation.
   * This parameter is required.
   */
  omega: number;

  /**
   * Oscillatory weight function type.
   * - 'cos': w(x) = cos(omega * x)
   * - 'sin': w(x) = sin(omega * x)
   * @default 'cos'
   */
  weight?: 'cos' | 'sin';

  /**
   * Maximum number of Chebyshev moments to compute.
   * Higher values improve accuracy for highly oscillatory functions.
   * @default 21
   */
  maxp1?: number;
}

/**
 * Options for Fourier integration over semi-infinite intervals (quad_fourier).
 *
 * Uses DQAWF internally. Computes integrals of the form:
 * ∫[a,∞) f(x) * cos(ω*x) dx  or  ∫[a,∞) f(x) * sin(ω*x) dx
 */
export interface QuadFourierOptions extends QuadOptionsBase {
  /**
   * Angular frequency of the oscillation.
   * This parameter is required.
   */
  omega: number;

  /**
   * Oscillatory weight function type.
   * - 'cos': w(x) = cos(omega * x)
   * - 'sin': w(x) = sin(omega * x)
   * @default 'cos'
   */
  weight?: 'cos' | 'sin';

  /**
   * Maximum number of cycles (periods of the oscillation).
   * @default 50
   */
  limlst?: number;

  /**
   * Maximum number of Chebyshev moments per cycle.
   * @default 21
   */
  maxp1?: number;
}

/**
 * Options for integration with algebraic-logarithmic singularities (quad_singular).
 *
 * Uses DQAWS internally. Computes integrals of the form:
 * ∫[a,b] f(x) * (x-a)^α * (b-x)^β * w(x) dx
 *
 * where w(x) can include logarithmic factors.
 */
export interface QuadSingularOptions extends QuadOptionsBase {
  /**
   * Exponent of (x-a) factor at the left endpoint.
   * Must be greater than -1.
   */
  alfa: number;

  /**
   * Exponent of (b-x) factor at the right endpoint.
   * Must be greater than -1.
   */
  beta: number;

  /**
   * Logarithmic weight function type:
   * - 1: w(x) = 1 (no logarithmic factor)
   * - 2: w(x) = log(x-a)
   * - 3: w(x) = log(b-x)
   * - 4: w(x) = log(x-a) * log(b-x)
   * @default 1
   */
  wgtfunc?: 1 | 2 | 3 | 4;
}

/**
 * Options for Cauchy principal value integration (quad_cauchy).
 *
 * Uses DQAWC internally. Computes the principal value:
 * P.V. ∫[a,b] f(x) / (x-c) dx
 *
 * where c is a point in (a, b).
 */
export interface QuadCauchyOptions extends QuadOptionsBase {
  /**
   * The singular point in the interval (a, b).
   * Must satisfy a < c < b.
   */
  c: number;
}

/**
 * Options for integration with user-specified break points (quad_break).
 *
 * Uses DQAGP internally. Useful when the integrand has known
 * discontinuities or singularities at specific points.
 */
export interface QuadBreakOptions extends QuadOptionsBase {
  /**
   * Array of break points within the interval (a, b).
   * These are points where the integrand may have difficulties
   * (discontinuities, peaks, etc.).
   */
  points: number[] | Float64Array;
}

/**
 * Options for non-adaptive integration (quad_ng).
 *
 * Uses DQNG internally. This is the simplest and fastest routine,
 * suitable only for smooth, well-behaved functions.
 */
export interface QuadNgOptions {
  /**
   * Absolute error tolerance.
   * @default 1e-8
   */
  epsabs?: number;

  /**
   * Relative error tolerance.
   * @default 1e-8
   */
  epsrel?: number;
}

// ============================================================
// RESULT INTERFACE
// ============================================================

/**
 * Result from a fixed-order Gauss-Kronrod rule (quad_fixed).
 *
 * Unlike QuadResult, this includes additional diagnostic information
 * about the integral approximation.
 */
export interface QuadFixedResult {
  /**
   * The computed integral approximation.
   */
  result: number;

  /**
   * Estimate of the absolute error.
   */
  abserr: number;

  /**
   * Approximation to the integral of |f(x)|.
   * Useful for detecting highly oscillatory behavior.
   */
  resabs: number;

  /**
   * Approximation to the integral of |f(x) - I/(b-a)|,
   * where I is the integral approximation.
   * Useful for detecting peaky integrands.
   */
  resasc: number;
}

/**
 * Result from a QUADPACK integration routine.
 */
export interface QuadResult {
  /**
   * The computed integral value.
   */
  result: number;

  /**
   * Estimate of the absolute error in the result.
   * The true error is likely smaller than this estimate.
   */
  abserr: number;

  /**
   * Number of integrand evaluations performed.
   */
  neval: number;

  /**
   * Number of subdivisions used in the adaptive algorithm.
   */
  subdivisions: number;

  /**
   * QUADPACK error code:
   * - 0: Normal termination, requested accuracy achieved
   * - 1: Maximum subdivisions reached
   * - 2: Roundoff error detected
   * - 3: Bad integrand behavior
   * - 4: Convergence failure
   * - 5: Integral is probably divergent
   * - 6: Invalid input
   */
  ier: number;

  /**
   * True if the integration was successful (ier === 0).
   */
  success: boolean;

  /**
   * Human-readable status message.
   */
  message: string;
}

// ============================================================
// UNIFIED API TYPES
// ============================================================

/**
 * Integration type for the unified integrate() function.
 */
export type IntegrationType =
  | 'general'
  | 'infinite'
  | 'oscillatory'
  | 'fourier'
  | 'singular'
  | 'cauchy'
  | 'breakpoints'
  | 'rule'
  | 'nonadaptive';

/**
 * Unified options for the integrate() function.
 *
 * Combines all type-specific options. The integration type is either
 * specified explicitly or auto-detected from the parameters.
 */
export interface IntegrateOptions extends QuadOptionsBase {
  /**
   * Explicit integration type. If not specified, the type is auto-detected
   * based on the bounds and other options provided.
   */
  type?: IntegrationType;

  // Oscillatory options
  /** Angular frequency for oscillatory/Fourier integration */
  omega?: number;
  /** Weight function: 'cos' or 'sin' */
  weight?: 'cos' | 'sin';
  /** Maximum Chebyshev moments */
  maxp1?: number;

  // Fourier-specific
  /** Maximum cycles for Fourier integration */
  limlst?: number;

  // Singular weight options
  /** Exponent at left endpoint for singular weight */
  alfa?: number;
  /** Exponent at right endpoint for singular weight */
  beta?: number;
  /** Logarithmic weight type (1-4) */
  wgtfunc?: 1 | 2 | 3 | 4;

  // Cauchy principal value
  /** Singular point for Cauchy principal value */
  c?: number;

  // Break points
  /** User-specified break points */
  points?: number[] | Float64Array;

  // Rule selection
  /** Gauss-Kronrod rule (1-6) */
  rule?: 1 | 2 | 3 | 4 | 5 | 6;
}
