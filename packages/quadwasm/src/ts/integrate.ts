/**
 * INTEGRATE - Unified Integration API
 *
 * Provides a single entry point that dispatches to the appropriate
 * QUADPACK routine based on the integration bounds and options.
 */

import { quad } from './quad.js';
import { quad_inf } from './quad_inf.js';
import { quad_osc, quad_fourier } from './quad_osc.js';
import { quad_singular, quad_cauchy } from './quad_singular.js';
import { quad_break } from './quad_break.js';
import { quad_rule } from './quad_rule.js';
import { quad_ng } from './quad_ng.js';
import type {
  IntegrandFunction,
  IntegrateOptions,
  QuadResult,
} from './high-level-types.js';

/**
 * Compute the definite integral of f(x).
 *
 * This is the unified integration API that automatically selects the
 * most appropriate QUADPACK routine based on the integration bounds
 * and options provided.
 *
 * **Auto-detection logic:**
 * 1. If bounds include ±Infinity:
 *    - If omega specified: quad_fourier (Fourier transform)
 *    - Otherwise: quad_inf (infinite interval)
 * 2. If omega specified: quad_osc (oscillatory)
 * 3. If c specified: quad_cauchy (Cauchy principal value)
 * 4. If alfa or beta specified: quad_singular (algebraic singularities)
 * 5. If points specified: quad_break (user break points)
 * 6. If rule specified: quad_rule (explicit GK rule)
 * 7. Otherwise: quad (general adaptive)
 *
 * @param f - The integrand function f(x)
 * @param a - Lower limit of integration (use -Infinity for lower infinite)
 * @param b - Upper limit of integration (use Infinity for upper infinite)
 * @param options - Integration options
 * @returns Integration result with value, error estimate, and diagnostics
 *
 * @example
 * ```ts
 * import { integrate } from 'quadwasm';
 *
 * // Simple integral (uses quad internally)
 * const r1 = await integrate(x => x * x, 0, 1);
 *
 * // Infinite interval (uses quad_inf internally)
 * const r2 = await integrate(x => Math.exp(-x * x), 0, Infinity);
 *
 * // Oscillatory (uses quad_osc internally)
 * const r3 = await integrate(x => x, 0, Math.PI, { omega: 10 });
 *
 * // Cauchy principal value (uses quad_cauchy internally)
 * const r4 = await integrate(x => 1, -1, 1, { c: 0 });
 *
 * // Singular weight (uses quad_singular internally)
 * const r5 = await integrate(Math.sin, 0, 1, { alfa: -0.5, beta: 0 });
 *
 * // With break points (uses quad_break internally)
 * const r6 = await integrate(x => x < 0.5 ? 0 : 1, 0, 1, { points: [0.5] });
 *
 * // Explicit GK rule (uses quad_rule internally)
 * const r7 = await integrate(x => x * x, 0, 1, { rule: 6 }); // GK61
 * ```
 */
export async function integrate(
  f: IntegrandFunction,
  a: number,
  b: number,
  options?: IntegrateOptions
): Promise<QuadResult> {
  const opts = options ?? {};

  // Check for explicit type override
  if (opts.type) {
    return dispatchByType(f, a, b, opts);
  }

  // Auto-detect based on bounds and options
  const aIsInf = !Number.isFinite(a);
  const bIsInf = !Number.isFinite(b);

  // Infinite interval cases
  if (aIsInf || bIsInf) {
    if (opts.omega !== undefined) {
      // Fourier transform over semi-infinite interval
      if (!aIsInf && bIsInf) {
        return quad_fourier(f, a, {
          omega: opts.omega,
          weight: opts.weight,
          epsabs: opts.epsabs,
          limlst: opts.limlst,
          limit: opts.limit,
          maxp1: opts.maxp1,
        });
      }
      throw new Error(
        'Fourier integration (quad_fourier) only supports [a, +∞) intervals'
      );
    }
    return quad_inf(f, a, b, {
      epsabs: opts.epsabs,
      epsrel: opts.epsrel,
      limit: opts.limit,
    });
  }

  // Finite interval cases
  if (opts.omega !== undefined) {
    return quad_osc(f, a, b, {
      omega: opts.omega,
      weight: opts.weight,
      epsabs: opts.epsabs,
      epsrel: opts.epsrel,
      limit: opts.limit,
      maxp1: opts.maxp1,
    });
  }

  if (opts.c !== undefined) {
    return quad_cauchy(f, a, b, {
      c: opts.c,
      epsabs: opts.epsabs,
      epsrel: opts.epsrel,
      limit: opts.limit,
    });
  }

  if (opts.alfa !== undefined || opts.beta !== undefined) {
    return quad_singular(f, a, b, {
      alfa: opts.alfa ?? 0,
      beta: opts.beta ?? 0,
      wgtfunc: opts.wgtfunc,
      epsabs: opts.epsabs,
      epsrel: opts.epsrel,
      limit: opts.limit,
    });
  }

  // Break points
  if (opts.points !== undefined) {
    return quad_break(f, a, b, {
      points: opts.points,
      epsabs: opts.epsabs,
      epsrel: opts.epsrel,
      limit: opts.limit,
    });
  }

  // Explicit rule selection
  if (opts.rule !== undefined) {
    return quad_rule(f, a, b, {
      rule: opts.rule,
      epsabs: opts.epsabs,
      epsrel: opts.epsrel,
      limit: opts.limit,
    });
  }

  // Default: general adaptive quadrature
  return quad(f, a, b, {
    epsabs: opts.epsabs,
    epsrel: opts.epsrel,
    limit: opts.limit,
  });
}

/**
 * Dispatch based on explicit type specification.
 */
async function dispatchByType(
  f: IntegrandFunction,
  a: number,
  b: number,
  opts: IntegrateOptions
): Promise<QuadResult> {
  switch (opts.type) {
    case 'general':
      return quad(f, a, b, {
        epsabs: opts.epsabs,
        epsrel: opts.epsrel,
        limit: opts.limit,
      });

    case 'infinite':
      return quad_inf(f, a, b, {
        epsabs: opts.epsabs,
        epsrel: opts.epsrel,
        limit: opts.limit,
      });

    case 'oscillatory':
      if (opts.omega === undefined) {
        throw new Error('omega is required for oscillatory integration');
      }
      return quad_osc(f, a, b, {
        omega: opts.omega,
        weight: opts.weight,
        epsabs: opts.epsabs,
        epsrel: opts.epsrel,
        limit: opts.limit,
        maxp1: opts.maxp1,
      });

    case 'fourier':
      if (opts.omega === undefined) {
        throw new Error('omega is required for Fourier integration');
      }
      return quad_fourier(f, a, {
        omega: opts.omega,
        weight: opts.weight,
        epsabs: opts.epsabs,
        limlst: opts.limlst,
        limit: opts.limit,
        maxp1: opts.maxp1,
      });

    case 'singular':
      return quad_singular(f, a, b, {
        alfa: opts.alfa ?? 0,
        beta: opts.beta ?? 0,
        wgtfunc: opts.wgtfunc,
        epsabs: opts.epsabs,
        epsrel: opts.epsrel,
        limit: opts.limit,
      });

    case 'cauchy':
      if (opts.c === undefined) {
        throw new Error('c is required for Cauchy principal value integration');
      }
      return quad_cauchy(f, a, b, {
        c: opts.c,
        epsabs: opts.epsabs,
        epsrel: opts.epsrel,
        limit: opts.limit,
      });

    case 'breakpoints':
      if (opts.points === undefined) {
        throw new Error('points is required for break point integration');
      }
      return quad_break(f, a, b, {
        points: opts.points,
        epsabs: opts.epsabs,
        epsrel: opts.epsrel,
        limit: opts.limit,
      });

    case 'rule':
      return quad_rule(f, a, b, {
        rule: opts.rule,
        epsabs: opts.epsabs,
        epsrel: opts.epsrel,
        limit: opts.limit,
      });

    case 'nonadaptive':
      return quad_ng(f, a, b, {
        epsabs: opts.epsabs,
        epsrel: opts.epsrel,
      });

    default:
      throw new Error(`Unknown integration type: ${opts.type}`);
  }
}
