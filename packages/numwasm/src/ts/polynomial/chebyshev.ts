/**
 * NumJS Chebyshev Polynomials
 *
 * Chebyshev polynomials of the first kind (T_n).
 * Defined by: T_n(x) = cos(n * arccos(x))
 *
 * Chebyshev polynomials have excellent numerical properties for
 * interpolation and approximation on [-1, 1].
 */

import { ABCPolyBase, companionEigenvalues } from './_polybase.js';
import { PolyError, trimseq, polymul as polyMul, polydiv as polyDiv } from './polyutils.js';

/**
 * Chebyshev polynomial of the first kind.
 *
 * Chebyshev polynomials are defined by:
 * T_n(x) = cos(n * arccos(x))
 *
 * Recurrence: T_{n+1}(x) = 2x*T_n(x) - T_{n-1}(x)
 * with T_0(x) = 1, T_1(x) = x
 *
 * @example
 * const p = new Chebyshev([1, 2, 3]); // 1*T_0 + 2*T_1 + 3*T_2
 */
export class Chebyshev extends ABCPolyBase {
  static defaultDomain: [number, number] = [-1, 1];
  static defaultWindow: [number, number] = [-1, 1];
  static basisName = 'T';

  /* ============ Basis Operations ============ */

  protected static _add(c1: number[], c2: number[]): number[] {
    const len = Math.max(c1.length, c2.length);
    const result = new Array(len).fill(0);

    for (let i = 0; i < c1.length; i++) result[i] += c1[i];
    for (let i = 0; i < c2.length; i++) result[i] += c2[i];

    return trimseq(result);
  }

  protected static _sub(c1: number[], c2: number[]): number[] {
    const len = Math.max(c1.length, c2.length);
    const result = new Array(len).fill(0);

    for (let i = 0; i < c1.length; i++) result[i] += c1[i];
    for (let i = 0; i < c2.length; i++) result[i] -= c2[i];

    return trimseq(result);
  }

  /**
   * Chebyshev multiplication using the identity:
   * T_m(x) * T_n(x) = 0.5 * (T_{m+n}(x) + T_{|m-n|}(x))
   */
  protected static _mul(c1: number[], c2: number[]): number[] {
    c1 = trimseq(c1);
    c2 = trimseq(c2);

    if (c1.length === 0 || c2.length === 0) {
      return [0];
    }

    const len = c1.length + c2.length - 1;
    const result = new Array(len).fill(0);

    for (let i = 0; i < c1.length; i++) {
      for (let j = 0; j < c2.length; j++) {
        const prod = c1[i] * c2[j] * 0.5;
        result[i + j] += prod;
        result[Math.abs(i - j)] += prod;
      }
    }

    return trimseq(result);
  }

  /**
   * Chebyshev division.
   * Convert to power series, divide, convert back.
   */
  protected static _div(c1: number[], c2: number[]): [number[], number[]] {
    c1 = trimseq(c1);
    c2 = trimseq(c2);

    if (c2.length === 1 && c2[0] === 0) {
      throw new PolyError('Division by zero polynomial');
    }

    if (c1.length < c2.length) {
      return [[0], c1];
    }

    // Convert to power series, divide, convert back
    const p1 = cheb2poly(c1);
    const p2 = cheb2poly(c2);

    const [quo, rem] = polyDiv(p1, p2);

    return [poly2cheb(quo), poly2cheb(rem)];
  }

  protected static _pow(c: number[], n: number, maxpow: number): number[] {
    if (n === 0) return [1];
    if (n === 1) return trimseq([...c]);

    // Check degree limit
    c = trimseq(c);
    if ((c.length - 1) * n > maxpow) {
      throw new PolyError(`Power would exceed maxpower (${maxpow})`);
    }

    // Binary exponentiation
    let result = [1];
    let base = [...c];
    let exp = n;

    while (exp > 0) {
      if (exp & 1) {
        result = Chebyshev._mul(result, base);
      }
      base = Chebyshev._mul(base, base);
      exp >>= 1;
    }

    return trimseq(result);
  }

  /**
   * Evaluate Chebyshev polynomial using Clenshaw's algorithm.
   *
   * Clenshaw's algorithm is numerically stable and O(n).
   */
  protected static _val(x: number | number[], c: number[]): number | number[] {
    if (typeof x === 'number') {
      return chebval(x, c);
    }
    return x.map(v => chebval(v, c));
  }

  /**
   * Chebyshev derivative.
   *
   * d/dx T_n(x) = n * U_{n-1}(x) where U is Chebyshev of second kind
   *
   * In the T basis, the derivative formula is:
   * c'_k = 2(k+1)*c_{k+1} + c'_{k+2} (for k >= 1)
   * c'_0 = c_1 + c'_2/2
   */
  protected static _der(c: number[], m: number, scl: number): number[] {
    let result = [...c];

    for (let cnt = 0; cnt < m; cnt++) {
      const n = result.length;
      if (n <= 1) return [0];

      const der = new Array(n - 1).fill(0);

      // Work backwards from highest degree
      for (let j = n - 1; j >= 2; j--) {
        der[j - 1] = 2 * j * result[j];
        der[j - 2] += der[j - 1] / (2 * (j - 1));
      }

      if (n >= 2) {
        der[0] += 2 * result[1];
      }

      // Apply scale
      result = der.map(v => v * scl);
    }

    return trimseq(result);
  }

  /**
   * Chebyshev integral.
   *
   * The integral of T_n(x) is:
   * - For n = 0: T_1(x)
   * - For n = 1: T_2(x)/4
   * - For n >= 2: T_{n+1}(x)/(2(n+1)) - T_{n-1}(x)/(2(n-1))
   */
  protected static _int(
    c: number[],
    m: number,
    k: number[],
    lbnd: number,
    scl: number
  ): number[] {
    let result = [...c];

    for (let cnt = 0; cnt < m; cnt++) {
      const n = result.length;
      const integ = new Array(n + 1).fill(0);

      // Integration formula
      for (let j = 0; j < n; j++) {
        if (j === 0) {
          // Integral of T_0 = 1 is T_1 = x
          integ[1] += result[0] * scl;
        } else if (j === 1) {
          // Integral of T_1 = x is T_0/2 + T_2/4
          integ[0] += result[1] * scl / 2;
          integ[2] += result[1] * scl / 4;
        } else {
          // Integral of T_n is T_{n+1}/(2(n+1)) - T_{n-1}/(2(n-1))
          integ[j + 1] += result[j] * scl / (2 * (j + 1));
          integ[j - 1] -= result[j] * scl / (2 * (j - 1));
        }
      }

      // Add integration constant by adjusting for lower bound
      const val = chebval(lbnd, integ);
      integ[0] += k[cnt] - val;

      result = integ;
    }

    return trimseq(result);
  }

  /**
   * Find Chebyshev polynomial roots using companion matrix.
   */
  protected static async _roots(c: number[]): Promise<number[]> {
    c = trimseq(c);

    if (c.length <= 1) return [];
    if (c.length === 2) {
      // Linear: c[0]*T_0 + c[1]*T_1 = c[0] + c[1]*x = 0
      return [-c[0] / c[1]];
    }

    // Use companion matrix
    const companion = Chebyshev._companion(c);
    const eigenvalues = await companionEigenvalues(companion);

    // Filter out roots outside [-1, 1] with some tolerance
    // (Chebyshev is only defined on [-1, 1], but companion matrix may give roots outside)
    return eigenvalues.sort((a, b) => a - b);
  }

  /**
   * Chebyshev Vandermonde matrix.
   *
   * Row i contains: [T_0(x[i]), T_1(x[i]), ..., T_deg(x[i])]
   */
  protected static _vander(x: number[], deg: number): number[][] {
    const result: number[][] = [];

    for (const xi of x) {
      const row: number[] = [1]; // T_0 = 1
      if (deg >= 1) {
        row.push(xi); // T_1 = x
      }
      // Recurrence: T_{n+1} = 2x*T_n - T_{n-1}
      for (let j = 2; j <= deg; j++) {
        row.push(2 * xi * row[j - 1] - row[j - 2]);
      }
      result.push(row);
    }

    return result;
  }

  /**
   * Chebyshev companion matrix.
   *
   * The companion matrix for Chebyshev polynomials has a symmetric
   * structure that provides better numerical accuracy for root finding.
   */
  protected static _companion(c: number[]): number[][] {
    c = trimseq(c);
    const n = c.length - 1;

    if (n < 1) {
      throw new PolyError('Companion matrix requires degree >= 1');
    }

    // Build Chebyshev companion matrix (symmetric tridiagonal form)
    const mat: number[][] = [];
    for (let i = 0; i < n; i++) {
      const row = new Array(n).fill(0);
      mat.push(row);
    }

    // First row
    mat[0][1] = 1;

    // Interior rows (tridiagonal structure)
    for (let i = 1; i < n - 1; i++) {
      mat[i][i - 1] = 0.5;
      mat[i][i + 1] = 0.5;
    }

    // Last row
    if (n >= 2) {
      mat[n - 1][n - 2] = 0.5;
    }

    // Modify last column for the polynomial coefficients
    const scale = -0.5 / c[n];
    for (let i = 0; i < n; i++) {
      mat[i][n - 1] += scale * c[i];
    }
    if (n >= 2) {
      mat[n - 2][n - 1] += 0.5;
    }

    return mat;
  }

  /**
   * Construct Chebyshev polynomial from roots.
   */
  protected static _fromroots(roots: number[]): number[] {
    if (roots.length === 0) return [1];

    // Build polynomial from roots in power basis, then convert
    let poly = [1];
    for (const r of roots) {
      poly = polyMul(poly, [-r, 1]);
    }

    return poly2cheb(poly);
  }
}

/* ============ Module-Level Functions ============ */

/**
 * Evaluate Chebyshev polynomial using Clenshaw's algorithm.
 *
 * @param x - Point at which to evaluate
 * @param c - Coefficient array [c0, c1, c2, ...] for c0*T_0 + c1*T_1 + c2*T_2 + ...
 * @returns Evaluated value
 *
 * @example
 * chebval(0.5, [1, 2, 3]) // 1*T_0(0.5) + 2*T_1(0.5) + 3*T_2(0.5)
 */
export function chebval(x: number, c: number[]): number {
  if (c.length === 0) return 0;
  if (c.length === 1) return c[0];
  if (c.length === 2) return c[0] + c[1] * x;

  // Clenshaw's algorithm
  const x2 = 2 * x;
  let c0 = c[c.length - 2];
  let c1 = c[c.length - 1];

  for (let i = c.length - 3; i >= 0; i--) {
    const tmp = c0;
    c0 = c[i] - c1;
    c1 = tmp + c1 * x2;
  }

  return c0 + c1 * x;
}

/**
 * Evaluate 2D Chebyshev polynomial.
 *
 * @param x - x coordinates
 * @param y - y coordinates
 * @param c - 2D coefficient array
 * @returns Array of evaluated values
 */
export function chebval2d(x: number[], y: number[], c: number[][]): number[] {
  if (x.length !== y.length) {
    throw new PolyError('x and y must have same length');
  }

  const result: number[] = [];

  for (let i = 0; i < x.length; i++) {
    // Evaluate as sum of c[j][k] * T_j(x) * T_k(y)
    let val = 0;
    for (let j = 0; j < c.length; j++) {
      const Tj = chebval(x[i], [0, ...new Array(j).fill(0), 1].slice(0, j + 1));
      for (let k = 0; k < c[j].length; k++) {
        const Tk = chebval(y[i], [0, ...new Array(k).fill(0), 1].slice(0, k + 1));
        val += c[j][k] * Tj * Tk;
      }
    }
    result.push(val);
  }

  return result;
}

/**
 * Convert power series coefficients to Chebyshev coefficients.
 *
 * @param pol - Power series coefficients
 * @returns Chebyshev coefficients
 */
export function poly2cheb(pol: number[]): number[] {
  pol = trimseq(pol);
  const n = pol.length;

  if (n === 0) return [0];
  if (n === 1) return [pol[0]];

  // Build up the result by adding each x^j term converted to Chebyshev basis
  const result = new Array(n).fill(0);
  for (let j = 0; j < n; j++) {
    // x^j in Chebyshev basis
    const xjCheb = _powerToCheb(j);
    for (let k = 0; k < xjCheb.length; k++) {
      if (k < result.length) {
        result[k] += pol[j] * xjCheb[k];
      }
    }
  }

  return trimseq(result);
}

/**
 * Convert Chebyshev coefficients to power series.
 *
 * @param c - Chebyshev coefficients
 * @returns Power series coefficients
 */
export function cheb2poly(c: number[]): number[] {
  c = trimseq(c);
  const n = c.length;

  if (n === 0) return [0];
  if (n === 1) return [c[0]];

  // Build power series from Chebyshev
  const result = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    // Add c[k] * T_k expressed in power basis
    const tk = _chebToPowerCoefs(k);
    for (let j = 0; j < tk.length; j++) {
      result[j] += c[k] * tk[j];
    }
  }

  return trimseq(result);
}

/**
 * Chebyshev Vandermonde matrix.
 *
 * @param x - Sample points
 * @param deg - Degree of polynomial
 * @returns Vandermonde matrix
 */
export function chebvander(x: number[], deg: number): number[][] {
  return Chebyshev['_vander'](x, deg);
}

/**
 * Chebyshev derivative.
 *
 * @param c - Coefficient array
 * @param m - Number of derivatives (default: 1)
 * @returns Derivative coefficients
 */
export function chebder(c: number[], m: number = 1): number[] {
  return Chebyshev['_der'](c, m, 1);
}

/**
 * Chebyshev integral.
 *
 * @param c - Coefficient array
 * @param m - Number of integrals (default: 1)
 * @param k - Integration constants (default: all zero)
 * @param lbnd - Lower bound (default: 0)
 * @returns Integral coefficients
 */
export function chebint(
  c: number[],
  m: number = 1,
  k: number[] = [],
  lbnd: number = 0
): number[] {
  const kPadded = [...k];
  while (kPadded.length < m) kPadded.push(0);
  return Chebyshev['_int'](c, m, kPadded, lbnd, 1);
}

/**
 * Chebyshev least squares fit.
 *
 * @param x - x values
 * @param y - y values
 * @param deg - Degree of fitting polynomial
 * @param rcond - Cutoff for small singular values
 * @param full - If true, return extra info
 * @param w - Weights
 * @returns Coefficient array (or [coefficients, info] if full=true)
 */
export async function chebfit(
  x: number[],
  y: number[],
  deg: number,
  rcond: number | null = null,
  full: boolean = false,
  w: number[] | null = null
): Promise<number[] | [number[], { residuals: number[]; rank: number; sv: number[]; rcond: number }]> {
  const result = await Chebyshev.fit(x, y, deg, null, rcond, full, w);

  if (full) {
    const [poly, info] = result as [Chebyshev, { residuals: number[]; rank: number; sv: number[]; rcond: number }];
    return [poly.coef, info];
  }

  return (result as Chebyshev).coef;
}

/**
 * Find Chebyshev polynomial roots.
 *
 * @param c - Coefficient array
 * @returns Array of roots
 */
export async function chebroots(c: number[]): Promise<number[]> {
  return Chebyshev['_roots'](c);
}

/**
 * Generate Chebyshev companion matrix.
 *
 * @param c - Coefficient array
 * @returns Companion matrix
 */
export function chebcompanion(c: number[]): number[][] {
  return Chebyshev['_companion'](c);
}

/**
 * Construct Chebyshev polynomial from roots.
 *
 * @param roots - Array of roots
 * @returns Coefficient array
 */
export function chebfromroots(roots: number[]): number[] {
  return Chebyshev['_fromroots'](roots);
}

/**
 * Interpolate a function using Chebyshev polynomials.
 *
 * Uses Chebyshev nodes for optimal interpolation on [-1, 1].
 *
 * @param func - Function to interpolate
 * @param deg - Degree of interpolating polynomial
 * @param domain - Domain for interpolation (default: [-1, 1])
 * @returns Chebyshev coefficients
 */
export async function chebinterpolate(
  func: (x: number) => number,
  deg: number,
  domain: [number, number] = [-1, 1]
): Promise<number[]> {
  // Use Chebyshev nodes: x_k = cos((2k+1)pi/(2(n+1))) for k = 0, ..., n
  const n = deg + 1;
  const nodes: number[] = [];
  for (let k = 0; k < n; k++) {
    nodes.push(Math.cos(((2 * k + 1) * Math.PI) / (2 * n)));
  }

  // Map nodes to domain
  const [d0, d1] = domain;
  const mappedNodes = nodes.map(x => ((d0 + d1) / 2) + ((d1 - d0) / 2) * x);

  // Evaluate function at mapped nodes
  const y = mappedNodes.map(func);

  // Fit Chebyshev polynomial (nodes are already in [-1, 1])
  return (await chebfit(nodes, y, deg)) as number[];
}

/**
 * Add two Chebyshev polynomials.
 */
export function chebadd(c1: number[], c2: number[]): number[] {
  return Chebyshev['_add'](c1, c2);
}

/**
 * Subtract two Chebyshev polynomials.
 */
export function chebsub(c1: number[], c2: number[]): number[] {
  return Chebyshev['_sub'](c1, c2);
}

/**
 * Multiply two Chebyshev polynomials.
 */
export function chebmul(c1: number[], c2: number[]): number[] {
  return Chebyshev['_mul'](c1, c2);
}

/**
 * Divide two Chebyshev polynomials.
 */
export function chebdiv(c1: number[], c2: number[]): [number[], number[]] {
  return Chebyshev['_div'](c1, c2);
}

/**
 * Raise Chebyshev polynomial to a power.
 */
export function chebpow(c: number[], pow: number, maxpower: number = 100): number[] {
  return Chebyshev['_pow'](c, pow, maxpower);
}

/* ============ Helper Functions ============ */

/**
 * Get T_k(x) expressed in power basis coefficients.
 *
 * Uses recurrence: T_{n+1} = 2x*T_n - T_{n-1}
 */
function _chebToPowerCoefs(k: number): number[] {
  if (k === 0) return [1];
  if (k === 1) return [0, 1];

  // Use recurrence: T_{n+1} = 2x*T_n - T_{n-1}
  let prev2 = [1]; // T_0
  let prev1 = [0, 1]; // T_1

  for (let n = 2; n <= k; n++) {
    // 2x * prev1
    const term1 = [0, ...prev1.map(c => 2 * c)];

    // - prev2
    const result = new Array(Math.max(term1.length, prev2.length)).fill(0);
    for (let i = 0; i < term1.length; i++) result[i] += term1[i];
    for (let i = 0; i < prev2.length; i++) result[i] -= prev2[i];

    prev2 = prev1;
    prev1 = trimseq(result);
  }

  return prev1;
}

/**
 * Get x^j expressed in Chebyshev basis.
 *
 * Uses the identity: x^n = 2^{1-n} * sum_{k} C(n, (n-k)/2) * T_k
 * where the sum is over k with same parity as n.
 */
function _powerToCheb(j: number): number[] {
  if (j === 0) return [1];
  if (j === 1) return [0, 1];

  // Build x^j by multiplying x^{j-1} by x
  // In Chebyshev: x * T_n = 0.5 * (T_{n+1} + T_{n-1}) for n >= 1
  //              x * T_0 = T_1
  let result = [0, 1]; // x = T_1

  for (let n = 2; n <= j; n++) {
    // Multiply current result by x (in Chebyshev basis)
    const newResult = new Array(result.length + 1).fill(0);

    for (let k = 0; k < result.length; k++) {
      if (k === 0) {
        // x * c_0 * T_0 = c_0 * T_1
        newResult[1] += result[0];
      } else {
        // x * c_k * T_k = 0.5 * c_k * (T_{k+1} + T_{k-1})
        newResult[k + 1] += 0.5 * result[k];
        newResult[k - 1] += 0.5 * result[k];
      }
    }

    result = trimseq(newResult);
  }

  return result;
}

