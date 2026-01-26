/**
 * NumJS Legendre Polynomials
 *
 * Legendre polynomials (P_n) are orthogonal polynomials on [-1, 1]
 * with weight function w(x) = 1.
 *
 * Recurrence: (n+1)*P_{n+1}(x) = (2n+1)*x*P_n(x) - n*P_{n-1}(x)
 */

import { ABCPolyBase, maxpower, companionEigenvalues } from './_polybase.js';
import { PolyError, trimseq, polymul as polyMul, polydiv as polyDiv } from './polyutils.js';

/**
 * Legendre polynomial.
 *
 * Legendre polynomials satisfy the recurrence:
 * P_{n+1}(x) = ((2n+1)*x*P_n(x) - n*P_{n-1}(x)) / (n+1)
 *
 * with P_0(x) = 1, P_1(x) = x
 *
 * @example
 * const p = new Legendre([1, 2, 3]); // 1*P_0 + 2*P_1 + 3*P_2
 */
export class Legendre extends ABCPolyBase {
  static defaultDomain: [number, number] = [-1, 1];
  static defaultWindow: [number, number] = [-1, 1];
  static basisName = 'P';

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
   * Legendre multiplication.
   *
   * Uses the linearization formula for Legendre polynomials.
   * For efficiency, converts to power series, multiplies, converts back.
   */
  protected static _mul(c1: number[], c2: number[]): number[] {
    c1 = trimseq(c1);
    c2 = trimseq(c2);

    if (c1.length === 0 || c2.length === 0) {
      return [0];
    }

    // Convert to power series, multiply, convert back
    const p1 = leg2poly(c1);
    const p2 = leg2poly(c2);
    const product = polyMul(p1, p2);

    return poly2leg(product);
  }

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
    const p1 = leg2poly(c1);
    const p2 = leg2poly(c2);
    const [quo, rem] = polyDiv(p1, p2);

    return [poly2leg(quo), poly2leg(rem)];
  }

  protected static _pow(c: number[], n: number, maxpow: number): number[] {
    if (n === 0) return [1];
    if (n === 1) return trimseq([...c]);

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
        result = Legendre._mul(result, base);
      }
      base = Legendre._mul(base, base);
      exp >>= 1;
    }

    return trimseq(result);
  }

  /**
   * Evaluate Legendre polynomial using Clenshaw-like recurrence.
   */
  protected static _val(x: number | number[], c: number[]): number | number[] {
    if (typeof x === 'number') {
      return legval(x, c);
    }
    return x.map(v => legval(v, c));
  }

  /**
   * Legendre derivative.
   *
   * d/dx P_n(x) can be computed using the recurrence:
   * (2n+1)*P_n = d/dx(P_{n+1} - P_{n-1})
   */
  protected static _der(c: number[], m: number, scl: number): number[] {
    let result = [...c];

    for (let cnt = 0; cnt < m; cnt++) {
      const n = result.length;
      if (n <= 1) return [0];

      const der = new Array(n - 1).fill(0);

      // Use the derivative relationship
      for (let j = n - 1; j >= 1; j--) {
        der[j - 1] = (2 * j - 1) * result[j];
        if (j >= 2) {
          der[j - 2] += der[j - 1] ? 0 : 0;
        }
      }

      // Accumulate using the recurrence
      for (let j = n - 3; j >= 0; j--) {
        der[j] += der[j + 2];
      }

      // Apply scale
      result = der.map(v => v * scl);
    }

    return trimseq(result);
  }

  /**
   * Legendre integral.
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

      // Integration formula based on:
      // integral(P_n) = (P_{n+1} - P_{n-1}) / (2n+1)
      for (let j = 0; j < n; j++) {
        if (j === 0) {
          // Integral of P_0 = 1 is x = P_1
          integ[1] += result[0] * scl;
        } else {
          // Integral of P_j
          integ[j + 1] += result[j] * scl / (2 * j + 1);
          if (j >= 1) {
            integ[j - 1] -= result[j] * scl / (2 * j + 1);
          }
        }
      }

      // Add integration constant
      const val = legval(lbnd, integ);
      integ[0] += k[cnt] - val;

      result = integ;
    }

    return trimseq(result);
  }

  protected static async _roots(c: number[]): Promise<number[]> {
    c = trimseq(c);

    if (c.length <= 1) return [];
    if (c.length === 2) {
      return [-c[0] / c[1]];
    }

    const companion = Legendre._companion(c);
    const eigenvalues = await companionEigenvalues(companion);
    return eigenvalues.sort((a, b) => a - b);
  }

  /**
   * Legendre Vandermonde matrix.
   */
  protected static _vander(x: number[], deg: number): number[][] {
    const result: number[][] = [];

    for (const xi of x) {
      const row: number[] = [1]; // P_0 = 1
      if (deg >= 1) {
        row.push(xi); // P_1 = x
      }
      // Recurrence: (n+1)*P_{n+1} = (2n+1)*x*P_n - n*P_{n-1}
      for (let j = 2; j <= deg; j++) {
        const pn = ((2 * j - 1) * xi * row[j - 1] - (j - 1) * row[j - 2]) / j;
        row.push(pn);
      }
      result.push(row);
    }

    return result;
  }

  /**
   * Legendre companion matrix.
   */
  protected static _companion(c: number[]): number[][] {
    c = trimseq(c);
    const n = c.length - 1;

    if (n < 1) {
      throw new PolyError('Companion matrix requires degree >= 1');
    }

    // Build companion matrix for Legendre polynomials
    // Similar structure to Chebyshev but with different coefficients
    const mat: number[][] = [];
    for (let i = 0; i < n; i++) {
      const row = new Array(n).fill(0);
      mat.push(row);
    }

    // Fill tridiagonal structure based on recurrence
    for (let i = 0; i < n - 1; i++) {
      // Subdiagonal: sqrt(i*(i+1)) / (2*i+1) related terms
      const k = i + 1;
      mat[i + 1][i] = k / Math.sqrt((2 * k - 1) * (2 * k + 1));
      mat[i][i + 1] = mat[i + 1][i];
    }

    // Modify for polynomial coefficients
    const scale = -1 / c[n];
    for (let i = 0; i < n; i++) {
      mat[i][n - 1] += scale * c[i] * Math.sqrt((2 * (n - 1) + 1) / (2 * i + 1));
    }

    return mat;
  }

  protected static _fromroots(roots: number[]): number[] {
    if (roots.length === 0) return [1];

    let poly = [1];
    for (const r of roots) {
      poly = polyMul(poly, [-r, 1]);
    }

    return poly2leg(poly);
  }
}

/* ============ Module-Level Functions ============ */

/**
 * Evaluate Legendre polynomial using Clenshaw-like algorithm.
 *
 * @param x - Point at which to evaluate
 * @param c - Coefficient array
 * @returns Evaluated value
 */
export function legval(x: number, c: number[]): number {
  if (c.length === 0) return 0;
  if (c.length === 1) return c[0];
  if (c.length === 2) return c[0] + c[1] * x;

  // Use Clenshaw's algorithm adapted for Legendre
  let c0 = c[c.length - 2];
  let c1 = c[c.length - 1];

  for (let i = c.length - 3; i >= 0; i--) {
    const n = i + 1;
    const tmp = c0;
    c0 = c[i] - c1 * n / (n + 1);
    c1 = tmp + c1 * x * (2 * n + 1) / (n + 1);
  }

  return c0 + c1 * x;
}

/**
 * Convert power series to Legendre coefficients.
 */
export function poly2leg(pol: number[]): number[] {
  pol = trimseq(pol);
  const n = pol.length;

  if (n === 0) return [0];
  if (n === 1) return [pol[0]];

  const result = new Array(n).fill(0);

  for (let j = 0; j < n; j++) {
    const xjLeg = _powerToLeg(j);
    for (let k = 0; k < xjLeg.length; k++) {
      if (k < result.length) {
        result[k] += pol[j] * xjLeg[k];
      }
    }
  }

  return trimseq(result);
}

/**
 * Convert Legendre coefficients to power series.
 */
export function leg2poly(c: number[]): number[] {
  c = trimseq(c);
  const n = c.length;

  if (n === 0) return [0];
  if (n === 1) return [c[0]];

  const result = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    const pk = _legToPowerCoefs(k);
    for (let j = 0; j < pk.length; j++) {
      result[j] += c[k] * pk[j];
    }
  }

  return trimseq(result);
}

/**
 * Legendre Vandermonde matrix.
 */
export function legvander(x: number[], deg: number): number[][] {
  return Legendre['_vander'](x, deg);
}

/**
 * Legendre derivative.
 */
export function legder(c: number[], m: number = 1): number[] {
  return Legendre['_der'](c, m, 1);
}

/**
 * Legendre integral.
 */
export function legint(
  c: number[],
  m: number = 1,
  k: number[] = [],
  lbnd: number = 0
): number[] {
  const kPadded = [...k];
  while (kPadded.length < m) kPadded.push(0);
  return Legendre['_int'](c, m, kPadded, lbnd, 1);
}

/**
 * Legendre least squares fit.
 */
export async function legfit(
  x: number[],
  y: number[],
  deg: number,
  rcond: number | null = null,
  full: boolean = false,
  w: number[] | null = null
): Promise<number[] | [number[], { residuals: number[]; rank: number; sv: number[]; rcond: number }]> {
  const result = await Legendre.fit(x, y, deg, null, rcond, full, w);

  if (full) {
    const [poly, info] = result as [Legendre, { residuals: number[]; rank: number; sv: number[]; rcond: number }];
    return [poly.coef, info];
  }

  return (result as Legendre).coef;
}

/**
 * Find Legendre polynomial roots.
 */
export async function legroots(c: number[]): Promise<number[]> {
  return Legendre['_roots'](c);
}

/**
 * Generate Legendre companion matrix.
 */
export function legcompanion(c: number[]): number[][] {
  return Legendre['_companion'](c);
}

/**
 * Construct Legendre polynomial from roots.
 */
export function legfromroots(roots: number[]): number[] {
  return Legendre['_fromroots'](roots);
}

/**
 * Add two Legendre polynomials.
 */
export function legadd(c1: number[], c2: number[]): number[] {
  return Legendre['_add'](c1, c2);
}

/**
 * Subtract two Legendre polynomials.
 */
export function legsub(c1: number[], c2: number[]): number[] {
  return Legendre['_sub'](c1, c2);
}

/**
 * Multiply two Legendre polynomials.
 */
export function legmul(c1: number[], c2: number[]): number[] {
  return Legendre['_mul'](c1, c2);
}

/**
 * Divide two Legendre polynomials.
 */
export function legdiv(c1: number[], c2: number[]): [number[], number[]] {
  return Legendre['_div'](c1, c2);
}

/**
 * Raise Legendre polynomial to a power.
 */
export function legpow(c: number[], pow: number, maxpower: number = 100): number[] {
  return Legendre['_pow'](c, pow, maxpower);
}

/* ============ Helper Functions ============ */

/**
 * Get P_k(x) expressed in power basis.
 */
function _legToPowerCoefs(k: number): number[] {
  if (k === 0) return [1];
  if (k === 1) return [0, 1];

  // Use recurrence: (n+1)*P_{n+1} = (2n+1)*x*P_n - n*P_{n-1}
  let prev2 = [1]; // P_0
  let prev1 = [0, 1]; // P_1

  for (let n = 1; n < k; n++) {
    // (n+1)*P_{n+1} = (2n+1)*x*P_n - n*P_{n-1}
    const xTimesP = [0, ...prev1.map(c => ((2 * n + 1) * c) / (n + 1))];
    const minusPrev = prev2.map(c => (-n * c) / (n + 1));

    const result = new Array(Math.max(xTimesP.length, minusPrev.length)).fill(0);
    for (let i = 0; i < xTimesP.length; i++) result[i] += xTimesP[i];
    for (let i = 0; i < minusPrev.length; i++) result[i] += minusPrev[i];

    prev2 = prev1;
    prev1 = trimseq(result);
  }

  return prev1;
}

/**
 * Get x^j expressed in Legendre basis.
 */
function _powerToLeg(j: number): number[] {
  if (j === 0) return [1];
  if (j === 1) return [0, 1];

  // Build up by multiplying by x
  // x * P_n = ((n+1)*P_{n+1} + n*P_{n-1}) / (2n+1)
  let result = [0, 1]; // x = P_1

  for (let n = 2; n <= j; n++) {
    const newResult = new Array(result.length + 1).fill(0);

    for (let k = 0; k < result.length; k++) {
      if (k === 0) {
        // x * c_0 * P_0 = c_0 * P_1
        newResult[1] += result[0];
      } else {
        // x * c_k * P_k = c_k * ((k+1)*P_{k+1} + k*P_{k-1}) / (2k+1)
        newResult[k + 1] += (result[k] * (k + 1)) / (2 * k + 1);
        newResult[k - 1] += (result[k] * k) / (2 * k + 1);
      }
    }

    result = trimseq(newResult);
  }

  return result;
}
