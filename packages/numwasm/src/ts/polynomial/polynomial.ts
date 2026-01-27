/**
 * NumJS Power Series Polynomial
 *
 * Standard polynomial representation using power series:
 * p(x) = c[0] + c[1]*x + c[2]*x^2 + ... + c[n]*x^n
 */

import { ABCPolyBase, companionEigenvalues } from './_polybase.js';
import {
  PolyError,
  trimseq,
  polyadd,
  polysub,
  polymul,
  polydiv,
  polypow,
} from './polyutils.js';

/**
 * Power series polynomial.
 *
 * Represents polynomials in the standard form:
 * p(x) = c[0] + c[1]*x + c[2]*x^2 + ... + c[n]*x^n
 *
 * @example
 * const p = new Polynomial([1, 2, 3]); // 1 + 2x + 3x^2
 * p.call(2);  // 1 + 4 + 12 = 17
 */
export class Polynomial extends ABCPolyBase {
  static defaultDomain: [number, number] = [-1, 1];
  static defaultWindow: [number, number] = [-1, 1];
  static basisName = 'x';

  /* ============ Basis Operations ============ */

  protected static _add(c1: number[], c2: number[]): number[] {
    return polyadd(c1, c2);
  }

  protected static _sub(c1: number[], c2: number[]): number[] {
    return polysub(c1, c2);
  }

  protected static _mul(c1: number[], c2: number[]): number[] {
    return polymul(c1, c2);
  }

  protected static _div(c1: number[], c2: number[]): [number[], number[]] {
    return polydiv(c1, c2);
  }

  protected static _pow(c: number[], n: number, maxpow: number): number[] {
    return polypow(c, n, maxpow);
  }

  /**
   * Evaluate polynomial using Horner's method.
   *
   * Horner's method computes:
   * p(x) = c[0] + x*(c[1] + x*(c[2] + ... + x*c[n]))
   *
   * This is numerically stable and O(n).
   */
  protected static _val(x: number | number[], c: number[]): number | number[] {
    if (typeof x === 'number') {
      return polyval(x, c);
    }
    return x.map(v => polyval(v, c));
  }

  /**
   * Compute derivative coefficients.
   *
   * d/dx (c[0] + c[1]*x + c[2]*x^2 + ...) = c[1] + 2*c[2]*x + 3*c[3]*x^2 + ...
   */
  protected static _der(c: number[], m: number, scl: number): number[] {
    let result = [...c];

    for (let cnt = 0; cnt < m; cnt++) {
      if (result.length <= 1) {
        return [0];
      }

      const newCoef: number[] = [];
      for (let i = 1; i < result.length; i++) {
        newCoef.push(i * result[i] * scl);
      }
      result = newCoef;
    }

    return trimseq(result);
  }

  /**
   * Compute integral coefficients.
   *
   * Integral of c[0] + c[1]*x + c[2]*x^2 + ... is:
   * k + c[0]*x + (c[1]/2)*x^2 + (c[2]/3)*x^3 + ...
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
      const newCoef: number[] = [0];

      for (let i = 0; i < result.length; i++) {
        newCoef.push((result[i] * scl) / (i + 1));
      }

      // Adjust for lower bound and integration constant
      if (lbnd !== 0 || k[cnt] !== 0) {
        const val = polyval(lbnd, newCoef);
        newCoef[0] = k[cnt] - val;
      } else {
        newCoef[0] = k[cnt];
      }

      result = newCoef;
    }

    return trimseq(result);
  }

  /**
   * Find roots using companion matrix eigenvalues.
   */
  protected static async _roots(c: number[]): Promise<number[]> {
    c = trimseq(c);

    if (c.length <= 1) {
      return [];
    }

    if (c.length === 2) {
      // Linear: c[0] + c[1]*x = 0 => x = -c[0]/c[1]
      return [-c[0] / c[1]];
    }

    if (c.length === 3) {
      // Quadratic: c[0] + c[1]*x + c[2]*x^2 = 0
      const [a, b, cc] = [c[2], c[1], c[0]];
      const disc = b * b - 4 * a * cc;
      if (disc >= 0) {
        const sqrtDisc = Math.sqrt(disc);
        return [(-b + sqrtDisc) / (2 * a), (-b - sqrtDisc) / (2 * a)].sort(
          (x, y) => x - y
        );
      } else {
        // Complex roots - return empty for real polynomial
        return [];
      }
    }

    // Use companion matrix eigenvalues
    const companion = Polynomial._companion(c);
    const eigenvalues = await companionEigenvalues(companion);
    return eigenvalues.sort((a, b) => a - b);
  }

  /**
   * Generate Vandermonde matrix for power series.
   *
   * Row i contains: [1, x[i], x[i]^2, ..., x[i]^deg]
   */
  protected static _vander(x: number[], deg: number): number[][] {
    const result: number[][] = [];

    for (const xi of x) {
      const row: number[] = [1];
      for (let j = 1; j <= deg; j++) {
        row.push(row[j - 1] * xi);
      }
      result.push(row);
    }

    return result;
  }

  /**
   * Generate companion matrix.
   *
   * For polynomial c[0] + c[1]*x + ... + c[n]*x^n (monic form),
   * the companion matrix has eigenvalues equal to the roots.
   */
  protected static _companion(c: number[]): number[][] {
    c = trimseq(c);
    const n = c.length - 1;

    if (n < 1) {
      throw new PolyError('Companion matrix requires degree >= 1');
    }

    // Normalize so leading coefficient is 1 (monic)
    const monic = c.map(v => v / c[n]);

    // Build companion matrix
    // Standard form: subdiagonal of 1s, last column is -c[i]/c[n]
    const mat: number[][] = [];
    for (let i = 0; i < n; i++) {
      const row = new Array(n).fill(0);
      if (i < n - 1) {
        row[i + 1] = 1;
      }
      row[0] = -monic[i];
      mat.push(row);
    }

    // Transpose to get eigenvalues as roots
    // (companion matrix convention varies, this matches NumPy)
    const transposed: number[][] = [];
    for (let i = 0; i < n; i++) {
      transposed.push([]);
      for (let j = 0; j < n; j++) {
        transposed[i].push(mat[j][i]);
      }
    }

    return transposed;
  }

  /**
   * Construct polynomial from roots.
   *
   * Given roots r1, r2, ..., rn, constructs (x - r1)(x - r2)...(x - rn)
   */
  protected static _fromroots(roots: number[]): number[] {
    if (roots.length === 0) {
      return [1];
    }

    // Build polynomial from roots: prod(x - r)
    let result = [1];
    for (const r of roots) {
      result = polymul(result, [-r, 1]);
    }

    return result;
  }
}

/* ============ Module-Level Functions ============ */

/**
 * Evaluate polynomial at x using Horner's method.
 *
 * @param x - Point at which to evaluate
 * @param c - Coefficient array [c0, c1, c2, ...] for c0 + c1*x + c2*x^2 + ...
 * @returns Evaluated value
 *
 * @example
 * polyval(2, [1, 2, 3]) // 1 + 2*2 + 3*4 = 17
 */
export function polyval(x: number, c: number[]): number {
  if (c.length === 0) return 0;

  let result = c[c.length - 1];
  for (let i = c.length - 2; i >= 0; i--) {
    result = result * x + c[i];
  }

  return result;
}

/**
 * Evaluate 2D polynomial on x-y pairs.
 *
 * @param x - x coordinates
 * @param y - y coordinates
 * @param c - 2D coefficient array
 * @returns Array of evaluated values
 */
export function polyval2d(x: number[], y: number[], c: number[][]): number[] {
  if (x.length !== y.length) {
    throw new PolyError('x and y must have same length');
  }

  const result: number[] = [];

  for (let i = 0; i < x.length; i++) {
    let val = 0;
    for (let j = 0; j < c.length; j++) {
      for (let k = 0; k < c[j].length; k++) {
        val += c[j][k] * Math.pow(x[i], j) * Math.pow(y[i], k);
      }
    }
    result.push(val);
  }

  return result;
}

/**
 * Evaluate 3D polynomial on x-y-z triples.
 *
 * @param x - x coordinates
 * @param y - y coordinates
 * @param z - z coordinates
 * @param c - 3D coefficient array
 * @returns Array of evaluated values
 */
export function polyval3d(
  x: number[],
  y: number[],
  z: number[],
  c: number[][][]
): number[] {
  if (x.length !== y.length || y.length !== z.length) {
    throw new PolyError('x, y, and z must have same length');
  }

  const result: number[] = [];

  for (let i = 0; i < x.length; i++) {
    let val = 0;
    for (let j = 0; j < c.length; j++) {
      for (let k = 0; k < c[j].length; k++) {
        for (let l = 0; l < c[j][k].length; l++) {
          val += c[j][k][l] * Math.pow(x[i], j) * Math.pow(y[i], k) * Math.pow(z[i], l);
        }
      }
    }
    result.push(val);
  }

  return result;
}

/**
 * Generate Vandermonde matrix.
 *
 * @param x - Sample points
 * @param deg - Degree of polynomial
 * @returns Vandermonde matrix
 */
export function polyvander(x: number[], deg: number): number[][] {
  return Polynomial['_vander'](x, deg);
}

/**
 * Generate 2D Vandermonde matrix.
 *
 * @param x - x sample points
 * @param y - y sample points
 * @param deg - [degree in x, degree in y]
 * @returns 2D Vandermonde matrix
 */
export function polyvander2d(
  x: number[],
  y: number[],
  deg: [number, number]
): number[][] {
  if (x.length !== y.length) {
    throw new PolyError('x and y must have same length');
  }

  const [degX, degY] = deg;
  const result: number[][] = [];

  for (let i = 0; i < x.length; i++) {
    const row: number[] = [];
    for (let j = 0; j <= degX; j++) {
      for (let k = 0; k <= degY; k++) {
        row.push(Math.pow(x[i], j) * Math.pow(y[i], k));
      }
    }
    result.push(row);
  }

  return result;
}

/**
 * Polynomial derivative.
 *
 * @param c - Coefficient array
 * @param m - Number of derivatives (default: 1)
 * @returns Derivative coefficients
 */
export function polyder(c: number[], m: number = 1): number[] {
  return Polynomial['_der'](c, m, 1);
}

/**
 * Polynomial integral.
 *
 * @param c - Coefficient array
 * @param m - Number of integrals (default: 1)
 * @param k - Integration constants (default: all zero)
 * @param lbnd - Lower bound (default: 0)
 * @returns Integral coefficients
 */
export function polyint(
  c: number[],
  m: number = 1,
  k: number[] = [],
  lbnd: number = 0
): number[] {
  const kPadded = [...k];
  while (kPadded.length < m) kPadded.push(0);
  return Polynomial['_int'](c, m, kPadded, lbnd, 1);
}

/**
 * Least squares polynomial fit.
 *
 * @param x - x values
 * @param y - y values
 * @param deg - Degree of fitting polynomial
 * @param rcond - Cutoff for small singular values
 * @param full - If true, return extra info
 * @param w - Weights
 * @returns Coefficient array (or [coefficients, info] if full=true)
 */
export async function polyfit(
  x: number[],
  y: number[],
  deg: number,
  rcond: number | null = null,
  full: boolean = false,
  w: number[] | null = null
): Promise<number[] | [number[], { residuals: number[]; rank: number; sv: number[]; rcond: number }]> {
  const result = await Polynomial.fit(x, y, deg, null, rcond, full, w);

  if (full) {
    const [poly, info] = result as [Polynomial, { residuals: number[]; rank: number; sv: number[]; rcond: number }];
    return [poly.coef, info];
  }

  return (result as Polynomial).coef;
}

/**
 * Find polynomial roots.
 *
 * @param c - Coefficient array
 * @returns Array of roots
 */
export async function polyroots(c: number[]): Promise<number[]> {
  return Polynomial['_roots'](c);
}

/**
 * Generate companion matrix.
 *
 * @param c - Coefficient array
 * @returns Companion matrix
 */
export function polycompanion(c: number[]): number[][] {
  return Polynomial['_companion'](c);
}

/**
 * Construct polynomial from roots.
 *
 * @param roots - Array of roots
 * @returns Coefficient array
 */
export function polyfromroots(roots: number[]): number[] {
  return Polynomial['_fromroots'](roots);
}

/**
 * Add two polynomials.
 *
 * @param c1 - First coefficient array
 * @param c2 - Second coefficient array
 * @returns Sum coefficients
 */
export { polyadd } from './polyutils.js';

/**
 * Subtract two polynomials.
 *
 * @param c1 - First coefficient array
 * @param c2 - Second coefficient array
 * @returns Difference coefficients
 */
export { polysub } from './polyutils.js';

/**
 * Multiply two polynomials.
 *
 * @param c1 - First coefficient array
 * @param c2 - Second coefficient array
 * @returns Product coefficients
 */
export { polymul } from './polyutils.js';

/**
 * Divide two polynomials.
 *
 * @param c1 - Dividend coefficient array
 * @param c2 - Divisor coefficient array
 * @returns [quotient, remainder]
 */
export { polydiv } from './polyutils.js';

/**
 * Raise polynomial to a power.
 *
 * @param c - Coefficient array
 * @param pow - Non-negative integer power
 * @param maxpower - Maximum allowed power
 * @returns Power coefficients
 */
export { polypow } from './polyutils.js';
