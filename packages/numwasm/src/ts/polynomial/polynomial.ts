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

/* ============ NumPy Legacy Functions ============ */

/**
 * Find the roots of a polynomial with coefficients given in c.
 *
 * This is a NumPy legacy function. The coefficients are given in descending
 * powers order (highest power first), which is the opposite of the numpy.polynomial
 * module convention.
 *
 * NOTE: This is the legacy numpy.roots function. For the modern API, use
 * polyroots() from the polynomial module which uses ascending powers order.
 *
 * @param p - Coefficient array in descending powers order [c_n, c_{n-1}, ..., c_1, c_0]
 * @returns Array of roots
 *
 * @example
 * ```typescript
 * // Find roots of x^2 - 4 (coefficients: 1*x^2 + 0*x - 4)
 * const r = await roots([1, 0, -4]);  // [2, -2]
 *
 * // Find roots of x^2 + 2x + 1 = (x+1)^2
 * const r = await roots([1, 2, 1]);   // [-1, -1]
 * ```
 */
export async function roots(p: number[]): Promise<number[]> {
  // Convert from descending powers (NumPy legacy) to ascending powers (polynomial module)
  const ascending = [...p].reverse();
  return polyroots(ascending);
}

/**
 * Find the coefficients of a polynomial with the given sequence of roots.
 *
 * This is a NumPy legacy function. Returns coefficients in descending
 * powers order (highest power first), which is the opposite of the numpy.polynomial
 * module convention.
 *
 * NOTE: This is the legacy numpy.poly function. For the modern API, use
 * polyfromroots() from the polynomial module which uses ascending powers order.
 *
 * @param seq_of_zeros - A sequence of polynomial roots
 * @returns Polynomial coefficients in descending powers order [c_n, c_{n-1}, ..., c_1, c_0]
 *
 * @example
 * ```typescript
 * // Create polynomial with roots at 2 and -2
 * const c = poly([2, -2]);  // [1, 0, -4] representing x^2 - 4
 *
 * // Create polynomial with root at -1 (double root)
 * const c = poly([-1, -1]); // [1, 2, 1] representing x^2 + 2x + 1 = (x+1)^2
 * ```
 */
export function poly(seq_of_zeros: number[]): number[] {
  // Get coefficients in ascending powers order
  const ascending = polyfromroots(seq_of_zeros);
  // Convert to descending powers order (NumPy legacy convention)
  return ascending.reverse();
}

/**
 * A one-dimensional polynomial class.
 *
 * NOTE: This is a convenience class for use with the legacy numpy.poly1d API.
 * For new code, the Polynomial class is preferred.
 *
 * A poly1d instance represents a polynomial p(x):
 *   p(x) = c[0]*x^n + c[1]*x^(n-1) + ... + c[n-1]*x + c[n]
 *
 * where coefficients are stored in descending powers order.
 *
 * @example
 * ```typescript
 * // Create from coefficients (descending powers)
 * const p = new poly1d([1, 2, 3]);  // x^2 + 2x + 3
 *
 * // Evaluate at a point
 * p.call(2);  // 1*4 + 2*2 + 3 = 11
 *
 * // Create from roots
 * const p = new poly1d([1, -1], true);  // (x-1)(x+1) = x^2 - 1
 *
 * // Arithmetic operations
 * const sum = p.add(q);
 * const prod = p.mul(q);
 * ```
 */
export class poly1d {
  /** Coefficients in descending powers order */
  readonly c: number[];

  /** Degree of the polynomial */
  readonly order: number;

  /** Alias for c - coefficients */
  readonly coeffs: number[];

  /** Alias for c - coefficients */
  readonly coef: number[];

  /** Whether constructed from roots */
  readonly variable: string;

  /**
   * Create a polynomial.
   *
   * @param c_or_r - Coefficients (descending powers) or roots
   * @param r - If true, c_or_r are roots. Default false.
   * @param variable - Symbol for display. Default 'x'.
   */
  constructor(c_or_r: number[], r: boolean = false, variable: string = "x") {
    this.variable = variable;

    let coeffs: number[];
    if (r) {
      // Create from roots - use the poly function
      coeffs = poly(c_or_r);
    } else {
      // Direct coefficients
      coeffs = [...c_or_r];
    }

    // Trim leading zeros
    while (coeffs.length > 1 && coeffs[0] === 0) {
      coeffs.shift();
    }

    this.c = coeffs;
    this.coeffs = this.c;
    this.coef = this.c;
    this.order = Math.max(0, coeffs.length - 1);
  }

  /**
   * Evaluate the polynomial at x.
   */
  call(x: number): number {
    // Horner's method for descending powers
    let result = 0;
    for (const coef of this.c) {
      result = result * x + coef;
    }
    return result;
  }

  /**
   * Evaluate the polynomial at each element of x.
   */
  callArray(x: number[]): number[] {
    return x.map((v) => this.call(v));
  }

  /**
   * Return the roots of the polynomial.
   * Note: Unlike NumPy, this is async due to WASM dependency.
   */
  get r(): Promise<number[]> {
    return roots(this.c);
  }

  /**
   * Alias for r - roots.
   * Note: Unlike NumPy, this is async due to WASM dependency.
   */
  get roots(): Promise<number[]> {
    return this.r;
  }

  /**
   * Compute and return the roots of the polynomial (async method version).
   */
  async getRoots(): Promise<number[]> {
    return roots(this.c);
  }

  /**
   * Add another polynomial.
   */
  add(other: poly1d | number[]): poly1d {
    const otherCoeffs = other instanceof poly1d ? other.c : other;
    // Pad shorter array with zeros on the left
    const maxLen = Math.max(this.c.length, otherCoeffs.length);
    const a = new Array(maxLen - this.c.length).fill(0).concat(this.c);
    const b = new Array(maxLen - otherCoeffs.length).fill(0).concat(otherCoeffs);
    const result = a.map((v, i) => v + b[i]);
    return new poly1d(result, false, this.variable);
  }

  /**
   * Subtract another polynomial.
   */
  sub(other: poly1d | number[]): poly1d {
    const otherCoeffs = other instanceof poly1d ? other.c : other;
    const maxLen = Math.max(this.c.length, otherCoeffs.length);
    const a = new Array(maxLen - this.c.length).fill(0).concat(this.c);
    const b = new Array(maxLen - otherCoeffs.length).fill(0).concat(otherCoeffs);
    const result = a.map((v, i) => v - b[i]);
    return new poly1d(result, false, this.variable);
  }

  /**
   * Multiply by another polynomial.
   */
  mul(other: poly1d | number[] | number): poly1d {
    if (typeof other === "number") {
      return new poly1d(this.c.map((v) => v * other), false, this.variable);
    }
    const otherCoeffs = other instanceof poly1d ? other.c : other;
    // Convolution
    const m = this.c.length;
    const n = otherCoeffs.length;
    const result = new Array(m + n - 1).fill(0);
    for (let i = 0; i < m; i++) {
      for (let j = 0; j < n; j++) {
        result[i + j] += this.c[i] * otherCoeffs[j];
      }
    }
    return new poly1d(result, false, this.variable);
  }

  /**
   * Negate the polynomial.
   */
  neg(): poly1d {
    return new poly1d(this.c.map((v) => -v), false, this.variable);
  }

  /**
   * Raise to a non-negative integer power.
   */
  pow(n: number): poly1d {
    if (!Number.isInteger(n) || n < 0) {
      throw new Error("Power must be a non-negative integer");
    }
    if (n === 0) {
      return new poly1d([1], false, this.variable);
    }
    let result = new poly1d([1], false, this.variable);
    for (let i = 0; i < n; i++) {
      result = result.mul(this);
    }
    return result;
  }

  /**
   * Differentiate the polynomial.
   *
   * @param m - Order of differentiation (default 1)
   */
  deriv(m: number = 1): poly1d {
    let coeffs = [...this.c];
    for (let d = 0; d < m; d++) {
      if (coeffs.length <= 1) {
        coeffs = [0];
        break;
      }
      const deg = coeffs.length - 1;
      const newCoeffs: number[] = [];
      for (let i = 0; i < deg; i++) {
        newCoeffs.push(coeffs[i] * (deg - i));
      }
      coeffs = newCoeffs;
    }
    return new poly1d(coeffs, false, this.variable);
  }

  /**
   * Integrate the polynomial.
   *
   * @param m - Order of integration (default 1)
   * @param k - Integration constants (default all zeros)
   */
  integ(m: number = 1, k: number | number[] = 0): poly1d {
    const kArr = typeof k === "number" ? [k] : k;
    let coeffs = [...this.c];

    for (let d = 0; d < m; d++) {
      const deg = coeffs.length;
      const newCoeffs: number[] = [];
      for (let i = 0; i < deg; i++) {
        newCoeffs.push(coeffs[i] / (deg - i));
      }
      // Add integration constant
      const constant = kArr[d] !== undefined ? kArr[d] : 0;
      newCoeffs.push(constant);
      coeffs = newCoeffs;
    }
    return new poly1d(coeffs, false, this.variable);
  }

  /**
   * Return a string representation.
   */
  toString(): string {
    if (this.c.length === 0 || (this.c.length === 1 && this.c[0] === 0)) {
      return "0";
    }

    const x = this.variable;
    const terms: string[] = [];
    const deg = this.order;

    for (let i = 0; i < this.c.length; i++) {
      const coef = this.c[i];
      const power = deg - i;

      if (coef === 0) continue;

      let term = "";
      const absCoef = Math.abs(coef);

      if (power === 0) {
        term = absCoef.toString();
      } else if (power === 1) {
        term = absCoef === 1 ? x : `${absCoef}*${x}`;
      } else {
        term = absCoef === 1 ? `${x}^${power}` : `${absCoef}*${x}^${power}`;
      }

      if (terms.length === 0) {
        terms.push(coef < 0 ? `-${term}` : term);
      } else {
        terms.push(coef < 0 ? ` - ${term}` : ` + ${term}`);
      }
    }

    return terms.join("") || "0";
  }
}
