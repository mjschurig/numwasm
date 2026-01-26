/**
 * NumJS Hermite Polynomials (Physicist's)
 *
 * Physicist's Hermite polynomials (H_n) are orthogonal on (-inf, inf)
 * with weight function w(x) = exp(-x^2).
 *
 * Recurrence: H_{n+1}(x) = 2*x*H_n(x) - 2*n*H_{n-1}(x)
 */

import { ABCPolyBase, maxpower, companionEigenvalues } from './_polybase.js';
import { PolyError, trimseq, polymul as polyMul, polydiv as polyDiv } from './polyutils.js';

/**
 * Physicist's Hermite polynomial.
 *
 * Hermite polynomials satisfy the recurrence:
 * H_{n+1}(x) = 2*x*H_n(x) - 2*n*H_{n-1}(x)
 *
 * with H_0(x) = 1, H_1(x) = 2x
 *
 * These are orthogonal with weight exp(-x^2).
 *
 * @example
 * const p = new Hermite([1, 2, 3]); // 1*H_0 + 2*H_1 + 3*H_2
 */
export class Hermite extends ABCPolyBase {
  static defaultDomain: [number, number] = [-1, 1];
  static defaultWindow: [number, number] = [-1, 1];
  static basisName = 'H';

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

  protected static _mul(c1: number[], c2: number[]): number[] {
    c1 = trimseq(c1);
    c2 = trimseq(c2);

    if (c1.length === 0 || c2.length === 0) {
      return [0];
    }

    // Convert to power series, multiply, convert back
    const p1 = herm2poly(c1);
    const p2 = herm2poly(c2);
    const product = polyMul(p1, p2);

    return poly2herm(product);
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

    const p1 = herm2poly(c1);
    const p2 = herm2poly(c2);
    const [quo, rem] = polyDiv(p1, p2);

    return [poly2herm(quo), poly2herm(rem)];
  }

  protected static _pow(c: number[], n: number, maxpow: number): number[] {
    if (n === 0) return [1];
    if (n === 1) return trimseq([...c]);

    c = trimseq(c);
    if ((c.length - 1) * n > maxpow) {
      throw new PolyError(`Power would exceed maxpower (${maxpow})`);
    }

    let result = [1];
    let base = [...c];
    let exp = n;

    while (exp > 0) {
      if (exp & 1) {
        result = Hermite._mul(result, base);
      }
      base = Hermite._mul(base, base);
      exp >>= 1;
    }

    return trimseq(result);
  }

  protected static _val(x: number | number[], c: number[]): number | number[] {
    if (typeof x === 'number') {
      return hermval(x, c);
    }
    return x.map(v => hermval(v, c));
  }

  /**
   * Hermite derivative.
   *
   * d/dx H_n(x) = 2*n*H_{n-1}(x)
   */
  protected static _der(c: number[], m: number, scl: number): number[] {
    let result = [...c];

    for (let cnt = 0; cnt < m; cnt++) {
      const n = result.length;
      if (n <= 1) return [0];

      const der = new Array(n - 1).fill(0);

      for (let j = 1; j < n; j++) {
        der[j - 1] = 2 * j * result[j] * scl;
      }

      result = der;
    }

    return trimseq(result);
  }

  /**
   * Hermite integral.
   *
   * integral(H_n) = H_{n+1} / (2*(n+1))
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

      for (let j = 0; j < n; j++) {
        integ[j + 1] = result[j] * scl / (2 * (j + 1));
      }

      // Add integration constant
      const val = hermval(lbnd, integ);
      integ[0] = k[cnt] - val;

      result = integ;
    }

    return trimseq(result);
  }

  protected static async _roots(c: number[]): Promise<number[]> {
    c = trimseq(c);

    if (c.length <= 1) return [];
    if (c.length === 2) {
      // Linear: c[0] + c[1]*2x = 0 (H_1 = 2x)
      return [-c[0] / (2 * c[1])];
    }

    const companion = Hermite._companion(c);
    const eigenvalues = await companionEigenvalues(companion);
    return eigenvalues.sort((a, b) => a - b);
  }

  protected static _vander(x: number[], deg: number): number[][] {
    const result: number[][] = [];

    for (const xi of x) {
      const row: number[] = [1]; // H_0 = 1
      if (deg >= 1) {
        row.push(2 * xi); // H_1 = 2x
      }
      // Recurrence: H_{n+1} = 2x*H_n - 2n*H_{n-1}
      for (let j = 2; j <= deg; j++) {
        row.push(2 * xi * row[j - 1] - 2 * (j - 1) * row[j - 2]);
      }
      result.push(row);
    }

    return result;
  }

  protected static _companion(c: number[]): number[][] {
    c = trimseq(c);
    const n = c.length - 1;

    if (n < 1) {
      throw new PolyError('Companion matrix requires degree >= 1');
    }

    const mat: number[][] = [];
    for (let i = 0; i < n; i++) {
      const row = new Array(n).fill(0);
      mat.push(row);
    }

    // Build tridiagonal companion matrix
    for (let i = 0; i < n - 1; i++) {
      mat[i][i + 1] = 0.5;
      mat[i + 1][i] = i + 1;
    }

    // Modify for polynomial coefficients
    const scale = -0.5 / c[n];
    for (let i = 0; i < n; i++) {
      mat[i][n - 1] += scale * c[i];
    }

    return mat;
  }

  protected static _fromroots(roots: number[]): number[] {
    if (roots.length === 0) return [1];

    let poly = [1];
    for (const r of roots) {
      poly = polyMul(poly, [-r, 1]);
    }

    return poly2herm(poly);
  }
}

/* ============ Module-Level Functions ============ */

/**
 * Evaluate Hermite polynomial.
 */
export function hermval(x: number, c: number[]): number {
  if (c.length === 0) return 0;
  if (c.length === 1) return c[0];
  if (c.length === 2) return c[0] + c[1] * 2 * x;

  // Clenshaw-like algorithm for Hermite
  let c0 = c[c.length - 2];
  let c1 = c[c.length - 1];

  for (let i = c.length - 3; i >= 0; i--) {
    const n = i + 1;
    const tmp = c0;
    c0 = c[i] - c1 * 2 * n;
    c1 = tmp + c1 * 2 * x;
  }

  return c0 + c1 * 2 * x;
}

/**
 * Convert power series to Hermite coefficients.
 */
export function poly2herm(pol: number[]): number[] {
  pol = trimseq(pol);
  const n = pol.length;

  if (n === 0) return [0];
  if (n === 1) return [pol[0]];

  const result = new Array(n).fill(0);

  for (let j = 0; j < n; j++) {
    const xjHerm = _powerToHerm(j);
    for (let k = 0; k < xjHerm.length; k++) {
      if (k < result.length) {
        result[k] += pol[j] * xjHerm[k];
      }
    }
  }

  return trimseq(result);
}

/**
 * Convert Hermite coefficients to power series.
 */
export function herm2poly(c: number[]): number[] {
  c = trimseq(c);
  const n = c.length;

  if (n === 0) return [0];
  if (n === 1) return [c[0]];

  const result = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    const hk = _hermToPowerCoefs(k);
    for (let j = 0; j < hk.length; j++) {
      result[j] += c[k] * hk[j];
    }
  }

  return trimseq(result);
}

export function hermvander(x: number[], deg: number): number[][] {
  return Hermite['_vander'](x, deg);
}

export function hermder(c: number[], m: number = 1): number[] {
  return Hermite['_der'](c, m, 1);
}

export function hermint(
  c: number[],
  m: number = 1,
  k: number[] = [],
  lbnd: number = 0
): number[] {
  const kPadded = [...k];
  while (kPadded.length < m) kPadded.push(0);
  return Hermite['_int'](c, m, kPadded, lbnd, 1);
}

export async function hermfit(
  x: number[],
  y: number[],
  deg: number,
  rcond: number | null = null,
  full: boolean = false,
  w: number[] | null = null
): Promise<number[] | [number[], { residuals: number[]; rank: number; sv: number[]; rcond: number }]> {
  const result = await Hermite.fit(x, y, deg, null, rcond, full, w);

  if (full) {
    const [poly, info] = result as [Hermite, { residuals: number[]; rank: number; sv: number[]; rcond: number }];
    return [poly.coef, info];
  }

  return (result as Hermite).coef;
}

export async function hermroots(c: number[]): Promise<number[]> {
  return Hermite['_roots'](c);
}

export function hermcompanion(c: number[]): number[][] {
  return Hermite['_companion'](c);
}

export function hermfromroots(roots: number[]): number[] {
  return Hermite['_fromroots'](roots);
}

export function hermadd(c1: number[], c2: number[]): number[] {
  return Hermite['_add'](c1, c2);
}

export function hermsub(c1: number[], c2: number[]): number[] {
  return Hermite['_sub'](c1, c2);
}

export function hermmul(c1: number[], c2: number[]): number[] {
  return Hermite['_mul'](c1, c2);
}

export function hermdiv(c1: number[], c2: number[]): [number[], number[]] {
  return Hermite['_div'](c1, c2);
}

export function hermpow(c: number[], pow: number, maxpower: number = 100): number[] {
  return Hermite['_pow'](c, pow, maxpower);
}

/* ============ Helper Functions ============ */

/**
 * Get H_k(x) expressed in power basis.
 */
function _hermToPowerCoefs(k: number): number[] {
  if (k === 0) return [1];
  if (k === 1) return [0, 2];

  // Recurrence: H_{n+1} = 2x*H_n - 2n*H_{n-1}
  let prev2 = [1]; // H_0
  let prev1 = [0, 2]; // H_1

  for (let n = 1; n < k; n++) {
    // 2x * H_n
    const xTimesH = [0, ...prev1.map(c => 2 * c)];
    // -2n * H_{n-1}
    const minusPrev = prev2.map(c => -2 * n * c);

    const result = new Array(Math.max(xTimesH.length, minusPrev.length)).fill(0);
    for (let i = 0; i < xTimesH.length; i++) result[i] += xTimesH[i];
    for (let i = 0; i < minusPrev.length; i++) result[i] += minusPrev[i];

    prev2 = prev1;
    prev1 = trimseq(result);
  }

  return prev1;
}

/**
 * Get x^j expressed in Hermite basis.
 */
function _powerToHerm(j: number): number[] {
  if (j === 0) return [1];
  if (j === 1) return [0, 0.5]; // x = H_1/2

  // Build up by multiplying by x
  // x * H_n = 0.5 * H_{n+1} + n * H_{n-1}
  let result = [0, 0.5]; // x = H_1/2

  for (let n = 2; n <= j; n++) {
    const newResult = new Array(result.length + 1).fill(0);

    for (let k = 0; k < result.length; k++) {
      if (k === 0) {
        // x * c_0 * H_0 = c_0 * H_1 / 2
        newResult[1] += result[0] * 0.5;
      } else {
        // x * c_k * H_k = c_k * (0.5 * H_{k+1} + k * H_{k-1})
        newResult[k + 1] += result[k] * 0.5;
        newResult[k - 1] += result[k] * k;
      }
    }

    result = trimseq(newResult);
  }

  return result;
}
