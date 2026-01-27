/**
 * NumJS Laguerre Polynomials
 *
 * Laguerre polynomials (L_n) are orthogonal on [0, inf)
 * with weight function w(x) = exp(-x).
 *
 * Recurrence: (n+1)*L_{n+1}(x) = (2n+1-x)*L_n(x) - n*L_{n-1}(x)
 */

import { ABCPolyBase, companionEigenvalues } from './_polybase.js';
import { PolyError, trimseq, polymul as polyMul, polydiv as polyDiv } from './polyutils.js';

/**
 * Laguerre polynomial.
 *
 * Laguerre polynomials satisfy the recurrence:
 * (n+1)*L_{n+1}(x) = (2n+1-x)*L_n(x) - n*L_{n-1}(x)
 *
 * with L_0(x) = 1, L_1(x) = 1-x
 *
 * These are orthogonal on [0, inf) with weight exp(-x).
 *
 * @example
 * const p = new Laguerre([1, 2, 3]); // 1*L_0 + 2*L_1 + 3*L_2
 */
export class Laguerre extends ABCPolyBase {
  static defaultDomain: [number, number] = [0, 1];
  static defaultWindow: [number, number] = [0, 1];
  static basisName = 'L';

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

    const p1 = lag2poly(c1);
    const p2 = lag2poly(c2);
    const product = polyMul(p1, p2);

    return poly2lag(product);
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

    const p1 = lag2poly(c1);
    const p2 = lag2poly(c2);
    const [quo, rem] = polyDiv(p1, p2);

    return [poly2lag(quo), poly2lag(rem)];
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
        result = Laguerre._mul(result, base);
      }
      base = Laguerre._mul(base, base);
      exp >>= 1;
    }

    return trimseq(result);
  }

  protected static _val(x: number | number[], c: number[]): number | number[] {
    if (typeof x === 'number') {
      return lagval(x, c);
    }
    return x.map(v => lagval(v, c));
  }

  /**
   * Laguerre derivative.
   *
   * d/dx L_n(x) = -L'_{n-1}(x) where L'_n is the associated Laguerre polynomial
   * For standard Laguerre: d/dx L_n(x) = sum_{k=0}^{n-1} -L_k(x)
   */
  protected static _der(c: number[], m: number, scl: number): number[] {
    let result = [...c];

    for (let cnt = 0; cnt < m; cnt++) {
      const n = result.length;
      if (n <= 1) return [0];

      const der = new Array(n - 1).fill(0);

      // Derivative formula: c'_k = -sum_{j=k+1}^{n-1} c_j
      for (let k = n - 2; k >= 0; k--) {
        der[k] = -result[k + 1] * scl;
        if (k < n - 2) {
          der[k] += der[k + 1];
        }
      }

      result = der;
    }

    return trimseq(result);
  }

  /**
   * Laguerre integral.
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

      // Integration: integral(L_n) = L_n - L_{n+1}
      for (let j = 0; j < n; j++) {
        integ[j] += result[j] * scl;
        integ[j + 1] -= result[j] * scl;
      }

      // Add integration constant
      const val = lagval(lbnd, integ);
      integ[0] += k[cnt] - val;

      result = integ;
    }

    return trimseq(result);
  }

  protected static async _roots(c: number[]): Promise<number[]> {
    c = trimseq(c);

    if (c.length <= 1) return [];
    if (c.length === 2) {
      // Linear: c[0] + c[1]*(1-x) = 0
      // c[0] + c[1] - c[1]*x = 0
      // x = (c[0] + c[1]) / c[1]
      return [(c[0] + c[1]) / c[1]];
    }

    const companion = Laguerre._companion(c);
    const eigenvalues = await companionEigenvalues(companion);
    return eigenvalues.filter(r => r >= 0).sort((a, b) => a - b);
  }

  protected static _vander(x: number[], deg: number): number[][] {
    const result: number[][] = [];

    for (const xi of x) {
      const row: number[] = [1]; // L_0 = 1
      if (deg >= 1) {
        row.push(1 - xi); // L_1 = 1 - x
      }
      // Recurrence: (n+1)*L_{n+1} = (2n+1-x)*L_n - n*L_{n-1}
      for (let j = 2; j <= deg; j++) {
        const ln = ((2 * j - 1 - xi) * row[j - 1] - (j - 1) * row[j - 2]) / j;
        row.push(ln);
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

    // Build companion matrix based on recurrence
    // Diagonal entries: 2i + 1
    // Off-diagonal: -(i+1), -i
    for (let i = 0; i < n; i++) {
      mat[i][i] = 2 * i + 1;
      if (i < n - 1) {
        mat[i][i + 1] = -(i + 1);
        mat[i + 1][i] = -(i + 1);
      }
    }

    // Modify for polynomial coefficients
    const scale = 1 / c[n];
    for (let i = 0; i < n; i++) {
      mat[i][n - 1] += scale * c[i] * (n);
    }

    return mat;
  }

  protected static _fromroots(roots: number[]): number[] {
    if (roots.length === 0) return [1];

    let poly = [1];
    for (const r of roots) {
      poly = polyMul(poly, [-r, 1]);
    }

    return poly2lag(poly);
  }
}

/* ============ Module-Level Functions ============ */

/**
 * Evaluate Laguerre polynomial.
 */
export function lagval(x: number, c: number[]): number {
  if (c.length === 0) return 0;
  if (c.length === 1) return c[0];
  if (c.length === 2) return c[0] + c[1] * (1 - x);

  // Clenshaw-like algorithm for Laguerre
  let c0 = c[c.length - 2];
  let c1 = c[c.length - 1];

  for (let i = c.length - 3; i >= 0; i--) {
    const n = i + 1;
    const tmp = c0;
    c0 = c[i] - c1 * n / (n + 1);
    c1 = tmp + c1 * (2 * n + 1 - x) / (n + 1);
  }

  return c0 + c1 * (1 - x);
}

/**
 * Convert power series to Laguerre coefficients.
 */
export function poly2lag(pol: number[]): number[] {
  pol = trimseq(pol);
  const n = pol.length;

  if (n === 0) return [0];
  if (n === 1) return [pol[0]];

  const result = new Array(n).fill(0);

  for (let j = 0; j < n; j++) {
    const xjLag = _powerToLag(j);
    for (let k = 0; k < xjLag.length; k++) {
      if (k < result.length) {
        result[k] += pol[j] * xjLag[k];
      }
    }
  }

  return trimseq(result);
}

/**
 * Convert Laguerre coefficients to power series.
 */
export function lag2poly(c: number[]): number[] {
  c = trimseq(c);
  const n = c.length;

  if (n === 0) return [0];
  if (n === 1) return [c[0]];

  const result = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    const lk = _lagToPowerCoefs(k);
    for (let j = 0; j < lk.length; j++) {
      result[j] += c[k] * lk[j];
    }
  }

  return trimseq(result);
}

export function lagvander(x: number[], deg: number): number[][] {
  return Laguerre['_vander'](x, deg);
}

export function lagder(c: number[], m: number = 1): number[] {
  return Laguerre['_der'](c, m, 1);
}

export function lagint(
  c: number[],
  m: number = 1,
  k: number[] = [],
  lbnd: number = 0
): number[] {
  const kPadded = [...k];
  while (kPadded.length < m) kPadded.push(0);
  return Laguerre['_int'](c, m, kPadded, lbnd, 1);
}

export async function lagfit(
  x: number[],
  y: number[],
  deg: number,
  rcond: number | null = null,
  full: boolean = false,
  w: number[] | null = null
): Promise<number[] | [number[], { residuals: number[]; rank: number; sv: number[]; rcond: number }]> {
  const result = await Laguerre.fit(x, y, deg, null, rcond, full, w);

  if (full) {
    const [poly, info] = result as [Laguerre, { residuals: number[]; rank: number; sv: number[]; rcond: number }];
    return [poly.coef, info];
  }

  return (result as Laguerre).coef;
}

export async function lagroots(c: number[]): Promise<number[]> {
  return Laguerre['_roots'](c);
}

export function lagcompanion(c: number[]): number[][] {
  return Laguerre['_companion'](c);
}

export function lagfromroots(roots: number[]): number[] {
  return Laguerre['_fromroots'](roots);
}

export function lagadd(c1: number[], c2: number[]): number[] {
  return Laguerre['_add'](c1, c2);
}

export function lagsub(c1: number[], c2: number[]): number[] {
  return Laguerre['_sub'](c1, c2);
}

export function lagmul(c1: number[], c2: number[]): number[] {
  return Laguerre['_mul'](c1, c2);
}

export function lagdiv(c1: number[], c2: number[]): [number[], number[]] {
  return Laguerre['_div'](c1, c2);
}

export function lagpow(c: number[], pow: number, maxpower: number = 100): number[] {
  return Laguerre['_pow'](c, pow, maxpower);
}

/* ============ Helper Functions ============ */

/**
 * Get L_k(x) expressed in power basis.
 * L_n(x) = sum_{k=0}^n C(n,k) * (-x)^k / k!
 */
function _lagToPowerCoefs(k: number): number[] {
  if (k === 0) return [1];
  if (k === 1) return [1, -1];

  // Recurrence: (n+1)*L_{n+1} = (2n+1-x)*L_n - n*L_{n-1}
  let prev2 = [1]; // L_0
  let prev1 = [1, -1]; // L_1

  for (let n = 1; n < k; n++) {
    // (2n+1)*L_n
    const scaleL = prev1.map(c => ((2 * n + 1) * c) / (n + 1));
    // -x*L_n
    const xTimesL = [0, ...prev1.map(c => (-c) / (n + 1))];
    // -n*L_{n-1}
    const minusPrev = prev2.map(c => (-n * c) / (n + 1));

    const len = Math.max(scaleL.length, xTimesL.length, minusPrev.length);
    const result = new Array(len).fill(0);
    for (let i = 0; i < scaleL.length; i++) result[i] += scaleL[i];
    for (let i = 0; i < xTimesL.length; i++) result[i] += xTimesL[i];
    for (let i = 0; i < minusPrev.length; i++) result[i] += minusPrev[i];

    prev2 = prev1;
    prev1 = trimseq(result);
  }

  return prev1;
}

/**
 * Get x^j expressed in Laguerre basis.
 */
function _powerToLag(j: number): number[] {
  if (j === 0) return [1];
  if (j === 1) return [1, -1]; // x = L_0 - L_1 = 1 - (1-x) = x (wait, L_1 = 1-x, so x = 1 - L_1 = L_0 - L_1)

  // x = -L_1 + L_0 means coefficients are [1, -1]
  // Build up by using the recurrence
  // x*L_n = -(n+1)*L_{n+1} + (2n+1)*L_n - n*L_{n-1}
  let result = [1, -1]; // x

  for (let n = 2; n <= j; n++) {
    const newResult = new Array(result.length + 1).fill(0);

    for (let k = 0; k < result.length; k++) {
      // x * c_k * L_k = c_k * [-(k+1)*L_{k+1} + (2k+1)*L_k - k*L_{k-1}]
      newResult[k + 1] += result[k] * (-(k + 1));
      newResult[k] += result[k] * (2 * k + 1);
      if (k > 0) {
        newResult[k - 1] += result[k] * (-k);
      }
    }

    result = trimseq(newResult);
  }

  return result;
}
