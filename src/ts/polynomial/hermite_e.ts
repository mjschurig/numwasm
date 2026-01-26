/**
 * NumJS Hermite Polynomials (Probabilist's)
 *
 * Probabilist's Hermite polynomials (He_n) are orthogonal on (-inf, inf)
 * with weight function w(x) = exp(-x^2/2).
 *
 * Recurrence: He_{n+1}(x) = x*He_n(x) - n*He_{n-1}(x)
 */

import { ABCPolyBase, maxpower, companionEigenvalues } from './_polybase.js';
import { PolyError, trimseq, polymul as polyMul, polydiv as polyDiv } from './polyutils.js';

/**
 * Probabilist's Hermite polynomial (HermiteE).
 *
 * HermiteE polynomials satisfy the recurrence:
 * He_{n+1}(x) = x*He_n(x) - n*He_{n-1}(x)
 *
 * with He_0(x) = 1, He_1(x) = x
 *
 * These are orthogonal with weight exp(-x^2/2).
 *
 * @example
 * const p = new HermiteE([1, 2, 3]); // 1*He_0 + 2*He_1 + 3*He_2
 */
export class HermiteE extends ABCPolyBase {
  static defaultDomain: [number, number] = [-1, 1];
  static defaultWindow: [number, number] = [-1, 1];
  static basisName = 'He';

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

    const p1 = herme2poly(c1);
    const p2 = herme2poly(c2);
    const product = polyMul(p1, p2);

    return poly2herme(product);
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

    const p1 = herme2poly(c1);
    const p2 = herme2poly(c2);
    const [quo, rem] = polyDiv(p1, p2);

    return [poly2herme(quo), poly2herme(rem)];
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
        result = HermiteE._mul(result, base);
      }
      base = HermiteE._mul(base, base);
      exp >>= 1;
    }

    return trimseq(result);
  }

  protected static _val(x: number | number[], c: number[]): number | number[] {
    if (typeof x === 'number') {
      return hermeval(x, c);
    }
    return x.map(v => hermeval(v, c));
  }

  /**
   * HermiteE derivative.
   *
   * d/dx He_n(x) = n*He_{n-1}(x)
   */
  protected static _der(c: number[], m: number, scl: number): number[] {
    let result = [...c];

    for (let cnt = 0; cnt < m; cnt++) {
      const n = result.length;
      if (n <= 1) return [0];

      const der = new Array(n - 1).fill(0);

      for (let j = 1; j < n; j++) {
        der[j - 1] = j * result[j] * scl;
      }

      result = der;
    }

    return trimseq(result);
  }

  /**
   * HermiteE integral.
   *
   * integral(He_n) = He_{n+1} / (n+1)
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
        integ[j + 1] = result[j] * scl / (j + 1);
      }

      // Add integration constant
      const val = hermeval(lbnd, integ);
      integ[0] = k[cnt] - val;

      result = integ;
    }

    return trimseq(result);
  }

  protected static async _roots(c: number[]): Promise<number[]> {
    c = trimseq(c);

    if (c.length <= 1) return [];
    if (c.length === 2) {
      // Linear: c[0] + c[1]*x = 0
      return [-c[0] / c[1]];
    }

    const companion = HermiteE._companion(c);
    const eigenvalues = await companionEigenvalues(companion);
    return eigenvalues.sort((a, b) => a - b);
  }

  protected static _vander(x: number[], deg: number): number[][] {
    const result: number[][] = [];

    for (const xi of x) {
      const row: number[] = [1]; // He_0 = 1
      if (deg >= 1) {
        row.push(xi); // He_1 = x
      }
      // Recurrence: He_{n+1} = x*He_n - n*He_{n-1}
      for (let j = 2; j <= deg; j++) {
        row.push(xi * row[j - 1] - (j - 1) * row[j - 2]);
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
      mat[i][i + 1] = 1;
      mat[i + 1][i] = i + 1;
    }

    // Modify for polynomial coefficients
    const scale = -1 / c[n];
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

    return poly2herme(poly);
  }
}

/* ============ Module-Level Functions ============ */

/**
 * Evaluate HermiteE polynomial.
 */
export function hermeval(x: number, c: number[]): number {
  if (c.length === 0) return 0;
  if (c.length === 1) return c[0];
  if (c.length === 2) return c[0] + c[1] * x;

  // Clenshaw-like algorithm for HermiteE
  let c0 = c[c.length - 2];
  let c1 = c[c.length - 1];

  for (let i = c.length - 3; i >= 0; i--) {
    const n = i + 1;
    const tmp = c0;
    c0 = c[i] - c1 * n;
    c1 = tmp + c1 * x;
  }

  return c0 + c1 * x;
}

/**
 * Convert power series to HermiteE coefficients.
 */
export function poly2herme(pol: number[]): number[] {
  pol = trimseq(pol);
  const n = pol.length;

  if (n === 0) return [0];
  if (n === 1) return [pol[0]];

  const result = new Array(n).fill(0);

  for (let j = 0; j < n; j++) {
    const xjHerme = _powerToHerme(j);
    for (let k = 0; k < xjHerme.length; k++) {
      if (k < result.length) {
        result[k] += pol[j] * xjHerme[k];
      }
    }
  }

  return trimseq(result);
}

/**
 * Convert HermiteE coefficients to power series.
 */
export function herme2poly(c: number[]): number[] {
  c = trimseq(c);
  const n = c.length;

  if (n === 0) return [0];
  if (n === 1) return [c[0]];

  const result = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    const hek = _hermeToPowerCoefs(k);
    for (let j = 0; j < hek.length; j++) {
      result[j] += c[k] * hek[j];
    }
  }

  return trimseq(result);
}

export function hermevander(x: number[], deg: number): number[][] {
  return HermiteE['_vander'](x, deg);
}

export function hermeder(c: number[], m: number = 1): number[] {
  return HermiteE['_der'](c, m, 1);
}

export function hermeint(
  c: number[],
  m: number = 1,
  k: number[] = [],
  lbnd: number = 0
): number[] {
  const kPadded = [...k];
  while (kPadded.length < m) kPadded.push(0);
  return HermiteE['_int'](c, m, kPadded, lbnd, 1);
}

export async function hermefit(
  x: number[],
  y: number[],
  deg: number,
  rcond: number | null = null,
  full: boolean = false,
  w: number[] | null = null
): Promise<number[] | [number[], { residuals: number[]; rank: number; sv: number[]; rcond: number }]> {
  const result = await HermiteE.fit(x, y, deg, null, rcond, full, w);

  if (full) {
    const [poly, info] = result as [HermiteE, { residuals: number[]; rank: number; sv: number[]; rcond: number }];
    return [poly.coef, info];
  }

  return (result as HermiteE).coef;
}

export async function hermeroots(c: number[]): Promise<number[]> {
  return HermiteE['_roots'](c);
}

export function hermecompanion(c: number[]): number[][] {
  return HermiteE['_companion'](c);
}

export function hermefromroots(roots: number[]): number[] {
  return HermiteE['_fromroots'](roots);
}

export function hermeadd(c1: number[], c2: number[]): number[] {
  return HermiteE['_add'](c1, c2);
}

export function hermesub(c1: number[], c2: number[]): number[] {
  return HermiteE['_sub'](c1, c2);
}

export function hermemul(c1: number[], c2: number[]): number[] {
  return HermiteE['_mul'](c1, c2);
}

export function hermediv(c1: number[], c2: number[]): [number[], number[]] {
  return HermiteE['_div'](c1, c2);
}

export function hermepow(c: number[], pow: number, maxpower: number = 100): number[] {
  return HermiteE['_pow'](c, pow, maxpower);
}

/* ============ Helper Functions ============ */

/**
 * Get He_k(x) expressed in power basis.
 */
function _hermeToPowerCoefs(k: number): number[] {
  if (k === 0) return [1];
  if (k === 1) return [0, 1];

  // Recurrence: He_{n+1} = x*He_n - n*He_{n-1}
  let prev2 = [1]; // He_0
  let prev1 = [0, 1]; // He_1

  for (let n = 1; n < k; n++) {
    // x * He_n
    const xTimesHe = [0, ...prev1];
    // -n * He_{n-1}
    const minusPrev = prev2.map(c => -n * c);

    const result = new Array(Math.max(xTimesHe.length, minusPrev.length)).fill(0);
    for (let i = 0; i < xTimesHe.length; i++) result[i] += xTimesHe[i];
    for (let i = 0; i < minusPrev.length; i++) result[i] += minusPrev[i];

    prev2 = prev1;
    prev1 = trimseq(result);
  }

  return prev1;
}

/**
 * Get x^j expressed in HermiteE basis.
 */
function _powerToHerme(j: number): number[] {
  if (j === 0) return [1];
  if (j === 1) return [0, 1]; // x = He_1

  // Build up by multiplying by x
  // x * He_n = He_{n+1} + n * He_{n-1}
  let result = [0, 1]; // x = He_1

  for (let n = 2; n <= j; n++) {
    const newResult = new Array(result.length + 1).fill(0);

    for (let k = 0; k < result.length; k++) {
      if (k === 0) {
        // x * c_0 * He_0 = c_0 * He_1
        newResult[1] += result[0];
      } else {
        // x * c_k * He_k = c_k * (He_{k+1} + k * He_{k-1})
        newResult[k + 1] += result[k];
        newResult[k - 1] += result[k] * k;
      }
    }

    result = trimseq(newResult);
  }

  return result;
}
