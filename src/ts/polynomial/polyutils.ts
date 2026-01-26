/**
 * NumJS Polynomial Utilities
 *
 * Utility functions for polynomial operations, shared across all polynomial bases.
 */

import { NDArray } from '../NDArray.js';

/**
 * Error class for polynomial operations.
 */
export class PolyError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'PolyError';
  }
}

/**
 * Warning class for polynomial domain issues.
 */
export class PolyDomainWarning extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'PolyDomainWarning';
  }
}

/**
 * Remove trailing zeros from a sequence.
 *
 * @param seq - Input sequence
 * @returns Sequence with trailing zeros removed (minimum length 1)
 *
 * @example
 * trimseq([1, 2, 0, 0]) // [1, 2]
 * trimseq([0, 0, 0]) // [0]
 */
export function trimseq(seq: number[]): number[] {
  if (seq.length === 0) return [0];

  let end = seq.length;
  while (end > 1 && seq[end - 1] === 0) {
    end--;
  }

  return seq.slice(0, end);
}

/**
 * Remove small trailing coefficients.
 *
 * Trailing elements with absolute value less than or equal to `tol` are removed.
 *
 * @param c - Coefficient array
 * @param tol - Tolerance for considering coefficients as zero (default: 0)
 * @returns Trimmed coefficient array
 *
 * @example
 * trimcoef([1, 2, 1e-10], 1e-9) // [1, 2]
 */
export function trimcoef(c: number[], tol: number = 0): number[] {
  if (c.length === 0) return [0];

  let end = c.length;
  while (end > 1 && Math.abs(c[end - 1]) <= tol) {
    end--;
  }

  return c.slice(0, end);
}

/**
 * Convert input to a list of 1D coefficient arrays.
 *
 * This is used to normalize multiple coefficient array inputs for
 * operations that work with several polynomials.
 *
 * @param alist - List of array-like coefficient sequences
 * @param trim - If true, remove trailing zeros (default: true)
 * @returns List of 1D arrays
 *
 * @example
 * as_series([[1, 2, 0], [3, 4]]) // [[1, 2], [3, 4]]
 */
export function as_series(
  alist: (number[] | NDArray)[],
  trim: boolean = true
): number[][] {
  const result: number[][] = [];

  for (const arr of alist) {
    let coef: number[];

    if (arr instanceof NDArray) {
      if (arr.ndim !== 1) {
        throw new PolyError('Coefficient arrays must be 1-D');
      }
      coef = arr.toArray() as number[];
    } else {
      coef = [...arr];
    }

    if (coef.length === 0) {
      coef = [0];
    }

    if (trim) {
      coef = trimseq(coef);
    }

    result.push(coef);
  }

  return result;
}

/**
 * Get the domain appropriate for given data points.
 *
 * Returns the minimum and maximum of the data as the domain.
 *
 * @param x - Data points
 * @returns [min, max] domain
 *
 * @example
 * getdomain([1, 2, 3, 4, 5]) // [1, 5]
 */
export function getdomain(x: number[] | NDArray): [number, number] {
  const arr = x instanceof NDArray ? (x.toArray() as number[]) : x;

  if (arr.length === 0) {
    return [-1, 1];
  }

  const min = Math.min(...arr);
  const max = Math.max(...arr);

  if (min === max) {
    // Handle single point by expanding domain
    return [min - 1, max + 1];
  }

  return [min, max];
}

/**
 * Get linear mapping parameters between domains.
 *
 * Returns [offset, scale] such that new = offset + scale * old
 *
 * @param old - Original domain [min, max]
 * @param new_ - Target domain [min, max]
 * @returns [offset, scale]
 *
 * @example
 * mapparms([0, 1], [-1, 1]) // [-1, 2]
 * // new = -1 + 2 * old
 */
export function mapparms(
  old: [number, number],
  new_: [number, number]
): [number, number] {
  const [oldMin, oldMax] = old;
  const [newMin, newMax] = new_;

  const oldLen = oldMax - oldMin;
  const newLen = newMax - newMin;

  if (oldLen === 0) {
    throw new PolyError('Old domain has zero length');
  }

  const scl = newLen / oldLen;
  const off = newMin - oldMin * scl;

  return [off, scl];
}

/**
 * Map x from old domain to new domain.
 *
 * The mapping is: new = off + scl * old
 * where off and scl are computed from the domains.
 *
 * @param x - Values to map (number, array, or NDArray)
 * @param old - Original domain [min, max]
 * @param new_ - Target domain [min, max]
 * @returns Mapped values
 *
 * @example
 * mapdomain(0.5, [0, 1], [-1, 1]) // 0
 * mapdomain([0, 0.5, 1], [0, 1], [-1, 1]) // [-1, 0, 1]
 */
export function mapdomain(
  x: number[] | NDArray | number,
  old: [number, number],
  new_: [number, number]
): number[] | number {
  const [off, scl] = mapparms(old, new_);

  if (typeof x === 'number') {
    return off + scl * x;
  }

  const arr = x instanceof NDArray ? (x.toArray() as number[]) : x;
  return arr.map(v => off + scl * v);
}

/* ============ Generic Polynomial Arithmetic ============ */

/**
 * Add two coefficient arrays (power series).
 *
 * @param c1 - First coefficient array
 * @param c2 - Second coefficient array
 * @returns Sum of coefficients
 */
export function polyadd(c1: number[], c2: number[]): number[] {
  const len = Math.max(c1.length, c2.length);
  const result = new Array(len).fill(0);

  for (let i = 0; i < c1.length; i++) {
    result[i] += c1[i];
  }
  for (let i = 0; i < c2.length; i++) {
    result[i] += c2[i];
  }

  return trimseq(result);
}

/**
 * Subtract two coefficient arrays (power series).
 *
 * @param c1 - First coefficient array
 * @param c2 - Second coefficient array
 * @returns Difference of coefficients (c1 - c2)
 */
export function polysub(c1: number[], c2: number[]): number[] {
  const len = Math.max(c1.length, c2.length);
  const result = new Array(len).fill(0);

  for (let i = 0; i < c1.length; i++) {
    result[i] += c1[i];
  }
  for (let i = 0; i < c2.length; i++) {
    result[i] -= c2[i];
  }

  return trimseq(result);
}

/**
 * Multiply two coefficient arrays (power series).
 * Uses convolution.
 *
 * @param c1 - First coefficient array
 * @param c2 - Second coefficient array
 * @returns Product of coefficients
 */
export function polymul(c1: number[], c2: number[]): number[] {
  if (c1.length === 0 || c2.length === 0) {
    return [0];
  }

  const len = c1.length + c2.length - 1;
  const result = new Array(len).fill(0);

  for (let i = 0; i < c1.length; i++) {
    for (let j = 0; j < c2.length; j++) {
      result[i + j] += c1[i] * c2[j];
    }
  }

  return trimseq(result);
}

/**
 * Divide two coefficient arrays (power series).
 *
 * @param c1 - Dividend coefficient array
 * @param c2 - Divisor coefficient array
 * @returns [quotient, remainder]
 */
export function polydiv(c1: number[], c2: number[]): [number[], number[]] {
  c1 = trimseq([...c1]);
  c2 = trimseq([...c2]);

  if (c2.length === 1 && c2[0] === 0) {
    throw new PolyError('Division by zero polynomial');
  }

  if (c1.length < c2.length) {
    return [[0], c1];
  }

  const quotient: number[] = [];
  const remainder = [...c1];

  while (remainder.length >= c2.length) {
    const coef = remainder[remainder.length - 1] / c2[c2.length - 1];
    quotient.unshift(coef);

    for (let i = 0; i < c2.length; i++) {
      remainder[remainder.length - c2.length + i] -= coef * c2[i];
    }
    remainder.pop();
  }

  return [trimseq(quotient), trimseq(remainder)];
}

/**
 * Raise polynomial to a non-negative integer power.
 *
 * @param c - Coefficient array
 * @param pow - Non-negative integer power
 * @param maxpower - Maximum allowed power (default: 100)
 * @returns Coefficients of c**pow
 */
export function polypow(c: number[], pow: number, maxpower: number = 100): number[] {
  if (!Number.isInteger(pow) || pow < 0) {
    throw new PolyError('Power must be a non-negative integer');
  }

  if (pow === 0) {
    return [1];
  }

  if (pow === 1) {
    return trimseq([...c]);
  }

  // Check degree limit
  c = trimseq(c);
  if ((c.length - 1) * pow > maxpower) {
    throw new PolyError(`Power would exceed maxpower (${maxpower})`);
  }

  // Binary exponentiation
  let result = [1];
  let base = [...c];
  let n = pow;

  while (n > 0) {
    if (n & 1) {
      result = polymul(result, base);
    }
    base = polymul(base, base);
    n >>= 1;
  }

  return trimseq(result);
}
