# Phase 16c: numpy.polynomial Implementation Plan

Complete implementation roadmap for the NumJS-WASM polynomial module, providing NumPy-compatible polynomial classes for multiple basis systems.

---

## ⚠️ Implementation Guideline

**The original NumPy source code is available in `/numpy`.**

Key reference files:
- `numpy/polynomial/_polybase.py` - Abstract base class (1,191 lines)
- `numpy/polynomial/polyutils.py` - Utility functions (759 lines)
- `numpy/polynomial/polynomial.py` - Power series (1,625 lines)
- `numpy/polynomial/chebyshev.py` - Chebyshev polynomials (2,001 lines)
- `numpy/polynomial/legendre.py` - Legendre polynomials (1,603 lines)
- `numpy/polynomial/hermite.py` - Hermite polynomials (1,738 lines)
- `numpy/polynomial/hermite_e.py` - Probabilist's Hermite (1,640 lines)
- `numpy/polynomial/laguerre.py` - Laguerre polynomials (1,673 lines)

Implementation should follow NumPy's algorithms and error handling for consistency.

---

## Current State (Pre-Phase 16c)

```
src/ts/
├── NDArray.ts         # Core array class
├── types.ts           # Type definitions
├── broadcast.ts       # Broadcasting functions
└── index.ts           # Public exports
```

**Dependencies Required:**
- NDArray core
- Broadcasting
- Linear algebra (Phase 13) for fitting and root finding

---

## Phase 16c Dependency Tree

```
PHASE 16c: NUMPY.POLYNOMIAL
│
├── 16c.1 Utility Functions (polyutils)
│   ├── as_series(alist) → normalize to 1D arrays
│   ├── trimseq(seq) → remove trailing zeros
│   ├── trimcoef(c, tol) → trim small coefficients
│   ├── getdomain(x) → get appropriate domain
│   ├── mapdomain(x, old, new) → map between domains
│   ├── mapparms(old, new) → get mapping parameters
│   ├── PolyError → exception class
│   └── PolyDomainWarning → warning class
│
├── 16c.2 Abstract Base Class (ABCPolyBase)
│   ├── Properties
│   │   ├── coef → coefficient array
│   │   ├── domain → polynomial domain
│   │   ├── window → evaluation window
│   │   ├── degree → polynomial degree
│   │   └── symbol → variable symbol
│   │
│   ├── Arithmetic
│   │   ├── __add__, __radd__
│   │   ├── __sub__, __rsub__
│   │   ├── __mul__, __rmul__
│   │   ├── __truediv__
│   │   ├── __floordiv__, __mod__, __divmod__
│   │   ├── __pow__
│   │   ├── __neg__, __pos__
│   │   └── __call__ (evaluation)
│   │
│   ├── Calculus
│   │   ├── deriv(m) → m-th derivative
│   │   └── integ(m, k, lbnd) → m-th integral
│   │
│   ├── Analysis
│   │   ├── roots() → polynomial roots
│   │   ├── linspace(n) → evaluate at n points
│   │   ├── trim(tol) → trim small coefficients
│   │   ├── truncate(size) → truncate to size
│   │   └── cutdeg(deg) → cut to degree
│   │
│   ├── Conversion
│   │   ├── convert(kind, domain, window)
│   │   ├── cast(kind)
│   │   └── mapparms() → domain mapping params
│   │
│   ├── Class Methods
│   │   ├── basis(deg) → basis polynomial
│   │   ├── identity() → p(x) = x
│   │   ├── fromroots(roots) → from roots
│   │   ├── fit(x, y, deg) → least squares fit
│   │   └── copy()
│   │
│   └── Abstract Methods (implemented per basis)
│       ├── _add(c1, c2)
│       ├── _sub(c1, c2)
│       ├── _mul(c1, c2)
│       ├── _div(c1, c2)
│       ├── _pow(c, n)
│       ├── _val(x, c)
│       ├── _der(c, m, scl)
│       ├── _int(c, m, k, lbnd, scl)
│       └── _roots(c)
│
├── 16c.3 Polynomial Class (Power Series)
│   ├── polyadd(c1, c2), polysub, polymul, polydiv
│   ├── polypow(c, pow)
│   ├── polyval(x, c), polyval2d, polyval3d
│   ├── polyvander(x, deg), polyvander2d, polyvander3d
│   ├── polyder(c, m), polyint(c, m, k, lbnd)
│   ├── polyfit(x, y, deg)
│   ├── polyroots(c), polycompanion(c)
│   └── polyfromroots(roots)
│
├── 16c.4 Chebyshev Class
│   ├── chebadd, chebsub, chebmul, chebdiv
│   ├── chebpow(c, pow)
│   ├── chebval(x, c), chebval2d, chebval3d
│   ├── chebvander(x, deg)
│   ├── chebder(c, m), chebint(c, m, k, lbnd)
│   ├── chebfit(x, y, deg)
│   ├── chebroots(c), chebcompanion(c)
│   ├── chebfromroots(roots)
│   ├── chebinterpolate(func, deg)
│   ├── poly2cheb(pol) → convert from power series
│   └── cheb2poly(c) → convert to power series
│
├── 16c.5 Legendre Class
│   └── (Same pattern as Chebyshev with leg* prefix)
│
├── 16c.6 Hermite Class (Physicist's)
│   └── (Same pattern as Chebyshev with herm* prefix)
│
├── 16c.7 HermiteE Class (Probabilist's)
│   └── (Same pattern as Chebyshev with herme* prefix)
│
└── 16c.8 Laguerre Class
    └── (Same pattern as Chebyshev with lag* prefix)

Dependencies: NDArray, Broadcasting, Linear Algebra (Phase 13)
```

---

## Detailed Implementation Specifications

### 16c.1 Utility Functions

**File:** `src/ts/polynomial/polyutils.ts`

```typescript
import { NDArray } from '../NDArray.js';
import { DType } from '../types.js';

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
 * Convert input to a list of 1D coefficient arrays.
 *
 * @param alist - List of array-like coefficient sequences
 * @param trim - If true, remove trailing zeros
 * @returns List of 1D arrays
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
 * Remove trailing zeros from a sequence.
 *
 * @param seq - Input sequence
 * @returns Sequence with trailing zeros removed (minimum length 1)
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
 * @param c - Coefficient array
 * @param tol - Tolerance for considering coefficients as zero
 * @returns Trimmed coefficient array
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
 * Get the domain appropriate for given data points.
 *
 * @param x - Data points
 * @returns [min, max] domain
 */
export function getdomain(x: number[] | NDArray): [number, number] {
  const arr = x instanceof NDArray ? (x.toArray() as number[]) : x;

  if (arr.length === 0) {
    return [-1, 1];
  }

  const min = Math.min(...arr);
  const max = Math.max(...arr);

  if (min === max) {
    // Handle single point
    return [min - 1, max + 1];
  }

  return [min, max];
}

/**
 * Map x from old domain to new domain.
 *
 * The mapping is: new = off + scl * old
 * where off and scl are computed from the domains.
 *
 * @param x - Values to map
 * @param old - Original domain [min, max]
 * @param new_ - Target domain [min, max]
 * @returns Mapped values
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

/**
 * Get linear mapping parameters between domains.
 *
 * Returns [offset, scale] such that new = offset + scale * old
 *
 * @param old - Original domain [min, max]
 * @param new_ - Target domain [min, max]
 * @returns [offset, scale]
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
 * Add two coefficient arrays.
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
 * Subtract two coefficient arrays.
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
 */
export function polymul(c1: number[], c2: number[]): number[] {
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
 * @returns [quotient, remainder]
 */
export function polydiv(c1: number[], c2: number[]): [number[], number[]] {
  c1 = trimseq(c1);
  c2 = trimseq(c2);

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
```

---

### 16c.2 Abstract Base Class

**File:** `src/ts/polynomial/_polybase.ts`

```typescript
import { NDArray } from '../NDArray.js';
import { DType } from '../types.js';
import { PolyError, trimcoef, mapparms, mapdomain } from './polyutils.js';

/**
 * Maximum allowed polynomial degree.
 */
export const maxpower = 100;

/**
 * Abstract base class for polynomial series.
 *
 * This class provides the common interface and operations for all
 * polynomial types (Polynomial, Chebyshev, Legendre, etc.).
 *
 * Subclasses must implement the abstract static methods for
 * basis-specific arithmetic.
 */
export abstract class ABCPolyBase {
  /** Coefficient array (lowest degree first) */
  protected _coef: number[];

  /** Domain for the polynomial */
  protected _domain: [number, number];

  /** Window for evaluation */
  protected _window: [number, number];

  /** Symbol used in string representation */
  protected _symbol: string;

  /* ============ Abstract Properties ============ */

  /** Default domain for this polynomial type */
  static defaultDomain: [number, number] = [-1, 1];

  /** Default window for this polynomial type */
  static defaultWindow: [number, number] = [-1, 1];

  /** Name of the basis (for display) */
  static basisName: string = 'P';

  /* ============ Abstract Methods ============ */

  /**
   * Add two coefficient arrays in this basis.
   */
  protected static _add(c1: number[], c2: number[]): number[] {
    throw new Error('Must be implemented by subclass');
  }

  /**
   * Subtract two coefficient arrays in this basis.
   */
  protected static _sub(c1: number[], c2: number[]): number[] {
    throw new Error('Must be implemented by subclass');
  }

  /**
   * Multiply two coefficient arrays in this basis.
   */
  protected static _mul(c1: number[], c2: number[]): number[] {
    throw new Error('Must be implemented by subclass');
  }

  /**
   * Divide two coefficient arrays in this basis.
   * @returns [quotient, remainder]
   */
  protected static _div(c1: number[], c2: number[]): [number[], number[]] {
    throw new Error('Must be implemented by subclass');
  }

  /**
   * Raise coefficient array to integer power.
   */
  protected static _pow(c: number[], n: number, maxpower: number): number[] {
    throw new Error('Must be implemented by subclass');
  }

  /**
   * Evaluate polynomial at x.
   */
  protected static _val(x: number | number[], c: number[]): number | number[] {
    throw new Error('Must be implemented by subclass');
  }

  /**
   * Compute derivative of coefficient array.
   */
  protected static _der(c: number[], m: number, scl: number): number[] {
    throw new Error('Must be implemented by subclass');
  }

  /**
   * Compute integral of coefficient array.
   */
  protected static _int(
    c: number[],
    m: number,
    k: number[],
    lbnd: number,
    scl: number
  ): number[] {
    throw new Error('Must be implemented by subclass');
  }

  /**
   * Find roots of polynomial.
   */
  protected static _roots(c: number[]): number[] {
    throw new Error('Must be implemented by subclass');
  }

  /**
   * Generate Vandermonde matrix.
   */
  protected static _vander(x: number[], deg: number): number[][] {
    throw new Error('Must be implemented by subclass');
  }

  /**
   * Generate companion matrix for root finding.
   */
  protected static _companion(c: number[]): number[][] {
    throw new Error('Must be implemented by subclass');
  }

  /**
   * Construct polynomial from roots.
   */
  protected static _fromroots(roots: number[]): number[] {
    throw new Error('Must be implemented by subclass');
  }

  /* ============ Constructor ============ */

  /**
   * Create a polynomial.
   *
   * @param coef - Coefficients (lowest degree first)
   * @param domain - Domain for the polynomial
   * @param window - Window for evaluation
   * @param symbol - Variable symbol for display
   */
  constructor(
    coef: number[] | NDArray,
    domain: [number, number] | null = null,
    window: [number, number] | null = null,
    symbol: string = 'x'
  ) {
    // Convert and validate coefficients
    if (coef instanceof NDArray) {
      if (coef.ndim !== 1) {
        throw new PolyError('Coefficient array must be 1-D');
      }
      this._coef = coef.toArray() as number[];
    } else {
      this._coef = [...coef];
    }

    // Trim trailing zeros
    this._coef = trimcoef(this._coef);

    if (this._coef.length === 0) {
      this._coef = [0];
    }

    // Check degree limit
    if (this._coef.length > maxpower + 1) {
      throw new PolyError(`Polynomial degree exceeds maxpower (${maxpower})`);
    }

    // Set domain and window
    const ctor = this.constructor as typeof ABCPolyBase;
    this._domain = domain ?? [...ctor.defaultDomain] as [number, number];
    this._window = window ?? [...ctor.defaultWindow] as [number, number];
    this._symbol = symbol;
  }

  /* ============ Properties ============ */

  /**
   * Coefficient array (copy).
   */
  get coef(): number[] {
    return [...this._coef];
  }

  /**
   * Domain of the polynomial.
   */
  get domain(): [number, number] {
    return [...this._domain] as [number, number];
  }

  /**
   * Window of the polynomial.
   */
  get window(): [number, number] {
    return [...this._window] as [number, number];
  }

  /**
   * Polynomial degree.
   */
  get degree(): number {
    return this._coef.length - 1;
  }

  /**
   * Variable symbol.
   */
  get symbol(): string {
    return this._symbol;
  }

  /* ============ Arithmetic Methods ============ */

  /**
   * Add another polynomial or scalar.
   */
  add(other: ABCPolyBase | number): this {
    const ctor = this.constructor as typeof ABCPolyBase;

    if (typeof other === 'number') {
      const newCoef = [...this._coef];
      newCoef[0] = (newCoef[0] || 0) + other;
      return this._createNew(newCoef);
    }

    this._checkCompatible(other);
    const newCoef = ctor._add(this._coef, other._coef);
    return this._createNew(newCoef);
  }

  /**
   * Subtract another polynomial or scalar.
   */
  sub(other: ABCPolyBase | number): this {
    const ctor = this.constructor as typeof ABCPolyBase;

    if (typeof other === 'number') {
      const newCoef = [...this._coef];
      newCoef[0] = (newCoef[0] || 0) - other;
      return this._createNew(newCoef);
    }

    this._checkCompatible(other);
    const newCoef = ctor._sub(this._coef, other._coef);
    return this._createNew(newCoef);
  }

  /**
   * Multiply by another polynomial or scalar.
   */
  mul(other: ABCPolyBase | number): this {
    const ctor = this.constructor as typeof ABCPolyBase;

    if (typeof other === 'number') {
      const newCoef = this._coef.map(c => c * other);
      return this._createNew(newCoef);
    }

    this._checkCompatible(other);
    const newCoef = ctor._mul(this._coef, other._coef);
    return this._createNew(newCoef);
  }

  /**
   * Divide by another polynomial.
   * @returns [quotient, remainder]
   */
  divmod(other: ABCPolyBase): [this, this] {
    const ctor = this.constructor as typeof ABCPolyBase;

    if (typeof other === 'number') {
      throw new PolyError('Cannot divmod by scalar');
    }

    this._checkCompatible(other);
    const [quo, rem] = ctor._div(this._coef, other._coef);
    return [this._createNew(quo), this._createNew(rem)];
  }

  /**
   * Floor division.
   */
  floordiv(other: ABCPolyBase): this {
    return this.divmod(other)[0];
  }

  /**
   * Modulo.
   */
  mod(other: ABCPolyBase): this {
    return this.divmod(other)[1];
  }

  /**
   * Raise to integer power.
   */
  pow(n: number): this {
    if (!Number.isInteger(n) || n < 0) {
      throw new PolyError('Power must be a non-negative integer');
    }

    const ctor = this.constructor as typeof ABCPolyBase;
    const newCoef = ctor._pow(this._coef, n, maxpower);
    return this._createNew(newCoef);
  }

  /**
   * Negate polynomial.
   */
  neg(): this {
    return this._createNew(this._coef.map(c => -c));
  }

  /* ============ Evaluation ============ */

  /**
   * Evaluate polynomial at x.
   *
   * Maps x from domain to window before evaluation.
   */
  call(x: number | number[] | NDArray): number | number[] {
    const ctor = this.constructor as typeof ABCPolyBase;

    // Map x from domain to window
    const xMapped = mapdomain(
      x instanceof NDArray ? (x.toArray() as number[]) : x,
      this._domain,
      this._window
    );

    return ctor._val(xMapped, this._coef);
  }

  /**
   * Evaluate at n equally-spaced points in domain.
   */
  linspace(n: number = 100): { x: number[]; y: number[] } {
    const [d0, d1] = this._domain;
    const x: number[] = [];
    const step = (d1 - d0) / (n - 1);

    for (let i = 0; i < n; i++) {
      x.push(d0 + i * step);
    }

    const y = this.call(x) as number[];
    return { x, y };
  }

  /* ============ Calculus ============ */

  /**
   * Differentiate polynomial m times.
   */
  deriv(m: number = 1): this {
    if (!Number.isInteger(m) || m < 0) {
      throw new PolyError('Derivative order must be a non-negative integer');
    }

    if (m === 0) {
      return this.copy();
    }

    const ctor = this.constructor as typeof ABCPolyBase;
    const [, scl] = mapparms(this._domain, this._window);
    const newCoef = ctor._der(this._coef, m, scl);

    return this._createNew(newCoef);
  }

  /**
   * Integrate polynomial m times.
   *
   * @param m - Integration order
   * @param k - Integration constants (one per integration)
   * @param lbnd - Lower bound for definite integral
   */
  integ(m: number = 1, k: number[] = [], lbnd: number | null = null): this {
    if (!Number.isInteger(m) || m < 0) {
      throw new PolyError('Integration order must be a non-negative integer');
    }

    if (m === 0) {
      return this.copy();
    }

    const ctor = this.constructor as typeof ABCPolyBase;
    const [, scl] = mapparms(this._domain, this._window);
    const effectiveLbnd = lbnd ?? this._domain[0];

    // Pad k with zeros
    const kPadded = [...k];
    while (kPadded.length < m) {
      kPadded.push(0);
    }

    const newCoef = ctor._int(this._coef, m, kPadded, effectiveLbnd, 1 / scl);

    return this._createNew(newCoef);
  }

  /* ============ Analysis ============ */

  /**
   * Find the roots of the polynomial.
   */
  roots(): number[] {
    const ctor = this.constructor as typeof ABCPolyBase;
    const windowRoots = ctor._roots(this._coef);

    // Map roots from window to domain
    return mapdomain(windowRoots, this._window, this._domain) as number[];
  }

  /**
   * Trim small trailing coefficients.
   */
  trim(tol: number = 0): this {
    return this._createNew(trimcoef(this._coef, tol));
  }

  /**
   * Truncate to number of coefficients.
   */
  truncate(size: number): this {
    if (size < 1) {
      throw new PolyError('Size must be at least 1');
    }
    return this._createNew(this._coef.slice(0, size));
  }

  /**
   * Truncate to maximum degree.
   */
  cutdeg(deg: number): this {
    return this.truncate(deg + 1);
  }

  /* ============ Conversion ============ */

  /**
   * Get mapping parameters from domain to window.
   */
  mapparms(): [number, number] {
    return mapparms(this._domain, this._window);
  }

  /**
   * Create a copy.
   */
  copy(): this {
    return this._createNew([...this._coef]);
  }

  /* ============ Class Methods ============ */

  /**
   * Create a basis polynomial of given degree.
   *
   * @example
   * Polynomial.basis(2)  // x^2 for power series
   */
  static basis<T extends ABCPolyBase>(
    this: new (coef: number[], domain?: [number, number], window?: [number, number]) => T,
    deg: number,
    domain: [number, number] | null = null,
    window: [number, number] | null = null
  ): T {
    if (!Number.isInteger(deg) || deg < 0) {
      throw new PolyError('Degree must be a non-negative integer');
    }

    const coef = new Array(deg + 1).fill(0);
    coef[deg] = 1;

    return new this(coef, domain ?? undefined, window ?? undefined);
  }

  /**
   * Create identity polynomial p(x) = x.
   */
  static identity<T extends ABCPolyBase>(
    this: new (coef: number[], domain?: [number, number], window?: [number, number]) => T,
    domain: [number, number] | null = null,
    window: [number, number] | null = null
  ): T {
    return new this([0, 1], domain ?? undefined, window ?? undefined);
  }

  /**
   * Create polynomial from roots.
   */
  static fromroots<T extends ABCPolyBase>(
    this: { new (coef: number[], domain?: [number, number], window?: [number, number]): T; _fromroots(roots: number[]): number[] },
    roots: number[],
    domain: [number, number] | null = null,
    window: [number, number] | null = null
  ): T {
    const coef = this._fromroots(roots);
    return new this(coef, domain ?? undefined, window ?? undefined);
  }

  /**
   * Least squares fit to data.
   *
   * @param x - x values
   * @param y - y values
   * @param deg - Degree of fit
   * @param domain - Domain for fit
   * @param rcond - Cutoff for small singular values
   * @param full - If true, return extra info
   * @param w - Weights
   */
  static fit<T extends ABCPolyBase>(
    this: {
      new (coef: number[], domain?: [number, number], window?: [number, number]): T;
      _vander(x: number[], deg: number): number[][];
      defaultDomain: [number, number];
      defaultWindow: [number, number];
    },
    x: number[] | NDArray,
    y: number[] | NDArray,
    deg: number,
    domain: [number, number] | null = null,
    rcond: number | null = null,
    full: boolean = false,
    w: number[] | null = null
  ): T | [T, { resid: number[]; rank: number; sv: number[]; rcond: number }] {
    const xArr = x instanceof NDArray ? (x.toArray() as number[]) : x;
    const yArr = y instanceof NDArray ? (y.toArray() as number[]) : y;

    if (xArr.length !== yArr.length) {
      throw new PolyError('x and y must have same length');
    }

    // Determine domain from data
    const effectiveDomain = domain ?? [Math.min(...xArr), Math.max(...xArr)] as [number, number];

    // Map x to window
    const xMapped = mapdomain(xArr, effectiveDomain, this.defaultWindow) as number[];

    // Build Vandermonde matrix
    const vander = this._vander(xMapped, deg);

    // Apply weights if provided
    let weightedVander = vander;
    let weightedY = yArr;
    if (w !== null) {
      weightedVander = vander.map((row, i) => row.map(v => v * w[i]));
      weightedY = yArr.map((v, i) => v * w[i]);
    }

    // Solve least squares (using simple method, should use linalg.lstsq)
    const coef = _solveLeastSquares(weightedVander, weightedY, rcond);

    const result = new this(coef, effectiveDomain, this.defaultWindow);

    if (full) {
      // Return additional info (simplified)
      return [result, { resid: [], rank: deg + 1, sv: [], rcond: rcond ?? 0 }];
    }

    return result;
  }

  /* ============ Comparison ============ */

  /**
   * Check if same type as other.
   */
  hassametype(other: ABCPolyBase): boolean {
    return this.constructor === other.constructor;
  }

  /**
   * Check if same domain as other.
   */
  hassamedomain(other: ABCPolyBase): boolean {
    return this._domain[0] === other._domain[0] &&
           this._domain[1] === other._domain[1];
  }

  /**
   * Check if same window as other.
   */
  hassamewindow(other: ABCPolyBase): boolean {
    return this._window[0] === other._window[0] &&
           this._window[1] === other._window[1];
  }

  /**
   * Check if same coefficients as other.
   */
  hassamecoef(other: ABCPolyBase): boolean {
    if (this._coef.length !== other._coef.length) return false;
    return this._coef.every((c, i) => c === other._coef[i]);
  }

  /* ============ String Representation ============ */

  /**
   * String representation.
   */
  toString(): string {
    const ctor = this.constructor as typeof ABCPolyBase;
    const terms: string[] = [];

    for (let i = 0; i < this._coef.length; i++) {
      const c = this._coef[i];
      if (c === 0) continue;

      let term: string;
      if (i === 0) {
        term = c.toFixed(4).replace(/\.?0+$/, '');
      } else {
        const coefStr = c === 1 ? '' : (c === -1 ? '-' : `${c.toFixed(4).replace(/\.?0+$/, '')}*`);
        const basisStr = i === 1
          ? this._symbol
          : `${ctor.basisName}_${i}(${this._symbol})`;
        term = `${coefStr}${basisStr}`;
      }
      terms.push(term);
    }

    return terms.length > 0
      ? terms.join(' + ').replace(/\+ -/g, '- ')
      : '0';
  }

  /* ============ Private Helpers ============ */

  /**
   * Create new instance of same type.
   */
  protected _createNew(coef: number[]): this {
    const ctor = this.constructor as new (
      coef: number[],
      domain?: [number, number],
      window?: [number, number],
      symbol?: string
    ) => this;

    return new ctor(coef, this._domain, this._window, this._symbol);
  }

  /**
   * Check compatibility with another polynomial.
   */
  protected _checkCompatible(other: ABCPolyBase): void {
    if (!this.hassametype(other)) {
      throw new PolyError('Polynomial types do not match');
    }
    if (!this.hassamedomain(other)) {
      throw new PolyError('Polynomial domains do not match');
    }
    if (!this.hassamewindow(other)) {
      throw new PolyError('Polynomial windows do not match');
    }
  }
}

/* ============ Helper Functions ============ */

/**
 * Simple least squares solver.
 * In production, should use linalg.lstsq from Phase 13.
 */
function _solveLeastSquares(
  A: number[][],
  b: number[],
  rcond: number | null
): number[] {
  // Normal equations: A^T A x = A^T b
  const m = A.length;
  const n = A[0].length;

  // A^T A
  const AtA: number[][] = [];
  for (let i = 0; i < n; i++) {
    AtA[i] = [];
    for (let j = 0; j < n; j++) {
      let sum = 0;
      for (let k = 0; k < m; k++) {
        sum += A[k][i] * A[k][j];
      }
      AtA[i][j] = sum;
    }
  }

  // A^T b
  const Atb: number[] = [];
  for (let i = 0; i < n; i++) {
    let sum = 0;
    for (let k = 0; k < m; k++) {
      sum += A[k][i] * b[k];
    }
    Atb[i] = sum;
  }

  // Solve using Gaussian elimination
  return _solveLinear(AtA, Atb);
}

/**
 * Simple linear system solver (Gaussian elimination).
 */
function _solveLinear(A: number[][], b: number[]): number[] {
  const n = A.length;
  const aug: number[][] = A.map((row, i) => [...row, b[i]]);

  // Forward elimination
  for (let i = 0; i < n; i++) {
    // Find pivot
    let maxRow = i;
    for (let k = i + 1; k < n; k++) {
      if (Math.abs(aug[k][i]) > Math.abs(aug[maxRow][i])) {
        maxRow = k;
      }
    }
    [aug[i], aug[maxRow]] = [aug[maxRow], aug[i]];

    // Eliminate
    for (let k = i + 1; k < n; k++) {
      const c = aug[k][i] / aug[i][i];
      for (let j = i; j <= n; j++) {
        aug[k][j] -= c * aug[i][j];
      }
    }
  }

  // Back substitution
  const x = new Array(n).fill(0);
  for (let i = n - 1; i >= 0; i--) {
    x[i] = aug[i][n];
    for (let j = i + 1; j < n; j++) {
      x[i] -= aug[i][j] * x[j];
    }
    x[i] /= aug[i][i];
  }

  return x;
}
```

---

### 16c.3 Polynomial Class (Power Series)

**File:** `src/ts/polynomial/polynomial.ts`

```typescript
import { ABCPolyBase, maxpower } from './_polybase.js';
import {
  PolyError,
  trimseq,
  polyadd,
  polysub,
  polymul,
  polydiv,
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
  static basisName = 'P';

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
    if (n === 0) return [1];
    if (n === 1) return [...c];

    // Check degree limit
    if ((c.length - 1) * n > maxpow) {
      throw new PolyError(`Power would exceed maxpower (${maxpow})`);
    }

    // Binary exponentiation
    let result = [1];
    let base = [...c];

    while (n > 0) {
      if (n & 1) {
        result = polymul(result, base);
      }
      base = polymul(base, base);
      n >>= 1;
    }

    return trimseq(result);
  }

  /**
   * Evaluate polynomial using Horner's method.
   */
  protected static _val(x: number | number[], c: number[]): number | number[] {
    if (typeof x === 'number') {
      return polyval(x, c);
    }
    return x.map(v => polyval(v, c));
  }

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

  protected static _int(
    c: number[],
    m: number,
    k: number[],
    lbnd: number,
    scl: number
  ): number[] {
    let result = [...c];

    for (let cnt = 0; cnt < m; cnt++) {
      const newCoef: number[] = [k[cnt]];

      for (let i = 0; i < result.length; i++) {
        newCoef.push(result[i] * scl / (i + 1));
      }

      // Adjust for lower bound
      if (lbnd !== 0) {
        const val = polyval(lbnd, newCoef);
        newCoef[0] -= val - k[cnt];
      }

      result = newCoef;
    }

    return trimseq(result);
  }

  protected static _roots(c: number[]): number[] {
    c = trimseq(c);

    if (c.length <= 1) {
      return [];
    }

    if (c.length === 2) {
      // Linear: c[0] + c[1]*x = 0 => x = -c[0]/c[1]
      return [-c[0] / c[1]];
    }

    // Use companion matrix eigenvalues
    const companion = Polynomial._companion(c);
    return _eigenvalues(companion);
  }

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

  protected static _companion(c: number[]): number[][] {
    c = trimseq(c);
    const n = c.length - 1;

    if (n < 1) {
      throw new PolyError('Companion matrix requires degree >= 1');
    }

    // Normalize so leading coefficient is 1
    const monic = c.map(v => v / c[n]);

    // Build companion matrix
    const mat: number[][] = [];
    for (let i = 0; i < n; i++) {
      const row = new Array(n).fill(0);
      if (i > 0) {
        row[i - 1] = 1;
      }
      row[n - 1] = -monic[i];
      mat.push(row);
    }

    return mat;
  }

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
 * Evaluate 2D polynomial.
 */
export function polyval2d(
  x: number[],
  y: number[],
  c: number[][]
): number[] {
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
 * Generate Vandermonde matrix.
 */
export function polyvander(x: number[], deg: number): number[][] {
  return Polynomial['_vander'](x, deg);
}

/**
 * Polynomial derivative.
 */
export function polyder(c: number[], m: number = 1): number[] {
  return Polynomial['_der'](c, m, 1);
}

/**
 * Polynomial integral.
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
 */
export function polyfit(
  x: number[],
  y: number[],
  deg: number,
  rcond: number | null = null,
  full: boolean = false,
  w: number[] | null = null
): number[] | [number[], any] {
  const result = Polynomial.fit(x, y, deg, null, rcond, full, w);

  if (full) {
    const [poly, info] = result as [Polynomial, any];
    return [poly.coef, info];
  }

  return (result as Polynomial).coef;
}

/**
 * Find polynomial roots.
 */
export function polyroots(c: number[]): number[] {
  return Polynomial['_roots'](c);
}

/**
 * Generate companion matrix.
 */
export function polycompanion(c: number[]): number[][] {
  return Polynomial['_companion'](c);
}

/**
 * Construct polynomial from roots.
 */
export function polyfromroots(roots: number[]): number[] {
  return Polynomial['_fromroots'](roots);
}

/* ============ Helper Functions ============ */

/**
 * Compute eigenvalues of a matrix.
 * Simplified implementation - in production use linalg.eigvals
 */
function _eigenvalues(mat: number[][]): number[] {
  const n = mat.length;

  if (n === 1) {
    return [mat[0][0]];
  }

  if (n === 2) {
    // Quadratic formula for 2x2
    const a = 1;
    const b = -(mat[0][0] + mat[1][1]);
    const c = mat[0][0] * mat[1][1] - mat[0][1] * mat[1][0];

    const disc = b * b - 4 * a * c;
    if (disc >= 0) {
      const sqrtDisc = Math.sqrt(disc);
      return [(-b + sqrtDisc) / (2 * a), (-b - sqrtDisc) / (2 * a)];
    } else {
      // Complex roots - return real parts only for now
      return [-b / (2 * a), -b / (2 * a)];
    }
  }

  // For larger matrices, use QR iteration (simplified)
  // In production, use linalg.eigvals
  throw new PolyError('Root finding for degree > 2 requires linalg module');
}
```

---

### 16c.4 Chebyshev Class

**File:** `src/ts/polynomial/chebyshev.ts`

```typescript
import { ABCPolyBase, maxpower } from './_polybase.js';
import { PolyError, trimseq } from './polyutils.js';

/**
 * Chebyshev polynomial of the first kind.
 *
 * Chebyshev polynomials are defined by:
 * T_n(x) = cos(n * arccos(x))
 *
 * They have excellent numerical properties for interpolation
 * and approximation on [-1, 1].
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
    const len = c1.length + c2.length - 1;
    const result = new Array(len).fill(0);

    for (let i = 0; i < c1.length; i++) {
      for (let j = 0; j < c2.length; j++) {
        const prod = c1[i] * c2[j] * 0.5;
        result[i + j] += prod;
        result[Math.abs(i - j)] += prod;
      }
    }

    // The constant term gets doubled in this process, so halve it
    result[0] *= 0.5;

    return trimseq(result);
  }

  protected static _div(c1: number[], c2: number[]): [number[], number[]] {
    // Chebyshev division is complex - use iterative approach
    c1 = trimseq(c1);
    c2 = trimseq(c2);

    if (c2.length === 1 && c2[0] === 0) {
      throw new PolyError('Division by zero polynomial');
    }

    if (c1.length < c2.length) {
      return [[0], c1];
    }

    // Convert to power series, divide, convert back
    // This is a simplified approach
    const p1 = cheb2poly(c1);
    const p2 = cheb2poly(c2);

    // Power series division
    const [quo, rem] = _polyDivide(p1, p2);

    return [poly2cheb(quo), poly2cheb(rem)];
  }

  protected static _pow(c: number[], n: number, maxpow: number): number[] {
    if (n === 0) return [1];
    if (n === 1) return [...c];

    // Check degree limit
    if ((c.length - 1) * n > maxpow) {
      throw new PolyError(`Power would exceed maxpower (${maxpow})`);
    }

    // Binary exponentiation
    let result = [1];
    let base = [...c];

    while (n > 0) {
      if (n & 1) {
        result = Chebyshev._mul(result, base);
      }
      base = Chebyshev._mul(base, base);
      n >>= 1;
    }

    return trimseq(result);
  }

  /**
   * Evaluate Chebyshev polynomial using Clenshaw's algorithm.
   */
  protected static _val(x: number | number[], c: number[]): number | number[] {
    if (typeof x === 'number') {
      return chebval(x, c);
    }
    return x.map(v => chebval(v, c));
  }

  /**
   * Chebyshev derivative.
   * d/dx T_n(x) = n * U_{n-1}(x) where U is Chebyshev of second kind
   * In T basis: d/dx T_n = n * sum_{k=n-1,n-3,...} 2*T_k / c_k
   */
  protected static _der(c: number[], m: number, scl: number): number[] {
    let result = [...c];

    for (let cnt = 0; cnt < m; cnt++) {
      const n = result.length;
      if (n <= 1) return [0];

      const der = new Array(n - 1).fill(0);

      for (let j = n - 1; j >= 1; j--) {
        der[j - 1] = 2 * j * result[j];
        if (j >= 2) {
          der[j - 2] += der[j - 1] ? 0 : 0; // Accumulate
        }
      }

      // First coefficient adjustment
      der[0] /= 2;

      // Apply scale
      result = der.map(v => v * scl);
    }

    return trimseq(result);
  }

  /**
   * Chebyshev integral.
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
          integ[1] = result[0] * scl;
        } else if (j === 1) {
          integ[0] += result[1] * scl / 2;
          integ[2] = result[1] * scl / 4;
        } else {
          integ[j - 1] -= result[j] * scl / (2 * (j - 1));
          integ[j + 1] += result[j] * scl / (2 * (j + 1));
        }
      }

      // Add integration constant
      const val = chebval(lbnd, integ);
      integ[0] += k[cnt] - val;

      result = integ;
    }

    return trimseq(result);
  }

  /**
   * Find Chebyshev polynomial roots using companion matrix.
   */
  protected static _roots(c: number[]): number[] {
    c = trimseq(c);

    if (c.length <= 1) return [];
    if (c.length === 2) {
      // Linear: c[0]*T_0 + c[1]*T_1 = c[0] + c[1]*x = 0
      return [-c[0] / c[1]];
    }

    // Use companion matrix
    const companion = Chebyshev._companion(c);

    // In production, use linalg.eigvals
    // For now, convert to power series and find roots there
    const poly = cheb2poly(c);
    return _polyRoots(poly);
  }

  /**
   * Chebyshev Vandermonde matrix.
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
   */
  protected static _companion(c: number[]): number[][] {
    c = trimseq(c);
    const n = c.length - 1;

    if (n < 1) {
      throw new PolyError('Companion matrix requires degree >= 1');
    }

    // Build Chebyshev companion matrix
    const mat: number[][] = [];
    for (let i = 0; i < n; i++) {
      const row = new Array(n).fill(0);
      mat.push(row);
    }

    // Subdiagonal
    mat[1][0] = 1;
    for (let i = 2; i < n; i++) {
      mat[i][i - 1] = 0.5;
    }

    // Superdiagonal
    for (let i = 0; i < n - 1; i++) {
      mat[i][i + 1] = 0.5;
    }

    // Last column
    const scale = -0.5 / c[n];
    for (let i = 0; i < n; i++) {
      mat[i][n - 1] += scale * c[i];
    }
    mat[n - 2][n - 1] += 0.5;

    return mat;
  }

  /**
   * Construct Chebyshev polynomial from roots.
   */
  protected static _fromroots(roots: number[]): number[] {
    if (roots.length === 0) return [1];

    // Build polynomial from roots: prod(x - r)
    // Then convert to Chebyshev
    let poly = [1];
    for (const r of roots) {
      poly = _polyMul(poly, [-r, 1]);
    }

    return poly2cheb(poly);
  }
}

/* ============ Module-Level Functions ============ */

/**
 * Evaluate Chebyshev polynomial using Clenshaw's algorithm.
 */
export function chebval(x: number, c: number[]): number {
  if (c.length === 0) return 0;
  if (c.length === 1) return c[0];
  if (c.length === 2) return c[0] + c[1] * x;

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
 * Convert power series to Chebyshev coefficients.
 */
export function poly2cheb(pol: number[]): number[] {
  pol = trimseq(pol);
  const n = pol.length;

  if (n === 0) return [0];
  if (n === 1) return [pol[0]];

  // Use the transformation matrix approach
  // c_k = sum_j T_{kj} * p_j where T is the transformation matrix
  const result = new Array(n).fill(0);

  for (let k = 0; k < n; k++) {
    for (let j = k; j < n; j++) {
      // Coefficient of x^j in T_k
      result[k] += pol[j] * _chebCoefInPowerBasis(k, j);
    }
  }

  return trimseq(result);
}

/**
 * Convert Chebyshev coefficients to power series.
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
 */
export function chebvander(x: number[], deg: number): number[][] {
  return Chebyshev['_vander'](x, deg);
}

/**
 * Chebyshev derivative.
 */
export function chebder(c: number[], m: number = 1): number[] {
  return Chebyshev['_der'](c, m, 1);
}

/**
 * Chebyshev integral.
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
 */
export function chebfit(
  x: number[],
  y: number[],
  deg: number,
  rcond: number | null = null,
  full: boolean = false,
  w: number[] | null = null
): number[] | [number[], any] {
  const result = Chebyshev.fit(x, y, deg, null, rcond, full, w);

  if (full) {
    const [poly, info] = result as [Chebyshev, any];
    return [poly.coef, info];
  }

  return (result as Chebyshev).coef;
}

/**
 * Find Chebyshev polynomial roots.
 */
export function chebroots(c: number[]): number[] {
  return Chebyshev['_roots'](c);
}

/**
 * Generate Chebyshev companion matrix.
 */
export function chebcompanion(c: number[]): number[][] {
  return Chebyshev['_companion'](c);
}

/**
 * Construct Chebyshev polynomial from roots.
 */
export function chebfromroots(roots: number[]): number[] {
  return Chebyshev['_fromroots'](roots);
}

/**
 * Interpolate a function using Chebyshev polynomials.
 */
export function chebinterpolate(
  func: (x: number) => number,
  deg: number,
  domain: [number, number] = [-1, 1]
): number[] {
  // Use Chebyshev nodes
  const nodes: number[] = [];
  for (let k = 0; k < deg + 1; k++) {
    nodes.push(Math.cos((2 * k + 1) * Math.PI / (2 * (deg + 1))));
  }

  // Map to domain
  const [d0, d1] = domain;
  const mappedNodes = nodes.map(x => (d0 + d1) / 2 + (d1 - d0) / 2 * x);

  // Evaluate function
  const y = mappedNodes.map(func);

  // Fit
  return chebfit(nodes, y, deg) as number[];
}

/* ============ Helper Functions ============ */

/**
 * Get coefficient of x^j in T_k(x).
 */
function _chebCoefInPowerBasis(k: number, j: number): number {
  if ((k - j) % 2 !== 0) return 0;
  if (j > k) return 0;

  // T_k has coefficients involving binomial coefficients
  // This is a simplified computation
  const coefs = _chebToPowerCoefs(k);
  return j < coefs.length ? coefs[j] : 0;
}

/**
 * Get T_k(x) expressed in power basis coefficients.
 */
function _chebToPowerCoefs(k: number): number[] {
  if (k === 0) return [1];
  if (k === 1) return [0, 1];

  // Use recurrence: T_{n+1} = 2x*T_n - T_{n-1}
  let prev2 = [1];      // T_0
  let prev1 = [0, 1];   // T_1

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

// Power series helpers (imported from polynomial or defined here)
function _polyMul(c1: number[], c2: number[]): number[] {
  const len = c1.length + c2.length - 1;
  const result = new Array(len).fill(0);
  for (let i = 0; i < c1.length; i++) {
    for (let j = 0; j < c2.length; j++) {
      result[i + j] += c1[i] * c2[j];
    }
  }
  return trimseq(result);
}

function _polyDivide(c1: number[], c2: number[]): [number[], number[]] {
  c1 = trimseq(c1);
  c2 = trimseq(c2);

  if (c1.length < c2.length) return [[0], c1];

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

function _polyRoots(c: number[]): number[] {
  // Simplified root finding - in production use linalg
  c = trimseq(c);
  if (c.length <= 1) return [];
  if (c.length === 2) return [-c[0] / c[1]];
  if (c.length === 3) {
    const [a, b, c_] = [c[2], c[1], c[0]];
    const disc = b * b - 4 * a * c_;
    if (disc >= 0) {
      const sqrtDisc = Math.sqrt(disc);
      return [(-b + sqrtDisc) / (2 * a), (-b - sqrtDisc) / (2 * a)];
    }
    return [];
  }
  throw new PolyError('Root finding for degree > 2 requires linalg module');
}
```

---

### Module Index

**File:** `src/ts/polynomial/index.ts`

```typescript
/**
 * NumJS Polynomial Module
 *
 * Provides polynomial classes for multiple basis systems,
 * compatible with NumPy's numpy.polynomial module.
 */

// Utilities
export {
  PolyError,
  PolyDomainWarning,
  as_series,
  trimseq,
  trimcoef,
  getdomain,
  mapdomain,
  mapparms,
} from './polyutils.js';

// Base class
export { ABCPolyBase, maxpower } from './_polybase.js';

// Polynomial (power series)
export {
  Polynomial,
  polyval,
  polyval2d,
  polyvander,
  polyder,
  polyint,
  polyfit,
  polyroots,
  polycompanion,
  polyfromroots,
} from './polynomial.js';

// Chebyshev
export {
  Chebyshev,
  chebval,
  chebvander,
  chebder,
  chebint,
  chebfit,
  chebroots,
  chebcompanion,
  chebfromroots,
  chebinterpolate,
  poly2cheb,
  cheb2poly,
} from './chebyshev.js';

// Legendre (similar structure)
export { Legendre } from './legendre.js';

// Hermite
export { Hermite } from './hermite.js';

// HermiteE
export { HermiteE } from './hermite_e.js';

// Laguerre
export { Laguerre } from './laguerre.js';
```

---

## Implementation Order

```
Weeks 4-6: numpy.polynomial Module

Week 4: Foundation
├── Day 1: polyutils.ts
│   ├── trimseq, trimcoef
│   ├── as_series
│   ├── getdomain, mapdomain, mapparms
│   └── PolyError, PolyDomainWarning
│
├── Day 2-3: _polybase.ts (ABCPolyBase)
│   ├── Constructor, properties
│   ├── Arithmetic methods
│   ├── Calculus methods (deriv, integ)
│   └── Analysis methods (roots, trim)
│
├── Day 4-5: polynomial.ts (Polynomial class)
│   ├── Power series operations
│   ├── Horner's method evaluation
│   ├── Module-level functions
│   └── Tests

Week 5: Chebyshev & Legendre
├── Day 1-2: chebyshev.ts
│   ├── Chebyshev multiplication
│   ├── Clenshaw evaluation
│   ├── Basis conversions
│   └── Interpolation
│
├── Day 3-4: legendre.ts
│   ├── Legendre-specific arithmetic
│   ├── Evaluation using recurrence
│   └── Fitting
│
└── Day 5: Tests for Chebyshev & Legendre

Week 6: Hermite & Laguerre
├── Day 1-2: hermite.ts (Physicist's Hermite)
│   └── H_n operations
│
├── Day 2-3: hermite_e.ts (Probabilist's Hermite)
│   └── He_n operations
│
├── Day 4: laguerre.ts
│   └── L_n operations
│
└── Day 5: Integration tests & documentation
```

---

## Verification Plan

```bash
# Build
npm run build

# Run tests
npm test -- --grep "polynomial"

# Test cases:

# Power Series
✓ Polynomial([1,2,3]).call(2) equals 17
✓ Polynomial.deriv([1,2,3]) equals [2,6]
✓ Polynomial.integ([2,6]) equals [0,2,3]
✓ Polynomial.fit(x, y, 2) returns correct coefficients
✓ Polynomial.roots([6,-5,1]) equals [2,3]

# Chebyshev
✓ chebval(0.5, [1,2,3]) equals correct value
✓ poly2cheb → cheb2poly roundtrip is identity
✓ chebfit matches NumPy output
✓ chebinterpolate approximates function

# Cross-basis
✓ Converting Polynomial to Chebyshev preserves evaluation
✓ Degree is preserved across conversions
```

---

## Estimated Lines of Code

| File | Lines |
|------|-------|
| polyutils.ts | ~200 |
| _polybase.ts | ~450 |
| polynomial.ts | ~350 |
| chebyshev.ts | ~400 |
| legendre.ts | ~350 |
| hermite.ts | ~350 |
| hermite_e.ts | ~350 |
| laguerre.ts | ~350 |
| index.ts | ~80 |
| **Total** | **~2,880** |
