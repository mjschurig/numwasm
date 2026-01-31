/**
 * NumJS Polynomial Base Class
 *
 * Abstract base class for polynomial series of all types.
 * Provides the common interface and operations for Polynomial, Chebyshev,
 * Legendre, Hermite, HermiteE, and Laguerre polynomial classes.
 */

import { NDArray } from "../_core/NDArray.js";
import { DType } from "../types.js";
import { lstsq, eigvals } from "../linalg.js";
import {
  PolyError,
  trimcoef,
  mapparms,
  mapdomain,
  getdomain,
} from "./polyutils.js";

/**
 * Maximum allowed polynomial degree.
 * Prevents excessive memory usage and computation time.
 */
export const maxpower = 100;

/**
 * Abstract base class for polynomial series.
 *
 * This class provides the common interface and operations for all
 * polynomial types (Polynomial, Chebyshev, Legendre, etc.).
 *
 * Subclasses must implement the static methods for basis-specific arithmetic.
 *
 * Coefficients are stored in order of increasing degree, so:
 * - `[1, 2, 3]` represents `1*P_0 + 2*P_1 + 3*P_2`
 *
 * Where `P_n` is the nth basis polynomial (x^n for power series, T_n for Chebyshev, etc.)
 */
export abstract class ABCPolyBase {
  /** Coefficient array (lowest degree first) */
  protected _coef: number[];

  /** Domain for the polynomial */
  protected _domain: [number, number];

  /** Window for evaluation (basis-specific natural domain) */
  protected _window: [number, number];

  /** Symbol used in string representation */
  protected _symbol: string;

  /* ============ Static Properties (Override in Subclasses) ============ */

  /** Default domain for this polynomial type */
  static defaultDomain: [number, number] = [-1, 1];

  /** Default window for this polynomial type */
  static defaultWindow: [number, number] = [-1, 1];

  /** Name of the basis (for display) */
  static basisName: string = "P";

  /* ============ Abstract Static Methods ============ */
  /* These must be implemented by each polynomial type */

  /**
   * Add two coefficient arrays in this basis.
   */
  protected static _add(_c1: number[], _c2: number[]): number[] {
    throw new Error("Must be implemented by subclass");
  }

  /**
   * Subtract two coefficient arrays in this basis.
   */
  protected static _sub(_c1: number[], _c2: number[]): number[] {
    throw new Error("Must be implemented by subclass");
  }

  /**
   * Multiply two coefficient arrays in this basis.
   */
  protected static _mul(_c1: number[], _c2: number[]): number[] {
    throw new Error("Must be implemented by subclass");
  }

  /**
   * Divide two coefficient arrays in this basis.
   * @returns [quotient, remainder]
   */
  protected static _div(_c1: number[], _c2: number[]): [number[], number[]] {
    throw new Error("Must be implemented by subclass");
  }

  /**
   * Raise coefficient array to integer power.
   */
  protected static _pow(_c: number[], _n: number, _maxpow: number): number[] {
    throw new Error("Must be implemented by subclass");
  }

  /**
   * Evaluate polynomial at x.
   */
  protected static _val(
    _x: number | number[],
    _c: number[],
  ): number | number[] {
    throw new Error("Must be implemented by subclass");
  }

  /**
   * Compute derivative of coefficient array.
   * @param c - Coefficient array
   * @param m - Number of derivatives to take
   * @param scl - Scale factor for the derivative
   */
  protected static _der(_c: number[], _m: number, _scl: number): number[] {
    throw new Error("Must be implemented by subclass");
  }

  /**
   * Compute integral of coefficient array.
   * @param c - Coefficient array
   * @param m - Number of integrals to compute
   * @param k - Integration constants (one per integration)
   * @param lbnd - Lower bound for definite integral
   * @param scl - Scale factor for the integral
   */
  protected static _int(
    _c: number[],
    _m: number,
    _k: number[],
    _lbnd: number,
    _scl: number,
  ): number[] {
    throw new Error("Must be implemented by subclass");
  }

  /**
   * Find roots of polynomial from its coefficients.
   * Default implementation uses companion matrix eigenvalues.
   */
  protected static async _roots(_c: number[]): Promise<number[]> {
    throw new Error("Must be implemented by subclass");
  }

  /**
   * Generate Vandermonde matrix for given x values and degree.
   */
  protected static _vander(_x: number[], _deg: number): number[][] {
    throw new Error("Must be implemented by subclass");
  }

  /**
   * Generate companion matrix for root finding.
   */
  protected static _companion(_c: number[]): number[][] {
    throw new Error("Must be implemented by subclass");
  }

  /**
   * Construct polynomial coefficients from roots.
   */
  protected static _fromroots(_roots: number[]): number[] {
    throw new Error("Must be implemented by subclass");
  }

  /* ============ Constructor ============ */

  /**
   * Create a polynomial.
   *
   * @param coef - Coefficients (lowest degree first)
   * @param domain - Domain for the polynomial (default: class default)
   * @param window - Window for evaluation (default: class default)
   * @param symbol - Variable symbol for display (default: 'x')
   */
  constructor(
    coef: number[] | NDArray,
    domain: [number, number] | null = null,
    window: [number, number] | null = null,
    symbol: string = "x",
  ) {
    // Convert and validate coefficients
    if (coef instanceof NDArray) {
      if (coef.ndim !== 1) {
        throw new PolyError("Coefficient array must be 1-D");
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

    // Set domain and window from class defaults
    const ctor = this.constructor as typeof ABCPolyBase;
    this._domain = domain ?? ([...ctor.defaultDomain] as [number, number]);
    this._window = window ?? ([...ctor.defaultWindow] as [number, number]);
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

    if (typeof other === "number") {
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

    if (typeof other === "number") {
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

    if (typeof other === "number") {
      const newCoef = this._coef.map((c) => c * other);
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

    this._checkCompatible(other);
    const [quo, rem] = ctor._div(this._coef, other._coef);
    return [this._createNew(quo), this._createNew(rem)];
  }

  /**
   * Floor division (quotient only).
   */
  floordiv(other: ABCPolyBase): this {
    return this.divmod(other)[0];
  }

  /**
   * Modulo (remainder only).
   */
  mod(other: ABCPolyBase): this {
    return this.divmod(other)[1];
  }

  /**
   * Raise to integer power.
   */
  pow(n: number): this {
    if (!Number.isInteger(n) || n < 0) {
      throw new PolyError("Power must be a non-negative integer");
    }

    const ctor = this.constructor as typeof ABCPolyBase;
    const newCoef = ctor._pow(this._coef, n, maxpower);
    return this._createNew(newCoef);
  }

  /**
   * Negate polynomial.
   */
  neg(): this {
    return this._createNew(this._coef.map((c) => -c));
  }

  /* ============ Evaluation ============ */

  /**
   * Evaluate polynomial at x.
   *
   * Maps x from domain to window before evaluation.
   *
   * @param x - Point(s) at which to evaluate
   * @returns Evaluated value(s)
   */
  call(x: number | number[] | NDArray): number | number[] {
    const ctor = this.constructor as typeof ABCPolyBase;

    // Map x from domain to window
    const xInput = x instanceof NDArray ? (x.toArray() as number[]) : x;
    const xMapped = mapdomain(xInput, this._domain, this._window);

    return ctor._val(xMapped, this._coef);
  }

  /**
   * Evaluate at n equally-spaced points in domain.
   *
   * @param n - Number of points (default: 100)
   * @returns Object with x and y arrays
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
   *
   * @param m - Number of derivatives (default: 1)
   * @returns Derivative polynomial
   */
  deriv(m: number = 1): this {
    if (!Number.isInteger(m) || m < 0) {
      throw new PolyError("Derivative order must be a non-negative integer");
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
   * @param m - Integration order (default: 1)
   * @param k - Integration constants (one per integration, default: 0s)
   * @param lbnd - Lower bound for definite integral (default: domain[0])
   * @returns Integral polynomial
   */
  integ(m: number = 1, k: number[] = [], lbnd: number | null = null): this {
    if (!Number.isInteger(m) || m < 0) {
      throw new PolyError("Integration order must be a non-negative integer");
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
   *
   * Uses companion matrix eigenvalues via the linalg module.
   *
   * @returns Array of roots
   */
  async roots(): Promise<number[]> {
    const ctor = this.constructor as typeof ABCPolyBase;
    const windowRoots = await ctor._roots(this._coef);

    // Map roots from window to domain
    return mapdomain(windowRoots, this._window, this._domain) as number[];
  }

  /**
   * Trim small trailing coefficients.
   *
   * @param tol - Tolerance for considering coefficients as zero
   * @returns Trimmed polynomial
   */
  trim(tol: number = 0): this {
    return this._createNew(trimcoef(this._coef, tol));
  }

  /**
   * Truncate to number of coefficients.
   *
   * @param size - Number of coefficients to keep
   * @returns Truncated polynomial
   */
  truncate(size: number): this {
    if (size < 1) {
      throw new PolyError("Size must be at least 1");
    }
    return this._createNew(this._coef.slice(0, size));
  }

  /**
   * Truncate to maximum degree.
   *
   * @param deg - Maximum degree
   * @returns Truncated polynomial
   */
  cutdeg(deg: number): this {
    return this.truncate(deg + 1);
  }

  /* ============ Conversion ============ */

  /**
   * Get mapping parameters from domain to window.
   *
   * @returns [offset, scale]
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
   * Returns a polynomial that is 1 in the specified basis element and 0 elsewhere.
   *
   * @param deg - Degree of the basis polynomial
   * @param domain - Domain for the polynomial
   * @param window - Window for the polynomial
   *
   * @example
   * Polynomial.basis(2)  // x^2 for power series
   * Chebyshev.basis(2)   // T_2 for Chebyshev
   */
  static basis<T extends ABCPolyBase>(
    this: new (
      coef: number[],
      domain?: [number, number] | null,
      window?: [number, number] | null,
    ) => T,
    deg: number,
    domain: [number, number] | null = null,
    window: [number, number] | null = null,
  ): T {
    if (!Number.isInteger(deg) || deg < 0) {
      throw new PolyError("Degree must be a non-negative integer");
    }

    const coef = new Array(deg + 1).fill(0);
    coef[deg] = 1;

    return new this(coef, domain, window);
  }

  /**
   * Create identity polynomial p(x) = x.
   *
   * @param domain - Domain for the polynomial
   * @param window - Window for the polynomial
   */
  static identity<T extends ABCPolyBase>(
    this: new (
      coef: number[],
      domain?: [number, number] | null,
      window?: [number, number] | null,
    ) => T,
    domain: [number, number] | null = null,
    window: [number, number] | null = null,
  ): T {
    return new this([0, 1], domain, window);
  }

  /**
   * Create polynomial from roots.
   *
   * @param roots - Array of roots
   * @param domain - Domain for the polynomial
   * @param window - Window for the polynomial
   */
  static fromroots<T extends ABCPolyBase>(
    this: {
      new (
        coef: number[],
        domain?: [number, number] | null,
        window?: [number, number] | null,
      ): T;
      _fromroots(roots: number[]): number[];
    },
    roots: number[],
    domain: [number, number] | null = null,
    window: [number, number] | null = null,
  ): T {
    const coef = this._fromroots(roots);
    return new this(coef, domain, window);
  }

  /**
   * Least squares fit to data.
   *
   * @param x - x values (independent variable)
   * @param y - y values (dependent variable)
   * @param deg - Degree of fitting polynomial
   * @param domain - Domain for the polynomial (default: from data)
   * @param rcond - Cutoff for small singular values (default: machine precision)
   * @param full - If true, return extra fitting info
   * @param w - Weights for weighted least squares
   * @returns Fitted polynomial (or [polynomial, info] if full=true)
   */
  static async fit<T extends ABCPolyBase>(
    this: (new (
      coef: number[],
      domain?: [number, number] | null,
      window?: [number, number] | null,
    ) => T) & {
      defaultDomain: [number, number];
      defaultWindow: [number, number];
    },
    x: number[] | NDArray,
    y: number[] | NDArray,
    deg: number,
    domain: [number, number] | null = null,
    rcond: number | null = null,
    full: boolean = false,
    w: number[] | null = null,
  ): Promise<
    T | [T, { residuals: number[]; rank: number; sv: number[]; rcond: number }]
  > {
    const xArr = x instanceof NDArray ? (x.toArray() as number[]) : x;
    const yArr = y instanceof NDArray ? (y.toArray() as number[]) : y;

    if (xArr.length !== yArr.length) {
      throw new PolyError("x and y must have same length");
    }

    if (xArr.length === 0) {
      throw new PolyError("x and y must have at least one element");
    }

    // Determine domain from data if not provided
    const effectiveDomain = domain ?? getdomain(xArr);

    // Map x to window
    const xMapped = mapdomain(
      xArr,
      effectiveDomain,
      this.defaultWindow,
    ) as number[];

    // Build Vandermonde matrix (access via the concrete class)
    const vander = (
      this as unknown as { _vander(x: number[], deg: number): number[][] }
    )._vander(xMapped, deg);

    // Apply weights if provided
    let weightedVander = vander;
    let weightedY = yArr;
    if (w !== null) {
      if (w.length !== xArr.length) {
        throw new PolyError("Weights must have same length as x");
      }
      weightedVander = vander.map((row, i) =>
        row.map((v) => v * Math.sqrt(w[i])),
      );
      weightedY = yArr.map((v, i) => v * Math.sqrt(w[i]));
    }

    // Convert to NDArrays for lstsq
    const m = weightedVander.length;
    const n = weightedVander[0].length;
    const flatVander = weightedVander.flat();
    const vanderArr = await NDArray.fromTypedArray(
      new Float64Array(flatVander),
      [m, n],
      DType.Float64,
    );
    const yNDArr = await NDArray.fromTypedArray(
      new Float64Array(weightedY),
      [m],
      DType.Float64,
    );

    // Solve least squares using linalg module
    const result = await lstsq(vanderArr, yNDArr, rcond);
    const coef = result.x.toArray() as number[];

    // Clean up NDArrays
    vanderArr.dispose();
    yNDArr.dispose();
    result.x.dispose();
    result.residuals.dispose();
    result.s.dispose();

    const poly = new this(coef, effectiveDomain, this.defaultWindow);

    if (full) {
      return [
        poly,
        {
          residuals: [],
          rank: result.rank,
          sv: [],
          rcond: rcond ?? 0,
        },
      ];
    }

    return poly;
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
    return (
      this._domain[0] === other._domain[0] &&
      this._domain[1] === other._domain[1]
    );
  }

  /**
   * Check if same window as other.
   */
  hassamewindow(other: ABCPolyBase): boolean {
    return (
      this._window[0] === other._window[0] &&
      this._window[1] === other._window[1]
    );
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
        term = formatCoef(c);
      } else {
        const absC = Math.abs(c);
        const sign = c < 0 ? "-" : "+";
        const coefStr = absC === 1 ? "" : formatCoef(absC) + "*";
        const basisStr =
          i === 1 ? this._symbol : `${ctor.basisName}_${i}(${this._symbol})`;
        term = `${sign} ${coefStr}${basisStr}`;
      }
      terms.push(term);
    }

    if (terms.length === 0) {
      return "0";
    }

    let result = terms[0];
    for (let i = 1; i < terms.length; i++) {
      result += " " + terms[i];
    }
    return result;
  }

  /* ============ Private Helpers ============ */

  /**
   * Create new instance of same type.
   */
  protected _createNew(coef: number[]): this {
    const ctor = this.constructor as new (
      coef: number[],
      domain?: [number, number] | null,
      window?: [number, number] | null,
      symbol?: string,
    ) => this;

    return new ctor(coef, this._domain, this._window, this._symbol);
  }

  /**
   * Check compatibility with another polynomial.
   */
  protected _checkCompatible(other: ABCPolyBase): void {
    if (!this.hassametype(other)) {
      throw new PolyError("Polynomial types do not match");
    }
    if (!this.hassamedomain(other)) {
      throw new PolyError("Polynomial domains do not match");
    }
    if (!this.hassamewindow(other)) {
      throw new PolyError("Polynomial windows do not match");
    }
  }
}

/* ============ Helper Functions ============ */

/**
 * Format a coefficient for display.
 */
function formatCoef(c: number): string {
  if (Number.isInteger(c)) {
    return String(c);
  }
  // Remove trailing zeros
  return c.toFixed(6).replace(/\.?0+$/, "");
}

/**
 * Find eigenvalues of a companion matrix.
 * Used for root finding in polynomial classes.
 */
export async function companionEigenvalues(
  companion: number[][],
): Promise<number[]> {
  const n = companion.length;
  if (n === 0) return [];

  // Convert to NDArray
  const flat = companion.flat();
  const mat = await NDArray.fromTypedArray(
    new Float64Array(flat),
    [n, n],
    DType.Float64,
  );

  // Get eigenvalues
  const eigenNDArray = await eigvals(mat);
  const eigenvalues = eigenNDArray.toArray() as number[];

  // Clean up
  mat.dispose();
  eigenNDArray.dispose();

  return eigenvalues;
}
