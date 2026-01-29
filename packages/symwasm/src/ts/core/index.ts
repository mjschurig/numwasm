/**
 * Core symbolic expression types.
 * @module core
 */

import { loadWasmModule, getWasmModule } from '../wasm-loader.js';
import {
  SymEngineObject,
  SymEngineSet,
  SymEngineVec,
  SymEngineMap,
  checkException,
  createBasic,
  withTempString,
} from '../wasm-memory.js';
import { SymEngineTypeID } from '../wasm-types.js';

// Export WASM initialization for users
export { loadWasmModule };

// Re-export type ID enum for users
export { SymEngineTypeID };

// Internal WASM access helper
export function getWasm() {
  return getWasmModule();
}

/**
 * Base class for all symbolic expressions.
 * All symbolic objects (symbols, numbers, operations) extend this class.
 *
 * Expressions are backed by SymEngine WASM objects. The _obj field holds
 * the pointer wrapper; it may be null for sentinel constants (pi, E, etc.)
 * that aren't yet fully backed by WASM.
 */
export abstract class Expr {
  /** @internal WASM pointer wrapper - null for sentinel constants */
  protected _obj: SymEngineObject | null = null;

  /**
   * @internal Get the underlying WASM pointer
   * @throws Error if not backed by WASM object
   */
  getWasmPtr(): number {
    if (!this._obj) {
      throw new Error('Expr not backed by WASM object');
    }
    return this._obj.getPtr();
  }

  /**
   * @internal Check if this expression is backed by a WASM object
   */
  protected hasWasmBacking(): boolean {
    return this._obj !== null && this._obj.isValid();
  }

  /**
   * String representation of this expression.
   * Uses WASM _basic_str if available, otherwise falls back to subclass implementation.
   */
  toString(): string {
    if (this._obj && this._obj.isValid()) {
      return this._obj.toString();
    }
    return this._fallbackString();
  }

  /**
   * @internal Fallback string representation for subclasses
   * Used when WASM object is not available
   */
  protected abstract _fallbackString(): string;

  /**
   * Substitute expressions.
   * @param old Expression to replace, OR a Map/object of substitutions
   * @param new_ Replacement expression (only when old is an Expr)
   * @returns A new expression with substitutions applied
   *
   * @example
   * // Single substitution
   * expr.subs(x, new Integer(1))
   *
   * // Multiple substitutions with Map
   * expr.subs(new Map([[x, one], [y, two]]))
   *
   * // Multiple substitutions with object (symbols by name)
   * expr.subs({ x: 1, y: 2 })
   */
  subs(old: Expr | Map<Expr, Expr> | Record<string, Expr | number>, new_?: Expr): Expr {
    if (!this._obj || !this._obj.isValid()) {
      throw new Error('Cannot substitute in expression not backed by WASM object');
    }

    const wasm = getWasmModule();

    // Case 1: Single substitution subs(old, new_)
    if (old instanceof Expr && new_ !== undefined) {
      const result = createBasic();
      try {
        const code = wasm._basic_subs2(
          result.getPtr(),
          this._obj.getPtr(),
          old.getWasmPtr(),
          new_.getWasmPtr()
        );
        checkException(code);
        return exprFromWasm(result);
      } catch (e) {
        result.free();
        throw e;
      }
    }

    // Case 2: Map-based substitution subs(Map<Expr, Expr>)
    if (old instanceof Map) {
      return this._subsWithMap(old);
    }

    // Case 3: Object-based substitution subs({ x: 1, y: 2 })
    if (typeof old === 'object' && !(old instanceof Expr)) {
      return this._subsWithObject(old as Record<string, Expr | number>);
    }

    throw new Error('Invalid arguments to subs()');
  }

  /** @internal Map-based substitution */
  private _subsWithMap(substitutions: Map<Expr, Expr>): Expr {
    const wasm = getWasmModule();
    const map = new SymEngineMap();
    const result = createBasic();

    try {
      for (const [key, value] of substitutions) {
        map.insert(key.getWasmPtr(), value.getWasmPtr());
      }

      const code = wasm._basic_subs(result.getPtr(), this._obj!.getPtr(), map.getPtr());
      checkException(code);
      return exprFromWasm(result);
    } catch (e) {
      result.free();
      throw e;
    } finally {
      map.free();
    }
  }

  /** @internal Object-based substitution (creates symbols by name) */
  private _subsWithObject(substitutions: Record<string, Expr | number>): Expr {
    const exprMap = new Map<Expr, Expr>();
    for (const [name, value] of Object.entries(substitutions)) {
      const sym = new Symbol(name);
      const expr = typeof value === 'number' ? new Integer(value) : value;
      exprMap.set(sym, expr);
    }
    return this._subsWithMap(exprMap);
  }

  /**
   * Evaluate the expression numerically.
   * @param precision Number of bits of precision (default: 53 for double)
   * @returns A RealDouble or ComplexDouble expression
   */
  evalf(precision: number = 53): Expr {
    if (!this._obj || !this._obj.isValid()) {
      throw new Error('Cannot evaluate expression not backed by WASM object');
    }

    const wasm = getWasmModule();
    const result = createBasic();

    try {
      // Use EvalfDomain.Symbolic (2) to let SymEngine decide real/complex
      const code = wasm._basic_evalf(
        result.getPtr(),
        this._obj.getPtr(),
        precision,
        2 // EvalfDomain.Symbolic
      );
      checkException(code);
      return exprFromWasm(result);
    } catch (e) {
      result.free();
      throw e;
    }
  }

  /**
   * Evaluate and extract as JavaScript number.
   * @param precision Number of bits of precision (default: 53)
   * @returns The numerical value as a JavaScript number
   * @throws Error if result is complex (use evalfComplex instead)
   */
  evalfNumber(precision: number = 53): number {
    const result = this.evalf(precision);
    const typeId = result.get_type();

    if (typeId === SymEngineTypeID.SYMENGINE_REAL_DOUBLE) {
      const wasm = getWasmModule();
      return wasm._real_double_get_d(result.getWasmPtr());
    } else if (typeId === SymEngineTypeID.SYMENGINE_INTEGER) {
      return parseInt(result.toString(), 10);
    } else if (typeId === SymEngineTypeID.SYMENGINE_RATIONAL) {
      const [num, den] = result.toString().split('/').map(Number);
      return num / (den || 1);
    } else if (typeId === SymEngineTypeID.SYMENGINE_COMPLEX_DOUBLE) {
      throw new Error('Result is complex. Use evalfComplex() instead.');
    }

    return parseFloat(result.toString());
  }

  /**
   * Evaluate and extract as complex number.
   * @param precision Number of bits of precision (default: 53)
   * @returns Object with real and imag properties
   */
  evalfComplex(precision: number = 53): { real: number; imag: number } {
    const result = this.evalf(precision);
    const typeId = result.get_type();
    const wasm = getWasmModule();

    if (typeId === SymEngineTypeID.SYMENGINE_COMPLEX_DOUBLE) {
      // Parse the string representation of complex double
      // Format is typically "real + imag*I" or "real - imag*I" or "imag*I"
      return parseComplexString(result.toString());
    } else if (typeId === SymEngineTypeID.SYMENGINE_REAL_DOUBLE) {
      return { real: wasm._real_double_get_d(result.getWasmPtr()), imag: 0 };
    }

    return { real: parseFloat(result.toString()), imag: 0 };
  }

  /**
   * Check structural equality with another expression.
   * Two expressions are equal if they have the same mathematical structure.
   */
  equals(other: Expr): boolean {
    // If both have WASM backing, use WASM comparison
    if (this._obj && other._obj) {
      return this._obj.equals(other._obj);
    }
    // For sentinel constants without WASM backing, use reference equality
    return this === other;
  }

  /**
   * Get the hash code for this expression.
   * Equal expressions will have equal hash codes.
   * @throws Error if not backed by WASM object
   */
  hash(): number {
    if (!this._obj) {
      throw new Error('Cannot hash expression not backed by WASM object');
    }
    return this._obj.hash();
  }

  /**
   * Get the SymEngine type ID for this expression.
   * @throws Error if not backed by WASM object
   */
  get_type(): SymEngineTypeID {
    if (!this._obj) {
      throw new Error('Cannot get type of expression not backed by WASM object');
    }
    return this._obj.getType() as SymEngineTypeID;
  }

  /**
   * Return free symbols in this expression.
   * A free symbol is a Symbol that appears in the expression.
   */
  free_symbols(): Symbol[] {
    // No WASM backing means no symbols (e.g., sentinel constants)
    if (!this._obj || !this._obj.isValid()) {
      return [];
    }

    const wasm = getWasmModule();
    const set = new SymEngineSet();
    try {
      const code = wasm._basic_free_symbols(this._obj.getPtr(), set.getPtr());
      checkException(code);

      const symbols: Symbol[] = [];
      const count = set.size();
      for (let i = 0; i < count; i++) {
        const obj = set.get(i);
        symbols.push(Symbol._fromWasm(obj));
      }
      return symbols;
    } finally {
      set.free();
    }
  }

  /**
   * Free the underlying WASM memory.
   * After calling this, the expression becomes invalid.
   */
  free(): void {
    if (this._obj) {
      this._obj.free();
      this._obj = null;
    }
  }

  /**
   * Get the arguments (sub-expressions) of this expression.
   * For atoms (symbols, numbers), returns empty array.
   * For compound expressions (Add, Mul, Pow), returns the operands.
   */
  get_args(): Expr[] {
    if (!this._obj || !this._obj.isValid()) {
      return [];
    }

    const wasm = getWasmModule();
    const vec = new SymEngineVec();
    try {
      const code = wasm._basic_get_args(this._obj.getPtr(), vec.getPtr());
      checkException(code);
      const result: Expr[] = [];
      const count = vec.size();
      for (let i = 0; i < count; i++) {
        result.push(exprFromWasm(vec.get(i)));
      }
      return result;
    } finally {
      vec.free();
    }
  }
}

/**
 * A named symbolic variable.
 * Mirrors sympy.Symbol.
 */
export class Symbol extends Expr {
  readonly name: string;

  /**
   * Create a new symbolic variable.
   * @param name The name of the symbol
   * @param _assumptions Optional assumptions (not yet implemented)
   */
  constructor(name: string, _assumptions?: Record<string, boolean>) {
    super();
    this.name = name;

    // Create WASM-backed symbol
    const wasm = getWasmModule();
    const obj = createBasic();
    const code = withTempString(name, (namePtr) =>
      wasm._symbol_set(obj.getPtr(), namePtr)
    );
    checkException(code);
    this._obj = obj;
  }

  protected _fallbackString(): string {
    return this.name;
  }

  /**
   * @internal Create a Symbol from an existing WASM object.
   * Used by free_symbols() and other methods that return symbols.
   */
  static _fromWasm(obj: SymEngineObject): Symbol {
    // Use Object.create to bypass the constructor
    const sym = Object.create(Symbol.prototype) as Symbol;
    (sym as any)._obj = obj;
    (sym as any).name = obj.toString();
    return sym;
  }
}

/**
 * Create multiple symbols at once.
 * Mirrors sympy.symbols.
 *
 * @param names Space-separated or comma-separated symbol names
 * @param assumptions Optional assumptions to apply to all symbols
 * @returns Array of Symbol instances
 *
 * @example
 * const [x, y, z] = symbols('x y z');
 * const [a, b] = symbols('a, b');
 */
export function symbols(names: string, assumptions?: Record<string, boolean>): Symbol[] {
  // Split by whitespace or commas, filter empty strings
  const nameList = names
    .split(/[\s,]+/)
    .map((n) => n.trim())
    .filter((n) => n.length > 0);

  return nameList.map((name) => new Symbol(name, assumptions));
}

/** Exact integer value. Mirrors sympy.Integer. */
export class Integer extends Expr {
  readonly value: number;

  /**
   * Create an exact integer value.
   * @param value The integer value (must fit in a signed 32-bit integer)
   */
  constructor(value: number) {
    super();
    this.value = Math.trunc(value); // Ensure it's an integer

    // Create WASM-backed integer
    const wasm = getWasmModule();
    const obj = createBasic();
    const code = wasm._integer_set_si(obj.getPtr(), this.value);
    checkException(code);
    this._obj = obj;
  }

  protected _fallbackString(): string {
    return String(this.value);
  }

  /**
   * @internal Create an Integer from an existing WASM object.
   */
  static _fromWasm(obj: SymEngineObject): Integer {
    const int = Object.create(Integer.prototype) as Integer;
    (int as any)._obj = obj;
    (int as any).value = parseInt(obj.toString(), 10);
    return int;
  }
}

/** Exact rational number p/q. Mirrors sympy.Rational. */
export class Rational extends Expr {
  readonly p: number;
  readonly q: number;

  /**
   * Create an exact rational number p/q.
   * @param p The numerator (must fit in a signed 32-bit integer)
   * @param q The denominator (must fit in a signed 32-bit integer, non-zero)
   */
  constructor(p: number, q: number = 1) {
    super();
    if (q === 0) {
      throw new Error('Rational denominator cannot be zero');
    }
    this.p = Math.trunc(p);
    this.q = Math.trunc(q);

    // Create WASM-backed rational
    const wasm = getWasmModule();
    const obj = createBasic();
    const code = wasm._rational_set_si(obj.getPtr(), this.p, this.q);
    checkException(code);
    this._obj = obj;
  }

  protected _fallbackString(): string {
    return this.p + "/" + this.q;
  }

  /**
   * @internal Create a Rational from an existing WASM object.
   */
  static _fromWasm(obj: SymEngineObject): Rational {
    const rat = Object.create(Rational.prototype) as Rational;
    (rat as any)._obj = obj;
    // Parse the string representation "p/q" or just "p" for integers
    const str = obj.toString();
    if (str.includes('/')) {
      const parts = str.split('/');
      (rat as any).p = parseInt(parts[0], 10);
      (rat as any).q = parseInt(parts[1], 10);
    } else {
      (rat as any).p = parseInt(str, 10);
      (rat as any).q = 1;
    }
    return rat;
  }
}

/**
 * Floating-point number. Mirrors sympy.Float.
 * Note: SymEngine uses double precision internally (no arbitrary precision without MPFR).
 */
export class Float extends Expr {
  readonly value: number;

  /**
   * Create a floating-point number.
   * @param value The numeric value
   * @param _precision Ignored (SymEngine uses double precision without MPFR)
   */
  constructor(value: number, _precision?: number) {
    super();
    this.value = value;

    // Create WASM-backed real double
    const wasm = getWasmModule();
    const obj = createBasic();
    const code = wasm._real_double_set_d(obj.getPtr(), this.value);
    checkException(code);
    this._obj = obj;
  }

  protected _fallbackString(): string {
    return String(this.value);
  }

  /**
   * @internal Create a Float from an existing WASM object.
   */
  static _fromWasm(obj: SymEngineObject): Float {
    const f = Object.create(Float.prototype) as Float;
    (f as any)._obj = obj;
    (f as any).value = parseFloat(obj.toString());
    return f;
  }
}

/**
 * Exact complex number re + im*I.
 * Uses SymEngine's Complex class which stores rational real and imaginary parts.
 */
export class Complex extends Expr {
  readonly re: Integer | Rational;
  readonly im: Integer | Rational;

  /**
   * Create an exact complex number re + im*I.
   * @param re Real part (number, Integer, or Rational)
   * @param im Imaginary part (number, Integer, or Rational)
   */
  constructor(re: number | Integer | Rational, im: number | Integer | Rational) {
    super();

    // Convert numbers to Integer
    const reExpr = typeof re === 'number' ? new Integer(re) : re;
    const imExpr = typeof im === 'number' ? new Integer(im) : im;

    this.re = reExpr;
    this.im = imExpr;

    // Create WASM-backed complex
    const wasm = getWasmModule();
    const obj = createBasic();
    const code = wasm._complex_set(
      obj.getPtr(),
      (reExpr as any)._obj.getPtr(),
      (imExpr as any)._obj.getPtr()
    );
    checkException(code);
    this._obj = obj;
  }

  protected _fallbackString(): string {
    return `${this.re} + ${this.im}*I`;
  }

  /**
   * @internal Create a Complex from an existing WASM object.
   * Note: re/im fields are placeholders - full parsing not implemented.
   */
  static _fromWasm(obj: SymEngineObject): Complex {
    const c = Object.create(Complex.prototype) as Complex;
    (c as any)._obj = obj;
    // Placeholder values - full parsing of "a + b*I" format not implemented
    (c as any).re = new Integer(0);
    (c as any).im = new Integer(0);
    return c;
  }
}

/**
 * Symbolic addition. Mirrors sympy.Add.
 * Note: Use the add() function to create Add expressions.
 */
export class Add extends Expr {
  private _args: Expr[] | null = null;

  /** Get the terms of this addition (lazily extracted from WASM) */
  get args(): Expr[] {
    if (!this._args) {
      this._args = this._extractArgs();
    }
    return this._args;
  }

  private _extractArgs(): Expr[] {
    if (!this._obj) return [];
    const wasm = getWasmModule();
    const vec = new SymEngineVec();
    try {
      const code = wasm._basic_get_args(this._obj.getPtr(), vec.getPtr());
      checkException(code);
      const result: Expr[] = [];
      const count = vec.size();
      for (let i = 0; i < count; i++) {
        result.push(exprFromWasm(vec.get(i)));
      }
      return result;
    } finally {
      vec.free();
    }
  }

  protected _fallbackString(): string {
    return this.args.map((a) => a.toString()).join(' + ');
  }

  /**
   * @internal Create an Add from an existing WASM object.
   */
  static _fromWasm(obj: SymEngineObject): Add {
    const add = Object.create(Add.prototype) as Add;
    (add as any)._obj = obj;
    (add as any)._args = null;
    return add;
  }
}

/**
 * Symbolic multiplication. Mirrors sympy.Mul.
 * Note: Use the mul() function to create Mul expressions.
 */
export class Mul extends Expr {
  private _args: Expr[] | null = null;

  /** Get the factors of this multiplication (lazily extracted from WASM) */
  get args(): Expr[] {
    if (!this._args) {
      this._args = this._extractArgs();
    }
    return this._args;
  }

  private _extractArgs(): Expr[] {
    if (!this._obj) return [];
    const wasm = getWasmModule();
    const vec = new SymEngineVec();
    try {
      const code = wasm._basic_get_args(this._obj.getPtr(), vec.getPtr());
      checkException(code);
      const result: Expr[] = [];
      const count = vec.size();
      for (let i = 0; i < count; i++) {
        result.push(exprFromWasm(vec.get(i)));
      }
      return result;
    } finally {
      vec.free();
    }
  }

  protected _fallbackString(): string {
    return this.args.map((a) => a.toString()).join('*');
  }

  /**
   * @internal Create a Mul from an existing WASM object.
   */
  static _fromWasm(obj: SymEngineObject): Mul {
    const mul = Object.create(Mul.prototype) as Mul;
    (mul as any)._obj = obj;
    (mul as any)._args = null;
    return mul;
  }
}

/**
 * Symbolic exponentiation. Mirrors sympy.Pow.
 * Note: Use the pow() function to create Pow expressions.
 */
export class Pow extends Expr {
  private _base: Expr | null = null;
  private _exp: Expr | null = null;

  /** Get the base of this power */
  get base(): Expr {
    if (!this._base) {
      this._extractBaseExp();
    }
    return this._base!;
  }

  /** Get the exponent of this power */
  get exp(): Expr {
    if (!this._exp) {
      this._extractBaseExp();
    }
    return this._exp!;
  }

  private _extractBaseExp(): void {
    if (!this._obj) {
      this._base = new Integer(0);
      this._exp = new Integer(0);
      return;
    }
    const wasm = getWasmModule();
    const vec = new SymEngineVec();
    try {
      const code = wasm._basic_get_args(this._obj.getPtr(), vec.getPtr());
      checkException(code);
      if (vec.size() >= 2) {
        this._base = exprFromWasm(vec.get(0));
        this._exp = exprFromWasm(vec.get(1));
      } else {
        this._base = new Integer(0);
        this._exp = new Integer(0);
      }
    } finally {
      vec.free();
    }
  }

  protected _fallbackString(): string {
    return this.base + '**' + this.exp;
  }

  /**
   * @internal Create a Pow from an existing WASM object.
   */
  static _fromWasm(obj: SymEngineObject): Pow {
    const p = Object.create(Pow.prototype) as Pow;
    (p as any)._obj = obj;
    (p as any)._base = null;
    (p as any)._exp = null;
    return p;
  }
}

// ============================================================================
// Constant Classes (WASM-backed)
// ============================================================================

/**
 * Symbolic mathematical constant (pi, E, EulerGamma, Catalan, GoldenRatio).
 * WASM-backed version that supports arithmetic operations.
 */
export class Constant extends Expr {
  readonly name: string;

  private constructor(name: string) {
    super();
    this.name = name;
  }

  protected _fallbackString(): string {
    return this.name;
  }

  /** Create pi constant */
  static pi(): Constant {
    const c = new Constant('pi');
    const obj = createBasic();
    getWasmModule()._basic_const_pi(obj.getPtr());
    (c as any)._obj = obj;
    return c;
  }

  /** Create E (Euler's number) constant */
  static E(): Constant {
    const c = new Constant('E');
    const obj = createBasic();
    getWasmModule()._basic_const_E(obj.getPtr());
    (c as any)._obj = obj;
    return c;
  }

  /** Create EulerGamma constant */
  static EulerGamma(): Constant {
    const c = new Constant('EulerGamma');
    const obj = createBasic();
    getWasmModule()._basic_const_EulerGamma(obj.getPtr());
    (c as any)._obj = obj;
    return c;
  }

  /** Create Catalan constant */
  static Catalan(): Constant {
    const c = new Constant('Catalan');
    const obj = createBasic();
    getWasmModule()._basic_const_Catalan(obj.getPtr());
    (c as any)._obj = obj;
    return c;
  }

  /** Create GoldenRatio constant */
  static GoldenRatio(): Constant {
    const c = new Constant('GoldenRatio');
    const obj = createBasic();
    getWasmModule()._basic_const_GoldenRatio(obj.getPtr());
    (c as any)._obj = obj;
    return c;
  }

  /** @internal Create from WASM object */
  static _fromWasm(obj: SymEngineObject): Constant {
    const c = Object.create(Constant.prototype) as Constant;
    (c as any)._obj = obj;
    (c as any).name = obj.toString();
    return c;
  }
}

/**
 * Imaginary unit (i = sqrt(-1)).
 * Special constant backed by SymEngine's I constant.
 */
export class ImaginaryUnit extends Expr {
  private constructor() {
    super();
  }

  protected _fallbackString(): string {
    return 'I';
  }

  /** Create the imaginary unit I */
  static create(): ImaginaryUnit {
    const i = new ImaginaryUnit();
    const obj = createBasic();
    getWasmModule()._basic_const_I(obj.getPtr());
    (i as any)._obj = obj;
    return i;
  }

  /** @internal Create from WASM object */
  static _fromWasm(obj: SymEngineObject): ImaginaryUnit {
    const i = Object.create(ImaginaryUnit.prototype) as ImaginaryUnit;
    (i as any)._obj = obj;
    return i;
  }
}

/**
 * Infinity (positive, negative, or complex).
 * WASM-backed version supporting arithmetic operations.
 */
export class Infinity_ extends Expr {
  readonly direction: 'positive' | 'negative' | 'complex';

  private constructor(direction: 'positive' | 'negative' | 'complex') {
    super();
    this.direction = direction;
  }

  protected _fallbackString(): string {
    switch (this.direction) {
      case 'positive':
        return 'oo';
      case 'negative':
        return '-oo';
      case 'complex':
        return 'zoo';
    }
  }

  /** Create positive infinity */
  static positive(): Infinity_ {
    const inf = new Infinity_('positive');
    const obj = createBasic();
    getWasmModule()._basic_const_infinity(obj.getPtr());
    (inf as any)._obj = obj;
    return inf;
  }

  /** Create negative infinity */
  static negative(): Infinity_ {
    const inf = new Infinity_('negative');
    const obj = createBasic();
    getWasmModule()._basic_const_neginfinity(obj.getPtr());
    (inf as any)._obj = obj;
    return inf;
  }

  /** Create complex infinity */
  static complex(): Infinity_ {
    const inf = new Infinity_('complex');
    const obj = createBasic();
    getWasmModule()._basic_const_complex_infinity(obj.getPtr());
    (inf as any)._obj = obj;
    return inf;
  }

  /** @internal Create from WASM object */
  static _fromWasm(obj: SymEngineObject): Infinity_ {
    const str = obj.toString();
    let direction: 'positive' | 'negative' | 'complex';
    if (str === '-oo') direction = 'negative';
    else if (str === 'zoo') direction = 'complex';
    else direction = 'positive';

    const inf = Object.create(Infinity_.prototype) as Infinity_;
    (inf as any)._obj = obj;
    (inf as any).direction = direction;
    return inf;
  }
}

/**
 * Not a Number (undefined/indeterminate result).
 * WASM-backed version.
 */
export class NaN_ extends Expr {
  private constructor() {
    super();
  }

  protected _fallbackString(): string {
    return 'nan';
  }

  /** Create NaN */
  static create(): NaN_ {
    const n = new NaN_();
    const obj = createBasic();
    getWasmModule()._basic_const_nan(obj.getPtr());
    (n as any)._obj = obj;
    return n;
  }

  /** @internal Create from WASM object */
  static _fromWasm(obj: SymEngineObject): NaN_ {
    const n = Object.create(NaN_.prototype) as NaN_;
    (n as any)._obj = obj;
    return n;
  }
}

/**
 * Generic expression wrapper for types not yet fully implemented.
 * @internal
 */
class GenericExpr extends Expr {
  protected _fallbackString(): string {
    return this._obj ? this._obj.toString() : '<unknown>';
  }

  static _fromWasm(obj: SymEngineObject): GenericExpr {
    const e = Object.create(GenericExpr.prototype) as GenericExpr;
    (e as any)._obj = obj;
    return e;
  }
}

/**
 * Create an Expr subclass instance from a WASM object based on its type.
 * Used internally and by the matrices module to convert WASM basic pointers to Expr.
 */
export function exprFromWasm(obj: SymEngineObject): Expr {
  const typeId = obj.getType();
  switch (typeId) {
    case SymEngineTypeID.SYMENGINE_INTEGER:
      return Integer._fromWasm(obj);
    case SymEngineTypeID.SYMENGINE_RATIONAL:
      return Rational._fromWasm(obj);
    case SymEngineTypeID.SYMENGINE_REAL_DOUBLE:
      return Float._fromWasm(obj);
    case SymEngineTypeID.SYMENGINE_COMPLEX_DOUBLE:
      return GenericExpr._fromWasm(obj); // ComplexDouble - handled via evalfComplex
    case SymEngineTypeID.SYMENGINE_COMPLEX:
      return Complex._fromWasm(obj);
    case SymEngineTypeID.SYMENGINE_SYMBOL:
      return Symbol._fromWasm(obj);
    case SymEngineTypeID.SYMENGINE_ADD:
      return Add._fromWasm(obj);
    case SymEngineTypeID.SYMENGINE_MUL:
      return Mul._fromWasm(obj);
    case SymEngineTypeID.SYMENGINE_POW:
      return Pow._fromWasm(obj);
    case SymEngineTypeID.SYMENGINE_CONSTANT:
      return Constant._fromWasm(obj);
    case SymEngineTypeID.SYMENGINE_INFTY:
      return Infinity_._fromWasm(obj);
    case SymEngineTypeID.SYMENGINE_NOT_A_NUMBER:
      return NaN_._fromWasm(obj);
    default:
      return GenericExpr._fromWasm(obj);
  }
}

// ============================================================================
// Helper Functions
// ============================================================================

/**
 * Parse a complex number string representation from SymEngine.
 * Handles formats like: "1.0 + 2.0*I", "1.5 - 2.3*I", "2.3*I", "-2.3*I", "1.5"
 * Note: SymEngine uses format with spaces around +/- operators
 * @internal
 */
function parseComplexString(str: string): { real: number; imag: number } {
  const s = str.trim();

  // Pure imaginary: just "I" or "-I"
  if (s === 'I') {
    return { real: 0, imag: 1 };
  }
  if (s === '-I') {
    return { real: 0, imag: -1 };
  }

  // Check for *I pattern (pure imaginary): "2.3*I" or "-2.3*I"
  const pureImagMatch = s.match(/^(-?[\d.eE+-]+)\*I$/);
  if (pureImagMatch) {
    return { real: 0, imag: parseFloat(pureImagMatch[1]) };
  }

  // Complex with both parts with spaces: "real + imag*I" or "real - imag*I"
  // SymEngine uses format like "1.0 + 2.0*I"
  const complexWithSpacesMatch = s.match(/^(-?[\d.eE+-]+)\s*([+-])\s*([\d.eE+-]+)\*I$/);
  if (complexWithSpacesMatch) {
    const real = parseFloat(complexWithSpacesMatch[1]);
    const sign = complexWithSpacesMatch[2] === '-' ? -1 : 1;
    const imag = sign * parseFloat(complexWithSpacesMatch[3]);
    return { real, imag };
  }

  // Handle case where imag coefficient is 1: "real + I" or "real - I"
  const complexOneMatch = s.match(/^(-?[\d.eE+-]+)\s*([+-])\s*I$/);
  if (complexOneMatch) {
    const real = parseFloat(complexOneMatch[1]);
    const imag = complexOneMatch[2] === '-' ? -1 : 1;
    return { real, imag };
  }

  // Pure real number
  return { real: parseFloat(s), imag: 0 };
}

// ============================================================================
// Arithmetic Helper Functions
// ============================================================================

/**
 * Add two expressions: a + b
 * @param a First operand
 * @param b Second operand
 * @returns The sum (may be simplified by SymEngine)
 */
export function add(a: Expr, b: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  const code = wasm._basic_add(obj.getPtr(), a.getWasmPtr(), b.getWasmPtr());
  checkException(code);
  return exprFromWasm(obj);
}

/**
 * Subtract two expressions: a - b
 * @param a First operand
 * @param b Second operand
 * @returns The difference (may be simplified by SymEngine)
 */
export function sub(a: Expr, b: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  const code = wasm._basic_sub(obj.getPtr(), a.getWasmPtr(), b.getWasmPtr());
  checkException(code);
  return exprFromWasm(obj);
}

/**
 * Multiply two expressions: a * b
 * @param a First operand
 * @param b Second operand
 * @returns The product (may be simplified by SymEngine)
 */
export function mul(a: Expr, b: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  const code = wasm._basic_mul(obj.getPtr(), a.getWasmPtr(), b.getWasmPtr());
  checkException(code);
  return exprFromWasm(obj);
}

/**
 * Divide two expressions: a / b
 * @param a Numerator
 * @param b Denominator
 * @returns The quotient (may be simplified by SymEngine)
 */
export function div(a: Expr, b: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  const code = wasm._basic_div(obj.getPtr(), a.getWasmPtr(), b.getWasmPtr());
  checkException(code);
  return exprFromWasm(obj);
}

/**
 * Raise to a power: base ** exp
 * @param base The base
 * @param exp The exponent
 * @returns The power (may be simplified by SymEngine)
 */
export function pow(base: Expr, exp: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  const code = wasm._basic_pow(obj.getPtr(), base.getWasmPtr(), exp.getWasmPtr());
  checkException(code);
  return exprFromWasm(obj);
}

/**
 * Negate an expression: -a
 * @param a The operand
 * @returns The negation
 */
export function neg(a: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  const code = wasm._basic_neg(obj.getPtr(), a.getWasmPtr());
  checkException(code);
  return exprFromWasm(obj);
}

// ============================================================================
// Elementary Function Factories
// ============================================================================

/**
 * Internal helper to create a 1-argument function wrapper.
 * @internal
 */
function makeOneArgFunc(
  wasmFnName: keyof ReturnType<typeof getWasmModule>
): (x: Expr) => Expr {
  return (x: Expr): Expr => {
    const wasm = getWasmModule();
    const obj = createBasic();
    try {
      const fn = wasm[wasmFnName] as (result: number, arg: number) => number;
      const code = fn.call(wasm, obj.getPtr(), x.getWasmPtr());
      checkException(code);
      return exprFromWasm(obj);
    } catch (e) {
      obj.free();
      throw e;
    }
  };
}

/**
 * Internal helper to create a 2-argument function wrapper.
 * @internal
 */
function makeTwoArgFunc(
  wasmFnName: keyof ReturnType<typeof getWasmModule>
): (a: Expr, b: Expr) => Expr {
  return (a: Expr, b: Expr): Expr => {
    const wasm = getWasmModule();
    const obj = createBasic();
    try {
      const fn = wasm[wasmFnName] as (
        result: number,
        arg1: number,
        arg2: number
      ) => number;
      const code = fn.call(wasm, obj.getPtr(), a.getWasmPtr(), b.getWasmPtr());
      checkException(code);
      return exprFromWasm(obj);
    } catch (e) {
      obj.free();
      throw e;
    }
  };
}

// ============================================================================
// Trigonometric Functions
// ============================================================================

/** Sine: sin(x) */
export const sin = makeOneArgFunc('_basic_sin');
/** Cosine: cos(x) */
export const cos = makeOneArgFunc('_basic_cos');
/** Tangent: tan(x) */
export const tan = makeOneArgFunc('_basic_tan');
/** Cotangent: cot(x) = 1/tan(x) */
export const cot = makeOneArgFunc('_basic_cot');
/** Secant: sec(x) = 1/cos(x) */
export const sec = makeOneArgFunc('_basic_sec');
/** Cosecant: csc(x) = 1/sin(x) */
export const csc = makeOneArgFunc('_basic_csc');

// ============================================================================
// Inverse Trigonometric Functions
// ============================================================================

/** Arcsine: asin(x) */
export const asin = makeOneArgFunc('_basic_asin');
/** Arccosine: acos(x) */
export const acos = makeOneArgFunc('_basic_acos');
/** Arctangent: atan(x) */
export const atan = makeOneArgFunc('_basic_atan');
/** Arc-cotangent: acot(x) */
export const acot = makeOneArgFunc('_basic_acot');
/** Arc-secant: asec(x) */
export const asec = makeOneArgFunc('_basic_asec');
/** Arc-cosecant: acsc(x) */
export const acsc = makeOneArgFunc('_basic_acsc');
/** Two-argument arctangent: atan2(y, x) */
export const atan2 = makeTwoArgFunc('_basic_atan2');

// ============================================================================
// Hyperbolic Functions
// ============================================================================

/** Hyperbolic sine: sinh(x) */
export const sinh = makeOneArgFunc('_basic_sinh');
/** Hyperbolic cosine: cosh(x) */
export const cosh = makeOneArgFunc('_basic_cosh');
/** Hyperbolic tangent: tanh(x) */
export const tanh = makeOneArgFunc('_basic_tanh');
/** Hyperbolic cotangent: coth(x) */
export const coth = makeOneArgFunc('_basic_coth');
/** Hyperbolic secant: sech(x) */
export const sech = makeOneArgFunc('_basic_sech');
/** Hyperbolic cosecant: csch(x) */
export const csch = makeOneArgFunc('_basic_csch');

// ============================================================================
// Inverse Hyperbolic Functions
// ============================================================================

/** Inverse hyperbolic sine: asinh(x) */
export const asinh = makeOneArgFunc('_basic_asinh');
/** Inverse hyperbolic cosine: acosh(x) */
export const acosh = makeOneArgFunc('_basic_acosh');
/** Inverse hyperbolic tangent: atanh(x) */
export const atanh = makeOneArgFunc('_basic_atanh');
/** Inverse hyperbolic cotangent: acoth(x) */
export const acoth = makeOneArgFunc('_basic_acoth');
/** Inverse hyperbolic secant: asech(x) */
export const asech = makeOneArgFunc('_basic_asech');
/** Inverse hyperbolic cosecant: acsch(x) */
export const acsch = makeOneArgFunc('_basic_acsch');

// ============================================================================
// Exponential & Logarithmic Functions
// ============================================================================

/** Exponential: exp(x) = e^x */
export const exp = makeOneArgFunc('_basic_exp');
/** Natural logarithm: log(x) = ln(x) */
export const log = makeOneArgFunc('_basic_log');
/** Square root: sqrt(x) */
export const sqrt = makeOneArgFunc('_basic_sqrt');
/** Cube root: cbrt(x) */
export const cbrt = makeOneArgFunc('_basic_cbrt');
/** Lambert W function: lambertw(x) (principal branch) */
export const lambertw = makeOneArgFunc('_basic_lambertw');

// ============================================================================
// Other Mathematical Functions
// ============================================================================

/** Absolute value: abs(x) */
export const abs = makeOneArgFunc('_basic_abs');
/** Sign function: sign(x) returns -1, 0, or 1 */
export const sign = makeOneArgFunc('_basic_sign');
/** Floor function: floor(x) = largest integer ≤ x */
export const floor = makeOneArgFunc('_basic_floor');
/** Ceiling function: ceiling(x) = smallest integer ≥ x */
export const ceiling = makeOneArgFunc('_basic_ceiling');

// ============================================================================
// Special Functions
// ============================================================================

/** Gamma function: gamma(x) = Γ(x) = (x-1)! for positive integers */
export const gamma = makeOneArgFunc('_basic_gamma');
/** Log-gamma function: loggamma(x) = log(Γ(x)) */
export const loggamma = makeOneArgFunc('_basic_loggamma');
/** Error function: erf(x) */
export const erf = makeOneArgFunc('_basic_erf');
/** Complementary error function: erfc(x) = 1 - erf(x) */
export const erfc = makeOneArgFunc('_basic_erfc');
/** Riemann zeta function: zeta(s) */
export const zeta = makeOneArgFunc('_basic_zeta');
/** Dirichlet eta function: dirichlet_eta(s) */
export const dirichlet_eta = makeOneArgFunc('_basic_dirichlet_eta');
/** Beta function: beta(a, b) = Γ(a)Γ(b)/Γ(a+b) */
export const beta = makeTwoArgFunc('_basic_beta');
/** Lower incomplete gamma: lowergamma(s, x) = γ(s, x) */
export const lowergamma = makeTwoArgFunc('_basic_lowergamma');
/** Upper incomplete gamma: uppergamma(s, x) = Γ(s, x) */
export const uppergamma = makeTwoArgFunc('_basic_uppergamma');
/** Polygamma function: polygamma(n, x) = ψ^(n)(x) */
export const polygamma = makeTwoArgFunc('_basic_polygamma');
/** Kronecker delta: kronecker_delta(i, j) = 1 if i=j, else 0 */
export const kronecker_delta = makeTwoArgFunc('_basic_kronecker_delta');

// ============================================================================
// Additional Functions (Phase 2.1b)
// ============================================================================

/** Digamma function: digamma(x) = ψ(x) = d/dx ln(Γ(x)) */
export const digamma = makeOneArgFunc('_basic_digamma');

/** Complex conjugate: conjugate(a + bi) = a - bi */
export const conjugate = makeOneArgFunc('_basic_conjugate');

/** Real part: re(x) extracts the real component of a complex expression */
export function re(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._complex_base_real_part(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/** Imaginary part: im(x) extracts the imaginary component of a complex expression */
export function im(x: Expr): Expr {
  const wasm = getWasmModule();
  const obj = createBasic();
  try {
    const code = wasm._complex_base_imaginary_part(obj.getPtr(), x.getWasmPtr());
    checkException(code);
    return exprFromWasm(obj);
  } catch (e) {
    obj.free();
    throw e;
  }
}

/** Argument (phase) of complex number: arg(x) = atan2(im(x), re(x)) */
export function arg(x: Expr): Expr {
  return atan2(im(x), re(x));
}

/** Maximum of multiple values: Max(a, b, c, ...) */
export function Max(...args: Expr[]): Expr {
  const wasm = getWasmModule();
  const vec = new SymEngineVec();
  const result = createBasic();

  try {
    for (const a of args) {
      wasm._vecbasic_push_back(vec.getPtr(), a.getWasmPtr());
    }
    const code = wasm._basic_max(result.getPtr(), vec.getPtr());
    checkException(code);
    return exprFromWasm(result);
  } catch (e) {
    result.free();
    throw e;
  } finally {
    vec.free();
  }
}

/** Minimum of multiple values: Min(a, b, c, ...) */
export function Min(...args: Expr[]): Expr {
  const wasm = getWasmModule();
  const vec = new SymEngineVec();
  const result = createBasic();

  try {
    for (const a of args) {
      wasm._vecbasic_push_back(vec.getPtr(), a.getWasmPtr());
    }
    const code = wasm._basic_min(result.getPtr(), vec.getPtr());
    checkException(code);
    return exprFromWasm(result);
  } catch (e) {
    result.free();
    throw e;
  } finally {
    vec.free();
  }
}

// ============================================================================
// Calculus Functions (Phase 2.2)
// ============================================================================

/**
 * Differentiate an expression with respect to one or more symbols.
 *
 * Supports multiple calling conventions:
 * - `diff(expr, x)` — First derivative with respect to x
 * - `diff(expr, x, 2)` — Second derivative with respect to x
 * - `diff(expr, x, y)` — Mixed partial: ∂²f/∂x∂y
 * - `diff(expr, x, y, z)` — Mixed partial: ∂³f/∂x∂y∂z
 * - `diff(expr, x, 2, y, 3)` — ∂⁵f/∂x²∂y³
 *
 * @param expr The expression to differentiate
 * @param args Symbols and optional derivative orders
 * @returns The derivative expression
 *
 * @example
 * const [x, y] = symbols('x y');
 * diff(pow(x, 3), x);           // → 3*x**2
 * diff(sin(x), x);              // → cos(x)
 * diff(pow(x, 3), x, 2);        // → 6*x (second derivative)
 * diff(mul(x, y), x, y);        // → 1 (∂²(xy)/∂x∂y)
 * diff(pow(x, 2), x, 2, y, 1);  // → 0 (no y dependence)
 */
export function diff(expr: Expr, ...args: (Expr | number)[]): Expr {
  if (args.length === 0) {
    throw new Error('diff() requires at least one symbol argument');
  }

  const wasm = getWasmModule();
  let current = expr;

  // Parse args: alternating [symbol, count?] pattern
  // e.g., [x], [x, 2], [x, y], [x, 2, y, 3]
  let i = 0;
  while (i < args.length) {
    const arg = args[i];

    // Must be an Expr (symbol)
    if (!(arg instanceof Expr)) {
      throw new Error(
        `diff() argument at position ${i + 1} must be a symbol, got number without preceding symbol`
      );
    }

    const symbol = arg;
    i++;

    // Check if next argument is a number (derivative order)
    let order = 1;
    if (i < args.length && typeof args[i] === 'number') {
      order = args[i] as number;
      i++;
    }

    if (order < 0) {
      throw new Error('Derivative order must be non-negative');
    }

    // Apply differentiation 'order' times with respect to 'symbol'
    for (let j = 0; j < order; j++) {
      const result = createBasic();
      try {
        const code = wasm._basic_diff(
          result.getPtr(),
          current.getWasmPtr(),
          symbol.getWasmPtr()
        );
        checkException(code);
        current = exprFromWasm(result);
      } catch (e) {
        result.free();
        throw e;
      }
    }
  }

  return current;
}

/**
 * Compute Taylor series expansion around x=0.
 *
 * @param expr The expression to expand
 * @param x The expansion variable
 * @param x0 The expansion point (currently only 0 is supported)
 * @param n Number of terms (precision = n, so terms up to x^(n-1))
 * @returns The series expansion as a polynomial (without O() term)
 *
 * @example
 * series(sin(x), x);           // x - x³/6 + x⁵/120
 * series(exp(x), x, 0, 4);     // 1 + x + x²/2 + x³/6
 */
export function series(
  expr: Expr,
  x: Expr,
  x0: number | Expr = 0,
  n: number = 6
): Expr {
  // Validate x0 = 0 (only supported case)
  const x0Val = typeof x0 === 'number' ? x0 : null;
  if (x0Val !== 0) {
    throw new Error('series() currently only supports expansion around x=0');
  }

  const wasm = getWasmModule();
  const result = createBasic();

  try {
    const code = wasm._basic_series(
      result.getPtr(),
      expr.getWasmPtr(),
      x.getWasmPtr(),
      n
    );
    checkException(code);
    return exprFromWasm(result);
  } catch (e) {
    result.free();
    throw e;
  }
}

// ============================================================================
// Simplification Functions (Phase 2.4)
// ============================================================================

/**
 * Expand an expression by distributing multiplication over addition.
 *
 * @param expr The expression to expand
 * @returns The expanded expression
 *
 * @example
 * expand((x + 1)**2);         // x**2 + 2*x + 1
 * expand((a + b) * (c + d));  // a*c + a*d + b*c + b*d
 */
export function expand(expr: Expr): Expr {
  const wasm = getWasmModule();
  const result = createBasic();

  try {
    const code = wasm._basic_expand(result.getPtr(), expr.getWasmPtr());
    checkException(code);
    return exprFromWasm(result);
  } catch (e) {
    result.free();
    throw e;
  }
}

/**
 * Simplify an expression using heuristics.
 *
 * The simplification includes:
 * - Converting csc(x)**(-1) to sin(x)
 * - Converting sec(x)**(-1) to cos(x)
 * - Converting cot(x)**(-1) to tan(x)
 *
 * @param expr The expression to simplify
 * @returns The simplified expression
 *
 * @example
 * simplify(csc(x)**(-1));  // sin(x)
 * simplify(sec(x)**(-1));  // cos(x)
 */
export function simplify(expr: Expr): Expr {
  const wasm = getWasmModule();
  const result = createBasic();

  try {
    const code = wasm._basic_simplify(result.getPtr(), expr.getWasmPtr());
    checkException(code);
    return exprFromWasm(result);
  } catch (e) {
    result.free();
    throw e;
  }
}

/**
 * Extract the numerator of an expression.
 *
 * @param expr The expression
 * @returns The numerator
 *
 * @example
 * numer(x / y);  // x
 * numer(3);      // 3
 * numer(Rational(3, 4));  // 3
 */
export function numer(expr: Expr): Expr {
  const wasm = getWasmModule();
  const numerResult = createBasic();
  const denomResult = createBasic();

  try {
    const code = wasm._basic_as_numer_denom(
      numerResult.getPtr(),
      denomResult.getPtr(),
      expr.getWasmPtr()
    );
    checkException(code);
    denomResult.free();
    return exprFromWasm(numerResult);
  } catch (e) {
    numerResult.free();
    denomResult.free();
    throw e;
  }
}

/**
 * Extract the denominator of an expression.
 *
 * @param expr The expression
 * @returns The denominator
 *
 * @example
 * denom(x / y);  // y
 * denom(3);      // 1
 * denom(Rational(3, 4));  // 4
 */
export function denom(expr: Expr): Expr {
  const wasm = getWasmModule();
  const numerResult = createBasic();
  const denomResult = createBasic();

  try {
    const code = wasm._basic_as_numer_denom(
      numerResult.getPtr(),
      denomResult.getPtr(),
      expr.getWasmPtr()
    );
    checkException(code);
    numerResult.free();
    return exprFromWasm(denomResult);
  } catch (e) {
    numerResult.free();
    denomResult.free();
    throw e;
  }
}

/**
 * Simplify trigonometric expressions.
 * Currently uses general simplify() internally.
 *
 * @param expr The expression to simplify
 * @returns The simplified expression
 */
export function trigsimp(expr: Expr): Expr {
  return simplify(expr);
}

/**
 * Simplify expressions with radicals.
 * Currently uses general simplify() internally.
 *
 * @param expr The expression to simplify
 * @returns The simplified expression
 */
export function radsimp(expr: Expr): Expr {
  return simplify(expr);
}

/**
 * Simplify combinatorial expressions with powers.
 * Currently uses general simplify() internally.
 *
 * @param expr The expression to simplify
 * @returns The simplified expression
 */
export function powsimp(expr: Expr): Expr {
  return simplify(expr);
}

/**
 * Rewrite trigonometric functions as exponentials using Euler's formula.
 * sin(x) → (e^(ix) - e^(-ix)) / (2i)
 * cos(x) → (e^(ix) + e^(-ix)) / 2
 *
 * @param expr The expression to rewrite
 * @returns Expression with trig functions as exponentials
 */
export function rewrite_as_exp(expr: Expr): Expr {
  const wasm = getWasmModule();
  const result = createBasic();
  checkException(wasm._basic_rewrite_as_exp(result.getPtr(), expr.getWasmPtr()));
  return exprFromWasm(result);
}

/**
 * Rewrite trigonometric functions in terms of sine.
 * cos(x) → sin(x + π/2)
 * tan(x) → sin(x) / sin(π/2 - x)
 *
 * @param expr The expression to rewrite
 * @returns Expression with trig functions in terms of sine
 */
export function rewrite_as_sin(expr: Expr): Expr {
  const wasm = getWasmModule();
  const result = createBasic();
  checkException(wasm._basic_rewrite_as_sin(result.getPtr(), expr.getWasmPtr()));
  return exprFromWasm(result);
}

/**
 * Rewrite trigonometric functions in terms of cosine.
 * sin(x) → cos(π/2 - x)
 * tan(x) → cos(π/2 - x) / cos(x)
 *
 * @param expr The expression to rewrite
 * @returns Expression with trig functions in terms of cosine
 */
export function rewrite_as_cos(expr: Expr): Expr {
  const wasm = getWasmModule();
  const result = createBasic();
  checkException(wasm._basic_rewrite_as_cos(result.getPtr(), expr.getWasmPtr()));
  return exprFromWasm(result);
}

/**
 * Extract real and imaginary parts of a complex expression.
 * For x + iy, returns { real: x, imag: y }
 *
 * @param expr The expression to analyze
 * @returns Object with real and imag Expr properties
 */
export function as_real_imag(expr: Expr): { real: Expr; imag: Expr } {
  const wasm = getWasmModule();
  const realResult = createBasic();
  const imagResult = createBasic();
  checkException(wasm._basic_as_real_imag(realResult.getPtr(), imagResult.getPtr(), expr.getWasmPtr()));
  return {
    real: exprFromWasm(realResult),
    imag: exprFromWasm(imagResult),
  };
}

/**
 * Expand trigonometric functions.
 * Applies rewrite_as_exp followed by expand.
 *
 * @param expr The expression to expand
 * @returns Expression with expanded trig functions
 */
export function expand_trig(expr: Expr): Expr {
  return expand(rewrite_as_exp(expr));
}

/**
 * Expand complex expressions, extracting real and imaginary parts.
 * Alias for as_real_imag.
 *
 * @param expr The expression to expand
 * @returns Object with real and imag Expr properties
 */
export function expand_complex(expr: Expr): { real: Expr; imag: Expr } {
  return as_real_imag(expr);
}

// Private cache for singleton constants (lazy initialization)
let _zero: Integer | null = null;
let _one: Integer | null = null;
let _negOne: Integer | null = null;
let _half: Rational | null = null;
let _infinity: Infinity_ | null = null;
let _negInfinity: Infinity_ | null = null;
let _complexInfinity: Infinity_ | null = null;
let _nan: NaN_ | null = null;
let _pi: Constant | null = null;
let _E: Constant | null = null;
let _I: ImaginaryUnit | null = null;
let _eulerGamma: Constant | null = null;
let _catalan: Constant | null = null;
let _goldenRatio: Constant | null = null;

/**
 * Singleton-like namespace for special symbolic constants.
 * Mirrors sympy.S.
 * Constants are lazily initialized on first access.
 */
export const S = {
  /** The number zero. */
  get Zero(): Integer {
    if (!_zero) _zero = new Integer(0);
    return _zero;
  },
  /** The number one. */
  get One(): Integer {
    if (!_one) _one = new Integer(1);
    return _one;
  },
  /** Negative one. */
  get NegativeOne(): Integer {
    if (!_negOne) _negOne = new Integer(-1);
    return _negOne;
  },
  /** One half. */
  get Half(): Rational {
    if (!_half) _half = new Rational(1, 2);
    return _half;
  },
  /** Positive infinity. */
  get Infinity(): Infinity_ {
    if (!_infinity) _infinity = Infinity_.positive();
    return _infinity;
  },
  /** Negative infinity. */
  get NegativeInfinity(): Infinity_ {
    if (!_negInfinity) _negInfinity = Infinity_.negative();
    return _negInfinity;
  },
  /** Complex infinity (undirected). */
  get ComplexInfinity(): Infinity_ {
    if (!_complexInfinity) _complexInfinity = Infinity_.complex();
    return _complexInfinity;
  },
  /** Not a Number (undefined/indeterminate). */
  get NaN(): NaN_ {
    if (!_nan) _nan = NaN_.create();
    return _nan;
  },
};

/**
 * Helper to create a lazy proxy for a constant.
 * The proxy delegates all property access to the lazily-initialized constant.
 */
function lazyConstantProxy<T extends Expr>(getter: () => T): T {
  return new Proxy({} as T, {
    get(_target, prop) {
      const instance = getter();
      const value = (instance as any)[prop];
      if (typeof value === 'function') {
        return value.bind(instance);
      }
      return value;
    },
  });
}

/** Positive infinity. */
export const oo: Infinity_ = lazyConstantProxy(() => {
  if (!_infinity) _infinity = Infinity_.positive();
  return _infinity;
});

/** Pi (π). */
export const pi: Constant = lazyConstantProxy(() => {
  if (!_pi) _pi = Constant.pi();
  return _pi;
});

/** Euler's number (e). */
export const E: Constant = lazyConstantProxy(() => {
  if (!_E) _E = Constant.E();
  return _E;
});

/** Imaginary unit (i). */
export const I: ImaginaryUnit = lazyConstantProxy(() => {
  if (!_I) _I = ImaginaryUnit.create();
  return _I;
});

/** Euler-Mascheroni constant (γ ≈ 0.5772). */
export const EulerGamma: Constant = lazyConstantProxy(() => {
  if (!_eulerGamma) _eulerGamma = Constant.EulerGamma();
  return _eulerGamma;
});

/** Catalan's constant (≈ 0.9159). */
export const Catalan: Constant = lazyConstantProxy(() => {
  if (!_catalan) _catalan = Constant.Catalan();
  return _catalan;
});

/** Golden ratio (φ ≈ 1.618). */
export const GoldenRatio: Constant = lazyConstantProxy(() => {
  if (!_goldenRatio) _goldenRatio = Constant.GoldenRatio();
  return _goldenRatio;
});
