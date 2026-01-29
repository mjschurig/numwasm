/**
 * Core symbolic expression types.
 * @module core
 */

import { NotImplementedError } from '../errors.js';
import { loadWasmModule, getWasmModule } from '../wasm-loader.js';
import {
  SymEngineObject,
  SymEngineSet,
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
  protected getWasmPtr(): number {
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

  /** Substitute a symbol with another expression. */
  subs(_old: Expr, _new: Expr): Expr {
    throw new NotImplementedError('symwasm.core.Expr.subs');
  }

  /** Evaluate the expression numerically. */
  evalf(_precision?: number): number {
    throw new NotImplementedError('symwasm.core.Expr.evalf');
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

  constructor(value: number) {
    super();
    this.value = value;
    throw new NotImplementedError('symwasm.core.Integer');
  }

  protected _fallbackString(): string {
    return String(this.value);
  }
}

/** Exact rational number p/q. Mirrors sympy.Rational. */
export class Rational extends Expr {
  readonly p: number;
  readonly q: number;

  constructor(p: number, q: number) {
    super();
    this.p = p;
    this.q = q;
    throw new NotImplementedError('symwasm.core.Rational');
  }

  protected _fallbackString(): string {
    return this.p + "/" + this.q;
  }
}

/** Floating-point number with arbitrary precision. Mirrors sympy.Float. */
export class Float extends Expr {
  readonly value: number;

  constructor(value: number, _precision?: number) {
    super();
    this.value = value;
    throw new NotImplementedError('symwasm.core.Float');
  }

  protected _fallbackString(): string {
    return String(this.value);
  }
}

/** Symbolic addition. Mirrors sympy.Add. */
export class Add extends Expr {
  readonly args: Expr[];

  constructor(args: Expr[]) {
    super();
    this.args = args;
    throw new NotImplementedError('symwasm.core.Add');
  }

  protected _fallbackString(): string {
    return this.args.map(a => a.toString()).join(' + ');
  }
}

/** Symbolic multiplication. Mirrors sympy.Mul. */
export class Mul extends Expr {
  readonly args: Expr[];

  constructor(args: Expr[]) {
    super();
    this.args = args;
    throw new NotImplementedError('symwasm.core.Mul');
  }

  protected _fallbackString(): string {
    return this.args.map(a => a.toString()).join('*');
  }
}

/** Symbolic exponentiation. Mirrors sympy.Pow. */
export class Pow extends Expr {
  readonly base: Expr;
  readonly exp: Expr;

  constructor(base: Expr, exp: Expr) {
    super();
    this.base = base;
    this.exp = exp;
    throw new NotImplementedError('symwasm.core.Pow');
  }

  protected _fallbackString(): string {
    return this.base + '**' + this.exp;
  }
}

/**
 * @internal A sentinel expression for constants like pi, E, I, oo.
 * These are placeholder expressions without full WASM backing.
 */
class SentinelExpr extends Expr {
  private readonly _name: string;

  constructor(name: string) {
    super();
    this._name = name;
  }

  protected _fallbackString(): string {
    return this._name;
  }

  // Override to return empty array (no free symbols in constants)
  free_symbols(): Symbol[] {
    return [];
  }

  // Override equals to use reference equality for sentinels
  equals(other: Expr): boolean {
    return this === other;
  }
}

/**
 * Singleton-like namespace for special symbolic constants.
 * Mirrors sympy.S.
 */
export const S = {
  /** The number zero. */
  get Zero(): Expr { throw new NotImplementedError('symwasm.core.S.Zero'); },
  /** The number one. */
  get One(): Expr { throw new NotImplementedError('symwasm.core.S.One'); },
  /** Negative one. */
  get NegativeOne(): Expr { throw new NotImplementedError('symwasm.core.S.NegativeOne'); },
  /** One half. */
  get Half(): Expr { throw new NotImplementedError('symwasm.core.S.Half'); },
  /** Positive infinity. */
  get Infinity(): Expr { throw new NotImplementedError('symwasm.core.S.Infinity'); },
  /** Negative infinity. */
  get NegativeInfinity(): Expr { throw new NotImplementedError('symwasm.core.S.NegativeInfinity'); },
};

/** Positive infinity. */
export const oo: Expr = new SentinelExpr('oo');

/** Pi (π). */
export const pi: Expr = new SentinelExpr('pi');

/** Euler's number (e). */
export const E: Expr = new SentinelExpr('E');

/** Imaginary unit (i). */
export const I: Expr = new SentinelExpr('I');
