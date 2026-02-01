/**
 * Base expression class for all symbolic expressions.
 * @module core/expr
 */

import { getWasmModule } from '../wasm-loader.js';
import {
  SymEngineObject,
  SymEngineSet,
  SymEngineVec,
  SymEngineMap,
  checkException,
  createBasic,
} from '../wasm-memory.js';
import { SymEngineTypeID } from '../wasm-types.js';
import { parseComplexString } from './helpers.js';

// Forward declarations for lazy imports to avoid circular dependencies
let _Symbol: typeof import('./classes/symbol.js').Symbol | null = null;
let _Integer: typeof import('./numbers/integer.js').Integer | null = null;
let _exprFromWasm: ((obj: SymEngineObject) => Expr) | null = null;

function getSymbol() {
  if (!_Symbol) {
    _Symbol = require('./classes/symbol.js').Symbol;
  }
  return _Symbol;
}

function getInteger() {
  if (!_Integer) {
    _Integer = require('./numbers/integer.js').Integer;
  }
  return _Integer;
}

function getExprFromWasm() {
  if (!_exprFromWasm) {
    _exprFromWasm = require('./expr-factory.js').exprFromWasm;
  }
  return _exprFromWasm;
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
   */
  subs(old: Expr | Map<Expr, Expr> | Record<string, Expr | number>, new_?: Expr): Expr {
    if (!this._obj || !this._obj.isValid()) {
      throw new Error('Cannot substitute in expression not backed by WASM object');
    }

    const wasm = getWasmModule();
    const exprFromWasm = getExprFromWasm();

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
    const exprFromWasm = getExprFromWasm();

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
    const Symbol = getSymbol();
    const Integer = getInteger();
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
    const exprFromWasm = getExprFromWasm();

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
  free_symbols(): Expr[] {
    // No WASM backing means no symbols (e.g., sentinel constants)
    if (!this._obj || !this._obj.isValid()) {
      return [];
    }

    const wasm = getWasmModule();
    const set = new SymEngineSet();
    const Symbol = getSymbol();
    try {
      const code = wasm._basic_free_symbols(this._obj.getPtr(), set.getPtr());
      checkException(code);

      const symbols: Expr[] = [];
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
    const exprFromWasm = getExprFromWasm();
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
