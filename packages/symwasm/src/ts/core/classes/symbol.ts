/**
 * Symbol class for symbolic variables.
 * @module core/classes/symbol
 */

import { getWasmModule } from '../../wasm-loader.js';
import { SymEngineObject, createBasic, checkException, withTempString } from '../../wasm-memory.js';
import { Expr } from '../expr.js';

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
