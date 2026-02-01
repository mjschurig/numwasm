/**
 * Rational number type.
 * @module core/numbers/rational
 */

import { getWasmModule } from '../../wasm-loader.js';
import { SymEngineObject, createBasic, checkException } from '../../wasm-memory.js';
import { Expr } from '../expr.js';

/**
 * Exact rational number p/q. Mirrors sympy.Rational.
 */
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
