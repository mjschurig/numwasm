/**
 * Constant class for mathematical constants.
 * @module core/classes/constant
 */

import { getWasmModule } from '../../wasm-loader.js';
import { SymEngineObject, createBasic } from '../../wasm-memory.js';
import { Expr } from '../expr.js';

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
