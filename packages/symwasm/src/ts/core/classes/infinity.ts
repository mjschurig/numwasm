/**
 * Infinity class for infinite values.
 * @module core/classes/infinity
 */

import { getWasmModule } from '../../wasm-loader.js';
import { SymEngineObject, createBasic } from '../../wasm-memory.js';
import { Expr } from '../expr.js';

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
