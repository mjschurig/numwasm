/**
 * Core symbolic expression types.
 * @module core
 */

// WASM initialization
import { loadWasmModule, getWasmModule } from '../wasm-loader.js';
export { loadWasmModule };
export { SymEngineTypeID } from '../wasm-types.js';

// Internal WASM access helper
export function getWasm() {
  return getWasmModule();
}

// Base class
export { Expr } from './expr.js';

// Expression factory
export { exprFromWasm } from './expr-factory.js';

// Classes
export { Symbol, symbols } from './classes/symbol.js';
export { Add } from './classes/add.js';
export { Mul } from './classes/mul.js';
export { Pow } from './classes/pow.js';
export { Constant } from './classes/constant.js';
export { ImaginaryUnit } from './classes/imaginary-unit.js';
export { Infinity_ } from './classes/infinity.js';
export { NaN_ } from './classes/nan.js';

// Number types
export { Integer } from './numbers/integer.js';
export { Rational } from './numbers/rational.js';
export { Float } from './numbers/float.js';
export { Complex } from './numbers/complex.js';

// Arithmetic operations
export { add } from './operations/add.js';
export { sub } from './operations/sub.js';
export { mul } from './operations/mul.js';
export { div } from './operations/div.js';
export { pow } from './operations/pow.js';
export { neg } from './operations/neg.js';

// Constants
export { S } from './constants/singletons.js';
export { pi } from './constants/pi.js';
export { E } from './constants/e.js';
export { I } from './constants/i.js';
export { oo } from './constants/infinity.js';
export { EulerGamma } from './constants/euler-gamma.js';
export { Catalan } from './constants/catalan.js';
export { GoldenRatio } from './constants/golden-ratio.js';

// Re-export all functions from the functions module for backwards compatibility
export * from '../functions/index.js';

// Re-export calculus functions for backwards compatibility
export { diff } from '../calculus/diff.js';
export { series } from '../calculus/series.js';

// Re-export simplification functions for backwards compatibility
export {
  expand,
  simplify,
  numer,
  denom,
  trigsimp,
  radsimp,
  powsimp,
  rewrite_as_exp,
  rewrite_as_sin,
  rewrite_as_cos,
  as_real_imag,
  expand_trig,
  expand_complex,
} from '../simplify/index.js';
