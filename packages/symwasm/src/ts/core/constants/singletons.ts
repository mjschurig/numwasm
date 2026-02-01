/**
 * Singleton constants (S namespace).
 * @module core/constants/singletons
 */

import { Integer } from '../numbers/integer.js';
import { Rational } from '../numbers/rational.js';
import { Infinity_ } from '../classes/infinity.js';
import { NaN_ } from '../classes/nan.js';

// Private cache for singleton constants (lazy initialization)
let _zero: Integer | null = null;
let _one: Integer | null = null;
let _negOne: Integer | null = null;
let _half: Rational | null = null;
let _infinity: Infinity_ | null = null;
let _negInfinity: Infinity_ | null = null;
let _complexInfinity: Infinity_ | null = null;
let _nan: NaN_ | null = null;

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
