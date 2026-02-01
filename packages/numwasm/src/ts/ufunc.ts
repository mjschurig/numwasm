/**
 * NumJS Universal Functions (Ufuncs)
 *
 * Re-exports from ufunc/ for backwards compatibility.
 */

// Arithmetic operations
export {
  negative,
  positive,
  absolute,
  abs,
  fabs,
  sign,
  add,
  subtract,
  multiply,
  divide,
  true_divide,
  floor_divide,
  remainder,
  mod,
  fmod,
  power,
  pow,
  float_power,
  sqrt,
  square,
  cbrt,
  reciprocal,
} from "./ufunc/arithmetic.js";

// Trigonometric operations
export {
  sin,
  cos,
  tan,
  arcsin,
  arccos,
  arctan,
  asin,
  acos,
  atan,
  arctan2,
  atan2,
  sinh,
  cosh,
  tanh,
  arcsinh,
  arccosh,
  arctanh,
  asinh,
  acosh,
  atanh,
  degrees,
  rad2deg,
  radians,
  deg2rad,
  hypot,
} from "./ufunc/trigonometric.js";

// Exponential and logarithmic operations
export {
  exp,
  exp2,
  expm1,
  log,
  log2,
  log10,
  log1p,
  logaddexp,
  logaddexp2,
} from "./ufunc/exponential.js";

// Comparison operations
export {
  equal,
  not_equal,
  less,
  less_equal,
  greater,
  greater_equal,
  maximum,
  minimum,
  fmax,
  fmin,
} from "./ufunc/comparison.js";

// Logical operations
export {
  logical_not,
  logical_and,
  logical_or,
  logical_xor,
} from "./ufunc/logical.js";

// Bitwise operations
export {
  invert,
  bitwise_not,
  bitwise_invert,
  bitwise_and,
  bitwise_or,
  bitwise_xor,
  left_shift,
  right_shift,
  bitwise_left_shift,
  bitwise_right_shift,
  bitwise_count,
} from "./ufunc/bitwise.js";

// Rounding operations
export { floor, ceil, trunc, rint, round, around, fix } from "./ufunc/rounding.js";

// Special functions
export {
  copysign,
  signbit,
  frexp,
  ldexp,
  nextafter,
  spacing,
  modf,
  gcd,
  lcm,
  sinc,
  heaviside,
  divmod,
  conjugate,
  conj,
} from "./ufunc/special.js";
