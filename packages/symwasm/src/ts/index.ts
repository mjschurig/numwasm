/**
 * symwasm - SymPy-inspired symbolic mathematics in TypeScript.
 *
 * A comprehensive symbolic computation library built on SymEngine (C++ library)
 * compiled to WebAssembly. Provides capabilities for symbolic algebra, calculus,
 * equation solving, matrix operations, and expression manipulation.
 *
 * @packageDocumentation
 * @module symwasm
 *
 * @remarks
 * ## Getting Started
 *
 * Before using any symwasm functions, you must initialize the WASM module:
 *
 * ```typescript
 * import { loadWasmModule, Symbol, sin, cos, diff, pi } from 'symwasm';
 *
 * async function main() {
 *   // Initialize the WASM module (required once at startup)
 *   await loadWasmModule();
 *
 *   // Now you can use symbolic math!
 *   const x = new Symbol('x');
 *
 *   // Create symbolic expressions
 *   const expr = sin(x);           // sin(x)
 *   const derivative = diff(sin(x), x);  // cos(x)
 *
 *   // Evaluate at specific values
 *   const value = sin(pi).evalf();  // 0
 * }
 * ```
 *
 * ## Core Concepts
 *
 * ### Expressions
 * All symbolic values in symwasm are expressions (`Expr`). This includes:
 * - Symbols: `new Symbol('x')` - Variables
 * - Numbers: `new Integer(5)`, `new Rational(1, 2)`, `new Float(3.14)`
 * - Constants: `pi`, `E`, `I` (imaginary unit)
 * - Composite expressions: `add(x, y)`, `mul(x, y)`, `sin(x)`
 *
 * ### Arithmetic Operations
 * Basic operations are provided as functions:
 * - `add(a, b)` - Addition
 * - `sub(a, b)` - Subtraction
 * - `mul(a, b)` - Multiplication
 * - `div(a, b)` - Division
 * - `pow(a, b)` - Power
 * - `neg(a)` - Negation
 *
 * ### Functions
 * Mathematical functions work on symbolic expressions:
 * - Trigonometric: `sin`, `cos`, `tan`, `asin`, `acos`, `atan`
 * - Hyperbolic: `sinh`, `cosh`, `tanh`, `asinh`, `acosh`, `atanh`
 * - Exponential/Log: `exp`, `log`, `sqrt`, `pow`
 * - Special: `gamma`, `erf`, `zeta`, `beta`
 *
 * ### Calculus
 * Calculus operations for symbolic differentiation and series:
 * - `diff(expr, x)` - Differentiation
 * - `series(expr, x)` - Taylor series expansion
 *
 * ### Simplification
 * Transform and simplify expressions:
 * - `expand(expr)` - Expand products
 * - `simplify(expr)` - Apply simplification rules
 *
 * ## Module Organization
 *
 * The library is organized into several modules:
 * - {@link core} - Core expression types and operations
 * - {@link functions} - Mathematical functions
 * - {@link calculus} - Calculus operations
 * - {@link simplify} - Expression simplification
 * - {@link solvers} - Equation solving (planned)
 * - {@link matrices} - Matrix operations
 * - {@link printing} - Output formatting (planned)
 *
 * @example
 * ### Symbolic Differentiation
 * ```typescript
 * import { loadWasmModule, Symbol, sin, cos, diff, mul, pow } from 'symwasm';
 *
 * await loadWasmModule();
 * const x = new Symbol('x');
 *
 * // d/dx sin(x) = cos(x)
 * diff(sin(x), x);
 *
 * // d/dx x³ = 3x²
 * diff(pow(x, 3), x);
 *
 * // d²/dx² sin(x) = -sin(x)
 * diff(sin(x), x, 2);
 * ```
 *
 * @example
 * ### Taylor Series
 * ```typescript
 * import { loadWasmModule, Symbol, sin, exp, series } from 'symwasm';
 *
 * await loadWasmModule();
 * const x = new Symbol('x');
 *
 * // sin(x) = x - x³/6 + x⁵/120 - ...
 * series(sin(x), x, 0, 6);
 *
 * // e^x = 1 + x + x²/2 + x³/6 + ...
 * series(exp(x), x, 0, 5);
 * ```
 *
 * @example
 * ### Matrix Operations
 * ```typescript
 * import { loadWasmModule, Matrix, Symbol, eye, det } from 'symwasm';
 *
 * await loadWasmModule();
 *
 * // Create a symbolic matrix
 * const x = new Symbol('x');
 * const m = new Matrix([[1, 2], [3, x]]);
 *
 * // Compute determinant
 * m.det();  // x - 6
 *
 * // Matrix operations
 * const identity = eye(3);
 * m.transpose();
 * m.inv();
 * ```
 */

// Namespace exports (for import * as core from 'symwasm')
export * as core from './core/index.js';
export * as functions from './functions/index.js';
export * as calculus from './calculus/index.js';
export * as simplification from './simplify/index.js';
export * as solvers from './solvers/index.js';
export * as matrices from './matrices/index.js';
export * as printing from './printing/index.js';

// Error types
export { NotImplementedError } from './errors.js';

// Direct exports for convenience (most commonly used)
export {
  // WASM initialization
  loadWasmModule,
  // Base class
  Expr,
  // Type ID enum
  SymEngineTypeID,
  // Expression factory
  exprFromWasm,
  // Classes
  Symbol,
  symbols,
  Add,
  Mul,
  Pow,
  Constant,
  ImaginaryUnit,
  Infinity_,
  NaN_,
  // Number types
  Integer,
  Rational,
  Float,
  Complex,
  // Arithmetic operations
  add,
  sub,
  mul,
  div,
  pow,
  neg,
  // Constants
  S,
  pi,
  E,
  I,
  oo,
  EulerGamma,
  Catalan,
  GoldenRatio,
} from './core/index.js';

// Direct function exports
export {
  // Trig
  sin,
  cos,
  tan,
  cot,
  sec,
  csc,
  asin,
  acos,
  atan,
  acot,
  asec,
  acsc,
  atan2,
  // Hyperbolic
  sinh,
  cosh,
  tanh,
  coth,
  sech,
  csch,
  asinh,
  acosh,
  atanh,
  acoth,
  asech,
  acsch,
  // Exp/Log
  exp,
  log,
  sqrt,
  cbrt,
  lambertw,
  // Other
  abs,
  sign,
  floor,
  ceiling,
  // Special
  gamma,
  loggamma,
  digamma,
  polygamma,
  beta,
  lowergamma,
  uppergamma,
  erf,
  erfc,
  zeta,
  dirichlet_eta,
  kronecker_delta,
  // Complex
  conjugate,
  re,
  im,
  arg,
  // Min/Max
  Max,
  Min,
} from './functions/index.js';

// Calculus exports
export { diff, series, integrate, limit, summation } from './calculus/index.js';

// Simplification exports
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
  cse,
  type CSEResult,
} from './simplify/index.js';

// Matrix exports
export { Matrix, SparseMatrix, eye, zeros, ones, diag, jacobian } from './matrices/index.js';
