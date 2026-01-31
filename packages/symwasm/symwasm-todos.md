# symwasm Implementation Todos

Inventory of all modules/functions needed to replace stubs and fully implement the symwasm package, based on the SymEngine C++ library in `/packages/symwasm/reference/symengine`.

**Strategy**: Copy SymEngine C++ kernels as-is, compile to WebAssembly, and write thin TypeScript wrappers using the SymPy-like API.

Legend: ✅ = implemented, 🔲 = stubbed (exists but throws NotImplementedError), ⬜ = not yet created

---

## Workflow

For each function/module, follow these steps in order:

### Step 1: Identify the C++ source in SymEngine reference
- Locate the corresponding C++ header and implementation files in `/packages/symwasm/reference/symengine/symengine/`
- Each section below lists the exact C++ files to use
- **Copy the C++ code as-is** — do not reimplement in TypeScript

### Step 2: Configure WASM build system
- Add the identified C++ files to the Emscripten build configuration
- Ensure all dependencies (headers, source files) are included
- Use the C API wrappers from `cwrapper.h` and `cwrapper.cpp` for JS/WASM interop
- Configure memory management and exception handling for WASM environment

### Step 3: Compile C++ to WebAssembly
- Run Emscripten build to compile SymEngine C++ to `.wasm`
- Expose the necessary C API entry points via EMSCRIPTEN_BINDINGS or cwrapper
- Test that WASM module loads and basic functions are callable from JavaScript
- Verify memory management (allocate, free, reference counting)

### Step 4: Write TypeScript wrappers (thin layer only)
- Create TypeScript classes/functions that call into the WASM module via cwrapper
- Follow SymPy-like API conventions (Symbol, Add, diff, expand, etc.) for familiarity
- Handle type conversions between TypeScript and WASM (strings, numbers, arrays)
- Implement proper resource cleanup (free WASM memory when objects are garbage collected)
- **Do NOT reimplement algorithms** — only create bindings

### Step 5: Port tests from SymEngine
- Find corresponding tests in `/packages/symwasm/reference/symengine/symengine/tests/`
- Translate C++ test cases to TypeScript/Vitest format
- Preserve numerical tolerances and mathematical rigor
- Add TypeScript-specific tests for API ergonomics and type safety

---

## Implementation Phases

### Phase 1: Foundation (Core Types & WASM Setup)
**Goal**: Establish WASM build system and basic symbolic types
**C++ Files**: `basic.h`, `symbol.h`, `number.h`, `integer.h`, `rational.h`, `real_double.h`, `complex_double.h`, `add.h`, `mul.h`, `pow.h`, `constants.h`, `cwrapper.h`

#### 1.1 Build System & C API Wrapper ✅ COMPLETED
**SymEngine Files**: `cwrapper.h`, `cwrapper.cpp`
- ✅ Set up Emscripten build configuration (`scripts/build-wasm.sh` with GMP support)
- ✅ Expose C API wrapper functions via `cwrapper.h` (compiled `cwrapper.cpp` to WASM)
- ✅ Implement memory management bridge (SymEngine's RCP with reference counting)
- ✅ Set up exception handling for WASM (`checkException` in `wasm-memory.ts`)
- ✅ Create TypeScript type definitions for WASM module (`wasm-types.ts`, `wasm-loader.ts`)
- **Note**: CMake not used - direct Emscripten compilation more appropriate for WASM
- **Build Output**: `dist/wasm/symwasm.wasm` (563KB), `symwasm.cjs` (19KB), `symwasm.mjs` (19KB)
- **Dependencies**: GMP library compiled to WASM at `.gmp-build/`
- **Verified**: Integer arithmetic (3+5=8, 3*5=15, 3^5=243) and symbolic math (x+x=2*x, x*x=x**2) working

#### 1.2 Core Base Classes ✅ COMPLETED
**SymEngine Files**: `basic.h`, `basic-inl.h`, `dict.h`, `type_codes.h`
- ✅ `Expr` — Base class for all symbolic expressions (maps to `Basic` in SymEngine)
- ✅ `Expr.equals(other)` — Structural equality (`basic_eq`)
- ✅ `Expr.free_symbols()` — Get free symbols (`basic_free_symbols`)
- ✅ `Expr.hash()` — Hash code (`basic_hash`)
- ✅ `Expr.get_type()` — Type identification (`basic_get_type`)
- ✅ `Expr.toString()` — String representation (`basic_str`)
- ✅ `Expr.free()` — Release WASM memory
- ✅ `SymEngineSet` — Wrapper for CSetBasic container
- ✅ `Symbol._fromWasm()` — Internal factory for creating Symbol from WASM object
- **Note**: Full functionality requires WASM-backed expressions (Phase 1.3+). Sentinel constants (pi, E, I, oo) work with limited functionality.

#### 1.3 Symbols and Variables ✅ COMPLETED
**SymEngine Files**: `symbol.h`, `symbol.cpp`
- ✅ `Symbol(name, assumptions?)` — Create symbolic variable (uses `_symbol_set`)
- ✅ `symbols(names, assumptions?)` — Create multiple symbols from space/comma-separated string
- ⬜ `Dummy(name?)` — Temporary symbol (like SymPy Dummy)
- ⬜ `Wild(name)` — Wildcard for pattern matching
- **Note**: Also fixed WASM memory management to use `_basic_new_heap`/`_basic_free_heap` instead of stack allocation

#### 1.4 Number Types ✅ COMPLETED
**SymEngine Files**: `number.h`, `integer.h`, `rational.h`, `real_double.h`, `complex_double.h`, `complex.h`
- ✅ `Integer(value)` — Exact integer (maps to `Integer` class, uses `_integer_set_si`)
- ✅ `Rational(p, q)` — Exact rational p/q (maps to `Rational` class, uses `_rational_set_si`)
- ✅ `Float(value, precision?)` — Machine precision float (maps to `RealDouble`, uses `_real_double_set_d`)
- ✅ `Complex(re, im)` — Exact complex number (maps to `Complex` class, uses `_complex_set`)
- ✅ `S.Zero`, `S.One`, `S.NegativeOne`, `S.Half` — Lazy-initialized singleton constants
- ⬜ `ComplexDouble(re, im)` — Machine precision complex (no cwrapper creation function)
- ⬜ `RealMPFR(value, precision)` — Arbitrary precision real (requires MPFR library)
- ⬜ `ComplexMPC(re, im, precision)` — Arbitrary precision complex (requires MPC library)
- **Note**: SymEngine simplifies rationals automatically (e.g., 4/8 → 1/2). Float precision parameter ignored without MPFR. ComplexDouble/RealMPFR/ComplexMPC not available due to missing cwrapper support or library dependencies.

#### 1.5 Core Arithmetic Operations ✅ COMPLETED
**SymEngine Files**: `add.h`, `mul.h`, `pow.h`
- ✅ `Add` — Symbolic addition class with lazy args extraction
- ✅ `Mul` — Symbolic multiplication class with lazy args extraction
- ✅ `Pow` — Symbolic exponentiation class with base/exp properties
- ✅ `add(a, b)` — Add two expressions (`_basic_add`)
- ✅ `sub(a, b)` — Subtract two expressions (`_basic_sub`)
- ✅ `mul(a, b)` — Multiply two expressions (`_basic_mul`)
- ✅ `div(a, b)` — Divide two expressions (`_basic_div`)
- ✅ `pow(base, exp)` — Raise to power (`_basic_pow`)
- ✅ `neg(a)` — Negate expression (`_basic_neg`)
- ✅ `Expr.get_args()` — Get sub-expressions (`_basic_get_args`)
- ✅ `SymEngineVec` — Wrapper for CVecBasic container
- ✅ `exprFromWasm()` — Factory function for type-based Expr creation
- ⬜ Operator overloading support in TypeScript wrappers
- **Note**: SymEngine auto-simplifies: `x + x` → `2*x`, `x * x` → `x**2`. Results dispatched to correct Expr subclass (Integer, Add, Mul, Pow, etc.)

#### 1.6 Constants ✅ COMPLETED
**SymEngine Files**: `constants.h`, `constants.cpp`
- ✅ `pi` — Pi (π) → `Constant` class (WASM-backed via `_basic_const_pi`)
- ✅ `E` — Euler's number (e) → `Constant` class (WASM-backed via `_basic_const_E`)
- ✅ `I` — Imaginary unit (i) → `ImaginaryUnit` class (WASM-backed via `_basic_const_I`)
- ✅ `oo` — Positive infinity → `Infinity_` class (WASM-backed via `_basic_const_infinity`)
- ✅ `S.Zero`, `S.One`, `S.NegativeOne`, `S.Half` — Numeric constants (implemented in Phase 1.4)
- ✅ `S.Infinity`, `S.NegativeInfinity`, `S.ComplexInfinity` — Infinity variants (WASM-backed)
- ✅ `S.NaN` — Not a number → `NaN_` class (WASM-backed via `_basic_const_nan`)
- ✅ `EulerGamma` — Euler-Mascheroni constant γ → `Constant` class (WASM-backed)
- ✅ `Catalan` — Catalan's constant → `Constant` class (WASM-backed)
- ✅ `GoldenRatio` — Golden ratio φ → `Constant` class (WASM-backed)
- **Note**: All constants use lazy initialization via Proxy to defer WASM calls until first use. Constants support full arithmetic operations, hashing, and type identification.

#### 1.7 Basic Expression Manipulation ✅ COMPLETED
**SymEngine Files**: `subs.h`, `subs.cpp`
- ✅ `Expr.subs(old, new)` — Single substitution (WASM-backed via `_basic_subs2`)
- ✅ `Expr.subs(Map<Expr, Expr>)` — Multiple simultaneous substitutions (WASM-backed via `_basic_subs` with `CMapBasicBasic`)
- ✅ `Expr.subs({ x: 1, y: 2 })` — Object notation for symbol substitution (convenience wrapper)
- ⬜ `Expr.xreplace(dict)` — Exact structural replacement (not exposed in C API, cannot implement)
- **Note**: Map-based substitution is atomic (all substitutions happen simultaneously), unlike chained single substitutions which apply sequentially. Object notation creates symbols by name and converts numbers to Integer automatically.

#### 1.8 Numerical Evaluation ✅ COMPLETED
**SymEngine Files**: `eval.h`, `eval_double.h`, `eval_mpfr.h`
- ✅ `Expr.evalf(precision?)` — Numerical evaluation (WASM-backed via `_basic_evalf`)
- ✅ `Expr.evalfNumber()` — Extract JavaScript number from evaluated expression
- ✅ `Expr.evalfComplex()` — Extract complex number { real, imag }
- ⬜ `N(expr, n)` — Numerical evaluation (alias) - deferred
- **Note**: Default precision is 53 bits (double precision). Higher precision requires MPFR which is not compiled into current WASM build. Complex results parsed from string representation.

---

**Phase 1 Complete**: All core foundation features implemented (197 tests passing).

---

### Phase 2: Essential Functions (Calculus & Simplification)
**Goal**: Enable symbolic calculus and expression manipulation
**C++ Files**: `functions.h`, `derivative.h`, `series.h`, `expand.h`, `subs.h`

**Phase 2.1 Complete**: 52 elementary functions implemented (330 tests passing).
- Phase 2.1a: 45 core functions (trig, hyperbolic, exp/log, special functions)
- Phase 2.1b: 7 additional functions (digamma, conjugate, re, im, arg, Max, Min)

**Phase 2.2 Complete**: Symbolic differentiation implemented (362 tests passing).
- `diff(expr, x)` — First derivative
- `diff(expr, x, n)` — nth derivative via iteration
- `diff(expr, x, y)` — Multi-variable partial derivatives

**Phase 2.3 Complete**: Taylor series expansion implemented (382 tests passing).
- `series(expr, x)` — Taylor series around x=0
- `series(expr, x, 0, n)` — Series with configurable number of terms

#### 2.1 Elementary Functions ✅ COMPLETED
**SymEngine Files**: `functions.h`, `functions.cpp`

##### Exponential & Logarithmic ✅
- ✅ `exp(x)` — Exponential e^x (WASM-backed via `_basic_exp`)
- ✅ `log(x)` — Natural log (WASM-backed via `_basic_log`)
- ✅ `sqrt(x)` — Square root (WASM-backed via `_basic_sqrt`)
- ✅ `cbrt(x)` — Cube root (WASM-backed via `_basic_cbrt`)
- ✅ `abs(x)` — Absolute value (WASM-backed via `_basic_abs`)
- ✅ `lambertw(x)` — Lambert W function (WASM-backed via `_basic_lambertw`)

##### Trigonometric Functions ✅
- ✅ `sin(x)` — Sine (WASM-backed via `_basic_sin`)
- ✅ `cos(x)` — Cosine (WASM-backed via `_basic_cos`)
- ✅ `tan(x)` — Tangent (WASM-backed via `_basic_tan`)
- ✅ `cot(x)` — Cotangent (WASM-backed via `_basic_cot`)
- ✅ `sec(x)` — Secant (WASM-backed via `_basic_sec`)
- ✅ `csc(x)` — Cosecant (WASM-backed via `_basic_csc`)

##### Inverse Trigonometric ✅
- ✅ `asin(x)` — Arcsine (WASM-backed via `_basic_asin`)
- ✅ `acos(x)` — Arccosine (WASM-backed via `_basic_acos`)
- ✅ `atan(x)` — Arctangent (WASM-backed via `_basic_atan`)
- ✅ `acot(x)` — Arccotangent (WASM-backed via `_basic_acot`)
- ✅ `asec(x)` — Arcsecant (WASM-backed via `_basic_asec`)
- ✅ `acsc(x)` — Arccosecant (WASM-backed via `_basic_acsc`)
- ✅ `atan2(y, x)` — Two-argument arctangent (WASM-backed via `_basic_atan2`)

##### Hyperbolic Functions ✅
- ✅ `sinh(x)` — Hyperbolic sine (WASM-backed via `_basic_sinh`)
- ✅ `cosh(x)` — Hyperbolic cosine (WASM-backed via `_basic_cosh`)
- ✅ `tanh(x)` — Hyperbolic tangent (WASM-backed via `_basic_tanh`)
- ✅ `coth(x)` — Hyperbolic cotangent (WASM-backed via `_basic_coth`)
- ✅ `sech(x)` — Hyperbolic secant (WASM-backed via `_basic_sech`)
- ✅ `csch(x)` — Hyperbolic cosecant (WASM-backed via `_basic_csch`)

##### Inverse Hyperbolic ✅
- ✅ `asinh(x)` — Inverse hyperbolic sine (WASM-backed via `_basic_asinh`)
- ✅ `acosh(x)` — Inverse hyperbolic cosine (WASM-backed via `_basic_acosh`)
- ✅ `atanh(x)` — Inverse hyperbolic tangent (WASM-backed via `_basic_atanh`)
- ✅ `acoth(x)` — Inverse hyperbolic cotangent (WASM-backed via `_basic_acoth`)
- ✅ `asech(x)` — Inverse hyperbolic secant (WASM-backed via `_basic_asech`)
- ✅ `acsch(x)` — Inverse hyperbolic cosecant (WASM-backed via `_basic_acsch`)

##### Special Functions ✅
**SymEngine Files**: `functions.h` (Gamma, Beta, Erf, etc. classes)
- ✅ `gamma(x)` — Gamma function (WASM-backed via `_basic_gamma`)
- ✅ `loggamma(x)` — Log-gamma (WASM-backed via `_basic_loggamma`)
- ✅ `polygamma(n, x)` — Polygamma (WASM-backed via `_basic_polygamma`)
- ✅ `beta(x, y)` — Beta function (WASM-backed via `_basic_beta`)
- ✅ `lowergamma(s, x)` — Lower incomplete gamma (WASM-backed via `_basic_lowergamma`)
- ✅ `uppergamma(s, x)` — Upper incomplete gamma (WASM-backed via `_basic_uppergamma`)
- ✅ `erf(x)` — Error function (WASM-backed via `_basic_erf`)
- ✅ `erfc(x)` — Complementary error function (WASM-backed via `_basic_erfc`)
- ✅ `zeta(s)` — Riemann zeta (WASM-backed via `_basic_zeta`)
- ✅ `dirichlet_eta(s)` — Dirichlet eta (WASM-backed via `_basic_dirichlet_eta`)
- ✅ `kronecker_delta(i, j)` — Kronecker delta (WASM-backed via `_basic_kronecker_delta`)
- ✅ `floor(x)` — Floor function (WASM-backed via `_basic_floor`)
- ✅ `ceiling(x)` — Ceiling function (WASM-backed via `_basic_ceiling`)
- ✅ `sign(x)` — Sign function (WASM-backed via `_basic_sign`)
- ✅ `digamma(x)` — Digamma ψ(x) (WASM-backed via `_basic_digamma`, added C wrapper)
- ⬜ `LeviCivita(*indices)` — Levi-Civita symbol (not exposed in C API)
- ✅ `Max(...args)` — Maximum (WASM-backed via `_basic_max` with CVecBasic)
- ✅ `Min(...args)` — Minimum (WASM-backed via `_basic_min` with CVecBasic)

##### Complex Number Functions ✅
- ✅ `conjugate(x)` — Complex conjugate (WASM-backed via `_basic_conjugate`, added C wrapper)
- ✅ `re(x)` — Real part (WASM-backed via `_complex_base_real_part`)
- ✅ `im(x)` — Imaginary part (WASM-backed via `_complex_base_imaginary_part`)
- ✅ `arg(x)` — Argument (phase) (derived via `atan2(im(x), re(x))`)

**Note**: 52 functions implemented (45 original + 7 Phase 2.1b). SymEngine auto-simplifies (e.g., `sin(0)` → `0`, `sqrt(4)` → `2`). All functions support symbolic inputs and numerical evaluation via `evalf()`. Note: `re()` and `im()` only work on ComplexBase types (not integers that SymEngine simplifies from Complex(x, 0)).

#### 2.2 Calculus — Differentiation ✅ COMPLETED
**SymEngine Files**: `derivative.h`, `derivative.cpp`
- ✅ `diff(expr, symbol, n?)` — Differentiate (WASM-backed via `_basic_diff`, nth via iteration)
- ✅ `diff(expr, symbol1, symbol2, ...)` — Multi-variable partial derivatives (chained `_basic_diff` calls)
- ✅ `diff(expr, x, 2, y, 3)` — Mixed higher-order partial derivatives
- ⬜ `Derivative(expr, *symbols)` — Unevaluated derivative class (C API doesn't expose)
- ⬜ `fdiff(expr, argindex)` — Derivative w.r.t. function argument (C API doesn't expose)

**Note**: Core differentiation implemented with nth derivative and multi-variable support. Chain rule, product rule, quotient rule all work automatically. 33 tests covering polynomials, trig, exp/log, higher-order, and mixed partial derivatives.

#### 2.3 Calculus — Series Expansion ✅ COMPLETED
**SymEngine Files**: `series.h`, `series_generic.h`, `series_visitor.h`
- ✅ `series(expr, symbol, x0?, n?)` — Taylor series expansion (WASM-backed via `_basic_series`)
- ⬜ `series(expr, x, x0, n, dir)` — Series with direction (C++ API only supports x=0)
- ⬜ `Order(expr)` — Order term O(x^n) (not exposed in C API, result is polynomial)

**Note**: Series expansion around x=0 implemented with configurable number of terms. Supports exp, sin, cos, tan, log, and composed functions. 20 tests covering basic expansions, trig, exp/log, and error handling. Expansion around non-zero points not supported by underlying C++ API.

#### 2.4 Simplification ✅ COMPLETED
**SymEngine Files**: `expand.h`, `simplify.h`, `rewrite.h`, `as_real_imag.h`

##### Implemented Functions
- ✅ `expand(expr)` — Expand expressions (WASM-backed via `_basic_expand`)
- ✅ `simplify(expr)` — Simplify using heuristics (WASM-backed via `_basic_simplify`)
- ✅ `trigsimp(expr)` — Simplify trigonometric expressions (delegates to simplify)
- ✅ `radsimp(expr)` — Simplify radicals (delegates to simplify)
- ✅ `powsimp(expr)` — Simplify powers (delegates to simplify)
- ✅ `numer(expr)` — Extract numerator (WASM-backed via `_basic_as_numer_denom`)
- ✅ `denom(expr)` — Extract denominator (WASM-backed via `_basic_as_numer_denom`)
- ✅ `rewrite_as_exp(expr)` — Rewrite trig as exponentials (WASM-backed via `_basic_rewrite_as_exp`)
- ✅ `rewrite_as_sin(expr)` — Rewrite trig in terms of sine (WASM-backed via `_basic_rewrite_as_sin`)
- ✅ `rewrite_as_cos(expr)` — Rewrite trig in terms of cosine (WASM-backed via `_basic_rewrite_as_cos`)
- ✅ `as_real_imag(expr)` — Extract real and imaginary parts (WASM-backed via `_basic_as_real_imag`)
- ✅ `expand_trig(expr)` — Expand trigonometric (via rewrite_as_exp + expand)
- ✅ `expand_complex(expr)` — Expand re + i*im (alias for as_real_imag)

##### Removed (no SymEngine support)
- ~~`factor(expr)`~~ — Removed (complex template function, no C API)
- ~~`collect(expr, syms)`~~ — Removed (no C API)
- ~~`cancel(expr)`~~ — Removed (complex template function, no C API)
- ~~`expand_mul(expr)`~~ — Use general expand()
- ~~`expand_log(expr)`~~ — Not in SymEngine
- ~~`expand_power_base(expr)`~~ — Not in SymEngine
- ~~`expand_power_exp(expr)`~~ — Not in SymEngine
- ~~`together(expr)`~~ — Not in SymEngine
- ~~`apart(expr, x)`~~ — Not in SymEngine

**Note**: Full simplification implemented with 13 functions. The simplify function handles csc^(-1)→sin, sec^(-1)→cos, cot^(-1)→tan transformations. Rewrite functions allow converting between trig representations. as_real_imag extracts complex number components. 50+ tests covering polynomial expansion, trig simplification, rewrite functions, and complex number decomposition.

---

### Phase 3: Advanced Mathematics (Matrices, Polynomials, Solvers)
**Goal**: Linear algebra, polynomial operations, equation solving
**C++ Files**: `matrix.h`, `polys/`, `solve.h`

#### 3.1 Matrix Operations ✅ PARTIALLY COMPLETED
**SymEngine Files**: `matrix.h` (C API via `cwrapper.h`)

##### Dense Matrix — Construction ✅ COMPLETED
- ✅ `Matrix(data)` — Create from nested array (`dense_matrix_new_vec`)
- ✅ `Matrix.fromFlat(flat, rows, cols)` — Create from flat array (`dense_matrix_new_vec`)
- ✅ `eye(n, m?, k?)` — Identity matrix (`dense_matrix_eye`)
- ✅ `zeros(rows, cols)` — Zero matrix (`dense_matrix_zeros`)
- ✅ `ones(rows, cols)` — Ones matrix (`dense_matrix_ones`)
- ✅ `diag(values, k?)` — Diagonal matrix (`dense_matrix_diag`)

##### Dense Matrix — Properties ✅ COMPLETED
- ✅ `Matrix.get(i, j)` — Get element (`dense_matrix_get_basic`)
- ✅ `Matrix.set(i, j, val)` — Set element (`dense_matrix_set_basic`)
- ✅ `Matrix.rows` — Number of rows (`dense_matrix_rows`)
- ✅ `Matrix.cols` — Number of columns (`dense_matrix_cols`)
- ✅ `Matrix.shape` — Tuple (rows, cols)
- ✅ `Matrix.toString()` — String representation (`dense_matrix_str`)
- ✅ `Matrix.equals(other)` — Equality test (`dense_matrix_eq`)
- ✅ `Matrix.free()` — Free WASM memory

**Note**: 35 tests passing. Matrix construction, factory functions, properties, and element access all implemented.

##### Dense Matrix — Basic Operations (Stubs)
- 🔲 `Matrix.det()` — Determinant (`dense_matrix_det`)
- 🔲 `Matrix.inv()` — Inverse (`dense_matrix_inv`)
- 🔲 `Matrix.transpose()` — Transpose (`dense_matrix_transpose`)
- ⬜ `Matrix.add(other)` — Matrix addition (`dense_matrix_add_matrix`)
- ⬜ `Matrix.mul(other)` — Matrix multiplication (`dense_matrix_mul_matrix`)
- ⬜ `Matrix.addScalar(k)` — Add scalar (`dense_matrix_add_scalar`)
- ⬜ `Matrix.mulScalar(k)` — Multiply by scalar (`dense_matrix_mul_scalar`)

##### Dense Matrix — Submatrix Operations
- ✅ `Matrix.submatrix(r1, c1, r2, c2)` — Extract submatrix (`dense_matrix_submatrix`)
- ✅ `Matrix.rowJoin(other)` — Horizontal stack (`dense_matrix_row_join`)
- ✅ `Matrix.colJoin(other)` — Vertical stack (`dense_matrix_col_join`)
- ✅ `Matrix.rowDel(k)` — Delete row (`dense_matrix_row_del`)
- ✅ `Matrix.colDel(k)` — Delete column (`dense_matrix_col_del`)

##### Dense Matrix — Factorizations
- ✅ `Matrix.lu()` — LU decomposition (`dense_matrix_LU`)
- ✅ `Matrix.ldl()` — LDL decomposition (`dense_matrix_LDL`)
- ✅ `Matrix.fflu()` — Fraction-free LU (`dense_matrix_FFLU`)
- ✅ `Matrix.ffldu()` — Fraction-free LDU (`dense_matrix_FFLDU`)
- ✅ `Matrix.luSolve(b)` — Solve Ax=b using LU (`dense_matrix_LU_solve`)

##### Dense Matrix — Calculus
- ⬜ `Matrix.diff(x)` — Differentiate elements (`dense_matrix_diff`)
- ⬜ `jacobian(A, x)` — Jacobian matrix (`dense_matrix_jacobian`)

##### Sparse Matrix (CSR format)
- ⬜ `SparseMatrix(rows, cols)` — Create sparse matrix (`sparse_matrix_new`)
- ⬜ `SparseMatrix.get(i, j)` — Get element (`sparse_matrix_get_basic`)
- ⬜ `SparseMatrix.set(i, j, val)` — Set element (`sparse_matrix_set_basic`)
- ⬜ `SparseMatrix.toString()` — String representation (`sparse_matrix_str`)
- ⬜ `SparseMatrix.equals(other)` — Equality test (`sparse_matrix_eq`)

##### NOT Available in C API (C++ only)
- ~~`eigenvals()`~~ — `eigen_values()` not in cwrapper
- ~~`eigenvects()`~~ — Not exposed
- ~~`rref()`~~ — `reduced_row_echelon_form()` not in cwrapper
- ~~`trace()`~~ — Method exists but not in cwrapper
- ~~`rank()`~~ — Not exposed
- ~~`QR()`~~ — Not exposed
- ~~`cholesky()`~~ — Not exposed
- ~~`conjugate()`~~ — Matrix conjugate not exposed
- ~~`char_poly()`~~ — Characteristic polynomial not exposed
- ~~Symbolic matrices~~ — `MatrixSymbol`, `Identity`, `ZeroMatrix`, `DiagonalMatrix`, `MatrixAdd`, `MatrixMul`, `HadamardProduct`, `Trace` classes not in C API

#### 3.2 Polynomial Operations
**SymEngine Files**: `polys/` subdirectory (11 headers)

##### Polynomial Classes
**SymEngine Files**: `polys/uexprpoly.h`, `polys/uintpoly.h`, `polys/uratpoly.h`
- ⬜ `Poly(expr, *gens)` — Polynomial class → `UExprPoly`, `UIntPoly`, `URatPoly`
- ⬜ Univariate polynomial backends: `UIntPoly`, `URatPoly`, `UExprPoly`
- ⬜ Multivariate polynomial backends: `MIntPoly`, `MExprPoly`
- ⬜ FLINT backend: `UIntPolyFlint`, `URatPolyFlint`
- ⬜ Piranha backend: `UIntPolyPiranha`, `URatPolyPiranha`

##### Polynomial Operations — Basic
**SymEngine Files**: Polynomial class methods
- ⬜ `degree(poly, gen?)` — Degree
- ⬜ `LC(poly)` — Leading coefficient
- ⬜ `coeffs(poly)` — List coefficients
- ⬜ `eval(poly, x, a)` — Evaluate at point

##### Polynomial Operations — Arithmetic
**SymEngine Files**: Polynomial class methods
- ⬜ `div(f, g)` — Division
- ⬜ `quo(f, g)` — Quotient
- ⬜ `rem(f, g)` — Remainder

##### Polynomial Operations — GCD/LCM
**SymEngine Files**: Polynomial GCD methods
- ⬜ `gcd(f, g)` — Greatest common divisor
- ⬜ `lcm(f, g)` — Least common multiple
- ⬜ `gcdex(f, g)` — Extended GCD
- ⬜ `resultant(f, g, x)` — Resultant
- ⬜ `discriminant(f, x)` — Discriminant

##### Polynomial Operations — Factorization
**SymEngine Files**: Polynomial factorization methods, `polys/cancel.h`
- 🔲 `factor(poly)` — Factor polynomial
- ⬜ `factor_list(poly)` — List of (factor, multiplicity)
- ⬜ `sqf(poly)` — Square-free factorization
- ⬜ `sqf_list(poly)` — Square-free factors list

##### Polynomial Operations — Roots
**SymEngine Files**: Polynomial root-finding methods
- ⬜ `roots(poly)` — Find all roots
- ⬜ `nroots(poly, n?)` — Numerical roots

#### 3.3 Equation Solving
**SymEngine Files**: `solve.h`, `solve.cpp`

##### Current Stubs
- 🔲 `solve(expr, symbols?)` — Solve equations → `solve()` function
- 🔲 `solveset(expr, symbol, domain?)` — Solve returning set
- 🔲 `linsolve(system, symbols)` — Linear system → `linsolve()` function
- 🔲 `nonlinsolve(system, symbols)` — Nonlinear system
- 🔲 `dsolve(eq, func?)` — ODE solver

##### Priority Additions
**SymEngine Files**: `solve.h`
- ⬜ `solve_poly(poly, x)` — Solve polynomial
- ⬜ `solve_rational(expr, x)` — Solve rational equation
- ⬜ `solve_trig(expr, x)` — Solve trigonometric equation
- ⬜ `solve_linear(eq, x)` — Solve single linear equation
- ⬜ `vecbasic_linsolve(eqs, syms)` — Matrix-based linear solver

---

### Phase 4: Specialized Mathematics (Number Theory, Sets, Logic)
**Goal**: Number theory, set operations, boolean logic
**C++ Files**: `ntheory.h`, `sets.h`, `logic.h`

#### 4.1 Number Theory
**SymEngine Files**: `ntheory.h`, `ntheory_funcs.h`, `ntheory.cpp`

##### Prime Numbers
- ⬜ `isprime(n)` — Primality test → `probab_prime_p()`
- ⬜ `nextprime(n)` — Next prime → `nextprime()`
- ⬜ `primepi(n)` — Prime counting function
- ⬜ `primorial(n)` — Product of first n primes → `primorial()`

##### Divisors
- ⬜ `divisors(n)` — List all divisors
- ⬜ `divisor_count(n)` — Count divisors
- ⬜ `totient(n)` — Euler's totient φ(n)

##### GCD and LCM
- ⬜ `gcd(*args)` — GCD → `gcd()` function
- ⬜ `lcm(*args)` — LCM → `lcm()` function
- ⬜ `gcdex(a, b)` — Extended GCD → `gcd_ext()`

##### Modular Arithmetic
- ⬜ `mod(a, m)` — Modulo → `mod()`
- ⬜ `mod_inverse(a, m)` — Modular inverse → `mod_inverse()`
- ⬜ `crt(m, v)` — Chinese Remainder Theorem → `crt()`

##### Sequences
- ⬜ `factorial(n)` — Factorial → `factorial()`
- ⬜ `binomial(n, k)` — Binomial coefficient → `binomial()`
- ⬜ `fibonacci(n)` — nth Fibonacci → `fibonacci()`
- ⬜ `lucas(n)` — nth Lucas number → `lucas()`

##### Diophantine Equations
**SymEngine Files**: `diophantine.h`, `diophantine.cpp`
- ⬜ `diophantine(eq)` — Solve Diophantine equations

#### 4.2 Sets and Intervals
**SymEngine Files**: `sets.h`, `sets.cpp`

##### Set Types
- ⬜ `FiniteSet(...elements)` — Finite set → `FiniteSet` class
- ⬜ `Interval(a, b, left_open?, right_open?)` — Real interval → `Interval` class
- ⬜ `Union(*sets)` — Union → `Union` class
- ⬜ `Intersection(*sets)` — Intersection
- ⬜ `Complement(A, B)` — Complement → `Complement` class
- ⬜ `ImageSet(lambda, base_set)` — Image set → `ImageSet` class
- ⬜ `ConditionSet(symbol, condition, base_set)` — Conditional set → `ConditionSet` class

##### Special Sets
- ⬜ `EmptySet` — Empty set → `EmptySet` class
- ⬜ `UniversalSet` — Universal set
- ⬜ `Naturals` — Natural numbers ℕ
- ⬜ `Naturals0` — ℕ ∪ {0}
- ⬜ `Integers` — Integers ℤ → `Integers` class
- ⬜ `Rationals` — Rationals ℚ → `Rationals` class
- ⬜ `Reals` — Reals ℝ → `Reals` class
- ⬜ `Complexes` — Complex numbers ℂ → `Complexes` class

##### Set Operations
- ⬜ `set_union(sets)` — Union → `set_union()`
- ⬜ `set_intersection(sets)` — Intersection → `set_intersection()`
- ⬜ `set_complement(A, B)` — Complement → `set_complement()`
- ⬜ `contains(set, elem)` — Membership test → `contains()`
- ⬜ `Set.boundary` — Boundary
- ⬜ `Set.interior` — Interior
- ⬜ `Set.closure` — Closure

#### 4.3 Boolean Logic
**SymEngine Files**: `logic.h`, `logic.cpp`

##### Logical Operators
- ⬜ `And(*args)` — Logical AND → `And` class
- ⬜ `Or(*args)` — Logical OR → `Or` class
- ⬜ `Not(expr)` — Logical NOT → `Not` class
- ⬜ `Xor(*args)` — Logical XOR → `Xor` class
- ⬜ `Implies(p, q)` — Implication
- ⬜ `true` — Boolean true → `BooleanAtom` (true)
- ⬜ `false` — Boolean false → `BooleanAtom` (false)

##### Relational Operators
**SymEngine Files**: `logic.h` (relational classes)
- ⬜ `Eq(a, b)` — Equality → `Equality` class
- ⬜ `Ne(a, b)` — Inequality → `Unequality` class
- ⬜ `Lt(a, b)` — Less than → `LessThan` class
- ⬜ `Le(a, b)` — Less than or equal → `StrictLessThan` class
- ⬜ `Gt(a, b)` — Greater than
- ⬜ `Ge(a, b)` — Greater than or equal

##### Piecewise Functions
**SymEngine Files**: `functions.h` (`Piecewise` class)
- ⬜ `Piecewise(...args)` — Piecewise function → `Piecewise` class

---

### Phase 5: I/O & Developer Tools (Printing, Parsing, Code Generation)
**Goal**: String representation, parsing, code generation
**C++ Files**: `printers/`, `parser/`, `lambda_double.h`, `llvm_double.h`

#### 5.1 Printing & String Representation
**SymEngine Files**: `printers/` subdirectory

##### Current Stubs
**Files**: `printers/strprinter.h`, `printers/latex.h`, `printers/mathml.h`, `printers/unicode.h`
- 🔲 `latex(expr, options?)` — LaTeX → `latex()` in `printers/latex.h`
- 🔲 `mathml(expr, printer?)` — MathML → `mathml()` in `printers/mathml.h`
- 🔲 `pretty(expr, options?)` — Unicode pretty-print → `unicode()` in `printers/unicode.h`
- 🔲 `sstr(expr)` — Simple string → `str()` in `printers/strprinter.h`

##### Priority Additions
**Files**: `printers/strprinter.h`, `printers/codegen.h`
- ⬜ `srepr(expr)` — Detailed representation
- ⬜ `tree(expr)` — Tree structure representation

#### 5.2 Code Generation
**SymEngine Files**: `printers/codegen.h`, `printers/codegen.cpp`

##### C/C++ Code Generation
- ⬜ `ccode(expr, assign_to?)` — Generate C code → `ccode()` function
- ⬜ `cxxcode(expr, assign_to?)` — Generate C++ code

##### JavaScript Code Generation
- ⬜ `jscode(expr, assign_to?)` — Generate JavaScript → `jscode()` function

##### Other Languages
- ⬜ `octave_code(expr)` — Octave/MATLAB code
- ⬜ `rust_code(expr)` — Rust code (if available)

#### 5.3 Parsing
**SymEngine Files**: `parser/parser.h`, `parser/tokenizer.h`

##### String to Expression
- ⬜ `parse_expr(s, transformations?)` — Parse string → `parse()` function in `parser.h`
- ⬜ `sympify(s)` — Convert to symbolic expression

##### Parsing Options
- ⬜ `parse_expr(s, {evaluate: false})` — Parse without evaluation
- ⬜ `parse_expr(s, {local_dict: {...}})` — Custom symbols

#### 5.4 Lambda & Numerical Compilation
**SymEngine Files**: `lambda_double.h`, `llvm_double.h`, `eval_double.h`

##### Lambdify — Convert to Callable Functions
- ⬜ `lambdify(args, expr)` — Convert to JS function → `LambdaRealDoubleVisitor`
- ⬜ LLVM compilation (optional): `CLLVMDoubleVisitor`, `CLLVMFloatVisitor`

##### Common Subexpression Elimination
**SymEngine Files**: CSE utilities
- ⬜ `cse(exprs, symbols?)` — CSE optimization → `basic_cse()`

---

## SymEngine C++ File Reference

### Core Files (Phase 1)
```
symengine/basic.h              → Expr base class
symengine/symbol.h             → Symbol, Dummy
symengine/number.h             → Number hierarchy
symengine/integer.h            → Integer class
symengine/rational.h           → Rational class
symengine/real_double.h        → RealDouble class
symengine/complex_double.h     → ComplexDouble class
symengine/complex.h            → Complex class
symengine/add.h                → Add class
symengine/mul.h                → Mul class
symengine/pow.h                → Pow class
symengine/constants.h          → Pi, E, I, EulerGamma, etc.
symengine/cwrapper.h           → C API for WASM bindings
symengine/subs.h               → Substitution
symengine/eval_double.h        → Numerical evaluation
```

### Functions & Calculus (Phase 2)
```
symengine/functions.h          → All elementary & special functions
symengine/derivative.h         → Differentiation
symengine/series.h             → Series expansion
symengine/series_visitor.h     → Series algorithms
symengine/expand.h             → Expression expansion
```

### Linear Algebra (Phase 3)
```
symengine/matrix.h                           → Dense matrix operations
symengine/matrices/matrix_symbol.h           → Symbolic matrices
symengine/matrices/identity_matrix.h         → Identity matrix
symengine/matrices/zero_matrix.h             → Zero matrix
symengine/matrices/diagonal_matrix.h         → Diagonal matrix
symengine/matrices/matrix_add.h              → Matrix addition
symengine/matrices/matrix_mul.h              → Matrix multiplication
symengine/matrices/hadamard_product.h        → Element-wise product
symengine/matrices/trace.h                   → Trace
symengine/matrices/transpose.h               → Transpose
symengine/matrices/conjugate_matrix.h        → Conjugate
```

### Polynomials (Phase 3)
```
symengine/polys/uintpoly.h              → Univariate integer polynomial
symengine/polys/uratpoly.h              → Univariate rational polynomial
symengine/polys/uexprpoly.h             → Univariate expression polynomial
symengine/polys/msymenginepoly.h        → Multivariate polynomial
symengine/polys/uintpoly_flint.h        → FLINT backend (optional)
symengine/polys/uintpoly_piranha.h      → Piranha backend (optional)
symengine/polys/cancel.h                → Polynomial cancellation
symengine/polys/basic_conversions.h     → Conversions
```

### Solving (Phase 3)
```
symengine/solve.h              → Equation solving
```

### Number Theory (Phase 4)
```
symengine/ntheory.h            → Number theory functions
symengine/ntheory_funcs.h      → Prime, GCD, LCM, modular arithmetic
symengine/diophantine.h        → Diophantine equations
```

### Sets & Logic (Phase 4)
```
symengine/sets.h               → Set theory
symengine/logic.h              → Boolean logic & relations
```

### I/O (Phase 5)
```
symengine/printers/strprinter.h    → String printer
symengine/printers/latex.h         → LaTeX printer
symengine/printers/mathml.h        → MathML printer
symengine/printers/unicode.h       → Unicode printer
symengine/printers/codegen.h       → Code generation (C, JS, etc.)
symengine/parser/parser.h          → Expression parser
symengine/parser/tokenizer.h       → Tokenizer
symengine/lambda_double.h          → Lambda compilation
symengine/llvm_double.h            → LLVM JIT compilation (optional)
```

---

## Summary

| Phase | Module | ✅ Done | 🔲 Stubbed | ⬜ To Create | Total | SymEngine Files |
|-------|--------|---------|-----------|-------------|-------|-----------------|
| **1** | **Foundation** | | | | | |
| 1.1 | Build System | ~6 | 0 | 0 | ~6 | cwrapper.h/cpp ✅ |
| 1.2 | Core Base | ~9 | 0 | 0 | ~9 | basic.h ✅ |
| 1.3 | Symbols | ~2 | 0 | ~2 | ~4 | symbol.h ✅ |
| 1.4 | Numbers | ~5 | 0 | ~3 | ~8 | number.h, integer.h, rational.h, complex*.h ✅ |
| 1.5 | Arithmetic | ~12 | 0 | ~1 | ~13 | add.h, mul.h, pow.h ✅ |
| 1.6 | Constants | ~11 | 0 | 0 | ~11 | constants.h ✅ |
| 1.7 | Substitution | ~3 | 0 | ~1 | ~4 | subs.h ✅ |
| 1.8 | Evaluation | ~3 | 0 | ~1 | ~4 | eval*.h ✅ |
| **2** | **Essential Functions** | | | | | |
| 2.1 | Functions | 52 | 0 | ~2 | ~54 | functions.h ✅ |
| 2.2 | Differentiation | 1 | 0 | ~3 | ~4 | derivative.h ✅ |
| 2.3 | Series | 0 | 1 | ~2 | ~3 | series.h |
| 2 | Simplification | 0 | 7 | ~10 | ~17 | expand.h, subs.h |
| **3** | **Advanced Math** | | | | | |
| 3 | Matrices | 0 | 12 | ~30 | ~42 | matrix.h, matrices/ |
| 3 | Polynomials | 0 | 1 | ~25 | ~26 | polys/ |
| 3 | Solvers | 0 | 5 | ~5 | ~10 | solve.h |
| **4** | **Specialized** | | | | | |
| 4 | Number Theory | 0 | 0 | ~20 | ~20 | ntheory*.h, diophantine.h |
| 4 | Sets | 0 | 0 | ~20 | ~20 | sets.h |
| 4 | Logic | 0 | 0 | ~12 | ~12 | logic.h |
| **5** | **I/O & Tools** | | | | | |
| 5 | Printing | 0 | 4 | ~2 | ~6 | printers/ |
| 5 | Codegen | 0 | 0 | ~5 | ~5 | printers/codegen.h |
| 5 | Parsing | 0 | 0 | ~3 | ~3 | parser/ |
| 5 | Lambda/CSE | 0 | 0 | ~3 | ~3 | lambda_double.h, llvm_double.h |
| **Total** | | **~98** | **29** | **~145** | **~272** | **65 main headers + subdirs** |

---

## Implementation Strategy

### Phase 1 Priority (Weeks 1-3)
1. Set up Emscripten build system for SymEngine
2. Compile basic types (Symbol, Integer, Rational, Add, Mul, Pow) to WASM
3. Expose cwrapper C API to JavaScript
4. Create TypeScript wrappers for core types
5. Implement basic substitution and evaluation

### Phase 2 Priority (Weeks 4-6)
1. Compile all elementary functions (trig, exp, log, etc.)
2. Add differentiation support
3. Implement series expansion
4. Add expression simplification/expansion

### Phase 3 Priority (Weeks 7-10)
1. Matrix operations (dense matrices first)
2. Polynomial arithmetic and factorization
3. Equation solving (linear systems, polynomial roots)

### Phase 4 Priority (Weeks 11-12)
1. Number theory functions
2. Set theory and logic

### Phase 5 Priority (Weeks 13-14)
1. Printers (LaTeX, MathML, etc.)
2. Parser for string input
3. Code generation
4. Lambda compilation

---

## Testing Strategy

### Test Files Location
SymEngine tests: `/packages/symwasm/reference/symengine/symengine/tests/`

### Key Test Files to Port
```
symengine/tests/basic/test_basic.cpp        → Core expression tests
symengine/tests/basic/test_number.cpp       → Number type tests
symengine/tests/basic/test_functions.cpp    → Function tests
symengine/tests/basic/test_series.cpp       → Series expansion tests
symengine/tests/basic/test_subs.cpp         → Substitution tests
symengine/tests/basic/test_solve.cpp        → Solving tests
symengine/tests/basic/test_matrix.cpp       → Matrix tests
symengine/tests/basic/test_polynomial.cpp   → Polynomial tests
symengine/tests/basic/test_ntheory.cpp      → Number theory tests
```

### Testing Approach
1. Port each C++ test case to TypeScript/Vitest
2. Preserve numerical tolerances and test coverage
3. Add TypeScript-specific tests for API ergonomics
4. Test WASM memory management (no leaks)
5. Verify performance vs pure JavaScript implementations

---

## Build Dependencies

### Required
- **SymEngine C++ library** (already in `/packages/symwasm/reference/symengine/`)
- **Emscripten** (for C++ → WASM compilation)
- **CMake** (SymEngine build system)
- **GMP** (GNU Multiple Precision library) — for arbitrary precision integers

### Optional (for enhanced features)
- **MPFR** (arbitrary precision floating point)
- **MPC** (arbitrary precision complex numbers)
- **FLINT** (Fast Library for Number Theory) — polynomial performance
- **LLVM** (for JIT compilation via `llvm_double.h`)

### Build Configuration
- Use SymEngine's existing CMakeLists.txt
- Add Emscripten toolchain file
- Configure optional dependencies based on needed features
- Export C API via cwrapper for JS interop
