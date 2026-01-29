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

#### 1.4 Number Types
**SymEngine Files**: `number.h`, `integer.h`, `rational.h`, `real_double.h`, `complex_double.h`, `complex.h`
- 🔲 `Integer(value)` — Exact integer (maps to `Integer` class)
- 🔲 `Rational(p, q)` — Exact rational p/q (maps to `Rational` class)
- 🔲 `Float(value, precision?)` — Arbitrary precision float (maps to `RealDouble` or `RealMPFR`)
- ⬜ `Complex(re, im)` — Exact complex number (maps to `Complex` class)
- ⬜ `RealDouble(value)` — Machine precision real (C++ `RealDouble`)
- ⬜ `ComplexDouble(re, im)` — Machine precision complex (C++ `ComplexDouble`)
- ⬜ `RealMPFR(value, precision)` — Arbitrary precision real (requires MPFR)
- ⬜ `ComplexMPC(re, im, precision)` — Arbitrary precision complex (requires MPC)

#### 1.5 Core Arithmetic Operations
**SymEngine Files**: `add.h`, `mul.h`, `pow.h`
- 🔲 `Add(args)` — Symbolic addition (C++ `Add` class)
- 🔲 `Mul(args)` — Symbolic multiplication (C++ `Mul` class)
- 🔲 `Pow(base, exp)` — Symbolic exponentiation (C++ `Pow` class)
- ⬜ Operator overloading support in TypeScript wrappers

#### 1.6 Constants
**SymEngine Files**: `constants.h`, `constants.cpp`
- 🔲 `pi` — Pi (π) → `Pi` class
- 🔲 `E` — Euler's number (e) → `E_` class
- 🔲 `I` — Imaginary unit (i) → `I_` class
- 🔲 `oo` — Infinity → `Infty` class
- 🔲 `S.Zero`, `S.One`, `S.NegativeOne`, `S.Half` — Numeric constants
- ⬜ `EulerGamma` — Euler-Mascheroni constant γ → `EulerGamma` class
- ⬜ `Catalan` — Catalan's constant → `Catalan` class
- ⬜ `GoldenRatio` — Golden ratio φ → `GoldenRatio` class
- ⬜ `NaN` — Not a number → `NaN_` class

#### 1.7 Basic Expression Manipulation
**SymEngine Files**: `subs.h`, `subs.cpp`
- 🔲 `Expr.subs(old, new)` — Substitution (`basic_subs`)
- ⬜ `Expr.subs(dict)` — Multiple substitutions (`basic_subs` with map)
- ⬜ `Expr.xreplace(dict)` — Exact structural replacement

#### 1.8 Numerical Evaluation
**SymEngine Files**: `eval.h`, `eval_double.h`, `eval_mpfr.h`
- 🔲 `Expr.evalf(precision?)` — Numerical evaluation (`eval_double`, `evalf`)
- ⬜ `N(expr, n)` — Numerical evaluation (alias)
- ⬜ `evalf_complex(expr, precision)` — Complex numerical evaluation

---

### Phase 2: Essential Functions (Calculus & Simplification)
**Goal**: Enable symbolic calculus and expression manipulation
**C++ Files**: `functions.h`, `derivative.h`, `series.h`, `expand.h`, `subs.h`

#### 2.1 Elementary Functions
**SymEngine Files**: `functions.h`, `functions.cpp`

##### Exponential & Logarithmic
- ⬜ `exp(x)` — Exponential e^x → `Exp` class
- ⬜ `log(x, base?)` — Natural log or log_base → `Log` class
- ⬜ `sqrt(x)` — Square root → `sqrt()` function
- ⬜ `cbrt(x)` — Cube root
- ⬜ `Abs(x)` — Absolute value → `Abs` class

##### Trigonometric Functions
- ⬜ `sin(x)` — Sine → `Sin` class
- ⬜ `cos(x)` — Cosine → `Cos` class
- ⬜ `tan(x)` — Tangent → `Tan` class
- ⬜ `cot(x)` — Cotangent → `Cot` class
- ⬜ `sec(x)` — Secant → `Sec` class
- ⬜ `csc(x)` — Cosecant → `Csc` class

##### Inverse Trigonometric
- ⬜ `asin(x)` — Arcsine → `ASin` class
- ⬜ `acos(x)` — Arccosine → `ACos` class
- ⬜ `atan(x)` — Arctangent → `ATan` class
- ⬜ `acot(x)` — Arccotangent → `ACot` class
- ⬜ `asec(x)` — Arcsecant → `ASec` class
- ⬜ `acsc(x)` — Arccosecant → `ACsc` class
- ⬜ `atan2(y, x)` — Two-argument arctangent → `ATan2` class

##### Hyperbolic Functions
- ⬜ `sinh(x)` — Hyperbolic sine → `Sinh` class
- ⬜ `cosh(x)` — Hyperbolic cosine → `Cosh` class
- ⬜ `tanh(x)` — Hyperbolic tangent → `Tanh` class
- ⬜ `coth(x)` — Hyperbolic cotangent → `Coth` class
- ⬜ `sech(x)` — Hyperbolic secant
- ⬜ `csch(x)` — Hyperbolic cosecant

##### Inverse Hyperbolic
- ⬜ `asinh(x)` — Inverse hyperbolic sine → `ASinh` class
- ⬜ `acosh(x)` — Inverse hyperbolic cosine → `ACosh` class
- ⬜ `atanh(x)` — Inverse hyperbolic tangent → `ATanh` class
- ⬜ `acoth(x)` — Inverse hyperbolic cotangent → `ACoth` class
- ⬜ `asech(x)` — Inverse hyperbolic secant
- ⬜ `acsch(x)` — Inverse hyperbolic cosecant

##### Special Functions
**SymEngine Files**: `functions.h` (Gamma, Beta, Erf, etc. classes)
- ⬜ `gamma(x)` — Gamma function → `Gamma` class
- ⬜ `loggamma(x)` — Log-gamma → `LogGamma` class
- ⬜ `digamma(x)` — Digamma ψ(x) → `Digamma` class
- ⬜ `polygamma(n, x)` — Polygamma ψ^(n)(x)
- ⬜ `beta(x, y)` — Beta function → `Beta` class
- ⬜ `lowergamma(s, x)` — Lower incomplete gamma → `LowerGamma` class
- ⬜ `uppergamma(s, x)` — Upper incomplete gamma → `UpperGamma` class
- ⬜ `erf(x)` — Error function → `Erf` class
- ⬜ `erfc(x)` — Complementary error function → `Erfc` class
- ⬜ `zeta(s)` — Riemann zeta → `Zeta` class
- ⬜ `dirichlet_eta(s)` — Dirichlet eta → `Dirichlet_eta` class
- ⬜ `lambertw(x)` — Lambert W → `LambertW` class
- ⬜ `KroneckerDelta(i, j)` — Kronecker delta → `KroneckerDelta` class
- ⬜ `LeviCivita(*indices)` — Levi-Civita symbol → `LeviCivita` class
- ⬜ `floor(x)` — Floor function → `Floor` class
- ⬜ `ceiling(x)` — Ceiling function → `Ceiling` class
- ⬜ `Max(...args)` — Maximum → `Max` class
- ⬜ `Min(...args)` — Minimum → `Min` class
- ⬜ `sign(x)` — Sign function → `Sign` class

##### Complex Number Functions
- ⬜ `conjugate(x)` — Complex conjugate → `Conjugate` class
- ⬜ `re(x)` — Real part
- ⬜ `im(x)` — Imaginary part
- ⬜ `arg(x)` — Argument (phase)

#### 2.2 Calculus — Differentiation
**SymEngine Files**: `derivative.h`, `derivative.cpp`
- 🔲 `diff(expr, symbol, n?)` — Differentiate → `Derivative` class and `diff()` function
- ⬜ `Derivative(expr, *symbols)` — Unevaluated derivative class
- ⬜ `diff(expr, symbol1, symbol2, ...)` — Multiple differentiation
- ⬜ `fdiff(expr, argindex)` — Derivative w.r.t. function argument

#### 2.3 Calculus — Series Expansion
**SymEngine Files**: `series.h`, `series_generic.h`, `series_visitor.h`
- 🔲 `series(expr, symbol, point?, n?)` — Power series → `series()` function
- ⬜ `series(expr, x, x0, n, dir)` — Series with direction
- ⬜ `Order(expr)` — Order term O(x^n)

#### 2.4 Simplification
**SymEngine Files**: `expand.h`, `subs.h`, and simplification utilities

##### Current Stubs
- 🔲 `expand(expr)` — Expand expressions → `expand()` function
- 🔲 `simplify(expr)` — Simplify using heuristics
- 🔲 `trigsimp(expr)` — Simplify trigonometric expressions
- 🔲 `radsimp(expr)` — Simplify radicals
- 🔲 `powsimp(expr)` — Simplify powers
- 🔲 `collect(expr, syms)` — Collect terms
- 🔲 `cancel(expr)` — Cancel rational functions

##### Priority Additions
- ⬜ `expand_mul(expr)` — Expand multiplication
- ⬜ `expand_trig(expr)` — Expand trigonometric
- ⬜ `expand_complex(expr)` — Expand re + i*im
- ⬜ `expand_log(expr)` — Expand logarithms
- ⬜ `expand_power_base(expr)` — Expand (a*b)**c
- ⬜ `expand_power_exp(expr)` — Expand a**(b+c)
- ⬜ `together(expr)` — Combine over common denominator
- ⬜ `apart(expr, x)` — Partial fractions
- ⬜ `numer(expr)` — Extract numerator
- ⬜ `denom(expr)` — Extract denominator

---

### Phase 3: Advanced Mathematics (Matrices, Polynomials, Solvers)
**Goal**: Linear algebra, polynomial operations, equation solving
**C++ Files**: `matrix.h`, `polys/`, `solve.h`

#### 3.1 Matrix Operations
**SymEngine Files**: `matrix.h`, `matrices/` subdirectory (24 headers)

##### Current Stubs — Dense Matrices
- 🔲 `Matrix(data)` — Dense matrix → `DenseMatrix` class
- 🔲 `Matrix.get(i, j)` — Get element
- 🔲 `Matrix.det()` — Determinant
- 🔲 `Matrix.inv()` — Inverse
- 🔲 `Matrix.transpose()` — Transpose
- 🔲 `Matrix.eigenvals()` — Eigenvalues
- 🔲 `Matrix.eigenvects()` — Eigenvectors
- 🔲 `Matrix.rref()` — Row echelon form
- 🔲 `eye(n)` — Identity matrix
- 🔲 `zeros(rows, cols)` — Zero matrix
- 🔲 `ones(rows, cols)` — Matrix of ones
- 🔲 `diag(...values)` — Diagonal matrix

##### Priority Additions — Matrix Construction
**SymEngine Files**: `matrix.h`
- ⬜ `DenseMatrix.from_list(list)` — From nested list
- ⬜ `DenseMatrix.from_flat(flat, rows, cols)` — From flat array

##### Priority Additions — Matrix Properties & Operations
**SymEngine Files**: `matrix.h`
- ⬜ `Matrix.rows` — Number of rows
- ⬜ `Matrix.cols` — Number of columns
- ⬜ `Matrix.shape` — (rows, cols)
- ⬜ `Matrix.add(other)` — Matrix addition
- ⬜ `Matrix.mul(other)` — Matrix multiplication
- ⬜ `Matrix.trace()` — Trace
- ⬜ `Matrix.conjugate()` — Conjugate
- ⬜ `Matrix.submatrix(i1, i2, j1, j2)` — Extract submatrix
- ⬜ `Matrix.row_join(other)` — Horizontal stack
- ⬜ `Matrix.col_join(other)` — Vertical stack
- ⬜ `Matrix.row_del(i)` — Delete row
- ⬜ `Matrix.col_del(j)` — Delete column

##### Priority Additions — Matrix Factorizations
**SymEngine Files**: `matrix.h` (LU, LDL, QR methods)
- ⬜ `Matrix.LU()` — LU decomposition
- ⬜ `Matrix.LDL()` — LDL decomposition
- ⬜ `Matrix.FFLU()` — Fraction-free LU
- ⬜ `Matrix.LU_solve(b)` — Solve using LU

##### Priority Additions — Matrix Calculus
**SymEngine Files**: `derivative.h`, `matrix.h`
- ⬜ `Matrix.diff(x)` — Differentiate each element
- ⬜ `jacobian(exprs, vars)` — Jacobian matrix

##### Symbolic Matrices
**SymEngine Files**: `matrices/matrix_symbol.h`, `matrices/identity_matrix.h`, etc.
- ⬜ `MatrixSymbol(name, n, m)` — Symbolic matrix
- ⬜ `Identity(n)` — Identity matrix symbol
- ⬜ `ZeroMatrix(n, m)` — Zero matrix symbol
- ⬜ `DiagonalMatrix(diag)` — Diagonal matrix symbol
- ⬜ `MatrixAdd` — Symbolic matrix addition
- ⬜ `MatrixMul` — Symbolic matrix multiplication
- ⬜ `HadamardProduct` — Element-wise product
- ⬜ `Trace` — Trace of symbolic matrix

##### Sparse Matrices
**SymEngine Files**: `matrix.h` (CSR/CSC support mentioned in cwrapper)
- ⬜ `SparseMatrix(rows, cols)` — Sparse matrix
- ⬜ `SparseMatrix.to_dense()` — Convert to dense

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
| 1.3 | Symbols | 0 | 2 | ~2 | ~4 | symbol.h |
| 1.4 | Numbers | 0 | 3 | ~5 | ~8 | number.h, integer.h, rational.h, complex*.h |
| 1.5 | Arithmetic | 0 | 3 | ~1 | ~4 | add.h, mul.h, pow.h |
| 1.6 | Constants | 0 | 6 | ~4 | ~10 | constants.h |
| 1.7 | Substitution | 0 | 1 | ~2 | ~3 | subs.h |
| 1.8 | Evaluation | 0 | 1 | ~2 | ~3 | eval*.h |
| **2** | **Essential Functions** | | | | | |
| 2 | Functions | 0 | 0 | ~50 | ~50 | functions.h |
| 2 | Calculus | 0 | 3 | ~7 | ~10 | derivative.h, series.h |
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
| **Total** | | **~15** | **48** | **~205** | **~268** | **65 main headers + subdirs** |

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
