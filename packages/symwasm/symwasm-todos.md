# symwasm Implementation Todos

Comprehensive API for a SymPy-inspired symbolic mathematics library in TypeScript, powered by SymEngine compiled to WebAssembly.

**Strategy**: Copy SymEngine C++ kernels as-is, compile to WebAssembly, and write thin TypeScript wrappers using the SymPy-like API.

Legend: ✅ = implemented, 🔲 = stubbed (exists but throws NotImplementedError), ⬜ = not yet created

---

## Table of Contents

1. [Core Foundation](#1-core-foundation)
2. [Elementary Functions](#2-elementary-functions)
3. [Calculus](#3-calculus)
4. [Simplification](#4-simplification)
5. [Matrices & Linear Algebra](#5-matrices--linear-algebra)
6. [Polynomials](#6-polynomials)
7. [Equation Solving](#7-equation-solving)
8. [Number Theory](#8-number-theory)
9. [Sets & Logic](#9-sets--logic)
10. [Printing & I/O](#10-printing--io)
11. [Lambda & Numerical](#11-lambda--numerical)
12. [Assumptions & Queries](#12-assumptions--queries)

---

## 1. Core Foundation

### 1.1 Symbols & Variables
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `Symbol(name, assumptions?)` | ✅ | Create symbolic variable | `core/symbol.ts` |
| `symbols(names, assumptions?)` | ✅ | Create multiple symbols from string | `core/symbols.ts` |
| `Dummy(name?)` | ⬜ | Temporary/unique symbol | `core/dummy.ts` |
| `Wild(name)` | ⬜ | Wildcard for pattern matching | `core/wild.ts` |

### 1.2 Number Types
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `Integer(value)` | ✅ | Exact integer | `core/numbers/integer.ts` |
| `Rational(p, q)` | ✅ | Exact rational p/q | `core/numbers/rational.ts` |
| `Float(value, precision?)` | ✅ | Floating-point number | `core/numbers/float.ts` |
| `Complex(re, im)` | ✅ | Complex number | `core/numbers/complex.ts` |
| `ComplexDouble(re, im)` | ⬜ | Machine-precision complex | `core/numbers/complex-double.ts` |
| `RealMPFR(value, precision)` | ⬜ | Arbitrary-precision real | `core/numbers/real-mpfr.ts` |
| `ComplexMPC(re, im, precision)` | ⬜ | Arbitrary-precision complex | `core/numbers/complex-mpc.ts` |

### 1.3 Arithmetic Operations
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `add(a, b)` | ✅ | Addition | `core/operations/add.ts` |
| `sub(a, b)` | ✅ | Subtraction | `core/operations/sub.ts` |
| `mul(a, b)` | ✅ | Multiplication | `core/operations/mul.ts` |
| `div(a, b)` | ✅ | Division | `core/operations/div.ts` |
| `pow(base, exp)` | ✅ | Exponentiation | `core/operations/pow.ts` |
| `neg(a)` | ✅ | Negation | `core/operations/neg.ts` |

### 1.4 Constants
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `pi` | ✅ | π = 3.14159... | `core/constants/pi.ts` |
| `E` | ✅ | Euler's number e | `core/constants/e.ts` |
| `I` | ✅ | Imaginary unit i | `core/constants/i.ts` |
| `oo` | ✅ | Positive infinity ∞ | `core/constants/infinity.ts` |
| `EulerGamma` | ✅ | Euler-Mascheroni γ | `core/constants/euler-gamma.ts` |
| `Catalan` | ✅ | Catalan's constant | `core/constants/catalan.ts` |
| `GoldenRatio` | ✅ | Golden ratio φ | `core/constants/golden-ratio.ts` |
| `S.Zero` | ✅ | Integer 0 | `core/constants/singletons.ts` |
| `S.One` | ✅ | Integer 1 | `core/constants/singletons.ts` |
| `S.Half` | ✅ | Rational 1/2 | `core/constants/singletons.ts` |
| `S.NegativeOne` | ✅ | Integer -1 | `core/constants/singletons.ts` |
| `S.Infinity` | ✅ | Positive infinity | `core/constants/singletons.ts` |
| `S.NegativeInfinity` | ✅ | Negative infinity | `core/constants/singletons.ts` |
| `S.ComplexInfinity` | ✅ | Complex infinity | `core/constants/singletons.ts` |
| `S.NaN` | ✅ | Not a number | `core/constants/singletons.ts` |

### 1.5 Expression Manipulation
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `Expr.subs(old, new)` | ✅ | Single substitution | `core/expr.ts` |
| `Expr.subs(Map)` | ✅ | Multiple simultaneous substitutions | `core/expr.ts` |
| `Expr.subs({...})` | ✅ | Object-notation substitution | `core/expr.ts` |
| `Expr.xreplace(dict)` | ⬜ | Exact structural replacement | `core/expr.ts` |
| `Expr.free_symbols()` | ✅ | Get free symbols | `core/expr.ts` |
| `Expr.get_args()` | ✅ | Get sub-expressions | `core/expr.ts` |
| `Expr.equals(other)` | ✅ | Structural equality | `core/expr.ts` |
| `Expr.hash()` | ✅ | Hash code | `core/expr.ts` |
| `Expr.get_type()` | ✅ | Type identification | `core/expr.ts` |

### 1.6 Numerical Evaluation
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `Expr.evalf(precision?)` | ✅ | Numerical evaluation | `core/expr.ts` |
| `Expr.evalfNumber()` | ✅ | Extract JS number | `core/expr.ts` |
| `Expr.evalfComplex()` | ✅ | Extract complex {real, imag} | `core/expr.ts` |
| `N(expr, n)` | ⬜ | Numerical evaluation alias | `core/numerical/n.ts` |

---

## 2. Elementary Functions

### 2.1 Exponential & Logarithmic
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `exp(x)` | ✅ | Exponential e^x | `functions/exp.ts` |
| `log(x)` | ✅ | Natural logarithm | `functions/log.ts` |
| `log(x, base)` | ⬜ | Logarithm with base | `functions/log.ts` |
| `sqrt(x)` | ✅ | Square root | `functions/sqrt.ts` |
| `cbrt(x)` | ✅ | Cube root | `functions/cbrt.ts` |
| `root(x, n)` | ⬜ | nth root | `functions/root.ts` |
| `abs(x)` | ✅ | Absolute value | `functions/abs.ts` |
| `sign(x)` | ✅ | Sign function | `functions/sign.ts` |
| `floor(x)` | ✅ | Floor | `functions/floor.ts` |
| `ceiling(x)` | ✅ | Ceiling | `functions/ceiling.ts` |
| `lambertw(x)` | ✅ | Lambert W function | `functions/lambertw.ts` |

### 2.2 Trigonometric
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `sin(x)` | ✅ | Sine | `functions/trig/sin.ts` |
| `cos(x)` | ✅ | Cosine | `functions/trig/cos.ts` |
| `tan(x)` | ✅ | Tangent | `functions/trig/tan.ts` |
| `cot(x)` | ✅ | Cotangent | `functions/trig/cot.ts` |
| `sec(x)` | ✅ | Secant | `functions/trig/sec.ts` |
| `csc(x)` | ✅ | Cosecant | `functions/trig/csc.ts` |
| `asin(x)` | ✅ | Arcsine | `functions/trig/asin.ts` |
| `acos(x)` | ✅ | Arccosine | `functions/trig/acos.ts` |
| `atan(x)` | ✅ | Arctangent | `functions/trig/atan.ts` |
| `atan2(y, x)` | ✅ | Two-argument arctangent | `functions/trig/atan2.ts` |
| `acot(x)` | ✅ | Arccotangent | `functions/trig/acot.ts` |
| `asec(x)` | ✅ | Arcsecant | `functions/trig/asec.ts` |
| `acsc(x)` | ✅ | Arccosecant | `functions/trig/acsc.ts` |

### 2.3 Hyperbolic
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `sinh(x)` | ✅ | Hyperbolic sine | `functions/hyperbolic/sinh.ts` |
| `cosh(x)` | ✅ | Hyperbolic cosine | `functions/hyperbolic/cosh.ts` |
| `tanh(x)` | ✅ | Hyperbolic tangent | `functions/hyperbolic/tanh.ts` |
| `coth(x)` | ✅ | Hyperbolic cotangent | `functions/hyperbolic/coth.ts` |
| `sech(x)` | ✅ | Hyperbolic secant | `functions/hyperbolic/sech.ts` |
| `csch(x)` | ✅ | Hyperbolic cosecant | `functions/hyperbolic/csch.ts` |
| `asinh(x)` | ✅ | Inverse hyperbolic sine | `functions/hyperbolic/asinh.ts` |
| `acosh(x)` | ✅ | Inverse hyperbolic cosine | `functions/hyperbolic/acosh.ts` |
| `atanh(x)` | ✅ | Inverse hyperbolic tangent | `functions/hyperbolic/atanh.ts` |
| `acoth(x)` | ✅ | Inverse hyperbolic cotangent | `functions/hyperbolic/acoth.ts` |
| `asech(x)` | ✅ | Inverse hyperbolic secant | `functions/hyperbolic/asech.ts` |
| `acsch(x)` | ✅ | Inverse hyperbolic cosecant | `functions/hyperbolic/acsch.ts` |

### 2.4 Special Functions
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `gamma(x)` | ✅ | Gamma function Γ(x) | `functions/special/gamma.ts` |
| `loggamma(x)` | ✅ | Log-gamma ln(Γ(x)) | `functions/special/loggamma.ts` |
| `digamma(x)` | ✅ | Digamma ψ(x) | `functions/special/digamma.ts` |
| `polygamma(n, x)` | ✅ | Polygamma ψ^(n)(x) | `functions/special/polygamma.ts` |
| `beta(x, y)` | ✅ | Beta function B(x,y) | `functions/special/beta.ts` |
| `lowergamma(s, x)` | ✅ | Lower incomplete gamma | `functions/special/lowergamma.ts` |
| `uppergamma(s, x)` | ✅ | Upper incomplete gamma | `functions/special/uppergamma.ts` |
| `erf(x)` | ✅ | Error function | `functions/special/erf.ts` |
| `erfc(x)` | ✅ | Complementary error function | `functions/special/erfc.ts` |
| `zeta(s)` | ✅ | Riemann zeta ζ(s) | `functions/special/index.ts` |
| `dirichlet_eta(s)` | ✅ | Dirichlet eta η(s) | `functions/special/index.ts` |
| `kronecker_delta(i, j)` | ✅ | Kronecker delta δ_ij | `functions/special/index.ts` |

### 2.5 Complex Functions
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `conjugate(x)` | ✅ | Complex conjugate | `functions/complex/conjugate.ts` |
| `re(x)` | ✅ | Real part | `functions/complex/re.ts` |
| `im(x)` | ✅ | Imaginary part | `functions/complex/im.ts` |
| `arg(x)` | ✅ | Argument (phase) | `functions/complex/arg.ts` |
| `Abs(x)` | ⬜ | Complex absolute value | `functions/complex/abs.ts` |
| `polar_lift(x)` | ⬜ | Polar representation | `functions/complex/polar-lift.ts` |

### 2.6 Min/Max/Comparison
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `Max(...args)` | ✅ | Maximum | `functions/minmax/max.ts` |
| `Min(...args)` | ✅ | Minimum | `functions/minmax/min.ts` |

---

## 3. Calculus

### 3.1 Differentiation
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `diff(expr, x)` | ✅ | First derivative | `calculus/diff.ts` |
| `diff(expr, x, n)` | ✅ | nth derivative | `calculus/diff.ts` |
| `diff(expr, x, y, ...)` | ✅ | Partial derivatives | `calculus/diff.ts` |
| `Derivative(expr, *symbols)` | ⬜ | Unevaluated derivative | `calculus/derivative.ts` |
| `fdiff(expr, argindex)` | ⬜ | Derivative w.r.t. argument | `calculus/fdiff.ts` |

### 3.2 Integration
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `integrate(expr, x)` | 🔲 | Indefinite integral | `calculus/integrate.ts` |
| `integrate(expr, (x, a, b))` | 🔲 | Definite integral | `calculus/integrate.ts` |
| `Integral(expr, *limits)` | ⬜ | Unevaluated integral | `calculus/integral.ts` |
| `line_integrate(field, curve, params)` | ⬜ | Line integral | `calculus/line-integrate.ts` |

### 3.3 Limits
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `limit(expr, x, x0)` | 🔲 | Limit as x→x0 | `calculus/limit.ts` |
| `limit(expr, x, x0, '+')` | ⬜ | Right-hand limit | `calculus/limit.ts` |
| `limit(expr, x, x0, '-')` | ⬜ | Left-hand limit | `calculus/limit.ts` |
| `Limit(expr, x, x0, dir)` | ⬜ | Unevaluated limit | `calculus/limit-class.ts` |

### 3.4 Series
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `series(expr, x, x0?, n?)` | ✅ | Taylor series expansion | `calculus/series.ts` |
| `Order(expr)` | ⬜ | Order term O(x^n) | `calculus/order.ts` |
| `fourier_series(f, (x, a, b))` | ⬜ | Fourier series | `calculus/fourier-series.ts` |

### 3.5 Summation & Products
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `summation(f, (i, a, b))` | 🔲 | Symbolic summation Σ | `calculus/summation.ts` |
| `product(f, (i, a, b))` | ⬜ | Symbolic product Π | `calculus/product.ts` |
| `Sum(f, limits)` | ⬜ | Unevaluated sum | `calculus/sum.ts` |
| `Product(f, limits)` | ⬜ | Unevaluated product | `calculus/product-class.ts` |

---

## 4. Simplification

| Function | Status | Description | File |
|----------|--------|-------------|------|
| `simplify(expr)` | ✅ | General simplification | `simplify/simplify.ts` |
| `expand(expr)` | ✅ | Expand products/powers | `simplify/expand.ts` |
| `trigsimp(expr)` | ✅ | Simplify trig expressions | `simplify/trigsimp.ts` |
| `expand_trig(expr)` | ✅ | Expand trig functions | `simplify/expand-trig.ts` |
| `radsimp(expr)` | ✅ | Simplify radicals | `simplify/radsimp.ts` |
| `powsimp(expr)` | ✅ | Simplify powers | `simplify/powsimp.ts` |
| `expand_complex(expr)` | ✅ | Expand complex to re+i*im | `simplify/expand-complex.ts` |
| `numer(expr)` | ✅ | Extract numerator | `simplify/numer.ts` |
| `denom(expr)` | ✅ | Extract denominator | `simplify/denom.ts` |
| `rewrite_as_exp(expr)` | ✅ | Rewrite as exponentials | `simplify/rewrite-as-exp.ts` |
| `rewrite_as_sin(expr)` | ✅ | Rewrite in terms of sine | `simplify/rewrite-as-sin.ts` |
| `rewrite_as_cos(expr)` | ✅ | Rewrite in terms of cosine | `simplify/rewrite-as-cos.ts` |
| `as_real_imag(expr)` | ✅ | Extract real/imag parts | `simplify/as-real-imag.ts` |
| `cse(exprs)` | ✅ | Common subexpression elimination | `simplify/cse.ts` |

---

## 5. Matrices & Linear Algebra

### 5.1 Matrix Construction
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `Matrix(data)` | ✅ | Create from nested array | `matrices/matrix.ts` |
| `Matrix.fromFlat(flat, rows, cols)` | ✅ | Create from flat array | `matrices/matrix.ts` |
| `eye(n, m?, k?)` | ✅ | Identity matrix | `matrices/eye.ts` |
| `zeros(rows, cols)` | ✅ | Zero matrix | `matrices/zeros.ts` |
| `ones(rows, cols)` | ✅ | Ones matrix | `matrices/ones.ts` |
| `diag(values, k?)` | ✅ | Diagonal matrix | `matrices/diag.ts` |

### 5.2 Matrix Properties
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `Matrix.get(i, j)` | ✅ | Get element | `matrices/matrix.ts` |
| `Matrix.set(i, j, val)` | ✅ | Set element | `matrices/matrix.ts` |
| `Matrix.rows` | ✅ | Row count | `matrices/matrix.ts` |
| `Matrix.cols` | ✅ | Column count | `matrices/matrix.ts` |
| `Matrix.shape` | ✅ | Dimensions tuple | `matrices/matrix.ts` |
| `Matrix.equals(other)` | ✅ | Equality check | `matrices/matrix.ts` |

### 5.3 Matrix Operations
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `Matrix.add(other)` | ✅ | Addition | `matrices/matrix.ts` |
| `Matrix.mul(other)` | ✅ | Multiplication | `matrices/matrix.ts` |
| `Matrix.addScalar(k)` | ✅ | Scalar addition | `matrices/matrix.ts` |
| `Matrix.mulScalar(k)` | ✅ | Scalar multiplication | `matrices/matrix.ts` |
| `Matrix.transpose()` | ✅ | Transpose | `matrices/matrix.ts` |
| `Matrix.det()` | ✅ | Determinant | `matrices/matrix.ts` |
| `Matrix.inv()` | ✅ | Inverse | `matrices/matrix.ts` |

### 5.4 Matrix Suboperations
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `Matrix.submatrix(r1, c1, r2, c2)` | ✅ | Extract submatrix | `matrices/matrix.ts` |
| `Matrix.rowJoin(other)` | ✅ | Horizontal concatenation | `matrices/matrix.ts` |
| `Matrix.colJoin(other)` | ✅ | Vertical concatenation | `matrices/matrix.ts` |
| `Matrix.rowDel(k)` | ✅ | Delete row | `matrices/matrix.ts` |
| `Matrix.colDel(k)` | ✅ | Delete column | `matrices/matrix.ts` |

### 5.5 Matrix Factorizations
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `Matrix.lu()` | ✅ | LU decomposition | `matrices/matrix.ts` |
| `Matrix.ldl()` | ✅ | LDL decomposition | `matrices/matrix.ts` |
| `Matrix.fflu()` | ✅ | Fraction-free LU | `matrices/matrix.ts` |
| `Matrix.ffldu()` | ✅ | Fraction-free LDU | `matrices/matrix.ts` |
| `Matrix.luSolve(b)` | ✅ | Solve Ax=b via LU | `matrices/matrix.ts` |

### 5.6 Matrix Calculus
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `Matrix.diff(x)` | ✅ | Differentiate elements | `matrices/matrix.ts` |
| `jacobian(funcs, vars)` | ✅ | Jacobian matrix | `matrices/index.ts` |

### 5.7 Sparse Matrices
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `SparseMatrix(rows, cols)` | ✅ | Create sparse matrix | `matrices/index.ts` |
| `SparseMatrix.get(i, j)` | ✅ | Get element | `matrices/index.ts` |
| `SparseMatrix.set(i, j, val)` | ✅ | Set element | `matrices/index.ts` |
| `SparseMatrix.toString()` | ✅ | String representation | `matrices/index.ts` |
| `SparseMatrix.equals(other)` | ✅ | Equality check | `matrices/index.ts` |

---

## 6. Polynomials

### 6.1 Polynomial Construction
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `Poly(expr, *gens)` | ⬜ | Create polynomial | `polys/poly.ts` |
| `Poly.from_list(coeffs, gen)` | ⬜ | From coefficient list | `polys/poly.ts` |
| `Poly.from_dict(terms, gen)` | ⬜ | From term dictionary | `polys/poly.ts` |

### 6.2 Polynomial Properties
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `degree(poly, gen?)` | ⬜ | Degree | `polys/degree.ts` |
| `degree_list(poly)` | ⬜ | Multi-variable degrees | `polys/degree-list.ts` |
| `LC(poly)` | ⬜ | Leading coefficient | `polys/lc.ts` |
| `LT(poly)` | ⬜ | Leading term | `polys/lt.ts` |
| `LM(poly)` | ⬜ | Leading monomial | `polys/lm.ts` |
| `TC(poly)` | ⬜ | Trailing coefficient | `polys/tc.ts` |
| `coeffs(poly)` | ⬜ | All coefficients | `polys/coeffs.ts` |
| `monoms(poly)` | ⬜ | All monomials | `polys/monoms.ts` |
| `terms(poly)` | ⬜ | All terms | `polys/terms.ts` |
| `nth(poly, n)` | ⬜ | nth coefficient | `polys/nth.ts` |

### 6.3 Polynomial Arithmetic
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `pdiv(f, g)` | ⬜ | Division with quotient & remainder | `polys/pdiv.ts` |
| `pquo(f, g)` | ⬜ | Quotient only | `polys/pquo.ts` |
| `prem(f, g)` | ⬜ | Remainder only | `polys/prem.ts` |
| `pexquo(f, g)` | ⬜ | Exact quotient | `polys/pexquo.ts` |

### 6.4 Polynomial GCD/LCM
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `gcd(f, g)` | ⬜ | Greatest common divisor | `polys/gcd.ts` |
| `lcm(f, g)` | ⬜ | Least common multiple | `polys/lcm.ts` |
| `gcdex(f, g)` | ⬜ | Extended GCD (Bezout) | `polys/gcdex.ts` |
| `resultant(f, g, x)` | ⬜ | Resultant | `polys/resultant.ts` |
| `discriminant(f, x)` | ⬜ | Discriminant | `polys/discriminant.ts` |
| `subresultants(f, g, x)` | ⬜ | Subresultant PRS | `polys/subresultants.ts` |

### 6.5 Polynomial Factorization
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `factor(poly)` | 🔲 | Factor polynomial | `polys/factor.ts` |
| `factor_list(poly)` | ⬜ | List of (factor, multiplicity) | `polys/factor-list.ts` |
| `sqf(poly)` | ⬜ | Square-free factorization | `polys/sqf.ts` |
| `sqf_list(poly)` | ⬜ | Square-free factor list | `polys/sqf-list.ts` |
| `decompose(f)` | ⬜ | Functional decomposition | `polys/decompose.ts` |

### 6.6 Polynomial Roots
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `roots(poly)` | ⬜ | All roots (symbolic) | `polys/roots.ts` |
| `nroots(poly, n?)` | ⬜ | Numerical roots | `polys/nroots.ts` |
| `real_roots(poly)` | ⬜ | Real roots only | `polys/real-roots.ts` |
| `complex_roots(poly)` | ⬜ | Complex roots | `polys/complex-roots.ts` |
| `RootOf(poly, index)` | ⬜ | Indexed root | `polys/root-of.ts` |
| `CRootOf(poly, index)` | ⬜ | Complex indexed root | `polys/croot-of.ts` |

### 6.7 Polynomial Evaluation
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `Poly.eval(x, a)` | ⬜ | Evaluate at point | `polys/poly.ts` |
| `Poly.all_roots()` | ⬜ | All roots | `polys/poly.ts` |
| `Poly.count_roots()` | ⬜ | Count roots in interval | `polys/poly.ts` |
| `Poly.intervals()` | ⬜ | Isolating intervals | `polys/poly.ts` |

### 6.8 Polynomial Conversion
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `Poly.as_expr()` | ⬜ | Convert to expression | `polys/poly.ts` |
| `Poly.content()` | ⬜ | Content (GCD of coeffs) | `polys/poly.ts` |
| `Poly.primitive()` | ⬜ | Primitive part | `polys/poly.ts` |
| `Poly.monic()` | ⬜ | Monic form | `polys/poly.ts` |

### 6.9 Orthogonal Polynomials
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `chebyshevt(n, x)` | ⬜ | Chebyshev T_n(x) | `polys/orthogonal/chebyshevt.ts` |
| `chebyshevu(n, x)` | ⬜ | Chebyshev U_n(x) | `polys/orthogonal/chebyshevu.ts` |
| `legendre(n, x)` | ⬜ | Legendre P_n(x) | `polys/orthogonal/legendre.ts` |
| `hermite(n, x)` | ⬜ | Hermite H_n(x) | `polys/orthogonal/hermite.ts` |
| `laguerre(n, x)` | ⬜ | Laguerre L_n(x) | `polys/orthogonal/laguerre.ts` |
| `jacobi(n, a, b, x)` | ⬜ | Jacobi P_n^(a,b)(x) | `polys/orthogonal/jacobi.ts` |
| `gegenbauer(n, a, x)` | ⬜ | Gegenbauer C_n^a(x) | `polys/orthogonal/gegenbauer.ts` |

---

## 7. Equation Solving

### 7.1 Algebraic Equations
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `solve(expr, symbols?)` | 🔲 | Solve equation(s) | `solvers/solve.ts` |
| `solveset(expr, symbol, domain?)` | 🔲 | Solve returning set | `solvers/solveset.ts` |
| `solve_poly(poly, x)` | ⬜ | Solve polynomial | `solvers/solve-poly.ts` |
| `solve_rational(expr, x)` | ⬜ | Solve rational equation | `solvers/solve-rational.ts` |
| `solve_trig(expr, x)` | ⬜ | Solve trigonometric | `solvers/solve-trig.ts` |

### 7.2 Systems of Equations
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `linsolve(system, symbols)` | 🔲 | Linear system solver | `solvers/linsolve.ts` |
| `nonlinsolve(system, symbols)` | 🔲 | Nonlinear system solver | `solvers/nonlinsolve.ts` |
| `solve_linear_system(eqs, syms)` | ⬜ | Matrix-based linear solve | `solvers/solve-linear-system.ts` |

### 7.3 Differential Equations
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `dsolve(eq, func?)` | 🔲 | ODE solver | `solvers/dsolve.ts` |
| `pdsolve(eq, func?)` | ⬜ | PDE solver | `solvers/pdsolve.ts` |
| `classify_ode(eq)` | ⬜ | Classify ODE type | `solvers/classify-ode.ts` |
| `checkodesol(eq, sol)` | ⬜ | Verify ODE solution | `solvers/checkodesol.ts` |

### 7.4 Recurrence Relations
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `rsolve(eq, func)` | ⬜ | Recurrence relation solver | `solvers/rsolve.ts` |
| `rsolve_poly(eq, func)` | ⬜ | Polynomial recurrence | `solvers/rsolve-poly.ts` |
| `rsolve_ratio(eq, func)` | ⬜ | Rational recurrence | `solvers/rsolve-ratio.ts` |
| `rsolve_hyper(eq, func)` | ⬜ | Hypergeometric recurrence | `solvers/rsolve-hyper.ts` |

### 7.5 Inequalities
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `reduce_inequalities(ineqs, syms)` | ⬜ | Solve inequalities | `solvers/reduce-inequalities.ts` |
| `solve_univariate_inequality(ineq, x)` | ⬜ | Single-variable inequality | `solvers/solve-univariate-inequality.ts` |
| `solve_rational_inequalities(ineqs, x)` | ⬜ | Rational inequalities | `solvers/solve-rational-inequalities.ts` |

---

## 8. Number Theory

### 8.1 Primality & Factorization
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `isprime(n)` | ⬜ | Primality test | `ntheory/isprime.ts` |
| `nextprime(n)` | ⬜ | Next prime | `ntheory/nextprime.ts` |
| `prevprime(n)` | ⬜ | Previous prime | `ntheory/prevprime.ts` |
| `primepi(n)` | ⬜ | Prime counting π(n) | `ntheory/primepi.ts` |
| `prime(n)` | ⬜ | nth prime | `ntheory/prime.ts` |
| `primorial(n)` | ⬜ | Product of primes ≤ n | `ntheory/primorial.ts` |
| `primerange(a, b)` | ⬜ | Primes in range | `ntheory/primerange.ts` |
| `factorint(n)` | ⬜ | Integer factorization | `ntheory/factorint.ts` |

### 8.2 Divisibility
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `divisors(n)` | ⬜ | List all divisors | `ntheory/divisors.ts` |
| `divisor_count(n)` | ⬜ | Count divisors τ(n) | `ntheory/divisor-count.ts` |
| `divisor_sigma(n, k)` | ⬜ | Sum of divisor powers σ_k(n) | `ntheory/divisor-sigma.ts` |
| `totient(n)` | ⬜ | Euler's totient φ(n) | `ntheory/totient.ts` |
| `mobius(n)` | ⬜ | Möbius function μ(n) | `ntheory/mobius.ts` |

### 8.3 GCD & LCM
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `gcd(*args)` | ⬜ | Greatest common divisor | `ntheory/gcd.ts` |
| `lcm(*args)` | ⬜ | Least common multiple | `ntheory/lcm.ts` |
| `gcdex(a, b)` | ⬜ | Extended GCD | `ntheory/gcdex.ts` |
| `ilcm(*args)` | ⬜ | Integer LCM | `ntheory/ilcm.ts` |
| `igcd(*args)` | ⬜ | Integer GCD | `ntheory/igcd.ts` |

### 8.4 Modular Arithmetic
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `mod(a, m)` | ⬜ | Modulo | `ntheory/mod.ts` |
| `Mod(a, m)` | ⬜ | Symbolic modulo | `ntheory/mod-class.ts` |
| `mod_inverse(a, m)` | ⬜ | Modular inverse | `ntheory/mod-inverse.ts` |
| `is_primitive_root(a, p)` | ⬜ | Primitive root check | `ntheory/is-primitive-root.ts` |
| `primitive_root(p)` | ⬜ | Find primitive root | `ntheory/primitive-root.ts` |
| `discrete_log(a, b, n)` | ⬜ | Discrete logarithm | `ntheory/discrete-log.ts` |
| `crt(m, v)` | ⬜ | Chinese Remainder Theorem | `ntheory/crt.ts` |
| `sqrt_mod(a, p)` | ⬜ | Modular square root | `ntheory/sqrt-mod.ts` |
| `nthroot_mod(a, n, p)` | ⬜ | Modular nth root | `ntheory/nthroot-mod.ts` |
| `quadratic_residue(a, p)` | ⬜ | Quadratic residue | `ntheory/quadratic-residue.ts` |
| `legendre_symbol(a, p)` | ⬜ | Legendre symbol | `ntheory/legendre-symbol.ts` |
| `jacobi_symbol(a, n)` | ⬜ | Jacobi symbol | `ntheory/jacobi-symbol.ts` |

### 8.5 Sequences
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `factorial(n)` | ⬜ | Factorial n! | `ntheory/factorial.ts` |
| `factorial2(n)` | ⬜ | Double factorial n!! | `ntheory/factorial2.ts` |
| `binomial(n, k)` | ⬜ | Binomial coefficient | `ntheory/binomial.ts` |
| `fibonacci(n)` | ⬜ | Fibonacci F_n | `ntheory/fibonacci.ts` |
| `lucas(n)` | ⬜ | Lucas L_n | `ntheory/lucas.ts` |
| `bell(n)` | ⬜ | Bell number B_n | `ntheory/bell.ts` |
| `bernoulli(n)` | ⬜ | Bernoulli B_n | `ntheory/bernoulli.ts` |
| `catalan(n)` | ⬜ | Catalan C_n | `ntheory/catalan.ts` |
| `euler(n)` | ⬜ | Euler number E_n | `ntheory/euler.ts` |
| `stirling(n, k, kind)` | ⬜ | Stirling numbers | `ntheory/stirling.ts` |

### 8.6 Partitions
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `npartitions(n)` | ⬜ | Partition count p(n) | `ntheory/npartitions.ts` |
| `partitions(n)` | ⬜ | Generate partitions | `ntheory/partitions.ts` |

### 8.7 Diophantine Equations
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `diophantine(eq)` | ⬜ | Solve Diophantine equation | `ntheory/diophantine.ts` |
| `diop_solve(eq, syms)` | ⬜ | General Diophantine solver | `ntheory/diop-solve.ts` |
| `diop_linear(eq)` | ⬜ | Linear Diophantine | `ntheory/diop-linear.ts` |
| `diop_quadratic(eq)` | ⬜ | Quadratic Diophantine | `ntheory/diop-quadratic.ts` |

---

## 9. Sets & Logic

### 9.1 Set Types
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `FiniteSet(...elements)` | ⬜ | Finite set | `sets/finite-set.ts` |
| `Interval(a, b, left_open?, right_open?)` | ⬜ | Real interval | `sets/interval.ts` |
| `Union(*sets)` | ⬜ | Union of sets | `sets/union.ts` |
| `Intersection(*sets)` | ⬜ | Intersection | `sets/intersection.ts` |
| `Complement(A, B)` | ⬜ | Set complement | `sets/complement.ts` |
| `SymmetricDifference(A, B)` | ⬜ | Symmetric difference | `sets/symmetric-difference.ts` |
| `ProductSet(*sets)` | ⬜ | Cartesian product | `sets/product-set.ts` |
| `ImageSet(lambda, base_set)` | ⬜ | Image set | `sets/image-set.ts` |
| `ConditionSet(sym, cond, base)` | ⬜ | Conditional set | `sets/condition-set.ts` |

### 9.2 Special Sets
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `EmptySet` | ⬜ | Empty set ∅ | `sets/empty-set.ts` |
| `UniversalSet` | ⬜ | Universal set | `sets/universal-set.ts` |
| `Naturals` | ⬜ | Natural numbers ℕ | `sets/naturals.ts` |
| `Naturals0` | ⬜ | ℕ ∪ {0} | `sets/naturals0.ts` |
| `Integers` | ⬜ | Integers ℤ | `sets/integers.ts` |
| `Rationals` | ⬜ | Rationals ℚ | `sets/rationals.ts` |
| `Reals` | ⬜ | Real numbers ℝ | `sets/reals.ts` |
| `Complexes` | ⬜ | Complex numbers ℂ | `sets/complexes.ts` |

### 9.3 Set Operations
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `set_union(*sets)` | ⬜ | Union | `sets/operations/set-union.ts` |
| `set_intersection(*sets)` | ⬜ | Intersection | `sets/operations/set-intersection.ts` |
| `set_complement(A, B)` | ⬜ | Complement | `sets/operations/set-complement.ts` |
| `contains(set, elem)` | ⬜ | Membership test | `sets/operations/contains.ts` |
| `is_subset(A, B)` | ⬜ | Subset check | `sets/operations/is-subset.ts` |
| `is_superset(A, B)` | ⬜ | Superset check | `sets/operations/is-superset.ts` |
| `is_proper_subset(A, B)` | ⬜ | Proper subset | `sets/operations/is-proper-subset.ts` |
| `Set.boundary` | ⬜ | Boundary | `sets/set.ts` |
| `Set.interior` | ⬜ | Interior | `sets/set.ts` |
| `Set.closure` | ⬜ | Closure | `sets/set.ts` |
| `Set.measure` | ⬜ | Measure (length/area) | `sets/set.ts` |

### 9.4 Boolean Logic
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `And(*args)` | ⬜ | Logical AND | `logic/and.ts` |
| `Or(*args)` | ⬜ | Logical OR | `logic/or.ts` |
| `Not(expr)` | ⬜ | Logical NOT | `logic/not.ts` |
| `Xor(*args)` | ⬜ | Logical XOR | `logic/xor.ts` |
| `Nand(*args)` | ⬜ | Logical NAND | `logic/nand.ts` |
| `Nor(*args)` | ⬜ | Logical NOR | `logic/nor.ts` |
| `Implies(p, q)` | ⬜ | Implication p → q | `logic/implies.ts` |
| `Equivalent(p, q)` | ⬜ | Equivalence p ↔ q | `logic/equivalent.ts` |
| `true` / `false` | ⬜ | Boolean atoms | `logic/boolean-atoms.ts` |

### 9.5 Relational Operators
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `Eq(a, b)` | ⬜ | Equality | `logic/eq.ts` |
| `Ne(a, b)` | ⬜ | Inequality | `logic/ne.ts` |
| `Lt(a, b)` | ⬜ | Less than | `logic/lt.ts` |
| `Le(a, b)` | ⬜ | Less than or equal | `logic/le.ts` |
| `Gt(a, b)` | ⬜ | Greater than | `logic/gt.ts` |
| `Ge(a, b)` | ⬜ | Greater than or equal | `logic/ge.ts` |

### 9.6 Piecewise
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `Piecewise((expr1, cond1), ...)` | ⬜ | Piecewise function | `logic/piecewise.ts` |

---

## 10. Printing & I/O

### 10.1 String Representation
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `sstr(expr)` | 🔲 | Simple string | `printing/sstr.ts` |
| `srepr(expr)` | ⬜ | Repr-style string | `printing/srepr.ts` |
| `pretty(expr, opts?)` | 🔲 | Unicode pretty-print | `printing/pretty.ts` |
| `pprint(expr)` | ⬜ | Print pretty | `printing/pprint.ts` |
| `tree(expr)` | ⬜ | Tree structure view | `printing/tree.ts` |

### 10.2 Export Formats
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `latex(expr, opts?)` | 🔲 | LaTeX output | `printing/latex.ts` |
| `mathml(expr, printer?)` | 🔲 | MathML output | `printing/mathml.ts` |
| `dotprint(expr)` | ⬜ | Graphviz DOT format | `printing/dotprint.ts` |

### 10.3 Code Generation
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `ccode(expr, assign?)` | ⬜ | C code | `codegen/ccode.ts` |
| `cxxcode(expr, assign?)` | ⬜ | C++ code | `codegen/cxxcode.ts` |
| `jscode(expr, assign?)` | ⬜ | JavaScript code | `codegen/jscode.ts` |
| `pythoncode(expr)` | ⬜ | Python code | `codegen/pythoncode.ts` |
| `octave_code(expr)` | ⬜ | MATLAB/Octave code | `codegen/octave-code.ts` |
| `rust_code(expr)` | ⬜ | Rust code | `codegen/rust-code.ts` |
| `julia_code(expr)` | ⬜ | Julia code | `codegen/julia-code.ts` |

### 10.4 Parsing
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `parse_expr(s, opts?)` | ⬜ | Parse string to expression | `parsing/parse-expr.ts` |
| `sympify(s)` | ⬜ | Convert to symbolic | `parsing/sympify.ts` |
| `S(s)` | ⬜ | Sympify alias | `parsing/sympify.ts` |

---

## 11. Lambda & Numerical

### 11.1 Lambdify
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `lambdify(args, expr, modules?)` | ⬜ | Convert to callable | `lambdify/lambdify.ts` |
| `lambdify([x,y], expr, 'math')` | ⬜ | Use JS Math module | `lambdify/lambdify.ts` |

### 11.2 Numerical Utilities
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `nsimplify(x, rational?)` | ⬜ | Numerical to symbolic | `numerical/nsimplify.ts` |
| `nsimplify(x, [pi, E])` | ⬜ | Recognize constants | `numerical/nsimplify.ts` |
| `Float(x, dps)` | ⬜ | Set decimal places | `numerical/float-dps.ts` |

### 11.3 CSE Optimization
| Function | Status | Description | File |
|----------|--------|-------------|------|
| `cse(exprs, symbols?)` | ⬜ | Common subexpression elimination | `numerical/cse.ts` |

---

## 12. Assumptions & Queries

| Function | Status | Description | File |
|----------|--------|-------------|------|
| `Symbol(name, {positive: true})` | ⬜ | Symbol with assumption | `assumptions/symbol.ts` |
| `ask(query, assumptions)` | ⬜ | Query system | `assumptions/ask.ts` |
| `refine(expr, assumptions)` | ⬜ | Simplify under assumptions | `assumptions/refine.ts` |
| `Q.positive(x)` | ⬜ | Positive query | `assumptions/queries/positive.ts` |
| `Q.negative(x)` | ⬜ | Negative query | `assumptions/queries/negative.ts` |
| `Q.real(x)` | ⬜ | Real query | `assumptions/queries/real.ts` |
| `Q.integer(x)` | ⬜ | Integer query | `assumptions/queries/integer.ts` |
| `Q.even(x)` | ⬜ | Even query | `assumptions/queries/even.ts` |
| `Q.odd(x)` | ⬜ | Odd query | `assumptions/queries/odd.ts` |
| `Q.prime(x)` | ⬜ | Prime query | `assumptions/queries/prime.ts` |
| `Q.composite(x)` | ⬜ | Composite query | `assumptions/queries/composite.ts` |
| `Q.zero(x)` | ⬜ | Zero query | `assumptions/queries/zero.ts` |
| `Q.nonzero(x)` | ⬜ | Nonzero query | `assumptions/queries/nonzero.ts` |
| `Q.finite(x)` | ⬜ | Finite query | `assumptions/queries/finite.ts` |
| `Q.infinite(x)` | ⬜ | Infinite query | `assumptions/queries/infinite.ts` |

---

## Summary Statistics

| Category | Total | ✅ Done | 🔲 Stub | ⬜ TODO |
|----------|-------|---------|---------|---------|
| Core Foundation | 45 | 35 | 0 | 10 |
| Elementary Functions | 52 | 52 | 0 | 0 |
| Calculus | 25 | 4 | 5 | 16 |
| Simplification | 14 | 14 | 0 | 0 |
| Matrices | 36 | 36 | 0 | 0 |
| Polynomials | 55 | 0 | 1 | 54 |
| Equation Solving | 25 | 0 | 5 | 20 |
| Number Theory | 55 | 0 | 0 | 55 |
| Sets & Logic | 50 | 0 | 0 | 50 |
| Printing & I/O | 25 | 0 | 4 | 21 |
| Lambda & Numerical | 10 | 0 | 0 | 10 |
| Assumptions | 20 | 0 | 0 | 20 |
| **TOTAL** | **412** | **141** | **15** | **256** |

---

## Workflow

For each function/module, follow these steps in order:

### Step 1: Identify the C++ source in SymEngine reference
- Locate the corresponding C++ header and implementation files in `/packages/symwasm/reference/symengine/symengine/`
- **Copy the C++ code as-is** — do not reimplement in TypeScript

### Step 2: Configure WASM build system
- Add the identified C++ files to the Emscripten build configuration
- Use the C API wrappers from `cwrapper.h` and `cwrapper.cpp` for JS/WASM interop

### Step 3: Compile C++ to WebAssembly
- Run Emscripten build to compile SymEngine C++ to `.wasm`
- Expose the necessary C API entry points via cwrapper
- Verify memory management (allocate, free, reference counting)

### Step 4: Write TypeScript wrappers (thin layer only)
- Create TypeScript classes/functions that call into the WASM module via cwrapper
- Follow SymPy-like API conventions for familiarity
- **Do NOT reimplement algorithms** — only create bindings

### Step 5: Port tests from SymEngine
- Find corresponding tests in `/packages/symwasm/reference/symengine/symengine/tests/`
- Translate C++ test cases to TypeScript/Vitest format

---

## SymEngine C++ File Reference

### Core Files
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

### Functions & Calculus
```
symengine/functions.h          → All elementary & special functions
symengine/derivative.h         → Differentiation
symengine/series.h             → Series expansion
symengine/series_visitor.h     → Series algorithms
symengine/expand.h             → Expression expansion
```

### Linear Algebra
```
symengine/matrix.h             → Dense matrix operations
symengine/matrices/            → Symbolic matrices (not in C API)
```

### Polynomials
```
symengine/polys/uintpoly.h     → Univariate integer polynomial
symengine/polys/uratpoly.h     → Univariate rational polynomial
symengine/polys/uexprpoly.h    → Univariate expression polynomial
```

### Solving
```
symengine/solve.h              → Equation solving
```

### Number Theory
```
symengine/ntheory.h            → Number theory functions
symengine/ntheory_funcs.h      → Prime, GCD, LCM, modular arithmetic
symengine/diophantine.h        → Diophantine equations
```

### Sets & Logic
```
symengine/sets.h               → Set theory
symengine/logic.h              → Boolean logic & relations
```

### I/O
```
symengine/printers/strprinter.h    → String printer
symengine/printers/latex.h         → LaTeX printer
symengine/printers/mathml.h        → MathML printer
symengine/printers/unicode.h       → Unicode printer
symengine/printers/codegen.h       → Code generation
symengine/parser/parser.h          → Expression parser
symengine/lambda_double.h          → Lambda compilation
```

---

## Build Dependencies

### Required
- **SymEngine C++ library** (in `/packages/symwasm/reference/symengine/`)
- **Emscripten** (for C++ → WASM compilation)
- **GMP** (GNU Multiple Precision library)

### Optional (for enhanced features)
- **MPFR** (arbitrary precision floating point)
- **MPC** (arbitrary precision complex numbers)
- **FLINT** (Fast Library for Number Theory)
- **LLVM** (for JIT compilation)
