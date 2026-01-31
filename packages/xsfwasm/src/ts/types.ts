/**
 * XSF (Special Functions) WASM module type definitions
 *
 * Provides special mathematical functions including:
 * - Gamma functions (gamma, gammaln, rgamma)
 * - Beta functions (beta, betaln)
 * - Error functions (erf, erfc, erfcx, erfi)
 * - Bessel functions (j0, j1, jv, y0, y1, yv, i0, i1, iv, k0, k1)
 * - Combinatorial functions (binomial coefficients, Pochhammer symbol, permutations)
 */

export interface XSFModule {
  // ============================================================================
  // GAMMA FUNCTIONS
  // ============================================================================

  /**
   * Computes the gamma function Γ(x).
   *
   * The gamma function is defined as:
   *   Γ(x) = ∫₀^∞ t^(x-1) * e^(-t) dt
   *
   * For positive integers: Γ(n) = (n-1)!
   *
   * Properties:
   * - Γ(1) = 1
   * - Γ(1/2) = √π
   * - Γ(x+1) = x * Γ(x)
   *
   * @param x - Input value (must not be a non-positive integer)
   * @returns Γ(x), the gamma function evaluated at x
   *
   * @note Returns NaN for non-positive integers (poles of the function)
   * @note May overflow for large positive x (x > ~171.6)
   */
  _wasm_gamma(x: number): number;

  /**
   * Computes the natural logarithm of the absolute value of the gamma function.
   *
   * gammaln(x) = ln(|Γ(x)|)
   *
   * This function is useful when Γ(x) would overflow but its logarithm is needed.
   * Commonly used in statistical distributions and combinatorics.
   *
   * @param x - Input value (must not be a non-positive integer)
   * @returns ln(|Γ(x)|), the log-gamma function
   *
   * @note More numerically stable than log(gamma(x)) for large x
   * @note Returns Inf for non-positive integers
   */
  _wasm_gammaln(x: number): number;

  /**
   * Computes the reciprocal gamma function 1/Γ(x).
   *
   * rgamma(x) = 1/Γ(x)
   *
   * This function is entire (analytic everywhere) unlike Γ(x) which has poles.
   * It equals zero at non-positive integers where Γ(x) has poles.
   *
   * @param x - Input value
   * @returns 1/Γ(x), the reciprocal of the gamma function
   *
   * @note Returns 0 for non-positive integers (zeros of rgamma)
   * @note More numerically stable than 1/gamma(x) near poles
   */
  _wasm_rgamma(x: number): number;

  // ============================================================================
  // BETA FUNCTIONS
  // ============================================================================

  /**
   * Computes the beta function B(a, b).
   *
   * The beta function is defined as:
   *   B(a, b) = ∫₀^1 t^(a-1) * (1-t)^(b-1) dt
   *
   * Equivalently:
   *   B(a, b) = Γ(a) * Γ(b) / Γ(a + b)
   *
   * Properties:
   * - B(a, b) = B(b, a) (symmetric)
   * - B(1, 1) = 1
   * - B(a, 1) = 1/a
   *
   * @param a - First parameter (must be positive)
   * @param b - Second parameter (must be positive)
   * @returns B(a, b), the beta function
   *
   * @note Used extensively in probability theory (Beta distribution)
   * @note Related to binomial coefficients: C(n,k) = 1/((n+1)*B(n-k+1, k+1))
   */
  _wasm_beta(a: number, b: number): number;

  /**
   * Computes the natural logarithm of the beta function.
   *
   * betaln(a, b) = ln(B(a, b)) = ln(Γ(a)) + ln(Γ(b)) - ln(Γ(a + b))
   *
   * This function is more numerically stable than log(beta(a, b))
   * for large values of a and/or b.
   *
   * @param a - First parameter (must be positive)
   * @param b - Second parameter (must be positive)
   * @returns ln(B(a, b)), the log-beta function
   *
   * @note Essential for computing binomial coefficients with large n
   * @note Avoids overflow issues present in direct beta computation
   */
  _wasm_betaln(a: number, b: number): number;

  // ============================================================================
  // ERROR FUNCTIONS
  // ============================================================================

  /**
   * Computes the error function erf(x).
   *
   * The error function is defined as:
   *   erf(x) = (2/√π) * ∫₀^x e^(-t²) dt
   *
   * Properties:
   * - erf(0) = 0
   * - erf(∞) = 1
   * - erf(-x) = -erf(x) (odd function)
   * - erf(1) ≈ 0.8427
   *
   * @param x - Input value
   * @returns erf(x), ranging from -1 to 1
   *
   * @note Fundamental in probability/statistics (normal distribution CDF)
   * @note Related to normal CDF: Φ(x) = (1 + erf(x/√2))/2
   */
  _wasm_erf(x: number): number;

  /**
   * Computes the complementary error function erfc(x).
   *
   * erfc(x) = 1 - erf(x) = (2/√π) * ∫_x^∞ e^(-t²) dt
   *
   * Properties:
   * - erfc(0) = 1
   * - erfc(∞) = 0
   * - erfc(-∞) = 2
   *
   * @param x - Input value
   * @returns erfc(x) = 1 - erf(x)
   *
   * @note More accurate than 1 - erf(x) for large x
   * @note Essential for computing tail probabilities
   */
  _wasm_erfc(x: number): number;

  /**
   * Computes the scaled complementary error function erfcx(x).
   *
   * erfcx(x) = e^(x²) * erfc(x)
   *
   * This scaled version avoids underflow for large positive x where
   * erfc(x) becomes extremely small.
   *
   * Properties:
   * - erfcx(0) = 1
   * - erfcx(x) → 1/(x√π) as x → ∞
   * - erfcx(x) is always positive
   *
   * @param x - Input value
   * @returns erfcx(x) = exp(x²) * erfc(x)
   *
   * @note Numerically stable for all x
   * @note Used in computing Voigt profiles and plasma physics
   */
  _wasm_erfcx(x: number): number;

  /**
   * Computes the imaginary error function erfi(x).
   *
   * erfi(x) = -i * erf(ix) = (2/√π) * ∫₀^x e^(t²) dt
   *
   * Unlike erf(x), erfi(x) grows without bound for large |x|.
   *
   * Properties:
   * - erfi(0) = 0
   * - erfi(-x) = -erfi(x) (odd function)
   * - d/dx erfi(x) = (2/√π) * e^(x²)
   *
   * @param x - Input value
   * @returns erfi(x), the imaginary error function
   *
   * @note Grows rapidly: erfi(x) ~ e^(x²)/(x√π) for large x
   * @note Appears in solutions to certain differential equations
   */
  _wasm_erfi(x: number): number;

  // ============================================================================
  // BESSEL FUNCTIONS OF THE FIRST KIND (J)
  // ============================================================================

  /**
   * Computes the Bessel function of the first kind of order 0.
   *
   * J₀(x) is the solution to Bessel's equation x²y'' + xy' + x²y = 0
   * that is finite at the origin.
   *
   * Properties:
   * - J₀(0) = 1
   * - J₀(x) is an even function
   * - Has infinitely many zeros: ~2.405, 5.520, 8.654, ...
   * - |J₀(x)| ≤ 1 for all real x
   *
   * @param x - Input value (any real number)
   * @returns J₀(x), Bessel function of order 0
   *
   * @note Appears in solutions to Laplace's equation in cylindrical coordinates
   * @note Used in signal processing (FM modulation), vibration analysis
   */
  _wasm_j0(x: number): number;

  /**
   * Computes the Bessel function of the first kind of order 1.
   *
   * J₁(x) is the solution to Bessel's equation for n=1.
   *
   * Properties:
   * - J₁(0) = 0
   * - J₁(x) is an odd function
   * - J₁(x) = -J₀'(x) (derivative relation)
   * - Has zeros at ~3.832, 7.016, 10.173, ...
   *
   * @param x - Input value (any real number)
   * @returns J₁(x), Bessel function of order 1
   *
   * @note Used in diffraction patterns (Airy disk)
   * @note Important in electromagnetic waveguides
   */
  _wasm_j1(x: number): number;

  /**
   * Computes the Bessel function of the first kind of order v.
   *
   * Jᵥ(x) generalizes J₀ and J₁ to arbitrary real order v.
   * It satisfies: x²y'' + xy' + (x² - v²)y = 0
   *
   * Properties:
   * - Jᵥ(0) = 0 for v > 0, J₀(0) = 1
   * - J₋ᵥ(x) = (-1)^v * Jᵥ(x) for integer v
   * - Recurrence: Jᵥ₊₁(x) = (2v/x)Jᵥ(x) - Jᵥ₋₁(x)
   *
   * @param v - Order of the Bessel function (any real number)
   * @param x - Argument (typically x ≥ 0 for non-integer v)
   * @returns Jᵥ(x), Bessel function of order v
   *
   * @note For non-integer v, Jᵥ(x) and J₋ᵥ(x) are linearly independent
   * @note Approaches zero as x → ∞ with oscillating behavior
   */
  _wasm_jv(v: number, x: number): number;

  // ============================================================================
  // BESSEL FUNCTIONS OF THE SECOND KIND (Y) - NEUMANN FUNCTIONS
  // ============================================================================

  /**
   * Computes the Bessel function of the second kind of order 0.
   *
   * Y₀(x), also called the Neumann function N₀(x), is the second
   * linearly independent solution to Bessel's equation for n=0.
   *
   * Properties:
   * - Y₀(x) → -∞ as x → 0⁺ (logarithmic singularity)
   * - Y₀(x) oscillates with decreasing amplitude as x → ∞
   * - Has zeros at ~0.894, 3.958, 7.086, ...
   *
   * @param x - Input value (must be positive, x > 0)
   * @returns Y₀(x), Neumann function of order 0
   *
   * @note Undefined for x ≤ 0
   * @note Forms complete solution with J₀: y = A*J₀(x) + B*Y₀(x)
   */
  _wasm_y0(x: number): number;

  /**
   * Computes the Bessel function of the second kind of order 1.
   *
   * Y₁(x) is the Neumann function of order 1.
   *
   * Properties:
   * - Y₁(x) → -∞ as x → 0⁺
   * - Y₁(x) is an odd function continuation (Y₁(-x) = -Y₁(x))
   * - Y₁(x) = -Y₀'(x) (derivative relation)
   *
   * @param x - Input value (must be positive, x > 0)
   * @returns Y₁(x), Neumann function of order 1
   *
   * @note Returns -Inf for x → 0⁺
   * @note Important in acoustic and electromagnetic wave problems
   */
  _wasm_y1(x: number): number;

  /**
   * Computes the Bessel function of the second kind of order v.
   *
   * Yᵥ(x) generalizes Y₀ and Y₁ to arbitrary real order v.
   * For integer n: Yₙ(x) = lim_{v→n} Yᵥ(x)
   *
   * Definition for non-integer v:
   *   Yᵥ(x) = [Jᵥ(x)cos(vπ) - J₋ᵥ(x)] / sin(vπ)
   *
   * Properties:
   * - Singular at x = 0 for all v
   * - Recurrence: Yᵥ₊₁(x) = (2v/x)Yᵥ(x) - Yᵥ₋₁(x)
   *
   * @param v - Order of the Bessel function (any real number)
   * @param x - Argument (must be positive, x > 0)
   * @returns Yᵥ(x), Neumann function of order v
   *
   * @note Required for complete solutions in cylindrical coordinates
   * @note Used in heat conduction, wave propagation problems
   */
  _wasm_yv(v: number, x: number): number;

  // ============================================================================
  // MODIFIED BESSEL FUNCTIONS OF THE FIRST KIND (I)
  // ============================================================================

  /**
   * Computes the modified Bessel function of the first kind of order 0.
   *
   * I₀(x) is the solution to the modified Bessel equation:
   *   x²y'' + xy' - x²y = 0
   * that is finite at the origin.
   *
   * Properties:
   * - I₀(0) = 1
   * - I₀(x) is an even function
   * - I₀(x) > 0 for all real x
   * - I₀(x) ~ e^x / √(2πx) as x → ∞
   *
   * @param x - Input value (any real number)
   * @returns I₀(x), modified Bessel function of order 0
   *
   * @note Exponentially growing for large |x|
   * @note Appears in statistics (von Mises distribution), heat conduction
   */
  _wasm_i0(x: number): number;

  /**
   * Computes the modified Bessel function of the first kind of order 1.
   *
   * I₁(x) is related to I₀ by: I₁(x) = I₀'(x)
   *
   * Properties:
   * - I₁(0) = 0
   * - I₁(x) is an odd function
   * - I₁(x) ~ e^x / √(2πx) as x → ∞
   *
   * @param x - Input value (any real number)
   * @returns I₁(x), modified Bessel function of order 1
   *
   * @note Used in cylindrical heat conduction problems
   * @note Appears in Bessel process (stochastic processes)
   */
  _wasm_i1(x: number): number;

  /**
   * Computes the modified Bessel function of the first kind of order v.
   *
   * Iᵥ(x) generalizes I₀ and I₁ to arbitrary real order v.
   * Related to Jᵥ: Iᵥ(x) = i^(-v) * Jᵥ(ix)
   *
   * Properties:
   * - Iᵥ(0) = 0 for v > 0, I₀(0) = 1
   * - I₋ᵥ(x) = Iᵥ(x) for integer v
   * - Recurrence: Iᵥ₊₁(x) = Iᵥ₋₁(x) - (2v/x)Iᵥ(x)
   *
   * @param v - Order of the Bessel function (any real number)
   * @param x - Argument (any real number for integer v; x ≥ 0 for non-integer v)
   * @returns Iᵥ(x), modified Bessel function of order v
   *
   * @note Monotonically increasing for v ≥ 0 and x > 0
   * @note Used in probability distributions on the circle
   */
  _wasm_iv(v: number, x: number): number;

  // ============================================================================
  // MODIFIED BESSEL FUNCTIONS OF THE SECOND KIND (K) - MACDONALD FUNCTIONS
  // ============================================================================

  /**
   * Computes the modified Bessel function of the second kind of order 0.
   *
   * K₀(x), also called the MacDonald function, is the second linearly
   * independent solution to the modified Bessel equation for n=0.
   *
   * Properties:
   * - K₀(x) → -ln(x/2) - γ as x → 0⁺ (logarithmic singularity)
   * - K₀(x) ~ √(π/(2x)) * e^(-x) as x → ∞
   * - K₀(x) > 0 for all x > 0
   * - K₀(x) is monotonically decreasing
   *
   * @param x - Input value (must be positive, x > 0)
   * @returns K₀(x), MacDonald function of order 0
   *
   * @note Exponentially decaying for large x
   * @note Used in potential theory, Green's functions for diffusion
   */
  _wasm_k0(x: number): number;

  /**
   * Computes the modified Bessel function of the second kind of order 1.
   *
   * K₁(x) is the MacDonald function of order 1.
   *
   * Properties:
   * - K₁(x) → 1/x as x → 0⁺
   * - K₁(x) ~ √(π/(2x)) * e^(-x) as x → ∞
   * - K₁(x) = -K₀'(x) (derivative relation)
   *
   * @param x - Input value (must be positive, x > 0)
   * @returns K₁(x), MacDonald function of order 1
   *
   * @note Returns +Inf for x → 0⁺
   * @note Important in cylindrical wave decay problems
   */
  _wasm_k1(x: number): number;

  // ============================================================================
  // COMBINATORIAL FUNCTIONS
  // ============================================================================

  /**
   * Computes the binomial coefficient C(n, k) = n! / (k! * (n-k)!).
   *
   * Also known as "n choose k", this counts the number of ways to choose
   * k items from n items without regard to order.
   *
   * Computed using:
   *   C(n, k) = Γ(n+1) / (Γ(k+1) * Γ(n-k+1))
   *
   * Properties:
   * - C(n, 0) = C(n, n) = 1
   * - C(n, k) = C(n, n-k) (symmetry)
   * - C(n, k) = C(n-1, k-1) + C(n-1, k) (Pascal's rule)
   *
   * @param n - Total number of items (can be any real number for generalized binomial)
   * @param k - Number of items to choose (can be any real number)
   * @returns C(n, k), the binomial coefficient
   *
   * @note For non-integer n, uses the gamma function generalization
   * @note May lose precision for very large n due to floating-point arithmetic
   */
  _wasm_binom(n: number, k: number): number;

  /**
   * Computes the exact binomial coefficient for integer arguments.
   *
   * Uses integer arithmetic where possible to maintain exact results
   * for small to moderate values of n and k.
   *
   * @param n - Total number of items (non-negative integer)
   * @param k - Number of items to choose (non-negative integer, k ≤ n)
   * @returns C(n, k) computed with maximum precision
   *
   * @note Returns 0 if k < 0 or k > n
   * @note May overflow for large n (n > ~170)
   * @note More accurate than binom() for integer arguments
   */
  _wasm_binom_exact(n: number, k: number): number;

  /**
   * Computes the Pochhammer symbol (rising factorial) (x)_m.
   *
   * Definition:
   *   (x)_m = x * (x+1) * (x+2) * ... * (x+m-1) = Γ(x+m) / Γ(x)
   *
   * Special cases:
   * - (x)_0 = 1
   * - (x)_1 = x
   * - (1)_m = m! (factorial)
   * - (n)_m for integer n: n!/(n-m)! if m ≤ n, else 0
   *
   * @param x - Base value (any real number)
   * @param m - Number of terms (typically a non-negative integer)
   * @returns (x)_m, the Pochhammer symbol / rising factorial
   *
   * @note Fundamental in hypergeometric functions
   * @note Also called "shifted factorial" or "rising factorial"
   */
  _wasm_poch(x: number, m: number): number;

  /**
   * Computes the exact number of permutations P(n, k) = n! / (n-k)!.
   *
   * This counts the number of ways to arrange k items from n items
   * where order matters.
   *
   * Properties:
   * - P(n, 0) = 1
   * - P(n, 1) = n
   * - P(n, n) = n!
   * - P(n, k) = n * P(n-1, k-1)
   *
   * @param n - Total number of items (non-negative integer)
   * @param k - Number of items to arrange (non-negative integer, k ≤ n)
   * @returns P(n, k) = n!/(n-k)!, the number of k-permutations of n
   *
   * @note Returns 0 if k > n or if either argument is negative
   * @note May overflow for large n
   */
  _wasm_perm_exact(n: number, k: number): number;

  // ============================================================================
  // AIRY FUNCTIONS
  // ============================================================================

  /**
   * Computes Airy functions Ai(x), Ai'(x), Bi(x), Bi'(x).
   * Results are written to the provided pointers.
   */
  _wasm_airy(x: number, ai: number, aip: number, bi: number, bip: number): void;

  // ============================================================================
  // DIGAMMA (PSI) FUNCTION
  // ============================================================================

  /**
   * Computes the digamma function ψ(x) = d/dx ln(Γ(x)).
   */
  _wasm_digamma(x: number): number;

  // ============================================================================
  // ELLIPTIC INTEGRALS
  // ============================================================================

  /** Complete elliptic integral of the first kind K(m). */
  _wasm_ellipk(m: number): number;

  /** Complete elliptic integral of the second kind E(m). */
  _wasm_ellipe(m: number): number;

  /** Incomplete elliptic integral of the first kind F(φ|m). */
  _wasm_ellipkinc(phi: number, m: number): number;

  /** Incomplete elliptic integral of the second kind E(φ|m). */
  _wasm_ellipeinc(phi: number, m: number): number;

  /** Jacobi elliptic functions sn, cn, dn, ph. */
  _wasm_ellipj(u: number, m: number, sn: number, cn: number, dn: number, ph: number): void;

  // ============================================================================
  // EXPONENTIAL INTEGRALS
  // ============================================================================

  /** Exponential integral E₁(x). */
  _wasm_exp1(x: number): number;

  /** Exponential integral Ei(x). */
  _wasm_expi(x: number): number;

  // ============================================================================
  // FRESNEL INTEGRALS
  // ============================================================================

  /** Fresnel integrals S(x) and C(x). */
  _wasm_fresnel(x: number, s: number, c: number): void;

  // ============================================================================
  // HYPERGEOMETRIC FUNCTION
  // ============================================================================

  /** Gauss hypergeometric function ₂F₁(a, b; c; x). */
  _wasm_hyp2f1(a: number, b: number, c: number, x: number): number;

  // ============================================================================
  // KELVIN FUNCTIONS
  // ============================================================================

  /** Kelvin function ber(x). */
  _wasm_ber(x: number): number;

  /** Kelvin function bei(x). */
  _wasm_bei(x: number): number;

  /** Kelvin function ker(x). */
  _wasm_ker(x: number): number;

  /** Kelvin function kei(x). */
  _wasm_kei(x: number): number;

  /** Derivative of Kelvin function ber'(x). */
  _wasm_berp(x: number): number;

  /** Derivative of Kelvin function bei'(x). */
  _wasm_beip(x: number): number;

  /** Derivative of Kelvin function ker'(x). */
  _wasm_kerp(x: number): number;

  /** Derivative of Kelvin function kei'(x). */
  _wasm_keip(x: number): number;

  // ============================================================================
  // LAMBERT W FUNCTION
  // ============================================================================

  /** Lambert W function W(x, k) for branch k. */
  _wasm_lambertw(x: number, k: number): number;

  // ============================================================================
  // SINE AND COSINE INTEGRALS
  // ============================================================================

  /** Sine and cosine integrals Si(x) and Ci(x). */
  _wasm_sici(x: number, si: number, ci: number): void;

  /** Hyperbolic sine and cosine integrals Shi(x) and Chi(x). */
  _wasm_shichi(x: number, shi: number, chi: number): void;

  // ============================================================================
  // SPHERICAL BESSEL FUNCTIONS
  // ============================================================================

  /** Spherical Bessel function of the first kind j_n(x). */
  _wasm_spherical_jn(n: number, x: number): number;

  /** Spherical Bessel function of the second kind y_n(x). */
  _wasm_spherical_yn(n: number, x: number): number;

  /** Modified spherical Bessel function of the first kind i_n(x). */
  _wasm_spherical_in(n: number, x: number): number;

  /** Modified spherical Bessel function of the second kind k_n(x). */
  _wasm_spherical_kn(n: number, x: number): number;

  // ============================================================================
  // STRUVE FUNCTIONS
  // ============================================================================

  /** Struve function H_v(x). */
  _wasm_struve_h(v: number, x: number): number;

  /** Modified Struve function L_v(x). */
  _wasm_struve_l(v: number, x: number): number;

  // ============================================================================
  // RIEMANN ZETA FUNCTION
  // ============================================================================

  /** Riemann zeta function ζ(x). */
  _wasm_zeta(x: number): number;

  /** Riemann zeta complement ζ(x) - 1. */
  _wasm_zetac(x: number): number;

  // ============================================================================
  // LEGENDRE POLYNOMIALS
  // ============================================================================

  /** Legendre polynomial P_n(x). */
  _wasm_legendre_p(n: number, x: number): number;

  // ============================================================================
  // STATISTICAL DISTRIBUTIONS
  // ============================================================================

  /** Normal distribution CDF. */
  _wasm_ndtr(x: number): number;

  /** Normal distribution inverse CDF (quantile function). */
  _wasm_ndtri(p: number): number;

  /** Log of normal distribution CDF. */
  _wasm_log_ndtr(x: number): number;

  /** Chi-square distribution CDF. */
  _wasm_chdtr(df: number, x: number): number;

  /** Chi-square distribution complement (survival function). */
  _wasm_chdtrc(df: number, x: number): number;

  /** Chi-square distribution inverse CDF. */
  _wasm_chdtri(df: number, p: number): number;

  /** F-distribution CDF. */
  _wasm_fdtr(dfn: number, dfd: number, x: number): number;

  /** F-distribution complement (survival function). */
  _wasm_fdtrc(dfn: number, dfd: number, x: number): number;

  /** F-distribution inverse CDF. */
  _wasm_fdtri(dfn: number, dfd: number, p: number): number;

  /** Gamma distribution CDF. */
  _wasm_gdtr(a: number, b: number, x: number): number;

  /** Gamma distribution complement (survival function). */
  _wasm_gdtrc(a: number, b: number, x: number): number;

  /** Poisson distribution CDF. */
  _wasm_pdtr(k: number, m: number): number;

  /** Poisson distribution complement. */
  _wasm_pdtrc(k: number, m: number): number;

  /** Poisson distribution inverse. */
  _wasm_pdtri(k: number, p: number): number;

  /** Binomial distribution CDF. */
  _wasm_bdtr(k: number, n: number, p: number): number;

  /** Binomial distribution complement. */
  _wasm_bdtrc(k: number, n: number, p: number): number;

  /** Binomial distribution inverse. */
  _wasm_bdtri(k: number, n: number, y: number): number;

  /** Negative binomial distribution CDF. */
  _wasm_nbdtr(k: number, n: number, p: number): number;

  /** Negative binomial distribution complement. */
  _wasm_nbdtrc(k: number, n: number, p: number): number;

  /** Negative binomial distribution inverse. */
  _wasm_nbdtri(k: number, n: number, p: number): number;

  /** Kolmogorov distribution CDF. */
  _wasm_kolmogorov(x: number): number;

  /** Kolmogorov distribution inverse CDF. */
  _wasm_kolmogi(p: number): number;

  /** Smirnov distribution CDF. */
  _wasm_smirnov(n: number, x: number): number;

  /** Smirnov distribution inverse CDF. */
  _wasm_smirnovi(n: number, p: number): number;

  /** Owen's T function. */
  _wasm_owens_t(h: number, a: number): number;

  // ============================================================================
  // MEMORY MANAGEMENT
  // ============================================================================

  /**
   * Allocates a block of memory in the WASM heap.
   *
   * The allocated memory is uninitialized (contains garbage values).
   * Must be freed with _free() when no longer needed to prevent memory leaks.
   *
   * @param size - Number of bytes to allocate
   * @returns Pointer (byte offset) to the allocated memory block, or 0 on failure
   *
   * @example
   * // Allocate space for 10 doubles
   * const ptr = module._malloc(10 * 8);
   * // ... use memory ...
   * module._free(ptr);
   */
  _malloc(size: number): number;

  /**
   * Frees a previously allocated block of memory.
   *
   * @param ptr - Pointer returned by _malloc()
   *
   * @warning Freeing the same pointer twice causes undefined behavior
   * @warning Freeing an invalid pointer causes undefined behavior
   */
  _free(ptr: number): void;

  // ============================================================================
  // EMSCRIPTEN RUNTIME METHODS
  // ============================================================================

  /**
   * Reads a value from the WASM memory at the specified pointer.
   *
   * @param ptr - Memory address (byte offset) to read from
   * @param type - Value type: "i8", "i16", "i32", "i64", "float", "double", or "*"
   * @returns The value at the specified address
   *
   * @example
   * const doubleValue = module.getValue(ptr, 'double');
   * const intValue = module.getValue(ptr, 'i32');
   */
  getValue(ptr: number, type: string): number;

  /**
   * Writes a value to the WASM memory at the specified pointer.
   *
   * @param ptr - Memory address (byte offset) to write to
   * @param value - The value to write
   * @param type - Value type: "i8", "i16", "i32", "i64", "float", "double", or "*"
   *
   * @example
   * module.setValue(ptr, 3.14159, 'double');
   * module.setValue(ptr, 42, 'i32');
   */
  setValue(ptr: number, value: number, type: string): void;

  // ============================================================================
  // HEAP VIEWS - Direct Access to WASM Memory
  // ============================================================================

  /**
   * Float64 (double precision) view of the WASM heap.
   * Each element is 8 bytes. Index = byteOffset / 8.
   *
   * @example
   * const ptr = module._malloc(80); // 10 doubles
   * module.HEAPF64[ptr / 8] = 1.5;
   * module.HEAPF64[ptr / 8 + 1] = 2.5;
   */
  HEAPF64: Float64Array;

  /**
   * Float32 (single precision) view of the WASM heap.
   * Each element is 4 bytes. Index = byteOffset / 4.
   */
  HEAPF32: Float32Array;

  /**
   * Int32 (signed 32-bit integer) view of the WASM heap.
   * Each element is 4 bytes. Index = byteOffset / 4.
   */
  HEAP32: Int32Array;

  /**
   * Int8 (signed 8-bit integer / byte) view of the WASM heap.
   * Index equals byte offset directly.
   */
  HEAP8: Int8Array;

  /**
   * Uint8 (unsigned 8-bit integer / byte) view of the WASM heap.
   * Index equals byte offset directly.
   */
  HEAPU8: Uint8Array;

  /**
   * Uint32 (unsigned 32-bit integer) view of the WASM heap.
   * Each element is 4 bytes. Index = byteOffset / 4.
   */
  HEAPU32: Uint32Array;
}

/**
 * Options for configuring the XSF WASM module loader.
 */
export interface XSFModuleOptions {
  /**
   * Custom function to locate WASM files.
   *
   * Override the default file location logic for loading the .wasm binary.
   * Useful when deploying to CDN or custom server configurations.
   *
   * @param path - The filename being requested (e.g., "xsf.wasm")
   * @param scriptDirectory - The directory of the JS loader script
   * @returns The full URL or path to the requested file
   *
   * @example
   * const module = await loadXSFModule({
   *   locateFile: (path) => `/wasm/${path}`
   * });
   */
  locateFile?: (path: string, scriptDirectory: string) => string;
}

/**
 * Factory function type for creating XSF WASM module instances.
 *
 * @param options - Optional configuration for the module loader
 * @returns Promise resolving to the initialized XSFModule
 *
 * @example
 * import createXSFModule from '@aspect/xsfwasm';
 * const xsf = await createXSFModule();
 * const result = xsf._wasm_gamma(5); // Returns 24 (4!)
 */
export type XSFModuleFactory = (
  options?: XSFModuleOptions
) => Promise<XSFModule>;
