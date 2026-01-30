/**
 * QUADPACK WASM Module Type Definitions
 *
 * QUADPACK is a Fortran library for numerical integration (quadrature),
 * developed at the Applied Mathematics Division of Argonne National Laboratory.
 * It is widely considered the gold standard for adaptive numerical integration.
 *
 * This module provides a comprehensive set of quadrature routines:
 *
 * **Adaptive Quadrature (Finite Intervals):**
 * - DQAGS/QAGS: General-purpose adaptive integration with extrapolation
 * - DQAG/QAG: Adaptive with user-selected Gauss-Kronrod rule
 * - DQAGP/QAGP: Adaptive with user-specified break points
 *
 * **Adaptive Quadrature (Infinite Intervals):**
 * - DQAGI/QAGI: Semi-infinite or doubly-infinite intervals
 *
 * **Non-Adaptive Quadrature:**
 * - DQNG/QNG: Simple non-adaptive Gauss-Kronrod
 *
 * **Special Integrands:**
 * - DQAWO/QAWO: Oscillatory integrands (cos/sin weights)
 * - DQAWF/QAWF: Fourier integrals over semi-infinite intervals
 * - DQAWS/QAWS: Algebraic-logarithmic endpoint singularities
 * - DQAWC/QAWC: Cauchy principal value integrals
 *
 * **Fixed-Order Rules:**
 * - DQK15, DQK21, DQK31, DQK41, DQK51, DQK61: Gauss-Kronrod rules
 *
 * Function naming convention (f2c):
 * - All functions have trailing underscore (Fortran calling convention)
 * - 'd' prefix = double precision, no prefix = single precision
 * - 'e' suffix = extended interface with subdivision arrays
 *
 * Parameter convention:
 * - All parameters are pointers (Fortran pass-by-reference)
 * - Function pointers use addFunction(fn, 'dd') for double→double callbacks
 *
 * Typical usage:
 * 1. Create JavaScript integrand function: f(x) → result
 * 2. Register with addFunction(f, 'dd') to get function pointer
 * 3. Allocate workspace arrays with _malloc
 * 4. Call quadrature routine
 * 5. Read result and error estimate
 * 6. Clean up with removeFunction and _free
 *
 * @see http://www.netlib.org/quadpack/
 */

export interface QUADPACKModule {
  // ============================================================
  // Double Precision - Adaptive Quadrature (Finite Interval)
  // ============================================================

  /**
   * DQAGS - Adaptive quadrature with singularities (simplified interface)
   * Computes integral of f(x) over [a,b] using adaptive subdivision
   * with extrapolation (epsilon algorithm).
   */
  _dqags_(
    f: number,        // D_fp function pointer
    a: number,        // doublereal* lower limit
    b: number,        // doublereal* upper limit
    epsabs: number,   // doublereal* absolute tolerance
    epsrel: number,   // doublereal* relative tolerance
    result: number,   // doublereal* computed result
    abserr: number,   // doublereal* error estimate
    neval: number,    // integer* number of function evaluations
    ier: number,      // integer* error code (0=success)
    limit: number,    // integer* max subdivisions
    lenw: number,     // integer* workspace size (4*limit)
    last: number,     // integer* actual subdivisions used
    iwork: number,    // integer* workspace (limit)
    work: number      // doublereal* workspace (lenw)
  ): void;

  /**
   * DQAGSE - Adaptive quadrature with singularities (extended interface)
   * Same as DQAGS but provides access to subdivision arrays.
   */
  _dqagse_(
    f: number,        // D_fp function pointer
    a: number,        // doublereal* lower limit
    b: number,        // doublereal* upper limit
    epsabs: number,   // doublereal* absolute tolerance
    epsrel: number,   // doublereal* relative tolerance
    limit: number,    // integer* max subdivisions
    result: number,   // doublereal* computed result
    abserr: number,   // doublereal* error estimate
    neval: number,    // integer* function evaluations
    ier: number,      // integer* error code
    alist: number,    // doublereal* left endpoints (limit)
    blist: number,    // doublereal* right endpoints (limit)
    rlist: number,    // doublereal* integral approximations (limit)
    elist: number,    // doublereal* error estimates (limit)
    iord: number,     // integer* subdivision order (limit)
    last: number      // integer* subdivisions used
  ): void;

  /**
   * DQAG - Adaptive quadrature with user-specified rule (simplified)
   */
  _dqag_(
    f: number,        // D_fp function pointer
    a: number,        // doublereal* lower limit
    b: number,        // doublereal* upper limit
    epsabs: number,   // doublereal* absolute tolerance
    epsrel: number,   // doublereal* relative tolerance
    key: number,      // integer* integration rule (1-6 for 15,21,31,41,51,61 point)
    result: number,   // doublereal* computed result
    abserr: number,   // doublereal* error estimate
    neval: number,    // integer* function evaluations
    ier: number,      // integer* error code
    limit: number,    // integer* max subdivisions
    lenw: number,     // integer* workspace size
    last: number,     // integer* subdivisions used
    iwork: number,    // integer* workspace
    work: number      // doublereal* workspace
  ): void;

  /**
   * DQAGE - Adaptive quadrature with user-specified rule (extended)
   */
  _dqage_(
    f: number,        // D_fp function pointer
    a: number,        // doublereal* lower limit
    b: number,        // doublereal* upper limit
    epsabs: number,   // doublereal* absolute tolerance
    epsrel: number,   // doublereal* relative tolerance
    key: number,      // integer* integration rule
    limit: number,    // integer* max subdivisions
    result: number,   // doublereal* computed result
    abserr: number,   // doublereal* error estimate
    neval: number,    // integer* function evaluations
    ier: number,      // integer* error code
    alist: number,    // doublereal* left endpoints
    blist: number,    // doublereal* right endpoints
    rlist: number,    // doublereal* integral approximations
    elist: number,    // doublereal* error estimates
    iord: number,     // integer* subdivision order
    last: number      // integer* subdivisions used
  ): void;

  // ============================================================
  // Double Precision - Adaptive Quadrature (Infinite Interval)
  // ============================================================

  /**
   * DQAGI - Adaptive quadrature for infinite intervals (simplified)
   * inf=1: (bound,+inf), inf=-1: (-inf,bound), inf=2: (-inf,+inf)
   */
  _dqagi_(
    f: number,        // D_fp function pointer
    bound: number,    // doublereal* finite bound
    inf: number,      // integer* interval type (1,-1,2)
    epsabs: number,   // doublereal* absolute tolerance
    epsrel: number,   // doublereal* relative tolerance
    result: number,   // doublereal* computed result
    abserr: number,   // doublereal* error estimate
    neval: number,    // integer* function evaluations
    ier: number,      // integer* error code
    limit: number,    // integer* max subdivisions
    lenw: number,     // integer* workspace size
    last: number,     // integer* subdivisions used
    iwork: number,    // integer* workspace
    work: number      // doublereal* workspace
  ): void;

  /**
   * DQAGIE - Adaptive quadrature for infinite intervals (extended)
   */
  _dqagie_(
    f: number,        // D_fp function pointer
    bound: number,    // doublereal* finite bound
    inf: number,      // integer* interval type
    epsabs: number,   // doublereal* absolute tolerance
    epsrel: number,   // doublereal* relative tolerance
    limit: number,    // integer* max subdivisions
    result: number,   // doublereal* computed result
    abserr: number,   // doublereal* error estimate
    neval: number,    // integer* function evaluations
    ier: number,      // integer* error code
    alist: number,    // doublereal* left endpoints
    blist: number,    // doublereal* right endpoints
    rlist: number,    // doublereal* integral approximations
    elist: number,    // doublereal* error estimates
    iord: number,     // integer* subdivision order
    last: number      // integer* subdivisions used
  ): void;

  // ============================================================
  // Double Precision - Non-Adaptive Quadrature
  // ============================================================

  /**
   * DQNG - Non-adaptive Gauss-Kronrod integration
   * Uses 10, 21, 43, and 87 point rules progressively.
   */
  _dqng_(
    f: number,        // D_fp function pointer
    a: number,        // doublereal* lower limit
    b: number,        // doublereal* upper limit
    epsabs: number,   // doublereal* absolute tolerance
    epsrel: number,   // doublereal* relative tolerance
    result: number,   // doublereal* computed result
    abserr: number,   // doublereal* error estimate
    neval: number,    // integer* function evaluations
    ier: number       // integer* error code
  ): void;

  // ============================================================
  // Double Precision - Oscillatory Integrals
  // ============================================================

  /**
   * DQAWO - Oscillatory integrals (simplified)
   * Computes integral of f(x)*w(x) where w(x)=cos(omega*x) or sin(omega*x)
   */
  _dqawo_(
    f: number,        // D_fp function pointer
    a: number,        // doublereal* lower limit
    b: number,        // doublereal* upper limit
    omega: number,    // doublereal* frequency parameter
    integr: number,   // integer* 1=cos, 2=sin
    epsabs: number,   // doublereal* absolute tolerance
    epsrel: number,   // doublereal* relative tolerance
    result: number,   // doublereal* computed result
    abserr: number,   // doublereal* error estimate
    neval: number,    // integer* function evaluations
    ier: number,      // integer* error code
    leniw: number,    // integer* integer workspace size
    maxp1: number,    // integer* max Chebyshev moments
    lenw: number,     // integer* workspace size
    last: number,     // integer* subdivisions used
    iwork: number,    // integer* workspace
    work: number      // doublereal* workspace
  ): void;

  /**
   * DQAWOE - Oscillatory integrals (extended)
   */
  _dqawoe_(
    f: number,        // D_fp function pointer
    a: number,        // doublereal* lower limit
    b: number,        // doublereal* upper limit
    omega: number,    // doublereal* frequency
    integr: number,   // integer* 1=cos, 2=sin
    epsabs: number,   // doublereal* absolute tolerance
    epsrel: number,   // doublereal* relative tolerance
    limit: number,    // integer* max subdivisions
    icall: number,    // integer* call indicator
    maxp1: number,    // integer* max Chebyshev moments
    result: number,   // doublereal* computed result
    abserr: number,   // doublereal* error estimate
    neval: number,    // integer* function evaluations
    ier: number,      // integer* error code
    last: number,     // integer* subdivisions used
    alist: number,    // doublereal* left endpoints
    blist: number,    // doublereal* right endpoints
    rlist: number,    // doublereal* integral approximations
    elist: number,    // doublereal* error estimates
    iord: number,     // integer* subdivision order
    nnlog: number,    // integer* log info
    momcom: number,   // integer* moment counter
    chebmo: number    // doublereal* Chebyshev moments
  ): void;

  /**
   * DQAWF - Fourier integrals (simplified)
   * Computes integral of f(x)*w(x) over [a,+inf)
   */
  _dqawf_(
    f: number,        // D_fp function pointer
    a: number,        // doublereal* lower limit
    omega: number,    // doublereal* frequency
    integr: number,   // integer* 1=cos, 2=sin
    epsabs: number,   // doublereal* absolute tolerance
    result: number,   // doublereal* computed result
    abserr: number,   // doublereal* error estimate
    neval: number,    // integer* function evaluations
    ier: number,      // integer* error code
    limlst: number,   // integer* cycle limit
    lst: number,      // integer* cycles used
    leniw: number,    // integer* integer workspace size
    maxp1: number,    // integer* max Chebyshev moments
    lenw: number,     // integer* workspace size
    iwork: number,    // integer* workspace
    work: number      // doublereal* workspace
  ): void;

  /**
   * DQAWFE - Fourier integrals (extended)
   */
  _dqawfe_(
    f: number,        // D_fp function pointer
    a: number,        // doublereal* lower limit
    omega: number,    // doublereal* frequency
    integr: number,   // integer* 1=cos, 2=sin
    epsabs: number,   // doublereal* absolute tolerance
    limlst: number,   // integer* cycle limit
    limit: number,    // integer* max subdivisions per cycle
    maxp1: number,    // integer* max Chebyshev moments
    result: number,   // doublereal* computed result
    abserr: number,   // doublereal* error estimate
    neval: number,    // integer* function evaluations
    ier: number,      // integer* error code
    rslst: number,    // doublereal* cycle results
    erlst: number,    // doublereal* cycle errors
    ierlst: number,   // integer* cycle error codes
    lst: number,      // integer* cycles used
    alist: number,    // doublereal* left endpoints
    blist: number,    // doublereal* right endpoints
    rlist: number,    // doublereal* integral approximations
    elist: number,    // doublereal* error estimates
    iord: number,     // integer* subdivision order
    nnlog: number,    // integer* log info
    chebmo: number    // doublereal* Chebyshev moments
  ): void;

  // ============================================================
  // Double Precision - Singular Weight Functions
  // ============================================================

  /**
   * DQAWS - Algebraic-logarithmic singularities (simplified)
   * Weight: (x-a)^alfa * (b-x)^beta * w(x)
   * integr: 1=1, 2=log(x-a), 3=log(b-x), 4=log(x-a)*log(b-x)
   */
  _dqaws_(
    f: number,        // D_fp function pointer
    a: number,        // doublereal* lower limit
    b: number,        // doublereal* upper limit
    alfa: number,     // doublereal* exponent at a (>-1)
    beta: number,     // doublereal* exponent at b (>-1)
    integr: number,   // integer* weight type (1-4)
    epsabs: number,   // doublereal* absolute tolerance
    epsrel: number,   // doublereal* relative tolerance
    result: number,   // doublereal* computed result
    abserr: number,   // doublereal* error estimate
    neval: number,    // integer* function evaluations
    ier: number,      // integer* error code
    limit: number,    // integer* max subdivisions
    lenw: number,     // integer* workspace size
    last: number,     // integer* subdivisions used
    iwork: number,    // integer* workspace
    work: number      // doublereal* workspace
  ): void;

  /**
   * DQAWSE - Algebraic-logarithmic singularities (extended)
   */
  _dqawse_(
    f: number,        // D_fp function pointer
    a: number,        // doublereal* lower limit
    b: number,        // doublereal* upper limit
    alfa: number,     // doublereal* exponent at a
    beta: number,     // doublereal* exponent at b
    integr: number,   // integer* weight type
    epsabs: number,   // doublereal* absolute tolerance
    epsrel: number,   // doublereal* relative tolerance
    limit: number,    // integer* max subdivisions
    result: number,   // doublereal* computed result
    abserr: number,   // doublereal* error estimate
    neval: number,    // integer* function evaluations
    ier: number,      // integer* error code
    alist: number,    // doublereal* left endpoints
    blist: number,    // doublereal* right endpoints
    rlist: number,    // doublereal* integral approximations
    elist: number,    // doublereal* error estimates
    iord: number,     // integer* subdivision order
    last: number      // integer* subdivisions used
  ): void;

  /**
   * DQAWC - Cauchy principal value (simplified)
   * Computes principal value of integral f(x)/(x-c) over [a,b]
   */
  _dqawc_(
    f: number,        // D_fp function pointer
    a: number,        // doublereal* lower limit
    b: number,        // doublereal* upper limit
    c: number,        // doublereal* singular point (a < c < b)
    epsabs: number,   // doublereal* absolute tolerance
    epsrel: number,   // doublereal* relative tolerance
    result: number,   // doublereal* computed result
    abserr: number,   // doublereal* error estimate
    neval: number,    // integer* function evaluations
    ier: number,      // integer* error code
    limit: number,    // integer* max subdivisions
    lenw: number,     // integer* workspace size
    last: number,     // integer* subdivisions used
    iwork: number,    // integer* workspace
    work: number      // doublereal* workspace
  ): void;

  /**
   * DQAWCE - Cauchy principal value (extended)
   */
  _dqawce_(
    f: number,        // D_fp function pointer
    a: number,        // doublereal* lower limit
    b: number,        // doublereal* upper limit
    c: number,        // doublereal* singular point
    epsabs: number,   // doublereal* absolute tolerance
    epsrel: number,   // doublereal* relative tolerance
    limit: number,    // integer* max subdivisions
    result: number,   // doublereal* computed result
    abserr: number,   // doublereal* error estimate
    neval: number,    // integer* function evaluations
    ier: number,      // integer* error code
    alist: number,    // doublereal* left endpoints
    blist: number,    // doublereal* right endpoints
    rlist: number,    // doublereal* integral approximations
    elist: number,    // doublereal* error estimates
    iord: number,     // integer* subdivision order
    last: number      // integer* subdivisions used
  ): void;

  // ============================================================
  // Double Precision - User-Supplied Break Points
  // ============================================================

  /**
   * DQAGP - Adaptive quadrature with break points (simplified)
   */
  _dqagp_(
    f: number,        // D_fp function pointer
    a: number,        // doublereal* lower limit
    b: number,        // doublereal* upper limit
    npts2: number,    // integer* number of break points + 2
    points: number,   // doublereal* break points array
    epsabs: number,   // doublereal* absolute tolerance
    epsrel: number,   // doublereal* relative tolerance
    result: number,   // doublereal* computed result
    abserr: number,   // doublereal* error estimate
    neval: number,    // integer* function evaluations
    ier: number,      // integer* error code
    leniw: number,    // integer* integer workspace size
    lenw: number,     // integer* workspace size
    last: number,     // integer* subdivisions used
    iwork: number,    // integer* workspace
    work: number      // doublereal* workspace
  ): void;

  /**
   * DQAGPE - Adaptive quadrature with break points (extended)
   */
  _dqagpe_(
    f: number,        // D_fp function pointer
    a: number,        // doublereal* lower limit
    b: number,        // doublereal* upper limit
    npts2: number,    // integer* number of break points + 2
    points: number,   // doublereal* break points array
    epsabs: number,   // doublereal* absolute tolerance
    epsrel: number,   // doublereal* relative tolerance
    limit: number,    // integer* max subdivisions
    result: number,   // doublereal* computed result
    abserr: number,   // doublereal* error estimate
    neval: number,    // integer* function evaluations
    ier: number,      // integer* error code
    alist: number,    // doublereal* left endpoints
    blist: number,    // doublereal* right endpoints
    rlist: number,    // doublereal* integral approximations
    elist: number,    // doublereal* error estimates
    pts: number,      // doublereal* sorted break points
    iord: number,     // integer* subdivision order
    level: number,    // integer* subdivision levels
    ndin: number,     // integer* break point info
    last: number      // integer* subdivisions used
  ): void;

  // ============================================================
  // Double Precision - Gauss-Kronrod Rules
  // ============================================================

  /** DQK15 - 15-point Gauss-Kronrod rule */
  _dqk15_(
    f: number,        // D_fp function pointer
    a: number,        // doublereal* lower limit
    b: number,        // doublereal* upper limit
    result: number,   // doublereal* integral approximation
    abserr: number,   // doublereal* error estimate
    resabs: number,   // doublereal* approximation to |f|
    resasc: number    // doublereal* approximation to |f-mean|
  ): void;

  /** DQK21 - 21-point Gauss-Kronrod rule */
  _dqk21_(
    f: number, a: number, b: number,
    result: number, abserr: number, resabs: number, resasc: number
  ): void;

  /** DQK31 - 31-point Gauss-Kronrod rule */
  _dqk31_(
    f: number, a: number, b: number,
    result: number, abserr: number, resabs: number, resasc: number
  ): void;

  /** DQK41 - 41-point Gauss-Kronrod rule */
  _dqk41_(
    f: number, a: number, b: number,
    result: number, abserr: number, resabs: number, resasc: number
  ): void;

  /** DQK51 - 51-point Gauss-Kronrod rule */
  _dqk51_(
    f: number, a: number, b: number,
    result: number, abserr: number, resabs: number, resasc: number
  ): void;

  /** DQK61 - 61-point Gauss-Kronrod rule */
  _dqk61_(
    f: number, a: number, b: number,
    result: number, abserr: number, resabs: number, resasc: number
  ): void;

  /** DQK15I - 15-point rule for infinite intervals */
  _dqk15i_(
    f: number,        // D_fp function pointer
    boun: number,     // doublereal* finite bound
    inf: number,      // integer* interval type
    a: number,        // doublereal* mapped lower limit
    b: number,        // doublereal* mapped upper limit
    result: number,   // doublereal* integral approximation
    abserr: number,   // doublereal* error estimate
    resabs: number,   // doublereal* approximation to |f|
    resasc: number    // doublereal* approximation to |f-mean|
  ): void;

  /** DQK15W - 15-point rule with weight function */
  _dqk15w_(
    f: number,        // D_fp function pointer
    w: number,        // D_fp weight function pointer
    p1: number,       // doublereal* weight parameter 1
    p2: number,       // doublereal* weight parameter 2
    p3: number,       // doublereal* weight parameter 3
    p4: number,       // doublereal* weight parameter 4
    kp: number,       // integer* weight function selector
    a: number,        // doublereal* lower limit
    b: number,        // doublereal* upper limit
    result: number,   // doublereal* integral approximation
    abserr: number,   // doublereal* error estimate
    resabs: number,   // doublereal* approximation to |f|
    resasc: number    // doublereal* approximation to |f-mean|
  ): void;

  // ============================================================
  // Single Precision - All functions (q prefix instead of dq)
  // ============================================================

  // Adaptive finite interval
  _qags_(f: number, a: number, b: number, epsabs: number, epsrel: number, result: number, abserr: number, neval: number, ier: number, limit: number, lenw: number, last: number, iwork: number, work: number): void;
  _qagse_(f: number, a: number, b: number, epsabs: number, epsrel: number, limit: number, result: number, abserr: number, neval: number, ier: number, alist: number, blist: number, rlist: number, elist: number, iord: number, last: number): void;
  _qag_(f: number, a: number, b: number, epsabs: number, epsrel: number, key: number, result: number, abserr: number, neval: number, ier: number, limit: number, lenw: number, last: number, iwork: number, work: number): void;
  _qage_(f: number, a: number, b: number, epsabs: number, epsrel: number, key: number, limit: number, result: number, abserr: number, neval: number, ier: number, alist: number, blist: number, rlist: number, elist: number, iord: number, last: number): void;

  // Adaptive infinite interval
  _qagi_(f: number, bound: number, inf: number, epsabs: number, epsrel: number, result: number, abserr: number, neval: number, ier: number, limit: number, lenw: number, last: number, iwork: number, work: number): void;
  _qagie_(f: number, bound: number, inf: number, epsabs: number, epsrel: number, limit: number, result: number, abserr: number, neval: number, ier: number, alist: number, blist: number, rlist: number, elist: number, iord: number, last: number): void;

  // Non-adaptive
  _qng_(f: number, a: number, b: number, epsabs: number, epsrel: number, result: number, abserr: number, neval: number, ier: number): void;

  // Oscillatory
  _qawo_(f: number, a: number, b: number, omega: number, integr: number, epsabs: number, epsrel: number, result: number, abserr: number, neval: number, ier: number, leniw: number, maxp1: number, lenw: number, last: number, iwork: number, work: number): void;
  _qawoe_(f: number, a: number, b: number, omega: number, integr: number, epsabs: number, epsrel: number, limit: number, icall: number, maxp1: number, result: number, abserr: number, neval: number, ier: number, last: number, alist: number, blist: number, rlist: number, elist: number, iord: number, nnlog: number, momcom: number, chebmo: number): void;
  _qawf_(f: number, a: number, omega: number, integr: number, epsabs: number, result: number, abserr: number, neval: number, ier: number, limlst: number, lst: number, leniw: number, maxp1: number, lenw: number, iwork: number, work: number): void;
  _qawfe_(f: number, a: number, omega: number, integr: number, epsabs: number, limlst: number, limit: number, maxp1: number, result: number, abserr: number, neval: number, ier: number, rslst: number, erlst: number, ierlst: number, lst: number, alist: number, blist: number, rlist: number, elist: number, iord: number, nnlog: number, chebmo: number): void;

  // Singular weights
  _qaws_(f: number, a: number, b: number, alfa: number, beta: number, integr: number, epsabs: number, epsrel: number, result: number, abserr: number, neval: number, ier: number, limit: number, lenw: number, last: number, iwork: number, work: number): void;
  _qawse_(f: number, a: number, b: number, alfa: number, beta: number, integr: number, epsabs: number, epsrel: number, limit: number, result: number, abserr: number, neval: number, ier: number, alist: number, blist: number, rlist: number, elist: number, iord: number, last: number): void;
  _qawc_(f: number, a: number, b: number, c: number, epsabs: number, epsrel: number, result: number, abserr: number, neval: number, ier: number, limit: number, lenw: number, last: number, iwork: number, work: number): void;
  _qawce_(f: number, a: number, b: number, c: number, epsabs: number, epsrel: number, limit: number, result: number, abserr: number, neval: number, ier: number, alist: number, blist: number, rlist: number, elist: number, iord: number, last: number): void;

  // Break points
  _qagp_(f: number, a: number, b: number, npts2: number, points: number, epsabs: number, epsrel: number, result: number, abserr: number, neval: number, ier: number, leniw: number, lenw: number, last: number, iwork: number, work: number): void;
  _qagpe_(f: number, a: number, b: number, npts2: number, points: number, epsabs: number, epsrel: number, limit: number, result: number, abserr: number, neval: number, ier: number, alist: number, blist: number, rlist: number, elist: number, pts: number, iord: number, level: number, ndin: number, last: number): void;

  // Gauss-Kronrod rules
  _qk15_(f: number, a: number, b: number, result: number, abserr: number, resabs: number, resasc: number): void;
  _qk21_(f: number, a: number, b: number, result: number, abserr: number, resabs: number, resasc: number): void;
  _qk31_(f: number, a: number, b: number, result: number, abserr: number, resabs: number, resasc: number): void;
  _qk41_(f: number, a: number, b: number, result: number, abserr: number, resabs: number, resasc: number): void;
  _qk51_(f: number, a: number, b: number, result: number, abserr: number, resabs: number, resasc: number): void;
  _qk61_(f: number, a: number, b: number, result: number, abserr: number, resabs: number, resasc: number): void;
  _qk15i_(f: number, boun: number, inf: number, a: number, b: number, result: number, abserr: number, resabs: number, resasc: number): void;
  _qk15w_(f: number, w: number, p1: number, p2: number, p3: number, p4: number, kp: number, a: number, b: number, result: number, abserr: number, resabs: number, resasc: number): void;

  // ============================================================
  // Memory Management
  // ============================================================

  _malloc(size: number): number;
  _free(ptr: number): void;

  // ============================================================
  // Emscripten Runtime Methods
  // ============================================================

  getValue(ptr: number, type: string): number;
  setValue(ptr: number, value: number, type: string): void;

  // ============================================================
  // Heap Views
  // ============================================================

  HEAPF64: Float64Array;
  HEAPF32: Float32Array;
  HEAP32: Int32Array;
  HEAP8: Int8Array;
  HEAPU8: Uint8Array;
  HEAPU32: Uint32Array;

  // ============================================================
  // Function Table Management (for JS→C callbacks)
  // ============================================================

  /**
   * Register a JS function for use as a callback from WASM.
   * For QUADPACK integrand functions, use signature 'dd' (double->double).
   */
  addFunction(func: Function, signature: string): number;

  /**
   * Remove a previously registered function.
   */
  removeFunction(ptr: number): void;
}

/**
 * Options for configuring the QUADPACK WASM module loader.
 */
export interface QUADPACKModuleOptions {
  /**
   * Custom function to locate WASM files.
   * @param path - Filename being requested (e.g., 'quadpack.wasm')
   * @param scriptDirectory - Directory of the JS loader script
   * @returns Full URL or path to the file
   */
  locateFile?: (path: string, scriptDirectory: string) => string;
}

/**
 * Factory function to create a QUADPACK WASM module instance.
 *
 * @param options - Optional configuration
 * @returns Promise resolving to initialized QUADPACK module
 *
 * @example
 * import createQUADPACKModule from 'quadwasm';
 *
 * const quad = await createQUADPACKModule();
 *
 * // Integrate f(x) = x^2 from 0 to 1
 * const f = (x: number) => x * x;
 * const fPtr = quad.addFunction(f, 'dd');
 *
 * // Allocate parameters and workspace
 * const limit = 50;
 * const lenw = 4 * limit;
 * const aPtr = quad._malloc(8);
 * const bPtr = quad._malloc(8);
 * // ... allocate other pointers ...
 *
 * quad.HEAPF64[aPtr / 8] = 0.0;  // lower limit
 * quad.HEAPF64[bPtr / 8] = 1.0;  // upper limit
 *
 * quad._dqags_(fPtr, aPtr, bPtr, epsabsPtr, epsrelPtr,
 *              resultPtr, abserrPtr, nevalPtr, ierPtr,
 *              limitPtr, lenwPtr, lastPtr, iworkPtr, workPtr);
 *
 * const result = quad.HEAPF64[resultPtr / 8];  // Should be ~0.333...
 *
 * quad.removeFunction(fPtr);
 * // ... free memory ...
 */
export type QUADPACKModuleFactory = (
  options?: QUADPACKModuleOptions
) => Promise<QUADPACKModule>;

// ============================================================
// QUADPACK ERROR CODES AND CONSTANTS
// ============================================================
// Error codes are returned in the 'ier' parameter.
// ier=0 indicates success; ier>0 indicates various issues.
// ============================================================

/**
 * QUADPACK error code descriptions.
 * Maps ier return codes to human-readable descriptions.
 */
export const QUADPACK_ERRORS: Record<number, string> = {
  0: 'Normal and reliable termination. Requested accuracy achieved.',
  1: 'Maximum number of subdivisions reached. May need to increase limit.',
  2: 'Roundoff error detected. Integral is probably divergent or slowly convergent.',
  3: 'Extremely bad integrand behavior at some points.',
  4: 'Algorithm does not converge. Roundoff error in extrapolation table.',
  5: 'Integral is probably divergent or slowly convergent.',
  6: 'Invalid input parameters.',
};

/**
 * QUADPACK integration rules for DQAG/QAG.
 *
 * Higher-order rules are more accurate for smooth functions but
 * require more function evaluations per subinterval.
 *
 * - GK15 (15-point): Good for most functions, low overhead
 * - GK21 (21-point): Slightly more accurate
 * - GK31 (31-point): Good balance of accuracy and efficiency
 * - GK41 (41-point): For smoother functions
 * - GK51 (51-point): Higher accuracy
 * - GK61 (61-point): Highest accuracy, most function evaluations
 */
export const enum QUADPACKRule {
  /** 15-point Gauss-Kronrod rule (7 Gauss + 8 Kronrod points) */
  GK15 = 1,
  /** 21-point Gauss-Kronrod rule (10 Gauss + 11 Kronrod points) */
  GK21 = 2,
  /** 31-point Gauss-Kronrod rule (15 Gauss + 16 Kronrod points) */
  GK31 = 3,
  /** 41-point Gauss-Kronrod rule (20 Gauss + 21 Kronrod points) */
  GK41 = 4,
  /** 51-point Gauss-Kronrod rule (25 Gauss + 26 Kronrod points) */
  GK51 = 5,
  /** 61-point Gauss-Kronrod rule (30 Gauss + 31 Kronrod points) */
  GK61 = 6,
}

/**
 * Infinite interval types for DQAGI/QAGI.
 *
 * Specifies which end(s) of the interval are infinite.
 * The 'bound' parameter specifies the finite endpoint.
 */
export const enum QUADPACKInf {
  /** Upper semi-infinite: integrate from bound to +∞ */
  UPPER = 1,
  /** Lower semi-infinite: integrate from -∞ to bound */
  LOWER = -1,
  /** Doubly infinite: integrate from -∞ to +∞ (bound is ignored) */
  BOTH = 2,
}

/**
 * Weight function types for DQAWS/QAWS.
 *
 * QAWS integrates f(x) * (x-a)^alfa * (b-x)^beta * w(x)
 * where w(x) is one of these weight functions.
 *
 * Used for integrands with algebraic-logarithmic endpoint singularities.
 */
export const enum QUADPACKWeight {
  /** w(x) = 1 (no logarithmic factor) */
  NONE = 1,
  /** w(x) = log(x-a) (logarithmic singularity at left endpoint) */
  LOG_LEFT = 2,
  /** w(x) = log(b-x) (logarithmic singularity at right endpoint) */
  LOG_RIGHT = 3,
  /** w(x) = log(x-a)*log(b-x) (logarithmic singularities at both endpoints) */
  LOG_BOTH = 4,
}

/**
 * Oscillatory weight types for DQAWO/QAWO and DQAWF/QAWF.
 *
 * QAWO integrates f(x) * w(x) where w(x) is cos(ω*x) or sin(ω*x).
 * QAWF integrates the same over semi-infinite intervals.
 */
export const enum QUADPACKOscillatory {
  /** w(x) = cos(omega*x) */
  COS = 1,
  /** w(x) = sin(omega*x) */
  SIN = 2,
}
