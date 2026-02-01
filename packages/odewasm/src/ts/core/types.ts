/**
 * ODE WASM Module Type Definitions
 *
 * High-quality ODE solvers from Ernst Hairer's collection, compiled to WebAssembly.
 * These are the same solvers used in SciPy's solve_ivp function.
 *
 * **Available Solvers:**
 *
 * - **DOPRI5 (RK45)**: Dormand-Prince 5(4) explicit Runge-Kutta method
 *   - Order 5 with order 4 error estimation
 *   - Good general-purpose solver for non-stiff problems
 *   - Supports dense output for interpolation
 *
 * - **DOP853**: Dormand-Prince 8(5,3) explicit Runge-Kutta method
 *   - Order 8 with orders 5 and 3 for error estimation
 *   - Higher accuracy than DOPRI5, fewer steps for smooth problems
 *   - Supports dense output
 *
 * - **RADAU5**: Implicit Runge-Kutta method (Radau IIA, order 5)
 *   - For stiff systems of ODEs
 *   - Uses Newton iteration with LU decomposition
 *   - Supports user-provided Jacobian for efficiency
 *
 * **Callback Architecture:**
 * The solvers use callbacks for the ODE function (fcn), Jacobian (jac),
 * and solution output (solout). Register JavaScript callbacks using
 * _wasm_set_*_callback() functions before calling the solver.
 *
 * **Typical Usage:**
 * 1. Set up callbacks with _wasm_set_fcn_callback() etc.
 * 2. Allocate work arrays with _wasm_*_work_size() helpers
 * 3. Call the solver (_wasm_dopri5, _wasm_dop853, or _wasm_radau5)
 * 4. Check idid for success/failure
 * 5. Use dense output functions (contd5, contd8, contr5) if needed
 *
 * @see http://www.unige.ch/~hairer/software.html
 */

export interface ODEModule {
  // ============================================================================
  // DOPRI5: Dormand-Prince 5(4) Explicit Runge-Kutta
  // ============================================================================

  /**
   * DOPRI5 - Dormand-Prince 5(4) explicit Runge-Kutta solver.
   *
   * Solves the initial value problem: dy/dt = f(t, y), y(x0) = y0
   *
   * This is equivalent to SciPy's RK45 method. Uses embedded error
   * estimation for adaptive step size control.
   *
   * @param n - Number of equations (dimension of y)
   * @param x - Pointer to initial time (input), final time reached (output)
   * @param y - Pointer to initial state vector (input), final state (output)
   * @param xend - Pointer to target end time
   * @param rtol - Pointer to relative tolerance (scalar or array)
   * @param atol - Pointer to absolute tolerance (scalar or array)
   * @param itol - Tolerance type: 0=scalar rtol/atol, 1=array rtol/atol
   * @param iout - Output mode: 0=no output, 1=call solout after each step
   * @param work - Pointer to double work array
   * @param lwork - Length of work array (use _wasm_dopri5_work_size)
   * @param iwork - Pointer to integer work array
   * @param liwork - Length of iwork array (use _wasm_dopri5_iwork_size)
   * @param idid - Pointer to return status (output): 1=success, negative=error
   */
  _wasm_dopri5(
    n: number,
    x: number,
    y: number,
    xend: number,
    rtol: number,
    atol: number,
    itol: number,
    iout: number,
    work: number,
    lwork: number,
    iwork: number,
    liwork: number,
    idid: number,
  ): void;

  // ============================================================================
  // DOP853: Dormand-Prince 8(5,3) Explicit Runge-Kutta
  // ============================================================================

  /**
   * DOP853 - Dormand-Prince 8(5,3) explicit Runge-Kutta solver.
   *
   * Higher-order version of DOPRI5 with order 8 accuracy. Uses embedded
   * methods of orders 5 and 3 for error estimation. More efficient than
   * DOPRI5 when high accuracy is required.
   *
   * Parameters are identical to DOPRI5.
   *
   * @param n - Number of equations
   * @param x - Pointer to current time (input/output)
   * @param y - Pointer to state vector (input/output)
   * @param xend - Pointer to target end time
   * @param rtol - Pointer to relative tolerance
   * @param atol - Pointer to absolute tolerance
   * @param itol - Tolerance type: 0=scalar, 1=array
   * @param iout - Output mode: 0=none, 1=call solout
   * @param work - Pointer to work array
   * @param lwork - Length of work (use _wasm_dop853_work_size)
   * @param iwork - Pointer to integer work array
   * @param liwork - Length of iwork (use _wasm_dop853_iwork_size)
   * @param idid - Pointer to return status
   */
  _wasm_dop853(
    n: number,
    x: number,
    y: number,
    xend: number,
    rtol: number,
    atol: number,
    itol: number,
    iout: number,
    work: number,
    lwork: number,
    iwork: number,
    liwork: number,
    idid: number,
  ): void;

  // ============================================================================
  // RADAU5: Implicit Runge-Kutta for Stiff Systems
  // ============================================================================

  /**
   * RADAU5 - Implicit Runge-Kutta solver for stiff ODEs.
   *
   * Uses the Radau IIA method of order 5. Suitable for stiff systems where
   * explicit methods would require impractically small step sizes.
   *
   * Can solve both explicit ODEs (dy/dt = f(t,y)) and implicit ODEs
   * (M * dy/dt = f(t,y)) where M is the mass matrix.
   *
   * @param n - Number of equations
   * @param x - Pointer to current time (input/output)
   * @param y - Pointer to state vector (input/output)
   * @param xend - Pointer to target end time
   * @param h - Pointer to initial step size (0 for automatic)
   * @param rtol - Pointer to relative tolerance
   * @param atol - Pointer to absolute tolerance
   * @param itol - Tolerance type: 0=scalar, 1=array
   * @param ijac - Jacobian: 0=internal numerical, 1=user-supplied
   * @param mljac - Jacobian bandwidth lower (n for full Jacobian)
   * @param mujac - Jacobian bandwidth upper (n for full Jacobian)
   * @param imas - Mass matrix: 0=identity (explicit ODE), 1=user-supplied
   * @param mlmas - Mass matrix bandwidth lower
   * @param mumas - Mass matrix bandwidth upper
   * @param iout - Output mode: 0=none, 1=call solout
   * @param work - Pointer to work array
   * @param lwork - Length of work (use _wasm_radau5_work_size)
   * @param iwork - Pointer to integer work array
   * @param liwork - Length of iwork (use _wasm_radau5_iwork_size)
   * @param idid - Pointer to return status
   */
  _wasm_radau5(
    n: number,
    x: number,
    y: number,
    xend: number,
    h: number,
    rtol: number,
    atol: number,
    itol: number,
    ijac: number,
    mljac: number,
    mujac: number,
    imas: number,
    mlmas: number,
    mumas: number,
    iout: number,
    work: number,
    lwork: number,
    iwork: number,
    liwork: number,
    idid: number,
  ): void;

  // ============================================================================
  // ODEX: GBS Extrapolation Method
  // ============================================================================

  /**
   * ODEX - GBS extrapolation solver based on explicit midpoint rule.
   *
   * An extrapolation algorithm with stepsize control, order selection,
   * and dense output. Based on the explicit midpoint rule, this method
   * is highly accurate for smooth problems and can achieve very high orders.
   *
   * @param n - Number of equations (dimension of y)
   * @param x - Pointer to initial time (input), final time reached (output)
   * @param y - Pointer to initial state vector (input), final state (output)
   * @param xend - Pointer to target end time
   * @param h - Pointer to initial step size (0 for automatic)
   * @param rtol - Pointer to relative tolerance (scalar or array)
   * @param atol - Pointer to absolute tolerance (scalar or array)
   * @param itol - Tolerance type: 0=scalar rtol/atol, 1=array rtol/atol
   * @param iout - Output mode: 0=no output, 1=after each step, 2=dense output
   * @param work - Pointer to double work array
   * @param lwork - Length of work array (use _wasm_odex_work_size)
   * @param iwork - Pointer to integer work array
   * @param liwork - Length of iwork array (use _wasm_odex_iwork_size)
   * @param idid - Pointer to return status (output): 1=success, negative=error
   */
  _wasm_odex(
    n: number,
    x: number,
    y: number,
    xend: number,
    h: number,
    rtol: number,
    atol: number,
    itol: number,
    iout: number,
    work: number,
    lwork: number,
    iwork: number,
    liwork: number,
    idid: number,
  ): void;

  // ============================================================================
  // DENSE OUTPUT INTERPOLATION
  // ============================================================================
  // These functions evaluate the solution at intermediate points using
  // the dense output computed during integration.
  // ============================================================================

  /**
   * CONTD5 - Dense output interpolation for DOPRI5.
   *
   * Evaluates the i-th component of the solution at time x using
   * the continuous extension from the last accepted step.
   *
   * @param ii - Component index (1-based)
   * @param x - Time at which to evaluate
   * @param con - Pointer to continuous output data from solver
   * @param icomp - Pointer to component selection array
   * @param nd - Number of dense output components
   * @returns The interpolated value y_ii(x)
   */
  _wasm_contd5(
    ii: number,
    x: number,
    con: number,
    icomp: number,
    nd: number,
  ): number;

  /**
   * CONTD8 - Dense output interpolation for DOP853.
   * Same parameters as CONTD5.
   */
  _wasm_contd8(
    ii: number,
    x: number,
    con: number,
    icomp: number,
    nd: number,
  ): number;

  /**
   * CONTR5 - Dense output interpolation for RADAU5.
   *
   * @param ii - Component index (1-based)
   * @param s - Scaled time within step
   * @param cont - Pointer to continuous output data
   * @param lrc - Length of cont array
   * @returns The interpolated value
   */
  _wasm_contr5(ii: number, s: number, cont: number, lrc: number): number;

  /**
   * CONTEX - Dense output interpolation for ODEX.
   *
   * Evaluates the i-th component of the solution at time x using
   * the continuous extension from the last accepted step.
   *
   * @param ii - Component index (1-based)
   * @param x - Time at which to evaluate
   * @param con - Pointer to continuous output data from solver
   * @param ncon - Length of continuous output array
   * @param icomp - Pointer to component selection array
   * @param nd - Number of dense output components
   * @returns The interpolated value y_ii(x)
   */
  _wasm_contex(
    ii: number,
    x: number,
    con: number,
    ncon: number,
    icomp: number,
    nd: number,
  ): number;

  // ============================================================================
  // CALLBACK SETTERS
  // ============================================================================
  // Register JavaScript callbacks before calling the solver.
  // Use addFunction() to create function pointers from JS functions.
  // ============================================================================

  /**
   * Set the ODE function callback.
   *
   * The callback signature is: fcn(n, x, y, f)
   * where n=dimension, x=time, y=state (input), f=derivatives (output)
   *
   * @param ptr - Function pointer from addFunction(fn, 'viiii')
   */
  _wasm_set_fcn_callback(ptr: number): void;

  /**
   * Set the solution output callback.
   *
   * Called after each accepted step if iout=1.
   * Callback signature: solout(nr, xold, x, y, n, con, icomp, nd, irtrn)
   *
   * @param ptr - Function pointer for solout callback
   */
  _wasm_set_solout_callback(ptr: number): void;

  /**
   * Set the Jacobian callback (for RADAU5 with ijac=1).
   *
   * Callback signature: jac(n, x, y, dfy, ldfy)
   * Computes df/dy at (x, y) and stores in dfy.
   *
   * @param ptr - Function pointer for Jacobian callback
   */
  _wasm_set_jac_callback(ptr: number): void;

  // ============================================================================
  // WORK ARRAY SIZE HELPERS
  // ============================================================================
  // Use these to determine required workspace sizes before calling solvers.
  // ============================================================================

  /**
   * Get required work array size for DOPRI5.
   * @param n - Number of equations
   * @param nrdens - Number of components for dense output (0 if not needed)
   * @returns Required size in doubles
   */
  _wasm_dopri5_work_size(n: number, nrdens: number): number;

  /**
   * Get required iwork array size for DOPRI5.
   * @param nrdens - Number of components for dense output
   * @returns Required size in integers
   */
  _wasm_dopri5_iwork_size(nrdens: number): number;

  /**
   * Get required work array size for DOP853.
   * @param n - Number of equations
   * @param nrdens - Number of components for dense output
   * @returns Required size in doubles
   */
  _wasm_dop853_work_size(n: number, nrdens: number): number;

  /**
   * Get required iwork array size for DOP853.
   * @param nrdens - Number of components for dense output
   * @returns Required size in integers
   */
  _wasm_dop853_iwork_size(nrdens: number): number;

  /**
   * Get required work array size for RADAU5.
   * @param n - Number of equations
   * @returns Required size in doubles
   */
  _wasm_radau5_work_size(n: number): number;

  /**
   * Get required iwork array size for RADAU5.
   * @param n - Number of equations
   * @returns Required size in integers
   */
  _wasm_radau5_iwork_size(n: number): number;

  /**
   * Get required work array size for ODEX.
   * @param n - Number of equations
   * @param nrdens - Number of components for dense output (0 if not used)
   * @param km - Maximum columns in extrapolation table (default 9)
   * @returns Required size in doubles
   */
  _wasm_odex_work_size(n: number, nrdens: number, km: number): number;

  /**
   * Get required iwork array size for ODEX.
   * @param nrdens - Number of components for dense output (0 if not used)
   * @param km - Maximum columns in extrapolation table (default 9)
   * @returns Required size in integers
   */
  _wasm_odex_iwork_size(nrdens: number, km: number): number;

  // ============================================================================
  // NETLIB SOLVERS
  // ============================================================================

  /** RKF45 solver entry point */
  _wasm_rkf45(
    neqn: number,
    y: number,
    t: number,
    tout: number,
    relerr: number,
    abserr: number,
    iflag: number,
    work: number,
    iwork: number,
  ): void;

  /** DVERK solver entry point */
  _wasm_dverk(
    n: number,
    x: number,
    y: number,
    xend: number,
    tol: number,
    ind: number,
    c: number,
    nw: number,
    w: number,
  ): void;

  /** ODE (Adams-Bashforth-Moulton) solver entry point */
  _wasm_ode(
    neqn: number,
    y: number,
    t: number,
    tout: number,
    relerr: number,
    abserr: number,
    iflag: number,
    work: number,
    iwork: number,
  ): void;

  /** VODE solver entry point */
  _wasm_vode(
    neq: number,
    y: number,
    t: number,
    tout: number,
    itol: number,
    rtol: number,
    atol: number,
    itask: number,
    istate: number,
    iopt: number,
    rwork: number,
    lrw: number,
    iwork: number,
    liw: number,
    mf: number,
  ): void;

  /** ZVODE solver entry point (complex-valued) */
  _wasm_zvode(
    neq: number,
    y: number,
    t: number,
    tout: number,
    itol: number,
    rtol: number,
    atol: number,
    itask: number,
    istate: number,
    iopt: number,
    zwork: number,
    lzw: number,
    rwork: number,
    lrw: number,
    iwork: number,
    liw: number,
    mf: number,
  ): void;

  /** VODPK solver entry point (Krylov methods) */
  _wasm_vodpk(
    neq: number,
    y: number,
    t: number,
    tout: number,
    itol: number,
    rtol: number,
    atol: number,
    itask: number,
    istate: number,
    iopt: number,
    rwork: number,
    lrw: number,
    iwork: number,
    liw: number,
    mf: number,
  ): void;

  /** RKSUITE setup */
  _wasm_rksuite_setup(
    neq: number,
    tstart: number,
    ystart: number,
    tend: number,
    tol: number,
    thres: number,
    method: number,
    task_code: number,
    errass: number,
    hstart: number,
    work: number,
    lenwrk: number,
    mesage: number,
  ): void;

  /** RKSUITE UT mode step */
  _wasm_rksuite_ut(
    twant: number,
    tgot: number,
    ygot: number,
    ypgot: number,
    ymax: number,
    work: number,
    uflag: number,
  ): void;

  /** RKSUITE CT mode step */
  _wasm_rksuite_ct(
    tnow: number,
    ynow: number,
    ypnow: number,
    work: number,
    cflag: number,
  ): void;

  /** RKSUITE statistics */
  _wasm_rksuite_stat(
    totfcn: number,
    stpcst: number,
    waste: number,
    stpsok: number,
    hnext: number,
  ): void;

  /** RKC solver entry point */
  _wasm_rkc(
    neqn: number,
    y: number,
    t: number,
    tend: number,
    rtol: number,
    atol: number,
    info: number,
    work: number,
    idid: number,
  ): void;

  /** RKC dense output interpolation */
  _wasm_rkc_int(t: number, yint: number, neqn: number, work: number): void;

  // ============================================================================
  // NETLIB WORK ARRAY SIZE HELPERS
  // ============================================================================

  /** RKF45 work array size */
  _wasm_rkf45_work_size(n: number): number;
  /** RKF45 iwork array size */
  _wasm_rkf45_iwork_size(): number;

  /** DVERK c array size */
  _wasm_dverk_c_size(n: number, error_control: number): number;
  /** DVERK work array size */
  _wasm_dverk_work_size(n: number): number;

  /** ODE work array size */
  _wasm_ode_work_size(n: number): number;
  /** ODE iwork array size */
  _wasm_ode_iwork_size(): number;

  /** VODE rwork array size */
  _wasm_vode_rwork_size(n: number, mf: number): number;
  /** VODE iwork array size */
  _wasm_vode_iwork_size(n: number, mf: number): number;

  /** ZVODE zwork array size */
  _wasm_zvode_zwork_size(n: number, mf: number): number;
  /** ZVODE rwork array size */
  _wasm_zvode_rwork_size(n: number, mf: number): number;

  /** VODPK rwork array size */
  _wasm_vodpk_rwork_size(n: number, maxl: number, maxp: number, lwp: number): number;
  /** VODPK iwork array size */
  _wasm_vodpk_iwork_size(n: number, liwp: number): number;

  /** RKSUITE work array size */
  _wasm_rksuite_work_size(n: number, method: number, errass: number): number;

  /** RKC work array size */
  _wasm_rkc_work_size(n: number, use_spcrad: number): number;

  // ============================================================================
  // NETLIB CALLBACK SETTERS
  // ============================================================================

  /** Set VODE/VODPK Jacobian callback */
  _wasm_set_jac_vode_callback(ptr: number): void;
  /** Set VODPK preconditioner solve callback */
  _wasm_set_psol_callback(ptr: number): void;
  /** Set RKC spectral radius callback */
  _wasm_set_spcrad_callback(ptr: number): void;

  // ============================================================================
  // MEMORY MANAGEMENT
  // ============================================================================

  /** Allocate memory (WASM wrapper) */
  _wasm_malloc(size: number): number;
  /** Free memory (WASM wrapper) */
  _wasm_free(ptr: number): void;
  /** Allocate memory (standard Emscripten) */
  _malloc(size: number): number;
  /** Free memory (standard Emscripten) */
  _free(ptr: number): void;

  // ============================================================================
  // EMSCRIPTEN RUNTIME
  // ============================================================================

  /** Read a value from WASM memory */
  getValue(ptr: number, type: string): number;
  /** Write a value to WASM memory */
  setValue(ptr: number, value: number, type: string): void;

  // ============================================================================
  // HEAP VIEWS
  // ============================================================================

  /** Float64 (double) view of WASM heap */
  HEAPF64: Float64Array;
  /** Float32 (float) view of WASM heap */
  HEAPF32: Float32Array;
  /** Int32 view of WASM heap */
  HEAP32: Int32Array;
  /** Int8 view of WASM heap */
  HEAP8: Int8Array;
  /** Uint8 view of WASM heap */
  HEAPU8: Uint8Array;
  /** Uint32 view of WASM heap */
  HEAPU32: Uint32Array;

  // ============================================================================
  // FUNCTION TABLE MANAGEMENT
  // ============================================================================

  /**
   * Register a JavaScript function for use as a WASM callback.
   *
   * @param func - JavaScript function to register
   * @param signature - Function signature (e.g., 'viiii' for void(int,int,int,int))
   * @returns Function pointer that can be passed to C functions
   */
  addFunction(func: Function, signature: string): number;

  /**
   * Remove a previously registered function.
   * @param ptr - Function pointer returned by addFunction
   */
  removeFunction(ptr: number): void;
}

/**
 * Options for configuring the ODE WASM module loader.
 */
export interface ODEModuleOptions {
  /**
   * Custom function to locate WASM files.
   * @param path - Filename being requested (e.g., 'ode.wasm')
   * @param scriptDirectory - Directory of the JS loader script
   * @returns Full URL or path to the file
   */
  locateFile?: (path: string, scriptDirectory: string) => string;
}

/**
 * Factory function to create an ODE WASM module instance.
 *
 * @param options - Optional configuration
 * @returns Promise resolving to initialized ODE module
 *
 * @example
 * import createODEModule from 'odewasm';
 *
 * const ode = await createODEModule();
 *
 * // Solve dy/dt = -y (exponential decay)
 * const n = 1;
 * const fcn = (n, x, y, f) => {
 *   ode.HEAPF64[f / 8] = -ode.HEAPF64[y / 8];
 * };
 * const fcnPtr = ode.addFunction(fcn, 'viiii');
 * ode._wasm_set_fcn_callback(fcnPtr);
 *
 * // ... set up arrays and call _wasm_dopri5 ...
 */
export type ODEModuleFactory = (
  options?: ODEModuleOptions,
) => Promise<ODEModule>;

// ============================================================================
// HIGH-LEVEL TYPESCRIPT TYPES
// ============================================================================
// These types are used by the high-level solve_ivp wrapper, not the raw WASM.
// ============================================================================

/**
 * ODE function type for the high-level interface.
 *
 * Defines the right-hand side of dy/dt = f(t, y).
 *
 * @param t - Current time
 * @param y - Current state vector
 * @returns The derivative dy/dt as an array
 *
 * @example
 * // Simple harmonic oscillator: y'' + y = 0
 * // Let y[0] = position, y[1] = velocity
 * const harmonic: ODEFunction = (t, y) => [y[1], -y[0]];
 */
export type ODEFunction = (t: number, y: number[]) => number[];

/**
 * Jacobian function type for stiff solvers.
 *
 * Computes the Jacobian matrix df/dy at (t, y).
 *
 * @param t - Current time
 * @param y - Current state vector
 * @returns The Jacobian matrix as row-major 2D array
 *
 * @example
 * // For dy/dt = -k*y (exponential decay)
 * const jac: JacobianFunction = (t, y) => [[-k]];
 */
export type JacobianFunction = (t: number, y: number[]) => number[][];

/**
 * Event function for detecting zero crossings.
 *
 * The solver can detect when this function crosses zero and optionally
 * terminate integration or record the event.
 *
 * @param t - Current time
 * @param y - Current state vector
 * @returns A scalar value; zero crossing triggers an event
 */
export type EventFunction = (t: number, y: number[]) => number;

/**
 * Solution output callback.
 *
 * Called after each accepted step to allow monitoring or early termination.
 *
 * @param t - Current time
 * @param y - Current state vector
 * @returns false to continue, true to terminate integration
 */
export type SoloutCallback = (t: number, y: number[]) => boolean;

/**
 * Available ODE integration methods.
 *
 * - 'RK45': DOPRI5 - explicit Runge-Kutta 5(4), good for non-stiff problems
 * - 'DOP853': Higher-order explicit RK 8(5,3), more accurate for smooth problems
 * - 'Radau': Implicit Radau IIA order 5, for stiff problems
 */
export type ODEMethod = ExplicitMethod | ImplicitMethod | StabilizedMethod;

/**
 * Options for the solve_ivp function.
 *
 * Similar to scipy.integrate.solve_ivp options.
 */
export interface SolveIVPOptions {
  /**
   * Integration method.
   * - 'RK45': Explicit Runge-Kutta 5(4) - good default for non-stiff problems
   * - 'DOP853': Higher-order RK 8(5,3) - more accurate, fewer steps for smooth problems
   * - 'Radau': Implicit Radau IIA - for stiff problems
   * @default 'RK45'
   */
  method?: ODEMethod;

  /**
   * Relative tolerance. Local error is estimated and kept below rtol * |y| + atol.
   * @default 1e-3
   */
  rtol?: number;

  /**
   * Absolute tolerance. Scalar or per-component array.
   * @default 1e-6
   */
  atol?: number | number[];

  /**
   * Maximum allowed step size. Use to prevent the solver from taking
   * steps that are too large for your problem.
   * @default Infinity
   */
  max_step?: number;

  /**
   * Initial step size. Set to 0 for automatic selection.
   * @default 0
   */
  first_step?: number;

  /**
   * Enable dense (continuous) output for interpolation between grid points.
   * @default false
   */
  dense_output?: boolean;

  /**
   * Times at which to store the computed solution.
   * If not provided, solution is stored at solver-chosen steps.
   */
  t_eval?: number[];

  /**
   * Jacobian function for stiff solvers.
   * Required for Radau with ijac=1 for efficiency.
   * If not provided, numerical differencing is used.
   */
  jac?: JacobianFunction;
}

/**
 * Dense output object for continuous solution interpolation.
 *
 * When dense_output=true, this object allows evaluating the solution
 * at any time within the integration interval, not just at grid points.
 */
export interface DenseOutput {
  /**
   * Interpolate the solution at time t.
   * @param t - Time at which to evaluate (must be within [t_min, t_max])
   * @returns Solution vector y(t)
   */
  (t: number): number[];

  /** Minimum valid interpolation time (start of integration) */
  t_min: number;

  /** Maximum valid interpolation time (end of integration) */
  t_max: number;
}

/**
 * Result from solve_ivp.
 *
 * Contains the computed solution, statistics, and success status.
 */
export interface SolveIVPResult {
  /** Time points at which solution was computed */
  t: Float64Array;

  /**
   * Solution values. y[i] is the i-th component of the state vector
   * as a Float64Array over all time points.
   *
   * To get y at time t[j]: [y[0][j], y[1][j], ...]
   */
  y: Float64Array[];

  /**
   * Dense output function for interpolation (only if dense_output=true).
   * Call sol(t) to evaluate solution at any time in [t[0], t[-1]].
   */
  sol?: DenseOutput;

  /** Number of function (f) evaluations */
  nfev: number;

  /** Number of Jacobian evaluations (for stiff solvers) */
  njev: number;

  /** Number of LU decompositions (for stiff solvers) */
  nlu: number;

  /**
   * Solver status code.
   * Positive values indicate success, negative values indicate errors.
   */
  status: number;

  /** Human-readable status message */
  message: string;

  /** True if integration completed successfully */
  success: boolean;
}

/**
 * Status messages for solver return codes (idid values).
 *
 * Maps the idid return value to human-readable descriptions.
 * Positive values indicate success, negative values indicate errors.
 */
export const STATUS_MESSAGES: Record<number, string> = {
  /** Normal completion, reached xend */
  1: "Integration successful.",
  /** Integration successful (stopped at xend) */
  2: "Integration successful (stopped at xend).",
  /** Input parameters are inconsistent */
  [-1]: "Input is not consistent.",
  /** Maximum number of steps exceeded (increase nmax in iwork) */
  [-2]: "Larger nmax is needed.",
  /** Step size became too small (likely stiff problem or singularity) */
  [-3]: "Step size too small.",
  /** Matrix is repeatedly singular (for implicit methods) */
  [-4]: "Matrix is repeatedly singular.",
};

export enum ExplicitMethod {
  /** DOPRI5 - Dormand-Prince 5(4) explicit Runge-Kutta */
  RK45 = "RK45",
  /** DOP853 - Dormand-Prince 8(5,3) high-order explicit Runge-Kutta */
  DOP853 = "DOP853",
  /** RKF45 - Runge-Kutta-Fehlberg 4(5) */
  RKF45 = "RKF45",
  /** DVERK - Verner 6(5) Runge-Kutta */
  DVERK = "DVERK",
  /** ODE - Adams-Bashforth-Moulton multistep */
  ODE = "ODE",
  /** RKSUITE - RK Suite (2,3) pair */
  RKSUITE23 = "RKSUITE23",
  /** RKSUITE - RK Suite (4,5) pair */
  RKSUITE45 = "RKSUITE45",
  /** RKSUITE - RK Suite (7,8) pair */
  RKSUITE78 = "RKSUITE78",
  /** ODEX - GBS extrapolation based on explicit midpoint rule */
  ODEX = "ODEX",
}

export enum ImplicitMethod {
  /** RADAU5 - Implicit Radau IIA order 5 for stiff systems */
  Radau = "Radau",
  /** VODE - Variable-coefficient ODE solver with BDF/Adams */
  VODE = "VODE",
  /** BDF - Backward Differentiation Formula via VODE */
  BDF = "BDF",
}

export enum StabilizedMethod {
  /** RKC - Runge-Kutta-Chebyshev for mildly stiff problems */
  RKC = "RKC",
}

// ============================================================================
// SOLVER-SPECIFIC OPTIONS AND RESULTS
// ============================================================================
// Each solver has its own options and result types for type safety.
// ============================================================================

/**
 * Base options shared by all ODE solvers.
 */
export interface ODESolverOptionsBase {
  /**
   * Relative tolerance. Local error is estimated and kept below rtol * |y| + atol.
   * @default 1e-3
   */
  rtol?: number;

  /**
   * Absolute tolerance. Scalar or per-component array.
   * @default 1e-6
   */
  atol?: number | number[];

  /**
   * Initial step size. Set to 0 for automatic selection.
   * @default 0
   */
  first_step?: number;

  /**
   * Enable dense (continuous) output for interpolation between grid points.
   * @default false
   */
  dense_output?: boolean;

  /**
   * Times at which to store the computed solution.
   * If not provided, solution is stored at solver-chosen steps.
   */
  t_eval?: number[];
}

/**
 * Options for the DOPRI5 solver (Dormand-Prince 5(4)).
 *
 * DOPRI5 is equivalent to SciPy's RK45. Good general-purpose solver for
 * non-stiff problems with moderate accuracy requirements.
 */
export interface Dopri5Options extends ODESolverOptionsBase {
  /**
   * Maximum allowed step size. Use to prevent the solver from taking
   * steps that are too large for your problem.
   * @default Infinity
   */
  max_step?: number;
}

/**
 * Options for the DOP853 solver (Dormand-Prince 8(5,3)).
 *
 * DOP853 is a higher-order explicit method. More efficient than DOPRI5
 * when high accuracy is required or the solution is smooth.
 */
export interface Dop853Options extends ODESolverOptionsBase {
  /**
   * Maximum allowed step size. Use to prevent the solver from taking
   * steps that are too large for your problem.
   * @default Infinity
   */
  max_step?: number;
}

/**
 * Options for the ODEX solver (GBS Extrapolation).
 *
 * ODEX is a GBS extrapolation method based on the explicit midpoint rule.
 * It features automatic order selection (up to very high orders) and is
 * highly efficient for smooth problems requiring high accuracy.
 */
export interface OdexOptions extends ODESolverOptionsBase {
  /**
   * Maximum allowed step size. Use to prevent the solver from taking
   * steps that are too large for your problem.
   * @default Infinity
   */
  max_step?: number;

  /**
   * Maximum number of columns in the extrapolation table.
   * Higher values allow higher-order approximations but require more work.
   * @default 9
   */
  max_columns?: number;

  /**
   * Step size sequence selector.
   * - 1: 2,4,6,8,10,12,14,16,... (default for iout≤1)
   * - 2: 2,4,8,12,16,20,24,28,...
   * - 3: 2,4,6,8,12,16,24,32,...
   * - 4: 2,6,10,14,18,22,26,30,... (default for iout≥2)
   * - 5: 4,8,12,16,20,24,28,32,...
   * @default undefined (automatic based on output mode)
   */
  step_sequence?: 1 | 2 | 3 | 4 | 5;
}

/**
 * Options for the RADAU5 solver (Implicit Radau IIA order 5).
 *
 * RADAU5 is for stiff systems where explicit methods would require
 * impractically small step sizes. Can use numerical or user-supplied Jacobian.
 */
export interface Radau5Options extends ODESolverOptionsBase {
  /**
   * Jacobian function for the ODE system.
   * If not provided, numerical differencing is used (less efficient).
   * Providing an analytical Jacobian improves performance for stiff systems.
   */
  jac?: JacobianFunction;
}

/**
 * Result from an ODE solver.
 *
 * Contains the computed solution, statistics, and success status.
 */
export interface ODESolverResult {
  /** Time points at which solution was computed */
  t: Float64Array;

  /**
   * Solution values. y[i] is the i-th component of the state vector
   * as a Float64Array over all time points.
   *
   * To get y at time t[j]: [y[0][j], y[1][j], ...]
   */
  y: Float64Array[];

  /**
   * Dense output function for interpolation (only if dense_output=true).
   * Call sol(t) to evaluate solution at any time in [t[0], t[-1]].
   */
  sol?: DenseOutput;

  /** Number of function (f) evaluations */
  nfev: number;

  /** Number of Jacobian evaluations (for stiff solvers, 0 otherwise) */
  njev: number;

  /** Number of LU decompositions (for stiff solvers, 0 otherwise) */
  nlu: number;

  /**
   * Solver status code.
   * Positive values indicate success, negative values indicate errors.
   */
  status: number;

  /** Human-readable status message */
  message: string;

  /** True if integration completed successfully */
  success: boolean;
}
