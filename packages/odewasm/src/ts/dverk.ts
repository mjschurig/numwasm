/**
 * DVERK - Verner 6(5) explicit Runge-Kutta solver
 *
 * Higher-order explicit RK method from Netlib. Uses Verner's 6(5) embedded
 * pair for adaptive step size control. More accurate than RKF45 for smooth problems.
 */

import { loadODEModule } from "./loader.js";
import type { ODEFunction, ODEModule, ODESolverResult } from "./types.js";

const DOUBLE_SIZE = 8;
const INT_SIZE = 4;

/**
 * Options for the DVERK solver.
 */
export interface DverkOptions {
  /**
   * Error tolerance. Controls local error per step.
   * @default 1e-6
   */
  tol?: number;

  /**
   * Error control mode (0-5):
   * - 0: No error control (first call uses tol)
   * - 1: Absolute error control
   * - 2: Relative error control
   * - 3: Mixed error control
   * - 4: Componentwise absolute error (requires n+30 size c array)
   * - 5: Componentwise relative error (requires n+30 size c array)
   * @default 2
   */
  error_control?: 0 | 1 | 2 | 3 | 4 | 5;

  /**
   * Minimum step size. Set to 0 for automatic.
   * @default 0
   */
  hmin?: number;

  /**
   * Initial step size. Set to 0 for automatic.
   * @default 0
   */
  hstart?: number;

  /**
   * Maximum step size.
   * @default Infinity
   */
  hmax?: number;

  /**
   * Maximum number of function evaluations.
   * @default 100000
   */
  max_nfcn?: number;
}

/**
 * Solve an ODE using the DVERK method (Verner 6(5) RK).
 *
 * DVERK is a higher-order explicit Runge-Kutta method with 6th order
 * accuracy and 5th order error estimation. More efficient than RKF45
 * for problems requiring higher accuracy.
 *
 * @param fun - ODE function: f(t, y) => dy/dt
 * @param t_span - Time span [t0, tf]
 * @param y0 - Initial state vector
 * @param options - Solver options
 * @returns Solution result
 *
 * @example
 * ```ts
 * import { dverk } from 'odewasm';
 *
 * // Harmonic oscillator: y'' + y = 0
 * const result = await dverk(
 *   (t, [y, v]) => [v, -y],
 *   [0, 10],
 *   [1, 0],
 *   { tol: 1e-10 }
 * );
 * ```
 */
export async function dverk(
  fun: ODEFunction,
  t_span: [number, number],
  y0: number[],
  options?: DverkOptions,
): Promise<ODESolverResult> {
  const Module = await loadODEModule();

  const {
    tol = 1e-6,
    error_control = 2,
    hmin = 0,
    hstart = 0,
    hmax = Infinity,
    max_nfcn = 100000,
  } = options ?? {};

  const n = y0.length;
  const [t0, tf] = t_span;

  const tValues: number[] = [t0];
  const yValues: number[][] = [y0.slice()];

  const { ind, nfev } = await solveDverkInternal(
    Module,
    n,
    t0,
    y0,
    tf,
    tol,
    error_control,
    hmin,
    hstart,
    hmax,
    max_nfcn,
    fun,
    tValues,
    yValues,
  );

  // Build result
  const success = ind === 3;
  const message = getDverkMessage(ind);

  const t = new Float64Array(tValues);
  const y: Float64Array[] = [];
  for (let i = 0; i < n; i++) {
    const yi = new Float64Array(tValues.length);
    for (let j = 0; j < tValues.length; j++) {
      yi[j] = yValues[j][i];
    }
    y.push(yi);
  }

  return {
    t,
    y,
    nfev,
    njev: 0,
    nlu: 0,
    status: ind,
    message,
    success,
  };
}

function getDverkMessage(ind: number): string {
  switch (ind) {
    case 3:
      return "Integration successful.";
    case 4:
      return "Computation stopped at x < xend.";
    case 5:
      return "Too many function evaluations.";
    case 6:
      return "Invalid input parameters.";
    default:
      return `Unknown status: ${ind}`;
  }
}

async function solveDverkInternal(
  Module: ODEModule,
  n: number,
  t0: number,
  y0: number[],
  tf: number,
  tol: number,
  error_control: number,
  hmin: number,
  hstart: number,
  hmax: number,
  max_nfcn: number,
  fun: ODEFunction,
  tValues: number[],
  yValues: number[][],
): Promise<{ ind: number; nfev: number }> {
  const cSize = Module._wasm_dverk_c_size(n, error_control);
  const nw = Module._wasm_dverk_work_size(n);

  // Allocate memory
  const xPtr = Module._malloc(DOUBLE_SIZE);
  const yPtr = Module._malloc(n * DOUBLE_SIZE);
  const indPtr = Module._malloc(INT_SIZE);
  const cPtr = Module._malloc(cSize * DOUBLE_SIZE);
  const wPtr = Module._malloc(nw * DOUBLE_SIZE);

  try {
    // Initialize memory
    Module.HEAPF64[xPtr >> 3] = t0;
    for (let i = 0; i < n; i++) {
      Module.HEAPF64[(yPtr >> 3) + i] = y0[i];
    }
    Module.HEAP32[indPtr >> 2] = 1; // First call

    // Initialize c array (options)
    for (let i = 0; i < cSize; i++) {
      Module.HEAPF64[(cPtr >> 3) + i] = 0;
    }
    // Set options in c array
    Module.HEAPF64[(cPtr >> 3) + 0] = error_control; // c[1] in Fortran
    Module.HEAPF64[(cPtr >> 3) + 1] = hmin; // c[2] = HMIN
    Module.HEAPF64[(cPtr >> 3) + 2] = hstart; // c[3] = HSTART
    Module.HEAPF64[(cPtr >> 3) + 3] = hmax === Infinity ? 0 : hmax; // c[4] = HMAX
    Module.HEAPF64[(cPtr >> 3) + 4] = max_nfcn; // c[5] = max function evaluations

    // Initialize work array
    for (let i = 0; i < nw; i++) {
      Module.HEAPF64[(wPtr >> 3) + i] = 0;
    }

    // Create FCN callback
    const fcnThunk = (
      _nVal: number,
      xVal: number,
      yVal: number,
      ypVal: number,
    ): void => {
      const yArr: number[] = [];
      for (let i = 0; i < n; i++) {
        yArr.push(Module.HEAPF64[(yVal >> 3) + i]);
      }
      const f = fun(xVal, yArr);
      for (let i = 0; i < n; i++) {
        Module.HEAPF64[(ypVal >> 3) + i] = f[i];
      }
    };

    const fcnPtr = Module.addFunction(fcnThunk, "vidii");
    Module._wasm_set_fcn_callback(fcnPtr);

    try {
      // Call DVERK
      Module._wasm_dverk(n, xPtr, yPtr, tf, tol, indPtr, cPtr, nw, wPtr);

      const ind = Module.HEAP32[indPtr >> 2];
      const xFinal = Module.HEAPF64[xPtr >> 3];
      const yFinal: number[] = [];
      for (let i = 0; i < n; i++) {
        yFinal.push(Module.HEAPF64[(yPtr >> 3) + i]);
      }
      tValues.push(xFinal);
      yValues.push(yFinal);

      // nfev is stored in c[24] (index 23)
      const nfev = Math.round(Module.HEAPF64[(cPtr >> 3) + 23]);

      return { ind, nfev };
    } finally {
      Module.removeFunction(fcnPtr);
    }
  } finally {
    Module._free(xPtr);
    Module._free(yPtr);
    Module._free(indPtr);
    Module._free(cPtr);
    Module._free(wPtr);
  }
}
