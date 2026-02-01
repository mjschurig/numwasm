/**
 * VODE - Variable-coefficient ODE solver
 *
 * VODE is the most versatile solver from the ODEPACK family. Supports both
 * Adams methods (non-stiff) and BDF methods (stiff). Can use numerical or
 * user-supplied Jacobians.
 */

import { loadODEModule } from "../core/loader.js";
import type {
  JacobianFunction,
  ODEFunction,
  ODEModule,
  ODESolverResult,
} from "../core/types.js";

const DOUBLE_SIZE = 8;
const INT_SIZE = 4;

/**
 * Method flag for VODE.
 * - 10: Adams method (non-stiff), no Jacobian needed
 * - 21: BDF method (stiff), full Jacobian (user-supplied)
 * - 22: BDF method (stiff), full Jacobian (internal numerical)
 * - 24: BDF method (stiff), banded Jacobian (user-supplied)
 * - 25: BDF method (stiff), banded Jacobian (internal numerical)
 */
export type VodeMethodFlag = 10 | 21 | 22 | 24 | 25;

/**
 * Options for the VODE solver.
 */
export interface VodeOptions {
  /**
   * Relative error tolerance.
   * @default 1e-3
   */
  rtol?: number;

  /**
   * Absolute error tolerance. Can be scalar or per-component array.
   * @default 1e-6
   */
  atol?: number | number[];

  /**
   * Method flag. See VodeMethodFlag type for options.
   * - 10: Adams (non-stiff)
   * - 22: BDF with internal Jacobian (stiff, recommended for most stiff problems)
   * @default 10
   */
  mf?: VodeMethodFlag;

  /**
   * Jacobian function (required for mf=21).
   * Returns the Jacobian matrix df/dy as a 2D array in column-major order.
   */
  jac?: JacobianFunction;

  /**
   * Lower bandwidth for banded Jacobian (mf=24,25).
   */
  ml?: number;

  /**
   * Upper bandwidth for banded Jacobian (mf=24,25).
   */
  mu?: number;

  /**
   * Maximum number of steps.
   * @default 500
   */
  max_steps?: number;

  /**
   * Maximum step size.
   * @default Infinity
   */
  max_step?: number;

  /**
   * Initial step size. Set to 0 for automatic.
   * @default 0
   */
  first_step?: number;
}

/**
 * Solve an ODE using VODE (Variable-coefficient ODE solver).
 *
 * VODE is the most versatile solver in the ODEPACK family. It automatically
 * selects between Adams methods (for non-stiff problems) and BDF methods
 * (for stiff problems) based on the mf parameter.
 *
 * @param fun - ODE function: f(t, y) => dy/dt
 * @param t_span - Time span [t0, tf]
 * @param y0 - Initial state vector
 * @param options - Solver options
 * @returns Solution result
 *
 * @example
 * ```ts
 * import { vode } from 'odewasm';
 *
 * // Non-stiff problem with Adams method
 * const result1 = await vode(
 *   (t, y) => [-y[0]],
 *   [0, 5],
 *   [1],
 *   { mf: 10 }
 * );
 *
 * // Stiff problem with BDF method
 * const result2 = await vode(
 *   (t, y) => [-1000 * y[0]],
 *   [0, 0.1],
 *   [1],
 *   { mf: 22, rtol: 1e-8 }
 * );
 * ```
 */
export async function vode(
  fun: ODEFunction,
  t_span: [number, number],
  y0: number[],
  options?: VodeOptions,
): Promise<ODESolverResult> {
  const Module = await loadODEModule();

  const {
    rtol = 1e-3,
    atol = 1e-6,
    mf = 10,
    jac,
    ml,
    mu,
    max_steps = 500,
    max_step = Infinity,
    first_step = 0,
  } = options ?? {};

  const n = y0.length;
  const [t0, tf] = t_span;

  // Determine tolerance mode
  const itol = typeof atol === "number" ? 1 : 2;
  const atolArray = typeof atol === "number" ? [atol] : atol;

  const tValues: number[] = [t0];
  const yValues: number[][] = [y0.slice()];

  const { istate, nfev, njev } = await solveVodeInternal(
    Module,
    n,
    t0,
    y0,
    tf,
    itol,
    rtol,
    atolArray,
    mf,
    jac,
    ml ?? n,
    mu ?? n,
    max_steps,
    max_step,
    first_step,
    fun,
    tValues,
    yValues,
  );

  // Build result
  const success = istate === 2;
  const message = getVodeMessage(istate);

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
    njev,
    nlu: 0,
    status: istate,
    message,
    success,
  };
}

function getVodeMessage(istate: number): string {
  switch (istate) {
    case 2:
      return "Integration successful.";
    case -1:
      return "Excess work done on this call (check max_steps).";
    case -2:
      return "Excess accuracy requested (tolerances too small).";
    case -3:
      return "Illegal input detected.";
    case -4:
      return "Repeated error test failures.";
    case -5:
      return "Repeated convergence failures.";
    case -6:
      return "Error weight became zero.";
    default:
      return `Unknown status: ${istate}`;
  }
}

async function solveVodeInternal(
  Module: ODEModule,
  n: number,
  t0: number,
  y0: number[],
  tf: number,
  itol: number,
  rtol: number,
  atolArray: number[],
  mf: VodeMethodFlag,
  jac: JacobianFunction | undefined,
  ml: number,
  mu: number,
  max_steps: number,
  max_step: number,
  first_step: number,
  fun: ODEFunction,
  tValues: number[],
  yValues: number[][],
): Promise<{ istate: number; nfev: number; njev: number }> {
  const lrw = Module._wasm_vode_rwork_size(n, mf);
  const liw = Module._wasm_vode_iwork_size(n, mf);

  // Allocate memory
  const tPtr = Module._malloc(DOUBLE_SIZE);
  const yPtr = Module._malloc(n * DOUBLE_SIZE);
  const rtolPtr = Module._malloc(DOUBLE_SIZE);
  const atolPtr = Module._malloc(atolArray.length * DOUBLE_SIZE);
  const istatePtr = Module._malloc(INT_SIZE);
  const rworkPtr = Module._malloc(lrw * DOUBLE_SIZE);
  const iworkPtr = Module._malloc(liw * INT_SIZE);

  try {
    // Initialize memory
    Module.HEAPF64[tPtr >> 3] = t0;
    for (let i = 0; i < n; i++) {
      Module.HEAPF64[(yPtr >> 3) + i] = y0[i];
    }
    Module.HEAPF64[rtolPtr >> 3] = rtol;
    for (let i = 0; i < atolArray.length; i++) {
      Module.HEAPF64[(atolPtr >> 3) + i] = atolArray[i];
    }
    Module.HEAP32[istatePtr >> 2] = 1; // First call

    // Initialize work arrays
    for (let i = 0; i < lrw; i++) {
      Module.HEAPF64[(rworkPtr >> 3) + i] = 0;
    }
    for (let i = 0; i < liw; i++) {
      Module.HEAP32[(iworkPtr >> 2) + i] = 0;
    }

    // Set options in iwork/rwork
    // RWORK(5) = H0 (initial step size)
    if (first_step > 0) {
      Module.HEAPF64[(rworkPtr >> 3) + 4] = first_step;
    }
    // RWORK(6) = HMAX (maximum step size)
    if (max_step !== Infinity) {
      Module.HEAPF64[(rworkPtr >> 3) + 5] = max_step;
    }
    // IWORK(6) = MXSTEP (maximum number of steps)
    Module.HEAP32[(iworkPtr >> 2) + 5] = max_steps;

    // Set bandwidth for banded Jacobian
    if (mf === 24 || mf === 25) {
      Module.HEAP32[(iworkPtr >> 2) + 0] = ml;
      Module.HEAP32[(iworkPtr >> 2) + 1] = mu;
    }

    // Create FCN callback
    const fcnThunk = (
      _neq: number,
      tVal: number,
      yVal: number,
      ydotVal: number,
    ): void => {
      const yArr: number[] = [];
      for (let i = 0; i < n; i++) {
        yArr.push(Module.HEAPF64[(yVal >> 3) + i]);
      }
      const f = fun(tVal, yArr);
      for (let i = 0; i < n; i++) {
        Module.HEAPF64[(ydotVal >> 3) + i] = f[i];
      }
    };

    const fcnPtr = Module.addFunction(fcnThunk, "vidii");
    Module._wasm_set_fcn_callback(fcnPtr);

    // Create JAC callback if provided
    let jacPtr = 0;
    if (jac && (mf === 21 || mf === 24)) {
      const jacThunk = (
        _neq: number,
        tVal: number,
        yVal: number,
        _mlVal: number,
        _muVal: number,
        pdVal: number,
        nrowpdVal: number,
      ): void => {
        const yArr: number[] = [];
        for (let i = 0; i < n; i++) {
          yArr.push(Module.HEAPF64[(yVal >> 3) + i]);
        }
        const J = jac(tVal, yArr);
        // Store in column-major order
        for (let j = 0; j < n; j++) {
          for (let i = 0; i < n; i++) {
            Module.HEAPF64[(pdVal >> 3) + i + j * nrowpdVal] = J[i][j];
          }
        }
      };

      jacPtr = Module.addFunction(jacThunk, "vidiiiiii");
      Module._wasm_set_jac_vode_callback(jacPtr);
    }

    try {
      const itask = 1; // Normal mode
      const iopt = first_step > 0 || max_step !== Infinity ? 1 : 0;

      Module._wasm_vode(
        n,
        yPtr,
        tPtr,
        tf,
        itol,
        rtolPtr,
        atolPtr,
        itask,
        istatePtr,
        iopt,
        rworkPtr,
        lrw,
        iworkPtr,
        liw,
        mf,
      );

      const istate = Module.HEAP32[istatePtr >> 2];
      const tFinal = Module.HEAPF64[tPtr >> 3];
      const yFinal: number[] = [];
      for (let i = 0; i < n; i++) {
        yFinal.push(Module.HEAPF64[(yPtr >> 3) + i]);
      }
      tValues.push(tFinal);
      yValues.push(yFinal);

      // Get statistics from iwork
      const nfev = Module.HEAP32[(iworkPtr >> 2) + 11]; // IWORK(12)
      const njev = Module.HEAP32[(iworkPtr >> 2) + 12]; // IWORK(13)

      return { istate, nfev, njev };
    } finally {
      Module.removeFunction(fcnPtr);
      if (jacPtr) {
        Module.removeFunction(jacPtr);
      }
    }
  } finally {
    Module._free(tPtr);
    Module._free(yPtr);
    Module._free(rtolPtr);
    Module._free(atolPtr);
    Module._free(istatePtr);
    Module._free(rworkPtr);
    Module._free(iworkPtr);
  }
}
