/**
 * RKSUITE - Runge-Kutta Suite
 *
 * RKSUITE provides three Runge-Kutta pairs of different orders:
 * - Method 1: (2,3) - Low order, efficient for loose tolerances
 * - Method 2: (4,5) - Medium order, good general purpose
 * - Method 3: (7,8) - High order, efficient for tight tolerances
 *
 * RKSUITE also provides error assessment and dense output.
 */

import { loadODEModule } from "../core/loader.js";
import type { ODEFunction, ODEModule, ODESolverResult } from "../core/types.js";

const DOUBLE_SIZE = 8;
const INT_SIZE = 4;

/**
 * RKSUITE method selection.
 * - 1: RK (2,3) - Low order, 3 stages
 * - 2: RK (4,5) - Medium order, 6 stages (default)
 * - 3: RK (7,8) - High order, 13 stages
 */
export type RksuiteMethod = 1 | 2 | 3;

/**
 * Options for the RKSUITE solver.
 */
export interface RksuiteOptions {
  /**
   * Error tolerance. Controls local error per step.
   * @default 1e-6
   */
  tol?: number;

  /**
   * Threshold array. Components of y below thres[i] are considered
   * effectively zero for error control purposes.
   * Can be scalar (applies to all) or array (per-component).
   * @default 1e-10
   */
  thres?: number | number[];

  /**
   * Runge-Kutta method to use.
   * - 1: (2,3) pair - 3 stages
   * - 2: (4,5) pair - 6 stages (recommended for most problems)
   * - 3: (7,8) pair - 13 stages
   * @default 2
   */
  method?: RksuiteMethod;

  /**
   * Enable global error assessment.
   * When true, provides an estimate of accumulated error.
   * @default false
   */
  errass?: boolean;

  /**
   * Initial step size. Set to 0 for automatic.
   * @default 0
   */
  first_step?: number;
}

/**
 * Solve an ODE using RKSUITE.
 *
 * RKSUITE is a well-tested suite of Runge-Kutta methods with three
 * pairs available: (2,3), (4,5), and (7,8). It includes sophisticated
 * error control and optional global error assessment.
 *
 * @param fun - ODE function: f(t, y) => dy/dt
 * @param t_span - Time span [t0, tf]
 * @param y0 - Initial state vector
 * @param options - Solver options
 * @returns Solution result
 *
 * @example
 * ```ts
 * import { rksuite } from 'odewasm';
 *
 * // Use high-order method for tight tolerance
 * const result = await rksuite(
 *   (t, y) => [-y[0]],
 *   [0, 10],
 *   [1],
 *   { method: 3, tol: 1e-12 }
 * );
 * ```
 */
export async function rksuite(
  fun: ODEFunction,
  t_span: [number, number],
  y0: number[],
  options?: RksuiteOptions,
): Promise<ODESolverResult> {
  const Module = await loadODEModule();

  const {
    tol = 1e-6,
    thres = 1e-10,
    method = 2,
    errass = false,
    first_step = 0,
  } = options ?? {};

  const n = y0.length;
  const [t0, tf] = t_span;

  // Prepare threshold array
  const thresArray = typeof thres === "number" ? new Array(n).fill(thres) : thres;

  const tValues: number[] = [t0];
  const yValues: number[][] = [y0.slice()];

  const { uflag, nfev } = await solveRksuiteInternal(
    Module,
    n,
    t0,
    y0,
    tf,
    tol,
    thresArray,
    method,
    errass,
    first_step,
    fun,
    tValues,
    yValues,
  );

  // Build result
  const success = uflag === 1;
  const message = getRksuiteMessage(uflag);

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
    status: uflag,
    message,
    success,
  };
}

function getRksuiteMessage(uflag: number): string {
  switch (uflag) {
    case 1:
      return "Integration successful.";
    case 2:
      return "Reached end of permitted range.";
    case 3:
      return "Tolerance too loose for method.";
    case 4:
      return "Too many steps, problem may be stiff.";
    case 5:
      return "Cannot start (bad input).";
    case 6:
      return "Catastrophic failure.";
    default:
      return `Unknown status: ${uflag}`;
  }
}

async function solveRksuiteInternal(
  Module: ODEModule,
  n: number,
  t0: number,
  y0: number[],
  tf: number,
  tol: number,
  thresArray: number[],
  method: RksuiteMethod,
  errass: boolean,
  first_step: number,
  fun: ODEFunction,
  tValues: number[],
  yValues: number[][],
): Promise<{ uflag: number; nfev: number }> {
  const lenwrk = Module._wasm_rksuite_work_size(n, method, errass ? 1 : 0);

  // Allocate memory
  const ystartPtr = Module._malloc(n * DOUBLE_SIZE);
  const thresPtr = Module._malloc(n * DOUBLE_SIZE);
  const workPtr = Module._malloc(lenwrk * DOUBLE_SIZE);
  const tgotPtr = Module._malloc(DOUBLE_SIZE);
  const ygotPtr = Module._malloc(n * DOUBLE_SIZE);
  const ypgotPtr = Module._malloc(n * DOUBLE_SIZE);
  const ymaxPtr = Module._malloc(n * DOUBLE_SIZE);
  const uflagPtr = Module._malloc(INT_SIZE);

  try {
    // Initialize ystart and thres
    for (let i = 0; i < n; i++) {
      Module.HEAPF64[(ystartPtr >> 3) + i] = y0[i];
      Module.HEAPF64[(thresPtr >> 3) + i] = thresArray[i];
    }

    // Initialize work array
    for (let i = 0; i < lenwrk; i++) {
      Module.HEAPF64[(workPtr >> 3) + i] = 0;
    }

    // Create FCN callback
    const fcnThunk = (
      _nVal: number,
      tVal: number,
      yVal: number,
      fVal: number,
    ): void => {
      const yArr: number[] = [];
      for (let i = 0; i < n; i++) {
        yArr.push(Module.HEAPF64[(yVal >> 3) + i]);
      }
      const f = fun(tVal, yArr);
      for (let i = 0; i < n; i++) {
        Module.HEAPF64[(fVal >> 3) + i] = f[i];
      }
    };

    const fcnPtr = Module.addFunction(fcnThunk, "vidii");
    Module._wasm_set_fcn_callback(fcnPtr);

    try {
      // Call SETUP
      Module._wasm_rksuite_setup(
        n,
        t0,
        ystartPtr,
        tf,
        tol,
        thresPtr,
        method,
        0, // task_code: 0 = 'U' (UT mode)
        errass ? 1 : 0,
        first_step,
        workPtr,
        lenwrk,
        0, // mesage: 0 = no messages
      );

      // Integrate using UT (usual task) mode
      // Call UT repeatedly until we reach tf
      let nfev = 0;
      let lastUflag = 1;

      while (true) {
        Module._wasm_rksuite_ut(
          tf,
          tgotPtr,
          ygotPtr,
          ypgotPtr,
          ymaxPtr,
          workPtr,
          uflagPtr,
        );

        const uflag = Module.HEAP32[uflagPtr >> 2];
        lastUflag = uflag;

        const tgot = Module.HEAPF64[tgotPtr >> 3];
        const ygot: number[] = [];
        for (let i = 0; i < n; i++) {
          ygot.push(Module.HEAPF64[(ygotPtr >> 3) + i]);
        }
        tValues.push(tgot);
        yValues.push(ygot);

        // Check completion or error
        if (uflag !== 1 || Math.abs(tgot - tf) < 1e-14) {
          break;
        }
      }

      // Get statistics
      const totfcnPtr = Module._malloc(INT_SIZE);
      const stpcstPtr = Module._malloc(INT_SIZE);
      const wastePtr = Module._malloc(DOUBLE_SIZE);
      const stpsokPtr = Module._malloc(INT_SIZE);
      const hnextPtr = Module._malloc(DOUBLE_SIZE);

      try {
        Module._wasm_rksuite_stat(
          totfcnPtr,
          stpcstPtr,
          wastePtr,
          stpsokPtr,
          hnextPtr,
        );
        nfev = Module.HEAP32[totfcnPtr >> 2];
      } finally {
        Module._free(totfcnPtr);
        Module._free(stpcstPtr);
        Module._free(wastePtr);
        Module._free(stpsokPtr);
        Module._free(hnextPtr);
      }

      return { uflag: lastUflag, nfev };
    } finally {
      Module.removeFunction(fcnPtr);
    }
  } finally {
    Module._free(ystartPtr);
    Module._free(thresPtr);
    Module._free(workPtr);
    Module._free(tgotPtr);
    Module._free(ygotPtr);
    Module._free(ypgotPtr);
    Module._free(ymaxPtr);
    Module._free(uflagPtr);
  }
}
