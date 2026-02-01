/**
 * VODPK - VODE with preconditioned Krylov methods
 *
 * VODPK is designed for large stiff systems where direct methods are
 * too expensive. Uses iterative Krylov methods (GMRES) with optional
 * user-supplied preconditioners.
 */

import { loadODEModule } from "../core/loader.js";
import type { ODEFunction, ODEModule, ODESolverResult } from "../core/types.js";

const DOUBLE_SIZE = 8;
const INT_SIZE = 4;

/**
 * Preconditioner solve function type.
 *
 * Solves the system P * x = b where P is the preconditioner.
 *
 * @param n - System dimension
 * @param t - Current time
 * @param y - Current state
 * @param b - Right-hand side (modified to solution on output)
 * @param lr - 1 for left preconditioning, 2 for right
 * @returns 0 for success, non-zero for failure
 */
export type PsolFunction = (
  n: number,
  t: number,
  y: number[],
  b: number[],
  lr: 1 | 2,
) => number;

/**
 * Method flag for VODPK.
 * - 10: Adams method (non-stiff)
 * - 21: BDF with left preconditioning
 * - 22: BDF with no preconditioning
 * - 23: BDF with right preconditioning
 * - 24: BDF with both left and right preconditioning
 */
export type VodpkMethodFlag = 10 | 21 | 22 | 23 | 24;

/**
 * Options for the VODPK solver.
 */
export interface VodpkOptions {
  /**
   * Relative error tolerance.
   * @default 1e-3
   */
  rtol?: number;

  /**
   * Absolute error tolerance.
   * @default 1e-6
   */
  atol?: number | number[];

  /**
   * Method flag.
   * - 22: BDF with no preconditioning (default, simplest)
   * @default 22
   */
  mf?: VodpkMethodFlag;

  /**
   * Preconditioner solve function.
   * Required for mf = 21, 23, or 24.
   */
  psol?: PsolFunction;

  /**
   * Maximum Krylov subspace dimension.
   * @default 5
   */
  maxl?: number;

  /**
   * Maximum number of restarts for GMRES.
   * @default 5
   */
  maxp?: number;

  /**
   * Maximum number of steps.
   * @default 500
   */
  max_steps?: number;
}

/**
 * Solve an ODE using VODPK (VODE with Krylov methods).
 *
 * VODPK is designed for large stiff systems where a full Jacobian would
 * be too expensive to compute and store. It uses GMRES with optional
 * preconditioning.
 *
 * @param fun - ODE function: f(t, y) => dy/dt
 * @param t_span - Time span [t0, tf]
 * @param y0 - Initial state vector
 * @param options - Solver options
 * @returns Solution result
 *
 * @example
 * ```ts
 * import { vodpk } from 'odewasm';
 *
 * // Large stiff system (e.g., discretized PDE)
 * const n = 100;
 * const y0 = new Array(n).fill(1);
 *
 * const result = await vodpk(
 *   (t, y) => y.map((yi, i) => -i * yi),  // Diagonal system
 *   [0, 1],
 *   y0,
 *   { mf: 22 }  // No preconditioning
 * );
 * ```
 */
export async function vodpk(
  fun: ODEFunction,
  t_span: [number, number],
  y0: number[],
  options?: VodpkOptions,
): Promise<ODESolverResult> {
  const Module = await loadODEModule();

  const {
    rtol = 1e-3,
    atol = 1e-6,
    mf = 22,
    psol,
    maxl = 5,
    maxp = 5,
    max_steps = 500,
  } = options ?? {};

  const n = y0.length;
  const [t0, tf] = t_span;

  // Determine tolerance mode
  const itol = typeof atol === "number" ? 1 : 2;
  const atolArray = typeof atol === "number" ? [atol] : atol;

  const tValues: number[] = [t0];
  const yValues: number[][] = [y0.slice()];

  const { istate, nfev } = await solveVodpkInternal(
    Module,
    n,
    t0,
    y0,
    tf,
    itol,
    rtol,
    atolArray,
    mf,
    psol,
    maxl,
    maxp,
    max_steps,
    fun,
    tValues,
    yValues,
  );

  // Build result
  const success = istate === 2;
  const message = getVodpkMessage(istate);

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
    status: istate,
    message,
    success,
  };
}

function getVodpkMessage(istate: number): string {
  switch (istate) {
    case 2:
      return "Integration successful.";
    case -1:
      return "Excess work done on this call.";
    case -2:
      return "Excess accuracy requested.";
    case -3:
      return "Illegal input detected.";
    case -4:
      return "Repeated error test failures.";
    case -5:
      return "Repeated convergence failures (Krylov iteration).";
    case -6:
      return "Error weight became zero.";
    case -7:
      return "Preconditioner solve failed.";
    default:
      return `Unknown status: ${istate}`;
  }
}

async function solveVodpkInternal(
  Module: ODEModule,
  n: number,
  t0: number,
  y0: number[],
  tf: number,
  itol: number,
  rtol: number,
  atolArray: number[],
  mf: VodpkMethodFlag,
  psol: PsolFunction | undefined,
  maxl: number,
  maxp: number,
  max_steps: number,
  fun: ODEFunction,
  tValues: number[],
  yValues: number[][],
): Promise<{ istate: number; nfev: number }> {
  // Work array sizes - use no preconditioner workspace for simplicity
  const lwp = 0;
  const liwp = 0;
  const lrw = Module._wasm_vodpk_rwork_size(n, maxl, maxp, lwp);
  const liw = Module._wasm_vodpk_iwork_size(n, liwp);

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

    // Set options in iwork
    Module.HEAP32[(iworkPtr >> 2) + 5] = max_steps; // MXSTEP
    Module.HEAP32[(iworkPtr >> 2) + 23] = maxl; // MAXL
    Module.HEAP32[(iworkPtr >> 2) + 24] = maxp; // KMP (restarts)

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

    // Create PSOL callback if provided
    let psolPtr = 0;
    if (psol && mf !== 22 && mf !== 10) {
      const psolThunk = (
        _neq: number,
        tVal: number,
        yVal: number,
        _savfVal: number,
        _wkVal: number,
        _hrl1: number,
        _wpVal: number,
        _iwpVal: number,
        bVal: number,
        lr: number,
      ): number => {
        const yArr: number[] = [];
        const bArr: number[] = [];
        for (let i = 0; i < n; i++) {
          yArr.push(Module.HEAPF64[(yVal >> 3) + i]);
          bArr.push(Module.HEAPF64[(bVal >> 3) + i]);
        }
        const ier = psol(n, tVal, yArr, bArr, lr as 1 | 2);
        // Write solution back to b
        for (let i = 0; i < n; i++) {
          Module.HEAPF64[(bVal >> 3) + i] = bArr[i];
        }
        return ier;
      };

      psolPtr = Module.addFunction(psolThunk, "iidiiiidiiii");
      Module._wasm_set_psol_callback(psolPtr);
    }

    try {
      const itask = 1; // Normal mode
      const iopt = 1; // Use optional inputs

      Module._wasm_vodpk(
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
      const nfev = Module.HEAP32[(iworkPtr >> 2) + 11];

      return { istate, nfev };
    } finally {
      Module.removeFunction(fcnPtr);
      if (psolPtr) {
        Module.removeFunction(psolPtr);
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
