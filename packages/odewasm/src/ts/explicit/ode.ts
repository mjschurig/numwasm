/**
 * ODE - Adams-Bashforth-Moulton multistep solver
 *
 * Classic multistep ODE solver from Netlib. Uses predictor-corrector
 * Adams-Bashforth-Moulton method. Good for non-stiff problems where
 * the solution is smooth.
 */

import { loadODEModule } from "../core/loader.js";
import type { ODEFunction, ODEModule, ODESolverResult } from "../core/types.js";

const DOUBLE_SIZE = 8;
const INT_SIZE = 4;

/**
 * Options for the ODE (Adams-Bashforth-Moulton) solver.
 */
export interface OdeOptions {
  /**
   * Relative error tolerance.
   * @default 1e-3
   */
  rtol?: number;

  /**
   * Absolute error tolerance.
   * @default 1e-6
   */
  atol?: number;

  /**
   * One-step mode: if true, return after each internal step.
   * @default false
   */
  one_step?: boolean;
}

/**
 * Solve an ODE using the Adams-Bashforth-Moulton method.
 *
 * This is a multistep predictor-corrector method. More efficient than
 * single-step methods for smooth problems because it reuses previous
 * function evaluations.
 *
 * @param fun - ODE function: f(t, y) => dy/dt
 * @param t_span - Time span [t0, tf]
 * @param y0 - Initial state vector
 * @param options - Solver options
 * @returns Solution result
 *
 * @example
 * ```ts
 * import { ode_abm } from 'odewasm';
 *
 * // Two-body problem
 * const result = await ode_abm(
 *   (t, [x, y, vx, vy]) => [vx, vy, -x/(x*x+y*y)**1.5, -y/(x*x+y*y)**1.5],
 *   [0, 10],
 *   [1, 0, 0, 1],
 *   { rtol: 1e-8 }
 * );
 * ```
 */
export async function ode_abm(
  fun: ODEFunction,
  t_span: [number, number],
  y0: number[],
  options?: OdeOptions,
): Promise<ODESolverResult> {
  const Module = await loadODEModule();

  const { rtol = 1e-3, atol = 1e-6, one_step = false } = options ?? {};

  const n = y0.length;
  const [t0, tf] = t_span;

  const tValues: number[] = [t0];
  const yValues: number[][] = [y0.slice()];

  const { iflag, nfev } = await solveOdeInternal(
    Module,
    n,
    t0,
    y0,
    tf,
    rtol,
    atol,
    one_step,
    fun,
    tValues,
    yValues,
  );

  // Build result
  const success = iflag === 2;
  const message = getOdeMessage(iflag);

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
    status: iflag,
    message,
    success,
  };
}

function getOdeMessage(iflag: number): string {
  switch (iflag) {
    case 2:
      return "Integration successful.";
    case 3:
      return "Integration successful (one-step mode).";
    case 4:
      return "Error tolerance could not be achieved.";
    case 5:
      return "Equations appear to be stiff.";
    case 6:
      return "Invalid input parameters.";
    default:
      return `Unknown status: ${iflag}`;
  }
}

async function solveOdeInternal(
  Module: ODEModule,
  n: number,
  t0: number,
  y0: number[],
  tf: number,
  rtol: number,
  atol: number,
  one_step: boolean,
  fun: ODEFunction,
  tValues: number[],
  yValues: number[][],
): Promise<{ iflag: number; nfev: number }> {
  const lwork = Module._wasm_ode_work_size(n);
  const liwork = Module._wasm_ode_iwork_size();

  // Allocate memory
  const tPtr = Module._malloc(DOUBLE_SIZE);
  const yPtr = Module._malloc(n * DOUBLE_SIZE);
  const relErrPtr = Module._malloc(DOUBLE_SIZE);
  const absErrPtr = Module._malloc(DOUBLE_SIZE);
  const iflagPtr = Module._malloc(INT_SIZE);
  const workPtr = Module._malloc(lwork * DOUBLE_SIZE);
  const iworkPtr = Module._malloc(liwork * INT_SIZE);

  try {
    // Initialize memory
    Module.HEAPF64[tPtr >> 3] = t0;
    for (let i = 0; i < n; i++) {
      Module.HEAPF64[(yPtr >> 3) + i] = y0[i];
    }
    Module.HEAPF64[relErrPtr >> 3] = rtol;
    Module.HEAPF64[absErrPtr >> 3] = atol;
    Module.HEAP32[iflagPtr >> 2] = 1; // First call

    // Initialize work arrays
    for (let i = 0; i < lwork; i++) {
      Module.HEAPF64[(workPtr >> 3) + i] = 0;
    }
    for (let i = 0; i < liwork; i++) {
      Module.HEAP32[(iworkPtr >> 2) + i] = 0;
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

    let nfev = 0;

    try {
      if (one_step) {
        // One-step mode: integrate step by step
        let t = t0;
        while (t < tf) {
          Module._wasm_ode(
            n,
            yPtr,
            tPtr,
            tf,
            relErrPtr,
            absErrPtr,
            iflagPtr,
            workPtr,
            iworkPtr,
          );

          const iflag = Module.HEAP32[iflagPtr >> 2];
          if (iflag !== 2 && iflag !== -2) {
            // Error
            return { iflag: Math.abs(iflag), nfev };
          }

          t = Module.HEAPF64[tPtr >> 3];
          const yNew: number[] = [];
          for (let i = 0; i < n; i++) {
            yNew.push(Module.HEAPF64[(yPtr >> 3) + i]);
          }
          tValues.push(t);
          yValues.push(yNew);

          // Set flag for continuation
          Module.HEAP32[iflagPtr >> 2] = -2;
        }
        return { iflag: 2, nfev };
      } else {
        // Normal mode: integrate to end
        Module._wasm_ode(
          n,
          yPtr,
          tPtr,
          tf,
          relErrPtr,
          absErrPtr,
          iflagPtr,
          workPtr,
          iworkPtr,
        );

        const iflag = Module.HEAP32[iflagPtr >> 2];
        const tFinal = Module.HEAPF64[tPtr >> 3];
        const yFinal: number[] = [];
        for (let i = 0; i < n; i++) {
          yFinal.push(Module.HEAPF64[(yPtr >> 3) + i]);
        }
        tValues.push(tFinal);
        yValues.push(yFinal);

        // Estimate nfev from work
        nfev = Module.HEAP32[iworkPtr >> 2];

        return { iflag: Math.abs(iflag), nfev };
      }
    } finally {
      Module.removeFunction(fcnPtr);
    }
  } finally {
    Module._free(tPtr);
    Module._free(yPtr);
    Module._free(relErrPtr);
    Module._free(absErrPtr);
    Module._free(iflagPtr);
    Module._free(workPtr);
    Module._free(iworkPtr);
  }
}
