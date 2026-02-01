/**
 * ODEX - GBS Extrapolation solver based on explicit midpoint rule
 *
 * Solves the initial value problem: dy/dt = f(t, y), y(t0) = y0
 *
 * ODEX is an extrapolation algorithm based on the explicit midpoint rule
 * with automatic stepsize control, order selection, and dense output.
 * It is highly efficient for smooth problems requiring high accuracy
 * as it can achieve very high orders of approximation.
 */

import { loadODEModule } from "../core/loader.js";
import type {
  OdexOptions,
  ODEFunction,
  ODEModule,
  ODESolverResult,
} from "../core/types.js";
import { STATUS_MESSAGES } from "../core/types.js";

const DOUBLE_SIZE = 8;
const INT_SIZE = 4;

/**
 * Solve an ODE using the ODEX method (GBS Extrapolation).
 *
 * ODEX is a GBS extrapolation method based on the explicit midpoint rule.
 * It features automatic order selection (up to very high orders) and is
 * highly efficient for smooth problems requiring high accuracy.
 *
 * @param fun - ODE function: f(t, y) => dy/dt
 * @param t_span - Time span [t0, tf]
 * @param y0 - Initial state vector
 * @param options - Solver options
 * @returns Solution result
 *
 * @example
 * ```ts
 * import { odex } from 'odewasm';
 *
 * // Simple harmonic oscillator: y'' + y = 0
 * // Let y[0] = position, y[1] = velocity
 * const result = await odex(
 *   (t, y) => [y[1], -y[0]],
 *   [0, 10],
 *   [1, 0],
 *   { rtol: 1e-10 }  // High accuracy
 * );
 *
 * console.log(result.t);  // time points
 * console.log(result.y);  // solution values
 * ```
 */
export async function odex(
  fun: ODEFunction,
  t_span: [number, number],
  y0: number[],
  options?: OdexOptions,
): Promise<ODESolverResult> {
  const Module = await loadODEModule();

  const {
    rtol = 1e-3,
    atol = 1e-6,
    max_step = Infinity,
    first_step = 0,
    dense_output = false,
    t_eval,
    max_columns = 9,
    step_sequence,
  } = options ?? {};

  const n = y0.length;
  const [t0, tf] = t_span;

  // Prepare tolerances
  const itol = typeof atol === "number" ? 0 : 1;
  const atolArray = typeof atol === "number" ? [atol] : atol;
  const rtolArray = itol === 1 ? new Array(n).fill(rtol) : [rtol];

  // Determine output mode
  let iout = 0;
  if (dense_output || t_eval) {
    iout = 2;
  }

  const tValues: number[] = [];
  const yValues: number[][] = [];

  const { idid, nfev } = await solveOdexInternal(
    Module,
    n,
    t0,
    y0,
    tf,
    rtolArray,
    atolArray,
    itol,
    iout,
    max_step,
    first_step,
    dense_output,
    t_eval,
    max_columns,
    step_sequence,
    fun,
    tValues,
    yValues,
  );

  // Build result
  const success = idid > 0;
  const message = STATUS_MESSAGES[idid] ?? `Unknown status: ${idid}`;

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
    status: idid,
    message,
    success,
  };
}

/**
 * Internal implementation of ODEX solver.
 */
async function solveOdexInternal(
  Module: ODEModule,
  n: number,
  t0: number,
  y0: number[],
  tf: number,
  rtolArray: number[],
  atolArray: number[],
  itol: number,
  iout: number,
  max_step: number,
  first_step: number,
  dense_output: boolean,
  t_eval: number[] | undefined,
  max_columns: number,
  step_sequence: number | undefined,
  fun: ODEFunction,
  tValues: number[],
  yValues: number[][],
): Promise<{ idid: number; nfev: number }> {
  // Work array sizes
  // ODEX formula: n*(km+5) + 5*km + 20 + (2*km*(km+2)+5)*nrdens
  const km = max_columns;
  const nrdens = dense_output || t_eval ? n : 0;
  const lwork = Module._wasm_odex_work_size(n, nrdens, km);
  const liwork = Module._wasm_odex_iwork_size(nrdens, km);

  // Allocate memory
  const xPtr = Module._malloc(DOUBLE_SIZE);
  const hPtr = Module._malloc(DOUBLE_SIZE);
  const yPtr = Module._malloc(n * DOUBLE_SIZE);
  const rtolPtr = Module._malloc(rtolArray.length * DOUBLE_SIZE);
  const atolPtr = Module._malloc(atolArray.length * DOUBLE_SIZE);
  const workPtr = Module._malloc(lwork * DOUBLE_SIZE);
  const iworkPtr = Module._malloc(liwork * INT_SIZE);
  const ididPtr = Module._malloc(INT_SIZE);

  let tEvalIdx = 0;

  try {
    // Initialize memory
    Module.HEAPF64[xPtr >> 3] = t0;
    Module.HEAPF64[hPtr >> 3] = first_step;
    for (let i = 0; i < n; i++) {
      Module.HEAPF64[(yPtr >> 3) + i] = y0[i];
    }
    for (let i = 0; i < rtolArray.length; i++) {
      Module.HEAPF64[(rtolPtr >> 3) + i] = rtolArray[i];
    }
    for (let i = 0; i < atolArray.length; i++) {
      Module.HEAPF64[(atolPtr >> 3) + i] = atolArray[i];
    }

    // Initialize work arrays to zero
    for (let i = 0; i < lwork; i++) {
      Module.HEAPF64[(workPtr >> 3) + i] = 0;
    }
    for (let i = 0; i < liwork; i++) {
      Module.HEAP32[(iworkPtr >> 2) + i] = 0;
    }

    // Set max step if specified (WORK(2) = HMAX)
    if (max_step !== Infinity) {
      Module.HEAPF64[(workPtr >> 3) + 1] = max_step;
    }

    // Set km if non-default (IWORK(2) = KM)
    if (max_columns !== 9) {
      Module.HEAP32[(iworkPtr >> 2) + 1] = max_columns;
    }

    // Set step sequence if specified (IWORK(3) = NSEQU)
    if (step_sequence !== undefined) {
      Module.HEAP32[(iworkPtr >> 2) + 2] = step_sequence;
    }

    // Set up dense output components (IWORK(8) = NRDENS)
    if (nrdens > 0) {
      Module.HEAP32[(iworkPtr >> 2) + 7] = nrdens;
      // IWORK(21),...,IWORK(NRDENS+20) = component indices (1-based)
      for (let i = 0; i < n; i++) {
        Module.HEAP32[(iworkPtr >> 2) + 20 + i] = i + 1;
      }
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

    // Create SOLOUT callback
    let soloutPtr = 0;
    if (iout > 0) {
      const soloutThunk = (
        _nr: number,
        told: number,
        tNew: number,
        yVal: number,
        _nVal: number,
        irtrn: number,
      ): void => {
        const yArr: number[] = [];
        for (let i = 0; i < n; i++) {
          yArr.push(Module.HEAPF64[(yVal >> 3) + i]);
        }

        if (t_eval) {
          while (tEvalIdx < t_eval.length && t_eval[tEvalIdx] <= tNew) {
            const tInterp = t_eval[tEvalIdx];
            if (tInterp >= told) {
              if (Math.abs(tInterp - tNew) < 1e-10) {
                tValues.push(tNew);
                yValues.push(yArr.slice());
              }
            }
            tEvalIdx++;
          }
        } else {
          tValues.push(tNew);
          yValues.push(yArr);
        }

        Module.HEAP32[irtrn >> 2] = 0;
      };

      soloutPtr = Module.addFunction(soloutThunk, "viddiii");
      Module._wasm_set_solout_callback(soloutPtr);
    }

    // Call ODEX
    Module._wasm_odex(
      n,
      xPtr,
      yPtr,
      tf,
      hPtr,
      rtolPtr,
      atolPtr,
      itol,
      iout,
      workPtr,
      lwork,
      iworkPtr,
      liwork,
      ididPtr,
    );

    const idid = Module.HEAP32[ididPtr >> 2];
    const nfev = Module.HEAP32[(iworkPtr >> 2) + 16]; // IWORK(17) = NFCN

    if (iout === 0) {
      tValues.push(t0);
      yValues.push(y0.slice());

      const xFinal = Module.HEAPF64[xPtr >> 3];
      const yFinal: number[] = [];
      for (let i = 0; i < n; i++) {
        yFinal.push(Module.HEAPF64[(yPtr >> 3) + i]);
      }
      tValues.push(xFinal);
      yValues.push(yFinal);
    }

    Module.removeFunction(fcnPtr);
    if (soloutPtr) {
      Module.removeFunction(soloutPtr);
    }

    return { idid, nfev };
  } finally {
    Module._free(xPtr);
    Module._free(hPtr);
    Module._free(yPtr);
    Module._free(rtolPtr);
    Module._free(atolPtr);
    Module._free(workPtr);
    Module._free(iworkPtr);
    Module._free(ididPtr);
  }
}
