/**
 * RADAU5 - Implicit Runge-Kutta solver for stiff ODEs
 *
 * Solves the initial value problem: dy/dt = f(t, y), y(t0) = y0
 *
 * Uses the Radau IIA method of order 5. Suitable for stiff systems where
 * explicit methods would require impractically small step sizes.
 * Can use numerical or user-supplied Jacobian.
 */

import { loadODEModule } from "../core/loader.js";
import type {
  JacobianFunction,
  ODEFunction,
  ODEModule,
  ODESolverResult,
  Radau5Options,
} from "../core/types.js";
import { STATUS_MESSAGES } from "../core/types.js";

const DOUBLE_SIZE = 8;
const INT_SIZE = 4;

/**
 * Solve an ODE using the RADAU5 method (Implicit Radau IIA order 5).
 *
 * RADAU5 is an implicit Runge-Kutta method suitable for stiff systems.
 * Uses Newton iteration with LU decomposition. Providing an analytical
 * Jacobian improves performance.
 *
 * @param fun - ODE function: f(t, y) => dy/dt
 * @param t_span - Time span [t0, tf]
 * @param y0 - Initial state vector
 * @param options - Solver options
 * @returns Solution result
 *
 * @example
 * ```ts
 * import { radau5 } from 'odewasm';
 *
 * // Stiff Van der Pol oscillator with mu = 1000
 * const mu = 1000;
 * const result = await radau5(
 *   (t, [y0, y1]) => [y1, mu * (1 - y0*y0) * y1 - y0],
 *   [0, 3000],
 *   [2, 0],
 *   {
 *     rtol: 1e-6,
 *     jac: (t, [y0, y1]) => [
 *       [0, 1],
 *       [-2 * mu * y0 * y1 - 1, mu * (1 - y0*y0)]
 *     ]
 *   }
 * );
 * ```
 */
export async function radau5(
  fun: ODEFunction,
  t_span: [number, number],
  y0: number[],
  options?: Radau5Options,
): Promise<ODESolverResult> {
  const Module = await loadODEModule();

  const {
    rtol = 1e-3,
    atol = 1e-6,
    first_step = 0,
    dense_output = false,
    t_eval,
    jac,
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
    iout = 1;
  }

  const tValues: number[] = [];
  const yValues: number[][] = [];

  const { idid, nfev, njev, nlu } = await solveRadau5Internal(
    Module,
    n,
    t0,
    y0,
    tf,
    rtolArray,
    atolArray,
    itol,
    iout,
    first_step,
    fun,
    jac,
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
    njev,
    nlu,
    status: idid,
    message,
    success,
  };
}

/**
 * Internal implementation of RADAU5 solver.
 */
async function solveRadau5Internal(
  Module: ODEModule,
  n: number,
  t0: number,
  y0: number[],
  tf: number,
  rtolArray: number[],
  atolArray: number[],
  itol: number,
  iout: number,
  first_step: number,
  fun: ODEFunction,
  jac: JacobianFunction | undefined,
  tValues: number[],
  yValues: number[][],
): Promise<{ idid: number; nfev: number; njev: number; nlu: number }> {
  // Work array sizes for full Jacobian
  const lwork = Module._wasm_radau5_work_size(n);
  const liwork = Module._wasm_radau5_iwork_size(n);

  // Jacobian settings
  const ijac = jac ? 1 : 0; // 0: numerical Jacobian, 1: user-supplied
  const mljac = n; // Full Jacobian
  const mujac = n;
  const imas = 0; // Identity mass matrix
  const mlmas = 0;
  const mumas = 0;

  // Allocate memory
  const xPtr = Module._malloc(DOUBLE_SIZE);
  const yPtr = Module._malloc(n * DOUBLE_SIZE);
  const rtolPtr = Module._malloc(rtolArray.length * DOUBLE_SIZE);
  const atolPtr = Module._malloc(atolArray.length * DOUBLE_SIZE);
  const workPtr = Module._malloc(lwork * DOUBLE_SIZE);
  const iworkPtr = Module._malloc(liwork * INT_SIZE);
  const ididPtr = Module._malloc(INT_SIZE);

  try {
    // Initialize memory
    Module.HEAPF64[xPtr >> 3] = t0;
    for (let i = 0; i < n; i++) {
      Module.HEAPF64[(yPtr >> 3) + i] = y0[i];
    }
    for (let i = 0; i < rtolArray.length; i++) {
      Module.HEAPF64[(rtolPtr >> 3) + i] = rtolArray[i];
    }
    for (let i = 0; i < atolArray.length; i++) {
      Module.HEAPF64[(atolPtr >> 3) + i] = atolArray[i];
    }

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

    // Create JAC callback if provided
    let jacPtr = 0;
    if (jac) {
      const jacThunk = (
        _nVal: number,
        tVal: number,
        yVal: number,
        dfyVal: number,
        ldfyVal: number,
      ): void => {
        const yArr: number[] = [];
        for (let i = 0; i < n; i++) {
          yArr.push(Module.HEAPF64[(yVal >> 3) + i]);
        }
        const J = jac(tVal, yArr);
        // Column-major storage (Fortran convention)
        for (let j = 0; j < n; j++) {
          for (let i = 0; i < n; i++) {
            Module.HEAPF64[(dfyVal >> 3) + i + j * ldfyVal] = J[i][j];
          }
        }
      };
      jacPtr = Module.addFunction(jacThunk, "vidiii");
      Module._wasm_set_jac_callback(jacPtr);
    }

    // Create SOLOUT callback
    let soloutPtr = 0;
    if (iout > 0) {
      const soloutThunk = (
        _nr: number,
        _told: number,
        tNew: number,
        yVal: number,
        _nVal: number,
        irtrn: number,
      ): void => {
        const yArr: number[] = [];
        for (let i = 0; i < n; i++) {
          yArr.push(Module.HEAPF64[(yVal >> 3) + i]);
        }
        tValues.push(tNew);
        yValues.push(yArr);
        Module.HEAP32[irtrn >> 2] = 0;
      };
      soloutPtr = Module.addFunction(soloutThunk, "viddiii");
      Module._wasm_set_solout_callback(soloutPtr);
    }

    // Call RADAU5
    Module._wasm_radau5(
      n,
      xPtr,
      yPtr,
      tf,
      first_step,
      rtolPtr,
      atolPtr,
      itol,
      ijac,
      mljac,
      mujac,
      imas,
      mlmas,
      mumas,
      iout,
      workPtr,
      lwork,
      iworkPtr,
      liwork,
      ididPtr,
    );

    const idid = Module.HEAP32[ididPtr >> 2];
    const nfev = Module.HEAP32[(iworkPtr >> 2) + 13]; // IWORK(14) = NFCN
    const njev = Module.HEAP32[(iworkPtr >> 2) + 14]; // IWORK(15) = NJAC
    const nlu = Module.HEAP32[(iworkPtr >> 2) + 18]; // IWORK(19) = NDEC

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
    if (jacPtr) Module.removeFunction(jacPtr);
    if (soloutPtr) Module.removeFunction(soloutPtr);

    return { idid, nfev, njev, nlu };
  } finally {
    Module._free(xPtr);
    Module._free(yPtr);
    Module._free(rtolPtr);
    Module._free(atolPtr);
    Module._free(workPtr);
    Module._free(iworkPtr);
    Module._free(ididPtr);
  }
}
