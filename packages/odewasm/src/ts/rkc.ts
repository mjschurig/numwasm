/**
 * RKC - Runge-Kutta-Chebyshev method
 *
 * RKC is an explicit stabilized method designed for mildly stiff problems
 * arising from parabolic PDEs. It automatically adjusts the number of
 * stages based on the spectral radius of the Jacobian.
 */

import { loadODEModule } from "./loader.js";
import type { ODEFunction, ODEModule, ODESolverResult } from "./types.js";

const DOUBLE_SIZE = 8;
const INT_SIZE = 4;

/**
 * Spectral radius estimation function.
 *
 * Returns an upper bound for the spectral radius of df/dy.
 * If not provided, RKC will estimate it internally (which adds
 * some computational cost).
 *
 * @param n - System dimension
 * @param t - Current time
 * @param y - Current state
 * @returns Estimated spectral radius (should be >= actual value)
 */
export type SpcradFunction = (n: number, t: number, y: number[]) => number;

/**
 * Options for the RKC solver.
 */
export interface RkcOptions {
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
   * Spectral radius estimation function.
   * If provided, RKC uses it to determine the number of stages.
   * If not provided, RKC estimates it internally (slower but automatic).
   */
  spcrad?: SpcradFunction;

  /**
   * Whether the Jacobian is constant (doesn't depend on y).
   * When true, the spectral radius is computed once and reused.
   * @default false
   */
  constant_jacobian?: boolean;
}

/**
 * Solve an ODE using RKC (Runge-Kutta-Chebyshev).
 *
 * RKC is designed for mildly stiff problems arising from spatial
 * discretization of parabolic PDEs (e.g., heat equation, diffusion).
 * It is explicit but uses a variable number of stages to achieve
 * stability for larger step sizes.
 *
 * @param fun - ODE function: f(t, y) => dy/dt
 * @param t_span - Time span [t0, tf]
 * @param y0 - Initial state vector
 * @param options - Solver options
 * @returns Solution result
 *
 * @example
 * ```ts
 * import { rkc } from 'odewasm';
 *
 * // Heat equation discretization
 * const n = 50;
 * const dx = 1 / (n + 1);
 * const result = await rkc(
 *   (t, y) => {
 *     const f = new Array(n).fill(0);
 *     for (let i = 0; i < n; i++) {
 *       const left = i === 0 ? 0 : y[i - 1];
 *       const right = i === n - 1 ? 0 : y[i + 1];
 *       f[i] = (left - 2 * y[i] + right) / (dx * dx);
 *     }
 *     return f;
 *   },
 *   [0, 0.1],
 *   new Array(n).fill(1),
 *   { rtol: 1e-4 }
 * );
 * ```
 */
export async function rkc(
  fun: ODEFunction,
  t_span: [number, number],
  y0: number[],
  options?: RkcOptions,
): Promise<ODESolverResult> {
  const Module = await loadODEModule();

  const {
    rtol = 1e-3,
    atol = 1e-6,
    spcrad,
    constant_jacobian = false,
  } = options ?? {};

  const n = y0.length;
  const [t0, tf] = t_span;

  // Determine tolerance mode
  const itol = typeof atol === "number" ? 0 : 1;
  const atolArray = typeof atol === "number" ? [atol] : atol;

  const tValues: number[] = [t0];
  const yValues: number[][] = [y0.slice()];

  const { idid, nfev } = await solveRkcInternal(
    Module,
    n,
    t0,
    y0,
    tf,
    rtol,
    atolArray,
    itol,
    spcrad,
    constant_jacobian,
    fun,
    tValues,
    yValues,
  );

  // Build result
  const success = idid === 1;
  const message = getRkcMessage(idid);

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

function getRkcMessage(idid: number): string {
  switch (idid) {
    case 1:
      return "Integration successful.";
    case 2:
      return "Interrupted by user (t is intermediate).";
    case -1:
      return "Invalid input parameters.";
    case -2:
      return "Spectral radius estimation failed.";
    case -3:
      return "Step size too small.";
    case -4:
      return "Too many function evaluations.";
    default:
      return `Unknown status: ${idid}`;
  }
}

async function solveRkcInternal(
  Module: ODEModule,
  n: number,
  t0: number,
  y0: number[],
  tf: number,
  rtol: number,
  atolArray: number[],
  itol: number,
  spcrad: SpcradFunction | undefined,
  constant_jacobian: boolean,
  fun: ODEFunction,
  tValues: number[],
  yValues: number[][],
): Promise<{ idid: number; nfev: number }> {
  // Work size depends on whether we provide spectral radius
  const use_internal_spcrad = spcrad === undefined ? 1 : 0;
  const lwork = Module._wasm_rkc_work_size(n, use_internal_spcrad);

  // Allocate memory
  const tPtr = Module._malloc(DOUBLE_SIZE);
  const yPtr = Module._malloc(n * DOUBLE_SIZE);
  const atolPtr = Module._malloc(atolArray.length * DOUBLE_SIZE);
  const infoPtr = Module._malloc(4 * INT_SIZE);
  const workPtr = Module._malloc(lwork * DOUBLE_SIZE);
  const ididPtr = Module._malloc(INT_SIZE);

  try {
    // Initialize memory
    Module.HEAPF64[tPtr >> 3] = t0;
    for (let i = 0; i < n; i++) {
      Module.HEAPF64[(yPtr >> 3) + i] = y0[i];
    }
    for (let i = 0; i < atolArray.length; i++) {
      Module.HEAPF64[(atolPtr >> 3) + i] = atolArray[i];
    }

    // Initialize info array
    // info[0] = 0 for first call
    // info[1] = 0 for scalar atol, 1 for array atol
    // info[2] = 0 if spectral radius provided, 1 if RKC should estimate
    // info[3] = 0 if Jacobian depends on y, 1 if constant
    Module.HEAP32[(infoPtr >> 2) + 0] = 0;
    Module.HEAP32[(infoPtr >> 2) + 1] = itol;
    Module.HEAP32[(infoPtr >> 2) + 2] = use_internal_spcrad;
    Module.HEAP32[(infoPtr >> 2) + 3] = constant_jacobian ? 1 : 0;

    // Initialize work array
    for (let i = 0; i < lwork; i++) {
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

    // Create SPCRAD callback if provided
    let spcradPtr = 0;
    if (spcrad) {
      const spcradThunk = (
        nVal: number,
        tVal: number,
        yVal: number,
      ): number => {
        const yArr: number[] = [];
        for (let i = 0; i < n; i++) {
          yArr.push(Module.HEAPF64[(yVal >> 3) + i]);
        }
        return spcrad(nVal, tVal, yArr);
      };

      spcradPtr = Module.addFunction(spcradThunk, "didi");
      Module._wasm_set_spcrad_callback(spcradPtr);
    }

    try {
      // Call RKC
      Module._wasm_rkc(
        n,
        yPtr,
        tPtr,
        tf,
        rtol,
        atolPtr,
        infoPtr,
        workPtr,
        ididPtr,
      );

      const idid = Module.HEAP32[ididPtr >> 2];
      const tFinal = Module.HEAPF64[tPtr >> 3];
      const yFinal: number[] = [];
      for (let i = 0; i < n; i++) {
        yFinal.push(Module.HEAPF64[(yPtr >> 3) + i]);
      }
      tValues.push(tFinal);
      yValues.push(yFinal);

      // nfev is stored in work[4] (0-indexed: work(5) in Fortran)
      const nfev = Math.round(Module.HEAPF64[(workPtr >> 3) + 4]);

      return { idid, nfev };
    } finally {
      Module.removeFunction(fcnPtr);
      if (spcradPtr) {
        Module.removeFunction(spcradPtr);
      }
    }
  } finally {
    Module._free(tPtr);
    Module._free(yPtr);
    Module._free(atolPtr);
    Module._free(infoPtr);
    Module._free(workPtr);
    Module._free(ididPtr);
  }
}
