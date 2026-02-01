/**
 * ZVODE - VODE for complex-valued ODEs
 *
 * ZVODE extends VODE to handle complex-valued differential equations.
 * Useful for quantum mechanics, wave equations, and other applications
 * involving complex state variables.
 */

import { loadODEModule } from "../core/loader.js";
import type { ODEModule, ODESolverResult } from "../core/types.js";

const DOUBLE_SIZE = 8;
const INT_SIZE = 4;

/**
 * Complex number type.
 */
export interface Complex {
  re: number;
  im: number;
}

/**
 * Complex ODE function type.
 */
export type ComplexODEFunction = (t: number, y: Complex[]) => Complex[];

/**
 * Method flag for ZVODE.
 * - 10: Adams method (non-stiff)
 * - 21: BDF method (stiff), full Jacobian (user-supplied)
 * - 22: BDF method (stiff), full Jacobian (internal numerical)
 */
export type ZvodeMethodFlag = 10 | 21 | 22;

/**
 * Options for the ZVODE solver.
 */
export interface ZvodeOptions {
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
   * Method flag.
   * - 10: Adams (non-stiff)
   * - 22: BDF with internal Jacobian (stiff)
   * @default 10
   */
  mf?: ZvodeMethodFlag;

  /**
   * Maximum number of steps.
   * @default 500
   */
  max_steps?: number;
}

/**
 * Result from ZVODE solver with complex values.
 */
export interface ZvodeSolverResult extends ODESolverResult {
  /**
   * Solution values as complex arrays.
   * y_complex[i] is the i-th component of the complex state vector.
   */
  y_complex: Complex[][];
}

/**
 * Solve a complex-valued ODE using ZVODE.
 *
 * ZVODE handles ODEs where both the state and derivatives are complex.
 * This is common in quantum mechanics and wave physics.
 *
 * @param fun - Complex ODE function: f(t, y) => dy/dt
 * @param t_span - Time span [t0, tf]
 * @param y0 - Initial state vector (complex)
 * @param options - Solver options
 * @returns Solution result with complex values
 *
 * @example
 * ```ts
 * import { zvode, Complex } from 'odewasm';
 *
 * // Schrodinger equation: i * dy/dt = H * y
 * // For a simple harmonic oscillator: dy/dt = -i * omega * y
 * const omega = 1;
 * const result = await zvode(
 *   (t, y) => [{ re: omega * y[0].im, im: -omega * y[0].re }],
 *   [0, 10],
 *   [{ re: 1, im: 0 }],
 *   { mf: 10 }
 * );
 * ```
 */
export async function zvode(
  fun: ComplexODEFunction,
  t_span: [number, number],
  y0: Complex[],
  options?: ZvodeOptions,
): Promise<ZvodeSolverResult> {
  const Module = await loadODEModule();

  const { rtol = 1e-3, atol = 1e-6, mf = 10, max_steps = 500 } = options ?? {};

  const n = y0.length;
  const [t0, tf] = t_span;

  const tValues: number[] = [t0];
  const yComplexValues: Complex[][] = [y0.slice()];

  const { istate, nfev, njev } = await solveZvodeInternal(
    Module,
    n,
    t0,
    y0,
    tf,
    rtol,
    atol,
    mf,
    max_steps,
    fun,
    tValues,
    yComplexValues,
  );

  // Build result
  const success = istate === 2;
  const message = getZvodeMessage(istate);

  const t = new Float64Array(tValues);

  // Create both Float64Array (real parts) and complex arrays
  const y: Float64Array[] = [];
  const y_complex: Complex[][] = [];

  for (let i = 0; i < n; i++) {
    const yi_real = new Float64Array(tValues.length);
    for (let j = 0; j < tValues.length; j++) {
      yi_real[j] = yComplexValues[j][i].re;
    }
    y.push(yi_real);
  }

  // Transpose yComplexValues from time-major to component-major
  for (let i = 0; i < n; i++) {
    const comp: Complex[] = [];
    for (let j = 0; j < tValues.length; j++) {
      comp.push(yComplexValues[j][i]);
    }
    y_complex.push(comp);
  }

  return {
    t,
    y,
    y_complex,
    nfev,
    njev,
    nlu: 0,
    status: istate,
    message,
    success,
  };
}

function getZvodeMessage(istate: number): string {
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
      return "Repeated convergence failures.";
    case -6:
      return "Error weight became zero.";
    default:
      return `Unknown status: ${istate}`;
  }
}

async function solveZvodeInternal(
  Module: ODEModule,
  n: number,
  t0: number,
  y0: Complex[],
  tf: number,
  rtol: number,
  atol: number,
  mf: ZvodeMethodFlag,
  max_steps: number,
  fun: ComplexODEFunction,
  tValues: number[],
  yComplexValues: Complex[][],
): Promise<{ istate: number; nfev: number; njev: number }> {
  const lzw = Module._wasm_zvode_zwork_size(n, mf);
  const lrw = Module._wasm_zvode_rwork_size(n, mf);
  const liw = Module._wasm_vode_iwork_size(n, mf);

  // Allocate memory
  const tPtr = Module._malloc(DOUBLE_SIZE);
  const yPtr = Module._malloc(2 * n * DOUBLE_SIZE); // Interleaved real/imag
  const rtolPtr = Module._malloc(DOUBLE_SIZE);
  const atolPtr = Module._malloc(DOUBLE_SIZE);
  const istatePtr = Module._malloc(INT_SIZE);
  const zworkPtr = Module._malloc(2 * lzw * DOUBLE_SIZE);
  const rworkPtr = Module._malloc(lrw * DOUBLE_SIZE);
  const iworkPtr = Module._malloc(liw * INT_SIZE);

  try {
    // Initialize memory
    Module.HEAPF64[tPtr >> 3] = t0;
    // Store y as interleaved real/imag pairs
    for (let i = 0; i < n; i++) {
      Module.HEAPF64[(yPtr >> 3) + 2 * i] = y0[i].re;
      Module.HEAPF64[(yPtr >> 3) + 2 * i + 1] = y0[i].im;
    }
    Module.HEAPF64[rtolPtr >> 3] = rtol;
    Module.HEAPF64[atolPtr >> 3] = atol;
    Module.HEAP32[istatePtr >> 2] = 1; // First call

    // Initialize work arrays
    for (let i = 0; i < 2 * lzw; i++) {
      Module.HEAPF64[(zworkPtr >> 3) + i] = 0;
    }
    for (let i = 0; i < lrw; i++) {
      Module.HEAPF64[(rworkPtr >> 3) + i] = 0;
    }
    for (let i = 0; i < liw; i++) {
      Module.HEAP32[(iworkPtr >> 2) + i] = 0;
    }

    // Set max steps in iwork
    Module.HEAP32[(iworkPtr >> 2) + 5] = max_steps;

    // Create FCN callback for complex values
    // The callback receives interleaved complex arrays
    const fcnThunk = (
      _neq: number,
      tVal: number,
      yVal: number,
      ydotVal: number,
    ): void => {
      // Read complex y
      const yArr: Complex[] = [];
      for (let i = 0; i < n; i++) {
        yArr.push({
          re: Module.HEAPF64[(yVal >> 3) + 2 * i],
          im: Module.HEAPF64[(yVal >> 3) + 2 * i + 1],
        });
      }
      // Call user function
      const f = fun(tVal, yArr);
      // Write complex ydot
      for (let i = 0; i < n; i++) {
        Module.HEAPF64[(ydotVal >> 3) + 2 * i] = f[i].re;
        Module.HEAPF64[(ydotVal >> 3) + 2 * i + 1] = f[i].im;
      }
    };

    const fcnPtr = Module.addFunction(fcnThunk, "vidii");
    Module._wasm_set_fcn_callback(fcnPtr);

    try {
      const itol = 1; // Scalar tolerances
      const itask = 1; // Normal mode
      const iopt = 0; // No optional inputs

      Module._wasm_zvode(
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
        zworkPtr,
        lzw,
        rworkPtr,
        lrw,
        iworkPtr,
        liw,
        mf,
      );

      const istate = Module.HEAP32[istatePtr >> 2];
      const tFinal = Module.HEAPF64[tPtr >> 3];

      // Read final complex y
      const yFinal: Complex[] = [];
      for (let i = 0; i < n; i++) {
        yFinal.push({
          re: Module.HEAPF64[(yPtr >> 3) + 2 * i],
          im: Module.HEAPF64[(yPtr >> 3) + 2 * i + 1],
        });
      }
      tValues.push(tFinal);
      yComplexValues.push(yFinal);

      // Get statistics from iwork
      const nfev = Module.HEAP32[(iworkPtr >> 2) + 11];
      const njev = Module.HEAP32[(iworkPtr >> 2) + 12];

      return { istate, nfev, njev };
    } finally {
      Module.removeFunction(fcnPtr);
    }
  } finally {
    Module._free(tPtr);
    Module._free(yPtr);
    Module._free(rtolPtr);
    Module._free(atolPtr);
    Module._free(istatePtr);
    Module._free(zworkPtr);
    Module._free(rworkPtr);
    Module._free(iworkPtr);
  }
}
