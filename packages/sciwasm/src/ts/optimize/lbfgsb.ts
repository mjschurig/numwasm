/**
 * L-BFGS-B optimization — TypeScript binding to WASM.
 *
 * Uses reverse communication: TypeScript calls setulb() in a loop,
 * evaluating f(x) and g(x) in JavaScript when requested by the C code.
 */

import { loadWasmModule } from '../wasm-loader.js';
import type { SciWasmModule } from '../wasm-types.js';
import type { ObjectiveFunction, OptimizeResult, LBFGSBOptions, Bounds, JacobianOption } from './types.js';

/** Task status codes from lbfgsb.c */
const TASK_START = 0;
const TASK_NEW_X = 1;
const TASK_FG = 3;
const TASK_CONVERGENCE = 4;
const TASK_STOP = 5;
// Reserved task codes for future use: TASK_WARNING=6, TASK_ERROR=7, TASK_ABNORMAL=8

const STATUS_MESSAGES: Record<number, string> = {
  0: 'Optimization terminated successfully.',
  1: 'Maximum number of iterations or function evaluations exceeded.',
  2: 'Abnormal termination or warning in L-BFGS-B.',
};

/**
 * Compute finite-difference gradient.
 */
function finiteDiffGradient(
  fun: ObjectiveFunction,
  x: number[],
  f0: number,
  eps: number
): number[] {
  const n = x.length;
  const grad = new Array(n);
  const xp = x.slice();
  for (let i = 0; i < n; i++) {
    const xi = xp[i];
    xp[i] = xi + eps;
    const fp = fun(xp);
    grad[i] = (fp - f0) / eps;
    xp[i] = xi;
  }
  return grad;
}

export async function minimizeLBFGSB(
  fun: ObjectiveFunction,
  x0: number[],
  options?: LBFGSBOptions,
  bounds?: Bounds,
  jac?: JacobianOption
): Promise<OptimizeResult> {
  const wasm: SciWasmModule = await loadWasmModule();
  const n = x0.length;

  const maxcor = options?.maxcor ?? 10;
  const m = maxcor;
  const ftol = options?.ftol ?? 2.2204460492503131e-09;
  const pgtol = options?.gtol ?? 1e-5;
  const eps = options?.eps ?? 1e-8;
  const maxfun = options?.maxfun ?? 15000;
  const maxiter = options?.maxiter ?? 15000;
  const maxls = options?.maxls ?? 20;

  // Convert ftol to factr (same as scipy)
  const machEps = 2.220446049250313e-16;
  const factr = ftol / machEps;

  // Determine if we have an analytical Jacobian
  const hasJac = typeof jac === 'function';
  const jacFn = hasJac ? jac as (x: number[]) => number[] : null;

  // Set up bounds arrays (nbd format)
  const nbd = new Int32Array(n);
  const lower = new Float64Array(n);
  const upper = new Float64Array(n);

  if (bounds) {
    for (let i = 0; i < n; i++) {
      const lb = bounds.lb[i];
      const ub = bounds.ub[i];
      const lInf = !isFinite(lb) || lb === -Infinity;
      const uInf = !isFinite(ub) || ub === Infinity;

      if (lInf && uInf) {
        nbd[i] = 0; // unbounded
      } else if (!lInf && uInf) {
        nbd[i] = 1; // lower bound only
        lower[i] = lb;
      } else if (!lInf && !uInf) {
        nbd[i] = 2; // both bounds
        lower[i] = lb;
        upper[i] = ub;
      } else {
        nbd[i] = 3; // upper bound only
        upper[i] = ub;
      }
    }
  }

  // Clip x0 to bounds
  const x0Clipped = x0.slice();
  if (bounds) {
    for (let i = 0; i < n; i++) {
      if (nbd[i] === 1 || nbd[i] === 2) {
        x0Clipped[i] = Math.max(x0Clipped[i], lower[i]);
      }
      if (nbd[i] === 2 || nbd[i] === 3) {
        x0Clipped[i] = Math.min(x0Clipped[i], upper[i]);
      }
    }
  }

  // Allocate WASM memory
  const waSize = 2 * m * n + 5 * n + 11 * m * m + 8 * m;
  const xPtr = wasm._malloc(n * 8);
  const lPtr = wasm._malloc(n * 8);
  const uPtr = wasm._malloc(n * 8);
  const nbdPtr = wasm._malloc(n * 4);
  const fPtr = wasm._malloc(8);
  const gPtr = wasm._malloc(n * 8);
  const waPtr = wasm._malloc(waSize * 8);
  const iwaPtr = wasm._malloc(3 * n * 4);
  const taskPtr = wasm._malloc(2 * 4);
  const lsavePtr = wasm._malloc(4 * 4);
  const isavePtr = wasm._malloc(44 * 4);
  const dsavePtr = wasm._malloc(29 * 8);
  const lnTaskPtr = wasm._malloc(2 * 4);

  try {
    // Write x0 to WASM heap
    for (let i = 0; i < n; i++) {
      wasm.HEAPF64[(xPtr >> 3) + i] = x0Clipped[i];
      wasm.HEAPF64[(lPtr >> 3) + i] = lower[i];
      wasm.HEAPF64[(uPtr >> 3) + i] = upper[i];
      wasm.HEAP32[(nbdPtr >> 2) + i] = nbd[i];
    }

    // Initialize work arrays to zero
    for (let i = 0; i < waSize; i++) {
      wasm.HEAPF64[(waPtr >> 3) + i] = 0.0;
    }
    for (let i = 0; i < 3 * n; i++) {
      wasm.HEAP32[(iwaPtr >> 2) + i] = 0;
    }

    // Set initial task to START (0)
    wasm.HEAP32[taskPtr >> 2] = TASK_START;
    wasm.HEAP32[(taskPtr >> 2) + 1] = 0;
    wasm.HEAP32[lnTaskPtr >> 2] = 0;
    wasm.HEAP32[(lnTaskPtr >> 2) + 1] = 0;

    // Initialize save arrays
    for (let i = 0; i < 4; i++) wasm.HEAP32[(lsavePtr >> 2) + i] = 0;
    for (let i = 0; i < 44; i++) wasm.HEAP32[(isavePtr >> 2) + i] = 0;
    for (let i = 0; i < 29; i++) wasm.HEAPF64[(dsavePtr >> 3) + i] = 0.0;

    // Initialize f to 0
    wasm.HEAPF64[fPtr >> 3] = 0.0;
    for (let i = 0; i < n; i++) {
      wasm.HEAPF64[(gPtr >> 3) + i] = 0.0;
    }

    let nfev = 0;
    let njev = 0;
    let nIterations = 0;
    let fval = 0.0;

    // Reverse communication loop
    const maxLoopIter = maxfun + maxiter + 1000; // safety limit
    for (let loop = 0; loop < maxLoopIter; loop++) {
      // Call setulb
      wasm._setulb(
        n, m, xPtr, lPtr, uPtr, nbdPtr, fPtr, gPtr,
        factr, pgtol, waPtr, iwaPtr, taskPtr, lsavePtr,
        isavePtr, dsavePtr, maxls, lnTaskPtr
      );

      const task0 = wasm.HEAP32[taskPtr >> 2];

      if (task0 === TASK_FG) {
        // C code wants f and g at current x
        const x: number[] = new Array(n);
        for (let i = 0; i < n; i++) {
          x[i] = wasm.HEAPF64[(xPtr >> 3) + i];
        }

        fval = fun(x);
        nfev++;
        wasm.HEAPF64[fPtr >> 3] = fval;

        let grad: number[];
        if (jacFn) {
          grad = jacFn(x);
          njev++;
        } else {
          grad = finiteDiffGradient(fun, x, fval, eps);
          nfev += n; // finite diff evaluates f n more times
          njev++;
        }

        for (let i = 0; i < n; i++) {
          wasm.HEAPF64[(gPtr >> 3) + i] = grad[i];
        }

      } else if (task0 === TASK_NEW_X) {
        // New iteration completed
        nIterations++;

        if (nIterations >= maxiter) {
          wasm.HEAP32[taskPtr >> 2] = TASK_STOP;
          wasm.HEAP32[(taskPtr >> 2) + 1] = 504;
        } else if (nfev > maxfun) {
          wasm.HEAP32[taskPtr >> 2] = TASK_STOP;
          wasm.HEAP32[(taskPtr >> 2) + 1] = 502;
        }

      } else {
        // Convergence, error, warning, or stop
        break;
      }
    }

    // Read final results
    const task0 = wasm.HEAP32[taskPtr >> 2];
    const x: number[] = new Array(n);
    for (let i = 0; i < n; i++) {
      x[i] = wasm.HEAPF64[(xPtr >> 3) + i];
    }
    fval = wasm.HEAPF64[fPtr >> 3];

    const jacResult: number[] = new Array(n);
    for (let i = 0; i < n; i++) {
      jacResult[i] = wasm.HEAPF64[(gPtr >> 3) + i];
    }

    let warnflag: number;
    if (task0 === TASK_CONVERGENCE) {
      warnflag = 0;
    } else if (nfev > maxfun || nIterations >= maxiter) {
      warnflag = 1;
    } else {
      warnflag = 2;
    }

    const message = STATUS_MESSAGES[warnflag] ?? `Unknown status: ${warnflag}`;

    return {
      x,
      success: warnflag === 0,
      status: warnflag,
      message,
      fun: fval,
      nfev,
      nit: nIterations,
      jac: jacResult,
      njev,
    };
  } finally {
    wasm._free(xPtr);
    wasm._free(lPtr);
    wasm._free(uPtr);
    wasm._free(nbdPtr);
    wasm._free(fPtr);
    wasm._free(gPtr);
    wasm._free(waPtr);
    wasm._free(iwaPtr);
    wasm._free(taskPtr);
    wasm._free(lsavePtr);
    wasm._free(isavePtr);
    wasm._free(dsavePtr);
    wasm._free(lnTaskPtr);
  }
}
