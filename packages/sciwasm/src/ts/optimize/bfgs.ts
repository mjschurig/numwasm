/**
 * BFGS quasi-Newton optimization — TypeScript binding to WASM.
 *
 * This is a thin interface layer. The numerical algorithm runs in C/WASM.
 */

import { loadWasmModule } from '../wasm-loader.js';
import type { SciWasmModule } from '../wasm-types.js';
import type { ObjectiveFunction, OptimizeResult, BFGSOptions, JacobianOption } from './types.js';

const STATUS_MESSAGES: Record<number, string> = {
  0: 'Optimization terminated successfully.',
  1: 'Maximum number of iterations has been exceeded.',
  2: 'Desired error not necessarily achieved due to precision loss.',
  3: 'NaN result encountered.',
};

/**
 * Size of the BFGSResult C struct in bytes:
 *   int nfev (4) + int njev (4) + int nit (4) + int status (4) + double fval (8) = 24
 */
const RESULT_STRUCT_SIZE = 24;

export async function minimizeBFGS(
  fun: ObjectiveFunction,
  x0: number[],
  options?: BFGSOptions,
  jac?: JacobianOption
): Promise<OptimizeResult> {
  const wasm: SciWasmModule = await loadWasmModule();
  const n = x0.length;

  const gtol = options?.gtol ?? 1e-5;
  const maxiter = options?.maxiter ?? 0;
  const c1 = options?.c1 ?? 1e-4;
  const c2 = options?.c2 ?? 0.9;
  const eps = options?.eps ?? 0;

  // Determine if we have an analytical Jacobian
  const hasJac = typeof jac === 'function' ? 1 : 0;
  const jacFn = typeof jac === 'function' ? jac : null;

  // Allocate WASM memory
  const x0Ptr = wasm._malloc(n * 8);
  const xOutPtr = wasm._malloc(n * 8);
  const hessInvOutPtr = wasm._malloc(n * n * 8);
  const jacOutPtr = wasm._malloc(n * 8);
  const resultPtr = wasm._malloc(RESULT_STRUCT_SIZE);

  // Register JS objective function as C function pointer
  const wrappedFn = wasm.addFunction(
    (xPtr: number, nn: number) => {
      const xArr: number[] = new Array(nn);
      for (let i = 0; i < nn; i++) {
        xArr[i] = wasm.HEAPF64[(xPtr >> 3) + i];
      }
      return fun(xArr);
    },
    'dii' // returns double, takes (int, int)
  );

  // Register JS gradient function as C function pointer
  // Signature: void grad(const double* x, double* grad_out, int n)
  const wrappedGrad = wasm.addFunction(
    (xPtr: number, gradOutPtr: number, nn: number) => {
      if (!jacFn) return; // Should not be called if has_jac=0

      const xArr: number[] = new Array(nn);
      for (let i = 0; i < nn; i++) {
        xArr[i] = wasm.HEAPF64[(xPtr >> 3) + i];
      }

      const grad = jacFn(xArr);
      for (let i = 0; i < nn; i++) {
        wasm.HEAPF64[(gradOutPtr >> 3) + i] = grad[i];
      }
    },
    'viii' // returns void, takes (int, int, int)
  );

  try {
    // Write x0 to WASM heap
    for (let i = 0; i < n; i++) {
      wasm.HEAPF64[(x0Ptr >> 3) + i] = x0[i];
    }

    // Call the C optimizer
    wasm._bfgs_minimize(
      n,
      x0Ptr,
      xOutPtr,
      hessInvOutPtr,
      jacOutPtr,
      resultPtr,
      gtol,
      maxiter,
      c1,
      c2,
      hasJac,
      eps,
      wrappedFn,
      wrappedGrad
    );

    // Read result struct from WASM memory
    // Struct layout: int nfev, int njev, int nit, int status, double fval
    const nfev = wasm.HEAP32[resultPtr >> 2];
    const njev = wasm.HEAP32[(resultPtr >> 2) + 1];
    const nit = wasm.HEAP32[(resultPtr >> 2) + 2];
    const status = wasm.HEAP32[(resultPtr >> 2) + 3];
    // fval is at offset 16 (4 ints = 16 bytes, aligned for double)
    const fval = wasm.HEAPF64[(resultPtr + 16) >> 3];

    // Read solution vector
    const x: number[] = new Array(n);
    for (let i = 0; i < n; i++) {
      x[i] = wasm.HEAPF64[(xOutPtr >> 3) + i];
    }

    // Read final Jacobian
    const jacResult: number[] = new Array(n);
    for (let i = 0; i < n; i++) {
      jacResult[i] = wasm.HEAPF64[(jacOutPtr >> 3) + i];
    }

    // Read inverse Hessian approximation
    const hessInv: number[][] = new Array(n);
    for (let i = 0; i < n; i++) {
      hessInv[i] = new Array(n);
      for (let j = 0; j < n; j++) {
        hessInv[i][j] = wasm.HEAPF64[(hessInvOutPtr >> 3) + i * n + j];
      }
    }

    const message = STATUS_MESSAGES[status] ?? `Unknown status: ${status}`;

    return {
      x,
      success: status === 0,
      status,
      message,
      fun: fval,
      nfev,
      nit,
      jac: jacResult,
      hess_inv: hessInv,
      njev,
    };
  } finally {
    // Free all WASM memory
    wasm.removeFunction(wrappedFn);
    wasm.removeFunction(wrappedGrad);
    wasm._free(x0Ptr);
    wasm._free(xOutPtr);
    wasm._free(hessInvOutPtr);
    wasm._free(jacOutPtr);
    wasm._free(resultPtr);
  }
}
