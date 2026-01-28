/**
 * Nelder-Mead simplex optimization — TypeScript binding to WASM.
 *
 * This is a thin interface layer. The numerical algorithm runs in C/WASM.
 */

import { loadWasmModule } from '../wasm-loader.js';
import type { SciWasmModule } from '../wasm-types.js';
import type { ObjectiveFunction, OptimizeResult, NelderMeadOptions, Bounds } from './types.js';

const STATUS_MESSAGES: Record<number, string> = {
  0: 'Optimization terminated successfully.',
  1: 'Maximum number of function evaluations has been exceeded.',
  2: 'Maximum number of iterations has been exceeded.',
};

/**
 * Size of the NelderMeadResult C struct in bytes:
 *   int nfev (4) + padding (4) + int nit (4) + int status (4) + double fval (8) = 24
 * But with alignment: nfev(4) + nit(4) + status(4) + pad(4) + fval(8) = 24
 */
const RESULT_STRUCT_SIZE = 24;

export async function minimizeNelderMead(
  fun: ObjectiveFunction,
  x0: number[],
  options?: NelderMeadOptions,
  bounds?: Bounds
): Promise<OptimizeResult> {
  const wasm: SciWasmModule = await loadWasmModule();
  const n = x0.length;
  const np1 = n + 1;

  const xatol = options?.xatol ?? 1e-4;
  const fatol = options?.fatol ?? 1e-4;
  const maxiter = options?.maxiter ?? 0;
  const maxfev = options?.maxfev ?? 0;
  const adaptive = options?.adaptive ? 1 : 0;
  const hasBounds = bounds ? 1 : 0;

  // Allocate WASM memory
  const x0Ptr = wasm._malloc(n * 8);
  const xOutPtr = wasm._malloc(n * 8);
  const simOutPtr = wasm._malloc(np1 * n * 8);
  const fsimOutPtr = wasm._malloc(np1 * 8);
  const resultPtr = wasm._malloc(RESULT_STRUCT_SIZE);
  const lowerPtr = bounds ? wasm._malloc(n * 8) : 0;
  const upperPtr = bounds ? wasm._malloc(n * 8) : 0;

  const hasInitSimplex = options?.initial_simplex != null;
  let initSimplexPtr = 0;
  if (hasInitSimplex) {
    initSimplexPtr = wasm._malloc(np1 * n * 8);
  }

  // Register JS objective function as C function pointer
  const wrappedFn = wasm.addFunction(
    (xPtr: number, nn: number) => {
      const xArr: number[] = new Array(nn);
      for (let i = 0; i < nn; i++) {
        xArr[i] = wasm.HEAPF64[(xPtr >> 3) + i];
      }
      return fun(xArr);
    },
    'dii' // returns double, takes (int, int) — pointer and length
  );

  try {
    // Write x0 to WASM heap
    for (let i = 0; i < n; i++) {
      wasm.HEAPF64[(x0Ptr >> 3) + i] = x0[i];
    }

    // Write bounds if provided
    if (bounds) {
      for (let i = 0; i < n; i++) {
        wasm.HEAPF64[(lowerPtr >> 3) + i] = bounds.lb[i];
        wasm.HEAPF64[(upperPtr >> 3) + i] = bounds.ub[i];
      }
    }

    // Write initial simplex if provided
    if (hasInitSimplex && options!.initial_simplex) {
      const simplex = options!.initial_simplex!;
      for (let i = 0; i < np1; i++) {
        for (let j = 0; j < n; j++) {
          wasm.HEAPF64[(initSimplexPtr >> 3) + i * n + j] = simplex[i][j];
        }
      }
    }

    // Call the C optimizer
    wasm._nelder_mead_minimize(
      n,
      x0Ptr,
      xOutPtr,
      simOutPtr,
      fsimOutPtr,
      resultPtr,
      xatol,
      fatol,
      maxiter,
      maxfev,
      adaptive,
      hasBounds,
      lowerPtr,
      upperPtr,
      initSimplexPtr,
      wrappedFn
    );

    // Read result struct from WASM memory
    // Struct layout: int nfev, int nit, int status, (pad), double fval
    const nfev = wasm.HEAP32[resultPtr >> 2];
    const nit = wasm.HEAP32[(resultPtr >> 2) + 1];
    const status = wasm.HEAP32[(resultPtr >> 2) + 2];
    // fval is at offset 16 (3 ints = 12 bytes, padded to 16 for double alignment)
    const fval = wasm.HEAPF64[(resultPtr + 16) >> 3];

    // Read solution vector
    const x: number[] = new Array(n);
    for (let i = 0; i < n; i++) {
      x[i] = wasm.HEAPF64[(xOutPtr >> 3) + i];
    }

    // Read final simplex
    const vertices: number[][] = new Array(np1);
    const values: number[] = new Array(np1);
    for (let i = 0; i < np1; i++) {
      vertices[i] = new Array(n);
      for (let j = 0; j < n; j++) {
        vertices[i][j] = wasm.HEAPF64[(simOutPtr >> 3) + i * n + j];
      }
      values[i] = wasm.HEAPF64[(fsimOutPtr >> 3) + i];
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
      final_simplex: { vertices, values },
    };
  } finally {
    // Free all WASM memory
    wasm.removeFunction(wrappedFn);
    wasm._free(x0Ptr);
    wasm._free(xOutPtr);
    wasm._free(simOutPtr);
    wasm._free(fsimOutPtr);
    wasm._free(resultPtr);
    if (lowerPtr) wasm._free(lowerPtr);
    if (upperPtr) wasm._free(upperPtr);
    if (initSimplexPtr) wasm._free(initSimplexPtr);
  }
}
