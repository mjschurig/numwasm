/**
 * Adaptive quadrature integration using QUADPACK via WebAssembly.
 * Mirrors scipy.integrate.quad.
 */

import { loadWasmModule } from '../wasm-loader.js';
import type { SciWasmModule } from '../wasm-types.js';
import type {
  IntegrandFunction,
  QuadOptions,
  QuadInfoDict,
  QuadResult,
  QuadFullResult,
  QuadFullResultWithMessage,
} from './types.js';

/**
 * Error messages corresponding to QUADPACK ier codes.
 */
const QUAD_MESSAGES: Record<number, string> = {
  1: "The maximum number of subdivisions has been achieved.\n  "
    + "If increasing the limit yields no improvement it is advised to "
    + "analyze\n  the integrand in order to determine the difficulties. "
    + "If the position of a\n  local difficulty can be determined "
    + "(singularity, discontinuity) one will\n  probably gain from "
    + "splitting up the interval and calling the integrator\n  on the "
    + "subranges. Perhaps a special-purpose integrator should be used.",
  2: "The occurrence of roundoff error is detected, which prevents\n  "
    + "the requested tolerance from being achieved. "
    + "The error may be\n  underestimated.",
  3: "Extremely bad integrand behavior occurs at some points of the\n  "
    + "integration interval.",
  4: "The algorithm does not converge. Roundoff error is detected\n  "
    + "in the extrapolation table. It is assumed that the requested "
    + "tolerance\n  cannot be achieved, and that the returned result "
    + "is the best which can be obtained.",
  5: "The integral is probably divergent, or slowly convergent.",
  6: "The input is invalid.",
  7: "Abnormal termination of the routine. The estimates for result\n  "
    + "and error are less reliable. It is assumed that the requested "
    + "accuracy\n  has not been achieved.",
};

/**
 * Compute a definite integral using adaptive quadrature.
 *
 * Integrates `func` from `a` to `b` (possibly infinite interval) using
 * the QUADPACK library compiled to WebAssembly.
 *
 * @param func - The integrand function f(x, ...args).
 * @param a - Lower limit of integration. Use -Infinity for -∞.
 * @param b - Upper limit of integration. Use Infinity for +∞.
 * @param options - Integration options.
 * @returns [result, abserr] or [result, abserr, infodict] if fullOutput.
 *
 * @example
 * ```ts
 * import { integrate } from 'sciwasm';
 *
 * // Simple integral of x^2 from 0 to 1
 * const [result, error] = await integrate.quad(x => x * x, 0, 1);
 * // result ≈ 0.333333...
 *
 * // Infinite interval: integral of exp(-x) from 0 to +∞
 * const [val, err] = await integrate.quad(x => Math.exp(-x), 0, Infinity);
 * // val ≈ 1.0
 *
 * // With extra arguments
 * const [v, e] = await integrate.quad((x, n, z) => Math.cos(n*x - z*Math.sin(x)) / Math.PI, 0, Math.PI, { args: [2, 1.8] });
 * ```
 */
export async function quad(
  func: IntegrandFunction,
  a: number,
  b: number,
  options?: QuadOptions,
): Promise<QuadResult | QuadFullResult | QuadFullResultWithMessage> {
  const {
    args: rawArgs,
    fullOutput = false,
    epsabs = 1.49e-8,
    epsrel = 1.49e-8,
    limit = 50,
    complexFunc = false,
  } = options ?? {};

  // Normalize args to array
  const args: number[] = rawArgs == null
    ? []
    : typeof rawArgs === 'number'
      ? [rawArgs]
      : rawArgs;

  // Shortcut for empty interval
  if (a === b) {
    if (!fullOutput) {
      return [0, 0] as QuadResult;
    }
    const infodict: QuadInfoDict = {
      neval: 0,
      last: 0,
      alist: new Float64Array(limit).fill(NaN),
      blist: new Float64Array(limit).fill(NaN),
      rlist: new Float64Array(limit),
      elist: new Float64Array(limit),
      iord: new Int32Array(limit),
    };
    return [0, 0, infodict] as QuadFullResult;
  }

  // Handle b < a by flipping and negating
  let flip = false;
  let lo = a;
  let hi = b;
  if (b < a) {
    flip = true;
    lo = b;
    hi = a;
  }

  // Complex function support: split into real/imag
  if (complexFunc) {
    const refunc: IntegrandFunction = (x, ...a) => {
      const val = func(x, ...a);
      // For real-valued returns, this is just the value
      return typeof val === 'object' ? (val as unknown as { re: number }).re : val;
    };
    const imfunc: IntegrandFunction = (x, ...a) => {
      const val = func(x, ...a);
      return typeof val === 'object' ? (val as unknown as { im: number }).im : 0;
    };

    const reResult = await quad(refunc, a, b, {
      ...options,
      complexFunc: false,
    });
    const imResult = await quad(imfunc, a, b, {
      ...options,
      complexFunc: false,
    });

    if (fullOutput) {
      return [
        reResult[0],
        imResult[0],
        {
          real: (reResult as QuadFullResult)[2],
          imag: (imResult as QuadFullResult)[2],
        },
      ] as unknown as QuadFullResult;
    }
    return [reResult[0], Math.max(reResult[1], imResult[1])] as QuadResult;
  }

  // Call the WASM QUADPACK implementation
  const result = await quadWasm(func, lo, hi, args, epsabs, epsrel, limit, fullOutput);

  // Apply flip
  if (flip) {
    result[0] = -result[0];
  }

  return result;
}

/**
 * Call QUADPACK via WASM.
 */
async function quadWasm(
  func: IntegrandFunction,
  a: number,
  b: number,
  args: number[],
  epsabs: number,
  epsrel: number,
  limit: number,
  fullOutput: boolean,
): Promise<QuadResult | QuadFullResult | QuadFullResultWithMessage> {
  const Module = await loadWasmModule();

  // Determine which routine to use based on bounds
  let infbounds = 0;
  let bound = 0;
  const aInf = a === -Infinity;
  const bInf = b === Infinity;

  if (!aInf && !bInf) {
    infbounds = 0; // finite: use dqagse
  } else if (!aInf && bInf) {
    infbounds = 1; // (a, +inf)
    bound = a;
  } else if (aInf && bInf) {
    infbounds = 2; // (-inf, +inf)
    bound = 0;
  } else if (aInf && !bInf) {
    infbounds = -1; // (-inf, b)
    bound = b;
  }

  // Create JS thunk that WASM will call via function pointer.
  // QUADPACK calls fcn(double* x) — receives pointer to x value.
  const thunk = (xPtr: number): number => {
    const x = Module.HEAPF64[xPtr >> 3];
    return func(x, ...args);
  };

  // Register the JS function in WASM function table.
  // Signature: 'di' = returns double, takes i32 (pointer)
  const funcPtr = Module.addFunction(thunk, 'di');

  // Allocate WASM heap memory for work arrays and output scalars
  const doubleSize = 8;
  const intSize = 4;

  // Output scalars (allocated as single-element arrays)
  const resultPtr = Module._malloc(doubleSize);
  const abserrPtr = Module._malloc(doubleSize);
  const nevalPtr = Module._malloc(intSize);
  const ierPtr = Module._malloc(intSize);
  const lastPtr = Module._malloc(intSize);

  // Work arrays
  const alistPtr = Module._malloc(limit * doubleSize);
  const blistPtr = Module._malloc(limit * doubleSize);
  const rlistPtr = Module._malloc(limit * doubleSize);
  const elistPtr = Module._malloc(limit * doubleSize);
  const iordPtr = Module._malloc(limit * intSize);

  try {
    if (infbounds === 0) {
      // Finite interval: dqagse
      Module._wasm_dqagse(
        funcPtr, a, b, epsabs, epsrel, limit,
        resultPtr, abserrPtr, nevalPtr, ierPtr,
        alistPtr, blistPtr, rlistPtr, elistPtr, iordPtr, lastPtr,
      );
    } else {
      // Infinite interval: dqagie
      Module._wasm_dqagie(
        funcPtr, bound, infbounds, epsabs, epsrel, limit,
        resultPtr, abserrPtr, nevalPtr, ierPtr,
        alistPtr, blistPtr, rlistPtr, elistPtr, iordPtr, lastPtr,
      );
    }

    // Read output values from WASM heap
    const result = Module.HEAPF64[resultPtr >> 3];
    const abserr = Module.HEAPF64[abserrPtr >> 3];
    const neval = Module.HEAP32[nevalPtr >> 2];
    const ier = Module.HEAP32[ierPtr >> 2];
    const last = Module.HEAP32[lastPtr >> 2];

    // Handle error codes
    if (ier !== 0) {
      const msg = QUAD_MESSAGES[ier] ?? "Unknown error.";

      if (ier === 6) {
        throw new ValueError(msg);
      }

      if (fullOutput) {
        const infodict = readInfoDict(Module, last, limit, alistPtr, blistPtr, rlistPtr, elistPtr, iordPtr, neval);
        return [result, abserr, infodict, msg] as QuadFullResultWithMessage;
      }

      // Emit warning but still return the result
      console.warn(`IntegrationWarning: ${msg}`);
    }

    if (fullOutput) {
      const infodict = readInfoDict(Module, last, limit, alistPtr, blistPtr, rlistPtr, elistPtr, iordPtr, neval);
      return [result, abserr, infodict] as QuadFullResult;
    }

    return [result, abserr] as QuadResult;
  } finally {
    // Clean up WASM resources
    Module.removeFunction(funcPtr);
    Module._free(resultPtr);
    Module._free(abserrPtr);
    Module._free(nevalPtr);
    Module._free(ierPtr);
    Module._free(lastPtr);
    Module._free(alistPtr);
    Module._free(blistPtr);
    Module._free(rlistPtr);
    Module._free(elistPtr);
    Module._free(iordPtr);
  }
}

/**
 * Read the info dictionary from WASM heap arrays.
 */
function readInfoDict(
  Module: SciWasmModule,
  last: number,
  limit: number,
  alistPtr: number,
  blistPtr: number,
  rlistPtr: number,
  elistPtr: number,
  iordPtr: number,
  neval: number,
): QuadInfoDict {
  // Copy data from WASM heap
  const alist = new Float64Array(limit);
  const blist = new Float64Array(limit);
  const rlist = new Float64Array(limit);
  const elist = new Float64Array(limit);
  const iord = new Int32Array(limit);

  for (let i = 0; i < limit; i++) {
    alist[i] = Module.HEAPF64[(alistPtr >> 3) + i];
    blist[i] = Module.HEAPF64[(blistPtr >> 3) + i];
    rlist[i] = Module.HEAPF64[(rlistPtr >> 3) + i];
    elist[i] = Module.HEAPF64[(elistPtr >> 3) + i];
    iord[i] = Module.HEAP32[(iordPtr >> 2) + i];
  }

  return { neval, last, alist, blist, rlist, elist, iord };
}

/**
 * ValueError for invalid input (ier=6).
 */
class ValueError extends Error {
  constructor(message: string) {
    super(message);
    this.name = 'ValueError';
  }
}
