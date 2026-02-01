/**
 * EIGN - Non-Symmetric Eigenvalue Solver
 *
 * Computes eigenvalues and eigenvectors of real non-symmetric matrices
 * using ARPACK's Implicitly Restarted Arnoldi Method.
 */

import { loadARPACKModule } from '../loader.js';
import type { ARPACKModule } from '../types.js';
import type { MatVecFunction, EignOptions, EignResult } from '../high-level-types.js';
import {
  dnaupdWorklSize,
  dneupdWorkevSize,
  defaultNcv,
  getDnaupdMessage,
  getDneupdMessage,
  allocateDoubles,
  allocateInts,
  allocateString,
  readDoubles,
  readInt,
  freeAll,
} from '../helpers.js';

/**
 * Compute eigenvalues and eigenvectors of a real non-symmetric matrix.
 *
 * Uses ARPACK's DNAUPD/DNEUPD routines with the Implicitly Restarted
 * Arnoldi Method. This handles general (non-symmetric) matrices where
 * eigenvalues may be complex.
 *
 * @param matvec - Function computing y = A*x
 * @param n - Dimension of the matrix (A is n×n)
 * @param nev - Number of eigenvalues to compute (must satisfy 0 < nev < n-1)
 * @param options - Solver options
 * @returns Eigenvalues (with real/imaginary parts) and eigenvectors
 *
 * @example
 * ```ts
 * import { eign } from 'arwasm';
 *
 * // Convection-diffusion matrix (non-symmetric)
 * const n = 100;
 * const matvec = (x: Float64Array): Float64Array => {
 *   const y = new Float64Array(n);
 *   for (let i = 0; i < n; i++) {
 *     y[i] = 2 * x[i];
 *     if (i > 0) y[i] -= x[i - 1];
 *     if (i < n - 1) y[i] -= 0.5 * x[i + 1]; // Asymmetric
 *   }
 *   return y;
 * };
 *
 * const result = await eign(matvec, n, 6, { which: 'LR' });
 * console.log('Largest real part eigenvalues:', result.eigenvaluesReal);
 * ```
 */
export async function eign(
  matvec: MatVecFunction,
  n: number,
  nev: number,
  options?: EignOptions
): Promise<EignResult> {
  // Validate inputs
  if (n < 1) {
    throw new Error('Matrix dimension n must be positive');
  }
  if (nev < 1 || nev >= n - 1) {
    throw new Error(
      `nev must satisfy 0 < nev < n-1 (got nev=${nev}, n=${n})`
    );
  }

  const Module = await loadARPACKModule();

  const {
    which = 'LM',
    tol = 0,
    ncv: userNcv,
    maxiter = n * 10,
    v0,
    return_eigenvectors = true,
    sigma,
    sigmai = 0,
    OPinv,
    Bmatvec,
    mode: userMode,
  } = options ?? {};

  // Calculate default ncv
  const ncv = userNcv ?? defaultNcv(n, nev, false);

  // Validate ncv: must satisfy nev+2 <= ncv <= n
  if (ncv < nev + 2 || ncv > n) {
    throw new Error(
      `ncv must satisfy nev+2 <= ncv <= n (got ncv=${ncv}, nev=${nev}, n=${n})`
    );
  }

  // Determine mode
  let mode = userMode ?? 1;
  if (sigma !== undefined && !userMode) {
    mode = 3; // shift-invert
  }

  const bmat = Bmatvec ? 'G' : 'I';

  let operator: MatVecFunction;
  if (mode === 1) {
    operator = matvec;
  } else if (mode === 3 || mode === 4) {
    if (!OPinv) {
      throw new Error('OPinv function required for shift-invert mode');
    }
    operator = OPinv;
  } else {
    operator = matvec;
  }

  return solveNonSymmetricInternal(
    Module,
    n,
    nev,
    which,
    bmat,
    tol,
    ncv,
    maxiter,
    v0 ?? null,
    return_eigenvectors,
    mode,
    sigma ?? 0,
    sigmai,
    operator,
    Bmatvec ?? null
  );
}

/**
 * Internal non-symmetric eigenvalue solver implementation.
 */
async function solveNonSymmetricInternal(
  Module: ARPACKModule,
  n: number,
  nev: number,
  which: string,
  bmat: string,
  tol: number,
  ncv: number,
  maxiter: number,
  v0: number[] | Float64Array | null,
  return_eigenvectors: boolean,
  mode: number,
  sigma: number,
  sigmai: number,
  operator: MatVecFunction,
  Bmatvec: MatVecFunction | null
): Promise<EignResult> {
  // Calculate work array sizes
  const lworkl = dnaupdWorklSize(ncv);

  // Allocate WASM memory
  const idoPtr = allocateInts(Module, [0], 1);
  const bmatPtr = allocateString(Module, bmat);
  const nPtr = allocateInts(Module, [n], 1);
  const whichPtr = allocateString(Module, which);
  const nevPtr = allocateInts(Module, [nev], 1);
  const tolPtr = allocateDoubles(Module, [tol], 1);
  const residPtr = allocateDoubles(Module, v0, n);
  const ncvPtr = allocateInts(Module, [ncv], 1);
  const vPtr = allocateDoubles(Module, null, n * ncv);
  const ldvPtr = allocateInts(Module, [n], 1);

  // iparam array (11 elements)
  const iparam = new Array(11).fill(0);
  iparam[0] = 1; // ishfts: exact shifts
  iparam[2] = maxiter; // mxiter
  iparam[3] = 1; // nb: block size (must be 1)
  iparam[6] = mode; // mode
  const iparamPtr = allocateInts(Module, iparam, 11);

  // ipntr has 14 elements for non-symmetric (vs 11 for symmetric)
  const ipntrPtr = allocateInts(Module, null, 14);
  const workdPtr = allocateDoubles(Module, null, 3 * n);
  const worklPtr = allocateDoubles(Module, null, lworkl);
  const lworklPtr = allocateInts(Module, [lworkl], 1);
  // info: 0 = random start, 1 = user-supplied residual
  const infoPtr = allocateInts(Module, [v0 ? 1 : 0], 1);

  // Pointers to track for cleanup
  const allPtrs = [
    idoPtr,
    bmatPtr,
    nPtr,
    whichPtr,
    nevPtr,
    tolPtr,
    residPtr,
    ncvPtr,
    vPtr,
    ldvPtr,
    iparamPtr,
    ipntrPtr,
    workdPtr,
    worklPtr,
    lworklPtr,
    infoPtr,
  ];

  try {
    let ido = 0;
    let nops = 0;

    // Reverse communication loop
    while (true) {
      Module._dnaupd_(
        idoPtr,
        bmatPtr,
        nPtr,
        whichPtr,
        nevPtr,
        tolPtr,
        residPtr,
        ncvPtr,
        vPtr,
        ldvPtr,
        iparamPtr,
        ipntrPtr,
        workdPtr,
        worklPtr,
        lworklPtr,
        infoPtr
      );

      ido = readInt(Module, idoPtr);

      if (ido === 99) {
        // Converged
        break;
      } else if (ido === -1 || ido === 1) {
        // Compute y = Op * x
        const ipntr1 = Module.HEAP32[(ipntrPtr >> 2) + 0] - 1; // x location
        const ipntr2 = Module.HEAP32[(ipntrPtr >> 2) + 1] - 1; // y location

        // Create Float64Array view of x from workd (zero-copy)
        const xOffset = (workdPtr >> 3) + ipntr1;
        const x = Module.HEAPF64.subarray(xOffset, xOffset + n);

        // Compute y = Op * x
        const y = operator(x);
        nops++;

        // Write y to workd
        const yOffset = (workdPtr >> 3) + ipntr2;
        for (let i = 0; i < n; i++) {
          Module.HEAPF64[yOffset + i] = y[i];
        }
      } else if (ido === 2) {
        // Compute y = B * x (for generalized problems)
        const ipntr1 = Module.HEAP32[(ipntrPtr >> 2) + 0] - 1;
        const ipntr2 = Module.HEAP32[(ipntrPtr >> 2) + 1] - 1;

        const xOffset = (workdPtr >> 3) + ipntr1;
        const x = Module.HEAPF64.subarray(xOffset, xOffset + n);

        const y = Bmatvec ? Bmatvec(x) : x;
        nops++;

        const yOffset = (workdPtr >> 3) + ipntr2;
        for (let i = 0; i < n; i++) {
          Module.HEAPF64[yOffset + i] = y[i];
        }
      } else if (ido === 3) {
        // Compute shifts - advanced usage, not typically needed
        // For now, just continue (ARPACK will use internal shifts)
        break;
      } else {
        // Unknown ido
        break;
      }
    }

    // Check convergence info
    const info = readInt(Module, infoPtr);
    const nconv = Module.HEAP32[(iparamPtr >> 2) + 4];
    const niter = Module.HEAP32[(iparamPtr >> 2) + 2];

    if (info < 0) {
      return {
        eigenvaluesReal: new Float64Array(0),
        eigenvaluesImag: new Float64Array(0),
        eigenvectors: return_eigenvectors ? [] : undefined,
        niter,
        nops,
        nconv: 0,
        info,
        success: false,
        message: getDnaupdMessage(info),
      };
    }

    // Extract eigenvalues and eigenvectors with DNEUPD
    const rvec = return_eigenvectors ? 1 : 0;
    const howmny = 'A';

    const rvecPtr = allocateInts(Module, [rvec], 1);
    const howmnyPtr = allocateString(Module, howmny);
    const selectPtr = allocateInts(Module, null, ncv);
    // nev+1 for possible complex conjugate pair
    const drPtr = allocateDoubles(Module, null, nev + 1);
    const diPtr = allocateDoubles(Module, null, nev + 1);
    const zPtr = allocateDoubles(Module, null, n * (nev + 1));
    const ldzPtr = allocateInts(Module, [n], 1);
    const sigmarPtr = allocateDoubles(Module, [sigma], 1);
    const sigmaiPtr = allocateDoubles(Module, [sigmai], 1);
    const workevPtr = allocateDoubles(Module, null, dneupdWorkevSize(ncv));

    const extractPtrs = [
      rvecPtr,
      howmnyPtr,
      selectPtr,
      drPtr,
      diPtr,
      zPtr,
      ldzPtr,
      sigmarPtr,
      sigmaiPtr,
      workevPtr,
    ];

    try {
      Module._dneupd_(
        rvecPtr,
        howmnyPtr,
        selectPtr,
        drPtr,
        diPtr,
        zPtr,
        ldzPtr,
        sigmarPtr,
        sigmaiPtr,
        workevPtr,
        bmatPtr,
        nPtr,
        whichPtr,
        nevPtr,
        tolPtr,
        residPtr,
        ncvPtr,
        vPtr,
        ldvPtr,
        iparamPtr,
        ipntrPtr,
        workdPtr,
        worklPtr,
        lworklPtr,
        infoPtr
      );

      const extractInfo = readInt(Module, infoPtr);

      if (extractInfo < 0) {
        return {
          eigenvaluesReal: new Float64Array(0),
          eigenvaluesImag: new Float64Array(0),
          eigenvectors: return_eigenvectors ? [] : undefined,
          niter,
          nops,
          nconv,
          info: extractInfo,
          success: false,
          message: getDneupdMessage(extractInfo),
        };
      }

      // Read eigenvalues (real and imaginary parts)
      const eigenvaluesReal = readDoubles(Module, drPtr, nconv);
      const eigenvaluesImag = readDoubles(Module, diPtr, nconv);

      // Read eigenvectors if requested
      let eigenvectors: Float64Array[] | undefined;
      if (return_eigenvectors && nconv > 0) {
        eigenvectors = [];
        for (let i = 0; i < nconv; i++) {
          const vec = new Float64Array(n);
          for (let j = 0; j < n; j++) {
            // Column-major storage: z(j,i) = z[j + i*n]
            vec[j] = Module.HEAPF64[(zPtr >> 3) + j + i * n];
          }
          eigenvectors.push(vec);
        }
      }

      return {
        eigenvaluesReal,
        eigenvaluesImag,
        eigenvectors,
        niter,
        nops,
        nconv,
        info: 0,
        success: true,
        message: 'Converged successfully',
      };
    } finally {
      freeAll(Module, extractPtrs);
    }
  } finally {
    freeAll(Module, allPtrs);
  }
}
