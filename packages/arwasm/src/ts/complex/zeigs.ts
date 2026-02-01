/**
 * ZEIGS - Complex Eigenvalue Solver
 *
 * Computes eigenvalues and eigenvectors of complex matrices
 * using ARPACK's Implicitly Restarted Arnoldi Method.
 */

import { loadARPACKModule } from '../loader.js';
import type { ARPACKModule } from '../types.js';
import type {
  ComplexMatVecFunction,
  ZeigsOptions,
  ZeigsResult,
} from '../high-level-types.js';
import type { ComplexArray } from '../helpers.js';
import {
  znaupdWorklSize,
  zneupdWorkevSize,
  znaupdRworkSize,
  defaultNcv,
  getZnaupdMessage,
  getZneupdMessage,
  allocateDoubles,
  allocateComplex,
  allocateInts,
  allocateString,
  readComplex,
  readInt,
  freeAll,
} from '../helpers.js';

/**
 * Compute eigenvalues and eigenvectors of a complex matrix.
 *
 * Uses ARPACK's ZNAUPD/ZNEUPD routines with the Implicitly Restarted
 * Arnoldi Method. This is efficient for large sparse complex matrices
 * where only a few eigenvalues are needed.
 *
 * @param matvec - Function computing y = A*x where A is a complex matrix
 * @param n - Dimension of the matrix (A is n×n)
 * @param nev - Number of eigenvalues to compute (must satisfy 0 < nev < n-1)
 * @param options - Solver options
 * @returns Complex eigenvalues and eigenvectors
 *
 * @example
 * ```ts
 * import { zeigs } from 'arwasm';
 *
 * // Complex matrix: diagonal 2+0.1i, off-diagonals -1
 * const n = 100;
 * const matvec = (x: ComplexArray): ComplexArray => {
 *   const re = new Float64Array(n);
 *   const im = new Float64Array(n);
 *   for (let i = 0; i < n; i++) {
 *     // Diagonal: (2 + 0.1i) * x[i]
 *     re[i] = 2 * x.re[i] - 0.1 * x.im[i];
 *     im[i] = 2 * x.im[i] + 0.1 * x.re[i];
 *     // Off-diagonals: -1 * x[i±1]
 *     if (i > 0) { re[i] -= x.re[i-1]; im[i] -= x.im[i-1]; }
 *     if (i < n-1) { re[i] -= x.re[i+1]; im[i] -= x.im[i+1]; }
 *   }
 *   return { re, im };
 * };
 *
 * const result = await zeigs(matvec, n, 6, { which: 'LM' });
 * console.log('Eigenvalues (real):', result.eigenvaluesReal);
 * console.log('Eigenvalues (imag):', result.eigenvaluesImag);
 * ```
 */
export async function zeigs(
  matvec: ComplexMatVecFunction,
  n: number,
  nev: number,
  options?: ZeigsOptions
): Promise<ZeigsResult> {
  // Validate inputs
  if (n < 1) {
    throw new Error('Matrix dimension n must be positive');
  }
  if (nev < 1 || nev >= n - 1) {
    throw new Error(`nev must satisfy 0 < nev < n-1 (got nev=${nev}, n=${n})`);
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
    OPinv,
    Bmatvec,
    mode: userMode,
  } = options ?? {};

  // Calculate default ncv (same as non-symmetric: nev+2 <= ncv <= n)
  const ncv = userNcv ?? defaultNcv(n, nev, false);

  // Validate ncv
  if (ncv <= nev + 1 || ncv > n) {
    throw new Error(
      `ncv must satisfy nev+2 <= ncv <= n (got ncv=${ncv}, nev=${nev}, n=${n})`
    );
  }

  // Determine mode
  let mode = userMode ?? 1;
  if (sigma !== undefined && !userMode) {
    mode = 3; // shift-invert
  }

  // Determine bmat
  const bmat = Bmatvec ? 'G' : 'I';

  // Select operator function based on mode
  let operator: ComplexMatVecFunction;
  if (mode === 1) {
    operator = matvec;
  } else if (mode === 3) {
    if (!OPinv) {
      throw new Error('OPinv function required for shift-invert mode');
    }
    operator = OPinv;
  } else {
    operator = matvec;
  }

  return solveComplexInternal(
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
    sigma ?? { re: 0, im: 0 },
    operator,
    Bmatvec ?? null
  );
}

/**
 * Internal complex eigenvalue solver implementation.
 */
async function solveComplexInternal(
  Module: ARPACKModule,
  n: number,
  nev: number,
  which: string,
  bmat: string,
  tol: number,
  ncv: number,
  maxiter: number,
  v0: ComplexArray | null,
  return_eigenvectors: boolean,
  mode: number,
  sigma: { re: number; im: number },
  operator: ComplexMatVecFunction,
  Bmatvec: ComplexMatVecFunction | null
): Promise<ZeigsResult> {
  // Calculate work array sizes
  const lworkl = znaupdWorklSize(ncv);
  const rworkSize = znaupdRworkSize(ncv);

  // Allocate WASM memory
  const idoPtr = allocateInts(Module, [0], 1);
  const bmatPtr = allocateString(Module, bmat);
  const nPtr = allocateInts(Module, [n], 1);
  const whichPtr = allocateString(Module, which);
  const nevPtr = allocateInts(Module, [nev], 1);
  const tolPtr = allocateDoubles(Module, [tol], 1);

  // residPtr: complex array of size n (2n doubles)
  const residPtr = allocateComplex(Module, v0, n);

  const ncvPtr = allocateInts(Module, [ncv], 1);

  // vPtr: complex array of size n*ncv (2*n*ncv doubles)
  const vPtr = allocateComplex(Module, null, n * ncv);

  const ldvPtr = allocateInts(Module, [n], 1);

  // iparam array (11 elements)
  const iparam = new Array(11).fill(0);
  iparam[0] = 1; // ishfts: exact shifts
  iparam[2] = maxiter; // mxiter
  iparam[3] = 1; // nb: block size (must be 1)
  iparam[6] = mode; // mode
  const iparamPtr = allocateInts(Module, iparam, 11);

  const ipntrPtr = allocateInts(Module, null, 14);

  // workdPtr: complex array of size 3*n (6*n doubles)
  const workdPtr = allocateComplex(Module, null, 3 * n);

  // worklPtr: complex array of size lworkl (2*lworkl doubles)
  const worklPtr = allocateComplex(Module, null, lworkl);

  const lworklPtr = allocateInts(Module, [lworkl], 1);

  // rworkPtr: real array of size ncv
  const rworkPtr = allocateDoubles(Module, null, rworkSize);

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
    rworkPtr,
    infoPtr,
  ];

  try {
    let ido = 0;
    let nops = 0;

    // Reverse communication loop
    while (true) {
      Module._znaupd_(
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
        rworkPtr,
        infoPtr
      );

      ido = readInt(Module, idoPtr);

      if (ido === 99) {
        // Converged
        break;
      } else if (ido === -1 || ido === 1) {
        // Compute y = Op * x
        // Read ipntr values (1-indexed in Fortran, convert to 0-indexed)
        const ipntr1 = Module.HEAP32[(ipntrPtr >> 2) + 0] - 1; // x location in workd (complex index)
        const ipntr2 = Module.HEAP32[(ipntrPtr >> 2) + 1] - 1; // y location in workd (complex index)

        // Create ComplexArray view of x from workd
        // workd is stored as interleaved complex: [re0, im0, re1, im1, ...]
        // ipntr values are 1-indexed complex positions
        const xOffset = (workdPtr >> 3) + ipntr1 * 2;
        const xRe = new Float64Array(n);
        const xIm = new Float64Array(n);
        for (let i = 0; i < n; i++) {
          xRe[i] = Module.HEAPF64[xOffset + 2 * i];
          xIm[i] = Module.HEAPF64[xOffset + 2 * i + 1];
        }

        // Compute y = Op * x
        const y = operator({ re: xRe, im: xIm });
        nops++;

        // Write y to workd
        const yOffset = (workdPtr >> 3) + ipntr2 * 2;
        for (let i = 0; i < n; i++) {
          Module.HEAPF64[yOffset + 2 * i] = y.re[i];
          Module.HEAPF64[yOffset + 2 * i + 1] = y.im[i];
        }
      } else if (ido === 2) {
        // Compute y = B * x (for generalized problems)
        const ipntr1 = Module.HEAP32[(ipntrPtr >> 2) + 0] - 1;
        const ipntr2 = Module.HEAP32[(ipntrPtr >> 2) + 1] - 1;

        const xOffset = (workdPtr >> 3) + ipntr1 * 2;
        const xRe = new Float64Array(n);
        const xIm = new Float64Array(n);
        for (let i = 0; i < n; i++) {
          xRe[i] = Module.HEAPF64[xOffset + 2 * i];
          xIm[i] = Module.HEAPF64[xOffset + 2 * i + 1];
        }

        const y = Bmatvec ? Bmatvec({ re: xRe, im: xIm }) : { re: xRe, im: xIm };
        nops++;

        const yOffset = (workdPtr >> 3) + ipntr2 * 2;
        for (let i = 0; i < n; i++) {
          Module.HEAPF64[yOffset + 2 * i] = y.re[i];
          Module.HEAPF64[yOffset + 2 * i + 1] = y.im[i];
        }
      } else {
        // Unknown ido, this shouldn't happen
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
        message: getZnaupdMessage(info),
      };
    }

    // Extract eigenvalues and eigenvectors with ZNEUPD
    const rvec = return_eigenvectors ? 1 : 0;
    const howmny = 'A';

    const rvecPtr = allocateInts(Module, [rvec], 1);
    const howmnyPtr = allocateString(Module, howmny);
    const selectPtr = allocateInts(Module, null, ncv);

    // dPtr: complex eigenvalues, size nev+1
    const dPtr = allocateComplex(Module, null, nev + 1);

    // zPtr: complex eigenvectors, size n*(nev+1)
    const zPtr = allocateComplex(Module, null, n * (nev + 1));

    const ldzPtr = allocateInts(Module, [n], 1);

    // sigmaPtr: complex shift, size 1 (2 doubles)
    const sigmaPtr = allocateDoubles(Module, [sigma.re, sigma.im], 2);

    // workevPtr: complex work array, size 2*ncv
    const workevSize = zneupdWorkevSize(ncv);
    const workevPtr = allocateComplex(Module, null, workevSize);

    const extractPtrs = [
      rvecPtr,
      howmnyPtr,
      selectPtr,
      dPtr,
      zPtr,
      ldzPtr,
      sigmaPtr,
      workevPtr,
    ];

    try {
      Module._zneupd_(
        rvecPtr,
        howmnyPtr,
        selectPtr,
        dPtr,
        zPtr,
        ldzPtr,
        sigmaPtr,
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
        rworkPtr,
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
          message: getZneupdMessage(extractInfo),
        };
      }

      // Read complex eigenvalues
      const eigenvalues = readComplex(Module, dPtr, nconv);

      // Read complex eigenvectors if requested
      let eigenvectors: ComplexArray[] | undefined;
      if (return_eigenvectors && nconv > 0) {
        eigenvectors = [];
        for (let i = 0; i < nconv; i++) {
          const vecRe = new Float64Array(n);
          const vecIm = new Float64Array(n);
          // Column-major storage: z(j,i) = z[j + i*n]
          // Complex interleaved: z[j + i*n] is at doubles [2*(j + i*n), 2*(j + i*n) + 1]
          const baseOffset = (zPtr >> 3) + i * n * 2;
          for (let j = 0; j < n; j++) {
            vecRe[j] = Module.HEAPF64[baseOffset + 2 * j];
            vecIm[j] = Module.HEAPF64[baseOffset + 2 * j + 1];
          }
          eigenvectors.push({ re: vecRe, im: vecIm });
        }
      }

      return {
        eigenvaluesReal: eigenvalues.re,
        eigenvaluesImag: eigenvalues.im,
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
