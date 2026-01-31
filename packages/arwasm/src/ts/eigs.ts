/**
 * EIGS - Symmetric Eigenvalue Solver
 *
 * Computes eigenvalues and eigenvectors of real symmetric matrices
 * using ARPACK's Implicitly Restarted Lanczos Method.
 */

import { loadARPACKModule } from './loader.js';
import type { ARPACKModule } from './types.js';
import type { MatVecFunction, EigsOptions, EigsResult } from './high-level-types.js';
import {
  dsaupdWorklSize,
  defaultNcv,
  getDsaupdMessage,
  getDseupdMessage,
  allocateDoubles,
  allocateInts,
  allocateString,
  readDoubles,
  readInt,
  freeAll,
} from './helpers.js';

/**
 * Compute eigenvalues and eigenvectors of a real symmetric matrix.
 *
 * Uses ARPACK's DSAUPD/DSEUPD routines with the Implicitly Restarted
 * Lanczos Method. This is efficient for large sparse symmetric matrices
 * where only a few eigenvalues are needed.
 *
 * @param matvec - Function computing y = A*x where A is symmetric
 * @param n - Dimension of the matrix (A is n×n)
 * @param nev - Number of eigenvalues to compute (must satisfy 0 < nev < n)
 * @param options - Solver options
 * @returns Eigenvalues and eigenvectors
 *
 * @example
 * ```ts
 * import { eigs } from 'arwasm';
 *
 * // 1D Laplacian: tridiagonal with 2 on diagonal, -1 on off-diagonals
 * const n = 100;
 * const matvec = (x: Float64Array): Float64Array => {
 *   const y = new Float64Array(n);
 *   for (let i = 0; i < n; i++) {
 *     y[i] = 2 * x[i];
 *     if (i > 0) y[i] -= x[i - 1];
 *     if (i < n - 1) y[i] -= x[i + 1];
 *   }
 *   return y;
 * };
 *
 * const result = await eigs(matvec, n, 6, { which: 'SM' });
 * console.log('Smallest eigenvalues:', result.eigenvalues);
 * ```
 */
export async function eigs(
  matvec: MatVecFunction,
  n: number,
  nev: number,
  options?: EigsOptions
): Promise<EigsResult> {
  // Validate inputs
  if (n < 1) {
    throw new Error('Matrix dimension n must be positive');
  }
  if (nev < 1 || nev >= n) {
    throw new Error(`nev must satisfy 0 < nev < n (got nev=${nev}, n=${n})`);
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

  // Calculate default ncv
  const ncv = userNcv ?? defaultNcv(n, nev, true);

  // Validate ncv
  if (ncv <= nev || ncv > n) {
    throw new Error(
      `ncv must satisfy nev < ncv <= n (got ncv=${ncv}, nev=${nev}, n=${n})`
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
  let operator: MatVecFunction;
  if (mode === 1) {
    operator = matvec;
  } else if (mode === 3 || mode === 4 || mode === 5) {
    if (!OPinv) {
      throw new Error('OPinv function required for shift-invert/buckling/cayley mode');
    }
    operator = OPinv;
  } else {
    operator = matvec;
  }

  return solveSymmetricInternal(
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
    operator,
    Bmatvec ?? null
  );
}

/**
 * Internal symmetric eigenvalue solver implementation.
 */
async function solveSymmetricInternal(
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
  operator: MatVecFunction,
  Bmatvec: MatVecFunction | null
): Promise<EigsResult> {
  // Calculate work array sizes
  const lworkl = dsaupdWorklSize(ncv);

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

  const ipntrPtr = allocateInts(Module, null, 11);
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
      Module._dsaupd_(
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
        // Read ipntr values (1-indexed in Fortran, convert to 0-indexed)
        const ipntr1 = Module.HEAP32[(ipntrPtr >> 2) + 0] - 1; // x location in workd
        const ipntr2 = Module.HEAP32[(ipntrPtr >> 2) + 1] - 1; // y location in workd

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

        const y = Bmatvec ? Bmatvec(x) : x; // If no B, assume identity
        nops++;

        const yOffset = (workdPtr >> 3) + ipntr2;
        for (let i = 0; i < n; i++) {
          Module.HEAPF64[yOffset + i] = y[i];
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
        eigenvalues: new Float64Array(0),
        eigenvectors: return_eigenvectors ? [] : undefined,
        niter,
        nops,
        nconv: 0,
        info,
        success: false,
        message: getDsaupdMessage(info),
      };
    }

    // Extract eigenvalues and eigenvectors with DSEUPD
    const rvec = return_eigenvectors ? 1 : 0;
    const howmny = 'A';

    const rvecPtr = allocateInts(Module, [rvec], 1);
    const howmnyPtr = allocateString(Module, howmny);
    const selectPtr = allocateInts(Module, null, ncv);
    const dPtr = allocateDoubles(Module, null, nev);
    const zPtr = allocateDoubles(Module, null, n * nev);
    const ldzPtr = allocateInts(Module, [n], 1);
    const sigmaPtr = allocateDoubles(Module, [sigma], 1);

    const extractPtrs = [
      rvecPtr,
      howmnyPtr,
      selectPtr,
      dPtr,
      zPtr,
      ldzPtr,
      sigmaPtr,
    ];

    try {
      Module._dseupd_(
        rvecPtr,
        howmnyPtr,
        selectPtr,
        dPtr,
        zPtr,
        ldzPtr,
        sigmaPtr,
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
          eigenvalues: new Float64Array(0),
          eigenvectors: return_eigenvectors ? [] : undefined,
          niter,
          nops,
          nconv,
          info: extractInfo,
          success: false,
          message: getDseupdMessage(extractInfo),
        };
      }

      // Read eigenvalues
      const eigenvalues = readDoubles(Module, dPtr, nconv);

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
        eigenvalues,
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
