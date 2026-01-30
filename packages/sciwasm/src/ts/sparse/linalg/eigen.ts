/**
 * Sparse eigenvalue solvers using ARPACK
 *
 * Provides eigenvalue and singular value decomposition for large sparse matrices:
 * - eigsh: Symmetric eigenvalues (Lanczos method)
 * - eigs: General eigenvalues (Arnoldi method)
 * - svds: Truncated SVD
 */

import type { SparseMatrix } from '../base.js';
import type {
  LinearOperator as ILinearOperator,
  LinearOperatorLike,
  EigshOptions,
  EigsOptions,
  EigsResult,
  SvdsOptions,
  SvdsResult,
} from './types.js';
import { aslinearoperator } from './interface.js';
import { loadARPACKModule, getARPACKModule } from 'arwasm';
import type { ARPACKModule } from 'arwasm';
import { DSAUPD_ERRORS, DSEUPD_ERRORS, DNAUPD_ERRORS, DNEUPD_ERRORS } from 'arwasm';
import { ARPACKError, ARPACKNoConvergence, DimensionMismatchError } from '../../errors.js';

// ============================================================
// WASM Memory Helpers
// ============================================================

function fromWasmF64(wasm: ARPACKModule, ptr: number, len: number): Float64Array {
  const result = new Float64Array(len);
  result.set(wasm.HEAPF64.subarray(ptr >> 3, (ptr >> 3) + len));
  return result;
}

function writeString(wasm: ARPACKModule, str: string, ptr: number): void {
  for (let i = 0; i < str.length; i++) {
    wasm.HEAP8[ptr + i] = str.charCodeAt(i);
  }
}

// ============================================================
// ARPACK Symmetric Eigenvalue Solver (eigsh)
// ============================================================

/**
 * Parameters manager for symmetric ARPACK solver
 */
class SymmetricARPACKParams {
  private wasm: ARPACKModule;
  private n: number;
  private nev: number;
  private ncv: number;
  private matvec: (x: Float64Array) => Float64Array;

  // WASM pointers (initialized in allocate())
  private idoPtr!: number;
  private bmatPtr!: number;
  private nPtr!: number;
  private whichPtr!: number;
  private nevPtr!: number;
  private tolPtr!: number;
  private residPtr!: number;
  private ncvPtr!: number;
  private vPtr!: number;
  private ldvPtr!: number;
  private iparamPtr!: number;
  private ipntrPtr!: number;
  private workdPtr!: number;
  private worklPtr!: number;
  private lworklPtr!: number;
  private infoPtr!: number;

  constructor(
    n: number,
    nev: number,
    matvec: (x: Float64Array) => Float64Array,
    options: EigshOptions = {}
  ) {
    this.wasm = getARPACKModule();
    this.n = n;
    this.nev = nev;
    this.ncv = options.ncv ?? Math.min(n, Math.max(2 * nev + 1, 20));
    this.matvec = matvec;

    this.allocate();
    this.initialize(options);
  }

  private allocate(): void {
    const wasm = this.wasm;
    const n = this.n;
    const ncv = this.ncv;
    const lworkl = ncv * (ncv + 8);

    // Integer scalars (4 bytes each)
    this.idoPtr = wasm._malloc(4);
    this.nPtr = wasm._malloc(4);
    this.nevPtr = wasm._malloc(4);
    this.ncvPtr = wasm._malloc(4);
    this.ldvPtr = wasm._malloc(4);
    this.lworklPtr = wasm._malloc(4);
    this.infoPtr = wasm._malloc(4);

    // Character/string parameters
    this.bmatPtr = wasm._malloc(4);
    this.whichPtr = wasm._malloc(4);

    // Double scalar
    this.tolPtr = wasm._malloc(8);

    // Arrays
    this.residPtr = wasm._malloc(n * 8);
    this.vPtr = wasm._malloc(n * ncv * 8);
    this.iparamPtr = wasm._malloc(11 * 4);
    this.ipntrPtr = wasm._malloc(11 * 4);
    this.workdPtr = wasm._malloc(3 * n * 8);
    this.worklPtr = wasm._malloc(lworkl * 8);
  }

  private initialize(options: EigshOptions): void {
    const wasm = this.wasm;
    const n = this.n;
    const ncv = this.ncv;

    // Set scalar values
    wasm.setValue(this.idoPtr, 0, 'i32');
    wasm.setValue(this.nPtr, n, 'i32');
    wasm.setValue(this.nevPtr, this.nev, 'i32');
    wasm.setValue(this.ncvPtr, ncv, 'i32');
    wasm.setValue(this.ldvPtr, n, 'i32');
    wasm.setValue(this.lworklPtr, ncv * (ncv + 8), 'i32');
    wasm.setValue(this.tolPtr, options.tol ?? 0, 'double');

    // BMAT = 'I' for standard eigenvalue problem
    writeString(wasm, 'I', this.bmatPtr);

    // WHICH string (2 characters)
    const which = options.which ?? 'LM';
    writeString(wasm, which, this.whichPtr);

    // INFO: 0 = random start, 1 = user supplied v0
    if (options.v0) {
      wasm.setValue(this.infoPtr, 1, 'i32');
      wasm.HEAPF64.set(options.v0, this.residPtr >> 3);
    } else {
      wasm.setValue(this.infoPtr, 0, 'i32');
    }

    // IPARAM initialization
    const iparam = new Int32Array(11);
    iparam[0] = 1;  // ISHIFT = 1 (exact shifts)
    iparam[2] = options.maxiter ?? n * 10;  // MXITER
    iparam[3] = 1;  // NB = 1 (block size)
    iparam[6] = 1;  // MODE = 1 (standard eigenvalue problem)
    wasm.HEAP32.set(iparam, this.iparamPtr >> 2);
  }

  iterate(): void {
    const wasm = this.wasm;
    const n = this.n;

    while (true) {
      // Call dsaupd
      wasm._dsaupd_(
        this.idoPtr, this.bmatPtr, this.nPtr, this.whichPtr,
        this.nevPtr, this.tolPtr, this.residPtr, this.ncvPtr,
        this.vPtr, this.ldvPtr, this.iparamPtr, this.ipntrPtr,
        this.workdPtr, this.worklPtr, this.lworklPtr, this.infoPtr
      );

      const ido = wasm.getValue(this.idoPtr, 'i32');

      if (ido === 1 || ido === -1) {
        // Compute y = A * x
        const ipntr = new Int32Array(11);
        for (let i = 0; i < 11; i++) {
          ipntr[i] = wasm.HEAP32[(this.ipntrPtr >> 2) + i];
        }

        // ARPACK uses 1-based Fortran indexing
        const xOffset = ipntr[0] - 1;
        const yOffset = ipntr[1] - 1;

        // Extract x from workd
        const x = new Float64Array(n);
        for (let i = 0; i < n; i++) {
          x[i] = wasm.HEAPF64[(this.workdPtr >> 3) + xOffset + i];
        }

        // Compute y = A * x
        const y = this.matvec(x);

        // Store y back to workd
        for (let i = 0; i < n; i++) {
          wasm.HEAPF64[(this.workdPtr >> 3) + yOffset + i] = y[i];
        }
      } else if (ido === 99) {
        // Converged
        break;
      } else {
        // Check for errors
        const info = wasm.getValue(this.infoPtr, 'i32');
        if (info !== 0) {
          const msg = DSAUPD_ERRORS[info] ?? `Unknown error code ${info}`;
          throw new ARPACKError(info, msg);
        }
        break;
      }
    }
  }

  extract(return_eigenvectors: boolean): EigsResult {
    const wasm = this.wasm;
    const n = this.n;
    const nev = this.nev;
    const ncv = this.ncv;

    // Allocate output arrays
    const dPtr = wasm._malloc(nev * 8);
    const zPtr = return_eigenvectors ? wasm._malloc(n * nev * 8) : 0;
    const ldzPtr = wasm._malloc(4);
    wasm.setValue(ldzPtr, n, 'i32');

    const sigmaPtr = wasm._malloc(8);
    wasm.setValue(sigmaPtr, 0, 'double');

    const rvecPtr = wasm._malloc(4);
    wasm.setValue(rvecPtr, return_eigenvectors ? 1 : 0, 'i32');

    const howmnyPtr = wasm._malloc(4);
    writeString(wasm, 'A', howmnyPtr);

    const selectPtr = wasm._malloc(ncv * 4);

    try {
      // Call dseupd
      wasm._dseupd_(
        rvecPtr, howmnyPtr, selectPtr, dPtr, zPtr, ldzPtr, sigmaPtr,
        this.bmatPtr, this.nPtr, this.whichPtr, this.nevPtr,
        this.tolPtr, this.residPtr, this.ncvPtr, this.vPtr, this.ldvPtr,
        this.iparamPtr, this.ipntrPtr, this.workdPtr, this.worklPtr,
        this.lworklPtr, this.infoPtr
      );

      const info = wasm.getValue(this.infoPtr, 'i32');
      if (info !== 0) {
        const msg = DSEUPD_ERRORS[info] ?? `Unknown error code ${info}`;
        throw new ARPACKError(info, msg);
      }

      // Read results
      const iparam = new Int32Array(11);
      for (let i = 0; i < 11; i++) {
        iparam[i] = wasm.HEAP32[(this.iparamPtr >> 2) + i];
      }
      const nconv = iparam[4];
      const iterations = iparam[2];

      // Copy eigenvalues
      const values = fromWasmF64(wasm, dPtr, nconv);

      // Copy eigenvectors (column-major)
      let vectors: Float64Array | undefined;
      if (return_eigenvectors && zPtr) {
        vectors = fromWasmF64(wasm, zPtr, n * nconv);
      }

      return { values, vectors, nconv, iterations };
    } finally {
      wasm._free(dPtr);
      if (zPtr) wasm._free(zPtr);
      wasm._free(ldzPtr);
      wasm._free(sigmaPtr);
      wasm._free(rvecPtr);
      wasm._free(howmnyPtr);
      wasm._free(selectPtr);
    }
  }

  free(): void {
    const wasm = this.wasm;
    wasm._free(this.idoPtr);
    wasm._free(this.bmatPtr);
    wasm._free(this.nPtr);
    wasm._free(this.whichPtr);
    wasm._free(this.nevPtr);
    wasm._free(this.tolPtr);
    wasm._free(this.residPtr);
    wasm._free(this.ncvPtr);
    wasm._free(this.vPtr);
    wasm._free(this.ldvPtr);
    wasm._free(this.iparamPtr);
    wasm._free(this.ipntrPtr);
    wasm._free(this.workdPtr);
    wasm._free(this.worklPtr);
    wasm._free(this.lworklPtr);
    wasm._free(this.infoPtr);
  }
}

// ============================================================
// ARPACK Non-Symmetric Eigenvalue Solver (eigs)
// ============================================================

/**
 * Parameters manager for non-symmetric ARPACK solver
 */
class NonSymmetricARPACKParams {
  private wasm: ARPACKModule;
  private n: number;
  private nev: number;
  private ncv: number;
  private matvec: (x: Float64Array) => Float64Array;

  // WASM pointers (initialized in allocate())
  private idoPtr!: number;
  private bmatPtr!: number;
  private nPtr!: number;
  private whichPtr!: number;
  private nevPtr!: number;
  private tolPtr!: number;
  private residPtr!: number;
  private ncvPtr!: number;
  private vPtr!: number;
  private ldvPtr!: number;
  private iparamPtr!: number;
  private ipntrPtr!: number;
  private workdPtr!: number;
  private worklPtr!: number;
  private lworklPtr!: number;
  private infoPtr!: number;

  constructor(
    n: number,
    nev: number,
    matvec: (x: Float64Array) => Float64Array,
    options: EigsOptions = {}
  ) {
    this.wasm = getARPACKModule();
    this.n = n;
    this.nev = nev;
    // For non-symmetric, ncv must satisfy nev+2 <= ncv <= n
    this.ncv = options.ncv ?? Math.min(n, Math.max(2 * nev + 1, 20));
    this.matvec = matvec;

    this.allocate();
    this.initialize(options);
  }

  private allocate(): void {
    const wasm = this.wasm;
    const n = this.n;
    const ncv = this.ncv;
    const lworkl = 3 * ncv * (ncv + 2);

    // Integer scalars
    this.idoPtr = wasm._malloc(4);
    this.nPtr = wasm._malloc(4);
    this.nevPtr = wasm._malloc(4);
    this.ncvPtr = wasm._malloc(4);
    this.ldvPtr = wasm._malloc(4);
    this.lworklPtr = wasm._malloc(4);
    this.infoPtr = wasm._malloc(4);

    // Character parameters
    this.bmatPtr = wasm._malloc(4);
    this.whichPtr = wasm._malloc(4);

    // Double scalar
    this.tolPtr = wasm._malloc(8);

    // Arrays
    this.residPtr = wasm._malloc(n * 8);
    this.vPtr = wasm._malloc(n * ncv * 8);
    this.iparamPtr = wasm._malloc(11 * 4);
    this.ipntrPtr = wasm._malloc(14 * 4);  // 14 for non-symmetric
    this.workdPtr = wasm._malloc(3 * n * 8);
    this.worklPtr = wasm._malloc(lworkl * 8);
  }

  private initialize(options: EigsOptions): void {
    const wasm = this.wasm;
    const n = this.n;
    const ncv = this.ncv;

    wasm.setValue(this.idoPtr, 0, 'i32');
    wasm.setValue(this.nPtr, n, 'i32');
    wasm.setValue(this.nevPtr, this.nev, 'i32');
    wasm.setValue(this.ncvPtr, ncv, 'i32');
    wasm.setValue(this.ldvPtr, n, 'i32');
    wasm.setValue(this.lworklPtr, 3 * ncv * (ncv + 2), 'i32');
    wasm.setValue(this.tolPtr, options.tol ?? 0, 'double');

    writeString(wasm, 'I', this.bmatPtr);
    const which = options.which ?? 'LM';
    writeString(wasm, which, this.whichPtr);

    if (options.v0) {
      wasm.setValue(this.infoPtr, 1, 'i32');
      wasm.HEAPF64.set(options.v0, this.residPtr >> 3);
    } else {
      wasm.setValue(this.infoPtr, 0, 'i32');
    }

    const iparam = new Int32Array(11);
    iparam[0] = 1;
    iparam[2] = options.maxiter ?? n * 10;
    iparam[3] = 1;
    iparam[6] = 1;
    wasm.HEAP32.set(iparam, this.iparamPtr >> 2);
  }

  iterate(): void {
    const wasm = this.wasm;
    const n = this.n;

    while (true) {
      wasm._dnaupd_(
        this.idoPtr, this.bmatPtr, this.nPtr, this.whichPtr,
        this.nevPtr, this.tolPtr, this.residPtr, this.ncvPtr,
        this.vPtr, this.ldvPtr, this.iparamPtr, this.ipntrPtr,
        this.workdPtr, this.worklPtr, this.lworklPtr, this.infoPtr
      );

      const ido = wasm.getValue(this.idoPtr, 'i32');

      if (ido === 1 || ido === -1) {
        const ipntr = new Int32Array(14);
        for (let i = 0; i < 14; i++) {
          ipntr[i] = wasm.HEAP32[(this.ipntrPtr >> 2) + i];
        }

        const xOffset = ipntr[0] - 1;
        const yOffset = ipntr[1] - 1;

        const x = new Float64Array(n);
        for (let i = 0; i < n; i++) {
          x[i] = wasm.HEAPF64[(this.workdPtr >> 3) + xOffset + i];
        }

        const y = this.matvec(x);

        for (let i = 0; i < n; i++) {
          wasm.HEAPF64[(this.workdPtr >> 3) + yOffset + i] = y[i];
        }
      } else if (ido === 99) {
        break;
      } else {
        const info = wasm.getValue(this.infoPtr, 'i32');
        if (info !== 0) {
          const msg = DNAUPD_ERRORS[info] ?? `Unknown error code ${info}`;
          throw new ARPACKError(info, msg);
        }
        break;
      }
    }
  }

  extract(return_eigenvectors: boolean): EigsResult {
    const wasm = this.wasm;
    const n = this.n;
    const nev = this.nev;
    const ncv = this.ncv;

    // Allocate output arrays
    const drPtr = wasm._malloc((nev + 1) * 8);  // Real parts
    const diPtr = wasm._malloc((nev + 1) * 8);  // Imaginary parts
    const zPtr = return_eigenvectors ? wasm._malloc(n * (nev + 1) * 8) : 0;
    const ldzPtr = wasm._malloc(4);
    wasm.setValue(ldzPtr, n, 'i32');

    const sigmarPtr = wasm._malloc(8);
    const sigmaiPtr = wasm._malloc(8);
    wasm.setValue(sigmarPtr, 0, 'double');
    wasm.setValue(sigmaiPtr, 0, 'double');

    const rvecPtr = wasm._malloc(4);
    wasm.setValue(rvecPtr, return_eigenvectors ? 1 : 0, 'i32');

    const howmnyPtr = wasm._malloc(4);
    writeString(wasm, 'A', howmnyPtr);

    const selectPtr = wasm._malloc(ncv * 4);
    const workevPtr = wasm._malloc(3 * ncv * 8);

    try {
      wasm._dneupd_(
        rvecPtr, howmnyPtr, selectPtr, drPtr, diPtr, zPtr, ldzPtr,
        sigmarPtr, sigmaiPtr, workevPtr,
        this.bmatPtr, this.nPtr, this.whichPtr, this.nevPtr,
        this.tolPtr, this.residPtr, this.ncvPtr, this.vPtr, this.ldvPtr,
        this.iparamPtr, this.ipntrPtr, this.workdPtr, this.worklPtr,
        this.lworklPtr, this.infoPtr
      );

      const info = wasm.getValue(this.infoPtr, 'i32');
      if (info !== 0) {
        const msg = DNEUPD_ERRORS[info] ?? `Unknown error code ${info}`;
        throw new ARPACKError(info, msg);
      }

      const iparam = new Int32Array(11);
      for (let i = 0; i < 11; i++) {
        iparam[i] = wasm.HEAP32[(this.iparamPtr >> 2) + i];
      }
      const nconv = iparam[4];
      const iterations = iparam[2];

      const values = fromWasmF64(wasm, drPtr, nconv);
      const valuesImag = fromWasmF64(wasm, diPtr, nconv);

      let vectors: Float64Array | undefined;
      if (return_eigenvectors && zPtr) {
        vectors = fromWasmF64(wasm, zPtr, n * nconv);
      }

      return { values, valuesImag, vectors, nconv, iterations };
    } finally {
      wasm._free(drPtr);
      wasm._free(diPtr);
      if (zPtr) wasm._free(zPtr);
      wasm._free(ldzPtr);
      wasm._free(sigmarPtr);
      wasm._free(sigmaiPtr);
      wasm._free(rvecPtr);
      wasm._free(howmnyPtr);
      wasm._free(selectPtr);
      wasm._free(workevPtr);
    }
  }

  free(): void {
    const wasm = this.wasm;
    wasm._free(this.idoPtr);
    wasm._free(this.bmatPtr);
    wasm._free(this.nPtr);
    wasm._free(this.whichPtr);
    wasm._free(this.nevPtr);
    wasm._free(this.tolPtr);
    wasm._free(this.residPtr);
    wasm._free(this.ncvPtr);
    wasm._free(this.vPtr);
    wasm._free(this.ldvPtr);
    wasm._free(this.iparamPtr);
    wasm._free(this.ipntrPtr);
    wasm._free(this.workdPtr);
    wasm._free(this.worklPtr);
    wasm._free(this.lworklPtr);
    wasm._free(this.infoPtr);
  }
}

// ============================================================
// Public API
// ============================================================

/**
 * Find k eigenvalues and eigenvectors of the real symmetric square matrix A.
 *
 * Uses ARPACK's Implicitly Restarted Lanczos Method.
 *
 * @param A - Sparse matrix or LinearOperator (must be symmetric)
 * @param k - Number of eigenvalues to compute (default: 6)
 * @param options - Solver options
 * @returns Eigenvalues and optionally eigenvectors
 *
 * @example
 * ```typescript
 * const A = csr_matrix([[4, 1], [1, 3]]);  // symmetric
 * const { values, vectors } = await eigsh(A, 2);
 * // values contains the 2 largest magnitude eigenvalues
 * ```
 */
export async function eigsh(
  A: LinearOperatorLike,
  k: number = 6,
  options: EigshOptions = {}
): Promise<EigsResult> {
  await loadARPACKModule();

  const Aop = aslinearoperator(A);
  const [m, n] = Aop.shape;

  if (m !== n) {
    throw new DimensionMismatchError(`Matrix must be square, got ${m}x${n}`);
  }

  if (k <= 0) {
    throw new Error(`k must be positive, got ${k}`);
  }

  if (k >= n) {
    throw new Error(`k must be less than n, got k=${k}, n=${n}`);
  }

  const return_eigenvectors = options.return_eigenvectors ?? true;

  const params = new SymmetricARPACKParams(n, k, Aop.matvec.bind(Aop), options);

  try {
    params.iterate();
    const result = params.extract(return_eigenvectors);

    if (result.nconv < k) {
      throw new ARPACKNoConvergence(result.nconv, k);
    }

    return result;
  } finally {
    params.free();
  }
}

/**
 * Find k eigenvalues and eigenvectors of a square matrix A.
 *
 * For non-symmetric matrices, eigenvalues may be complex.
 * Uses ARPACK's Implicitly Restarted Arnoldi Method.
 *
 * @param A - Sparse matrix or LinearOperator
 * @param k - Number of eigenvalues to compute (default: 6)
 * @param options - Solver options
 * @returns Eigenvalues (possibly complex) and optionally eigenvectors
 *
 * @example
 * ```typescript
 * const A = csr_matrix([[0, 1], [-1, 0]]);  // non-symmetric
 * const { values, valuesImag, vectors } = await eigs(A, 2);
 * // Eigenvalues are ±i (purely imaginary)
 * ```
 */
export async function eigs(
  A: LinearOperatorLike,
  k: number = 6,
  options: EigsOptions = {}
): Promise<EigsResult> {
  await loadARPACKModule();

  const Aop = aslinearoperator(A);
  const [m, n] = Aop.shape;

  if (m !== n) {
    throw new DimensionMismatchError(`Matrix must be square, got ${m}x${n}`);
  }

  if (k <= 0) {
    throw new Error(`k must be positive, got ${k}`);
  }

  // For non-symmetric, k must be < n - 1
  if (k >= n - 1) {
    throw new Error(`k must be less than n-1 for non-symmetric, got k=${k}, n=${n}`);
  }

  const return_eigenvectors = options.return_eigenvectors ?? true;

  const params = new NonSymmetricARPACKParams(n, k, Aop.matvec.bind(Aop), options);

  try {
    params.iterate();
    const result = params.extract(return_eigenvectors);

    if (result.nconv < k) {
      throw new ARPACKNoConvergence(result.nconv, k);
    }

    return result;
  } finally {
    params.free();
  }
}

/**
 * Compute the truncated SVD of a sparse matrix.
 *
 * Implemented via eigenvalue decomposition of A'A or AA'.
 *
 * @param A - Sparse matrix
 * @param k - Number of singular values to compute (default: 6)
 * @param options - Solver options
 * @returns Singular values and optionally singular vectors
 *
 * @example
 * ```typescript
 * const A = csr_matrix([[1, 2], [3, 4], [5, 6]]);
 * const { u, s, vt } = await svds(A, 2);
 * // A ≈ u @ diag(s) @ vt
 * ```
 */
export async function svds(
  A: SparseMatrix | LinearOperatorLike,
  k: number = 6,
  options: SvdsOptions = {}
): Promise<SvdsResult> {
  await loadARPACKModule();

  const Aop = aslinearoperator(A);
  const [m, n] = Aop.shape;

  if (k <= 0) {
    throw new Error(`k must be positive, got ${k}`);
  }

  // Choose whether to form A'A (m >= n) or AA' (m < n)
  const useATA = m >= n;
  const dim = useATA ? n : m;

  if (k >= dim) {
    throw new Error(`k must be less than min(m, n), got k=${k}, min(m,n)=${dim}`);
  }

  // Create the symmetric operator for A'A or AA'
  const symmetricOp: ILinearOperator = useATA
    ? {
        shape: [n, n] as [number, number],
        dtype: 'float64',
        matvec: (x: Float64Array) => Aop.rmatvec(Aop.matvec(x)),
        rmatvec: (x: Float64Array) => Aop.rmatvec(Aop.matvec(x)),
      }
    : {
        shape: [m, m] as [number, number],
        dtype: 'float64',
        matvec: (x: Float64Array) => Aop.matvec(Aop.rmatvec(x)),
        rmatvec: (x: Float64Array) => Aop.matvec(Aop.rmatvec(x)),
      };

  // Compute eigenvalues of the symmetric operator
  const eigResult = await eigsh(symmetricOp, k, {
    which: 'LM',  // Largest magnitude = largest singular values
    tol: options.tol,
    maxiter: options.maxiter,
    ncv: options.ncv,
    v0: options.v0,
    return_eigenvectors: options.return_singular_vectors !== false,
  });

  // Singular values are sqrt of eigenvalues
  const s = new Float64Array(eigResult.nconv);
  for (let i = 0; i < eigResult.nconv; i++) {
    s[i] = Math.sqrt(Math.max(0, eigResult.values[i]));
  }

  // Compute singular vectors if requested
  let u: Float64Array | undefined;
  let vt: Float64Array | undefined;

  if (options.return_singular_vectors !== false && eigResult.vectors) {
    const needU = options.return_singular_vectors === true || options.return_singular_vectors === 'u';
    const needVt = options.return_singular_vectors === true || options.return_singular_vectors === 'vh';
    const nconv = eigResult.nconv;

    if (useATA) {
      // Eigenvectors are v (right singular vectors)
      if (needVt) {
        // vt is k x n, stored as eigenvectors are n x k column-major
        // We need to transpose
        vt = new Float64Array(nconv * n);
        for (let j = 0; j < nconv; j++) {
          for (let i = 0; i < n; i++) {
            vt[j * n + i] = eigResult.vectors![j * n + i];
          }
        }
      }
      if (needU) {
        // u = A * v / sigma
        u = new Float64Array(m * nconv);
        for (let j = 0; j < nconv; j++) {
          const v_j = new Float64Array(n);
          for (let i = 0; i < n; i++) {
            v_j[i] = eigResult.vectors![j * n + i];
          }
          const u_j = Aop.matvec(v_j);
          if (s[j] > 1e-10) {
            for (let i = 0; i < m; i++) {
              u[j * m + i] = u_j[i] / s[j];
            }
          }
        }
      }
    } else {
      // Eigenvectors are u (left singular vectors)
      if (needU) {
        u = new Float64Array(m * nconv);
        for (let j = 0; j < nconv; j++) {
          for (let i = 0; i < m; i++) {
            u[j * m + i] = eigResult.vectors![j * m + i];
          }
        }
      }
      if (needVt) {
        // v = A' * u / sigma
        vt = new Float64Array(nconv * n);
        for (let j = 0; j < nconv; j++) {
          const u_j = new Float64Array(m);
          for (let i = 0; i < m; i++) {
            u_j[i] = eigResult.vectors![j * m + i];
          }
          const v_j = Aop.rmatvec(u_j);
          if (s[j] > 1e-10) {
            for (let i = 0; i < n; i++) {
              vt[j * n + i] = v_j[i] / s[j];
            }
          }
        }
      }
    }
  }

  return { u, s, vt };
}
