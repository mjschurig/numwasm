/**
 * ARPACK WASM module type definitions
 *
 * ARPACK is a collection of Fortran77 subroutines designed to solve
 * large scale eigenvalue problems using the Implicitly Restarted
 * Arnoldi Method (IRAM).
 *
 * This module supports double precision real eigenvalue problems:
 * - dsaupd/dseupd: Symmetric eigenvalues (Lanczos)
 * - dnaupd/dneupd: Non-symmetric eigenvalues (Arnoldi)
 */

/**
 * ARPACK WASM module interface
 */
export interface ARPACKModule {
  // ============================================================
  // Double precision symmetric eigenvalue solver (Lanczos)
  // ============================================================

  /**
   * Symmetric Arnoldi update iteration.
   * Reverse communication interface for computing eigenvalues
   * of a real symmetric operator.
   *
   * Parameters are passed as pointers to WASM memory.
   */
  _dsaupd_(
    idoPtr: number,      // Reverse communication flag (in/out)
    bmatPtr: number,     // 'I' for standard, 'G' for generalized
    nPtr: number,        // Dimension of the eigenproblem
    whichPtr: number,    // Which eigenvalues ('LA', 'SA', 'LM', 'SM', 'BE')
    nevPtr: number,      // Number of eigenvalues to compute
    tolPtr: number,      // Tolerance (0 = machine precision)
    residPtr: number,    // Residual vector (n)
    ncvPtr: number,      // Number of Lanczos vectors (ncv > nev)
    vPtr: number,        // Lanczos basis (n x ncv, column-major)
    ldvPtr: number,      // Leading dimension of v (usually n)
    iparamPtr: number,   // Integer parameter array (11)
    ipntrPtr: number,    // Pointer array (11)
    workdPtr: number,    // Work array (3*n)
    worklPtr: number,    // Work array (ncv*(ncv+8))
    lworklPtr: number,   // Length of workl
    infoPtr: number      // Error/info flag (in/out)
  ): void;

  /**
   * Symmetric eigenvalue post-processing.
   * Extracts eigenvalues and eigenvectors after dsaupd converges.
   */
  _dseupd_(
    rvecPtr: number,     // Compute eigenvectors? (logical as int)
    howmnyPtr: number,   // 'A' for all, 'S' for selected
    selectPtr: number,   // Selection array (ncv)
    dPtr: number,        // Eigenvalues output (nev)
    zPtr: number,        // Eigenvectors output (n x nev, column-major)
    ldzPtr: number,      // Leading dimension of z
    sigmaPtr: number,    // Shift (for shift-invert mode)
    bmatPtr: number,     // Same as dsaupd
    nPtr: number,
    whichPtr: number,
    nevPtr: number,
    tolPtr: number,
    residPtr: number,
    ncvPtr: number,
    vPtr: number,
    ldvPtr: number,
    iparamPtr: number,
    ipntrPtr: number,
    workdPtr: number,
    worklPtr: number,
    lworklPtr: number,
    infoPtr: number
  ): void;

  // ============================================================
  // Double precision non-symmetric eigenvalue solver (Arnoldi)
  // ============================================================

  /**
   * Non-symmetric Arnoldi update iteration.
   * Reverse communication interface for computing eigenvalues
   * of a general (non-symmetric) real matrix.
   */
  _dnaupd_(
    idoPtr: number,
    bmatPtr: number,
    nPtr: number,
    whichPtr: number,    // 'LM', 'SM', 'LR', 'SR', 'LI', 'SI'
    nevPtr: number,
    tolPtr: number,
    residPtr: number,
    ncvPtr: number,
    vPtr: number,
    ldvPtr: number,
    iparamPtr: number,
    ipntrPtr: number,    // (14 integers for non-symmetric)
    workdPtr: number,
    worklPtr: number,    // (3*ncv^2 + 6*ncv)
    lworklPtr: number,
    infoPtr: number
  ): void;

  /**
   * Non-symmetric eigenvalue post-processing.
   * Extracts eigenvalues and eigenvectors after dnaupd converges.
   */
  _dneupd_(
    rvecPtr: number,
    howmnyPtr: number,
    selectPtr: number,
    drPtr: number,       // Real part of eigenvalues (nev+1)
    diPtr: number,       // Imaginary part of eigenvalues (nev+1)
    zPtr: number,        // Eigenvectors
    ldzPtr: number,
    sigmarPtr: number,   // Real part of shift
    sigmaiPtr: number,   // Imaginary part of shift
    workevPtr: number,   // Work array (3*ncv)
    bmatPtr: number,
    nPtr: number,
    whichPtr: number,
    nevPtr: number,
    tolPtr: number,
    residPtr: number,
    ncvPtr: number,
    vPtr: number,
    ldvPtr: number,
    iparamPtr: number,
    ipntrPtr: number,
    workdPtr: number,
    worklPtr: number,
    lworklPtr: number,
    infoPtr: number
  ): void;

  // ============================================================
  // Memory management
  // ============================================================

  _malloc(size: number): number;
  _free(ptr: number): void;

  // ============================================================
  // Heap views
  // ============================================================

  HEAPF64: Float64Array;
  HEAPF32: Float32Array;
  HEAP32: Int32Array;
  HEAP8: Int8Array;
  HEAPU8: Uint8Array;

  // ============================================================
  // Emscripten utilities
  // ============================================================

  getValue(ptr: number, type: string): number;
  setValue(ptr: number, value: number, type: string): void;
}

export interface ARPACKModuleOptions {
  locateFile?: (path: string, scriptDirectory: string) => string;
}

export type ARPACKModuleFactory = (
  options?: ARPACKModuleOptions
) => Promise<ARPACKModule>;

// ============================================================
// ARPACK Error Messages
// ============================================================

/** Error messages for dsaupd */
export const DSAUPD_ERRORS: Record<number, string> = {
  0: 'Normal exit.',
  1: 'Maximum number of iterations taken. Check NCV, WHICH, TOL.',
  2: 'No longer an informational error. Deprecated starting with release 2 of ARPACK.',
  3: 'No shifts could be applied during a cycle of the Implicitly restarted Arnoldi iteration.',
  [-1]: 'N must be positive.',
  [-2]: 'NEV must be positive.',
  [-3]: 'NCV must be greater than NEV and less than or equal to N.',
  [-4]: 'The maximum number of Arnoldi update iterations allowed must be greater than zero.',
  [-5]: 'WHICH must be one of "LM", "SM", "LA", "SA" or "BE".',
  [-6]: 'BMAT must be one of "I" or "G".',
  [-7]: 'Length of private work array WORKL is not sufficient.',
  [-8]: 'Error return from LAPACK eigenvalue calculation.',
  [-9]: 'Starting vector is zero.',
  [-10]: 'IPARAM(7) must be 1, 2, 3, 4, or 5.',
  [-11]: 'IPARAM(7) = 1 and BMAT = "G" are incompatible.',
  [-12]: 'IPARAM(1) must be equal to 0 or 1.',
  [-13]: 'NEV and WHICH = "BE" are incompatible.',
  [-9999]: 'Could not build an Arnoldi factorization.',
};

/** Error messages for dnaupd */
export const DNAUPD_ERRORS: Record<number, string> = {
  0: 'Normal exit.',
  1: 'Maximum number of iterations taken.',
  2: 'No longer an informational error. Deprecated.',
  3: 'No shifts could be applied during implicit restart.',
  [-1]: 'N must be positive.',
  [-2]: 'NEV must be positive.',
  [-3]: 'NCV must satisfy NEV+2 <= NCV <= N.',
  [-4]: 'Maximum number of iterations must be > 0.',
  [-5]: 'WHICH must be one of "LM", "SM", "LR", "SR", "LI", "SI".',
  [-6]: 'BMAT must be "I" or "G".',
  [-7]: 'Length of WORKL not sufficient.',
  [-8]: 'Error return from LAPACK eigenvalue calculation.',
  [-9]: 'Starting vector is zero.',
  [-10]: 'IPARAM(7) must be 1, 2, 3, or 4.',
  [-11]: 'IPARAM(7) = 1 and BMAT = "G" are incompatible.',
  [-12]: 'IPARAM(1) must be 0 or 1.',
  [-9999]: 'Could not build an Arnoldi factorization.',
};

/** Error messages for dseupd */
export const DSEUPD_ERRORS: Record<number, string> = {
  0: 'Normal exit.',
  [-1]: 'N must be positive.',
  [-2]: 'NEV must be positive.',
  [-3]: 'NCV must be greater than NEV and less than or equal to N.',
  [-5]: 'WHICH must be one of "LM", "SM", "LA", "SA" or "BE".',
  [-6]: 'BMAT must be "I" or "G".',
  [-7]: 'Length of WORKL not sufficient.',
  [-8]: 'Error return from LAPACK eigenvalue calculation.',
  [-9]: 'Error return from computation of eigenvectors.',
  [-10]: 'IPARAM(7) must be 1, 2, 3, 4, or 5.',
  [-11]: 'IPARAM(7) = 1 and BMAT = "G" are incompatible.',
  [-12]: 'HOWMNY = "S" not yet implemented.',
  [-13]: 'HOWMNY must be "A" or "P" if RVEC = true.',
  [-14]: 'DSAUPD did not find any eigenvalues to sufficient accuracy.',
  [-15]: 'DSEUPD received an erroneous value of HOWMNY.',
  [-16]: 'HOWMNY = "A" and IPARAM(7) = 1 or 2 are incompatible.',
  [-17]: 'Problem in computing eigenvectors.',
};

/** Error messages for dneupd */
export const DNEUPD_ERRORS: Record<number, string> = {
  0: 'Normal exit.',
  1: 'The Schur form computed by LAPACK routine dlahqr could not be reordered.',
  [-1]: 'N must be positive.',
  [-2]: 'NEV must be positive.',
  [-3]: 'NCV must satisfy NEV+2 <= NCV <= N.',
  [-5]: 'WHICH must be one of "LM", "SM", "LR", "SR", "LI", "SI".',
  [-6]: 'BMAT must be "I" or "G".',
  [-7]: 'Length of WORKL not sufficient.',
  [-8]: 'Error return from LAPACK eigenvalue calculation.',
  [-9]: 'Error return from computation of eigenvectors.',
  [-10]: 'IPARAM(7) must be 1, 2, 3, or 4.',
  [-11]: 'IPARAM(7) = 1 and BMAT = "G" are incompatible.',
  [-12]: 'HOWMNY = "S" not yet implemented.',
  [-13]: 'HOWMNY must be "A" or "P" if RVEC = true.',
  [-14]: 'DNAUPD did not find any eigenvalues to sufficient accuracy.',
};

/** Valid 'which' values for symmetric eigensolvers */
export const WHICH_SYMMETRIC = ['LA', 'SA', 'LM', 'SM', 'BE'] as const;
export type WhichSymmetric = (typeof WHICH_SYMMETRIC)[number];

/** Valid 'which' values for non-symmetric eigensolvers */
export const WHICH_NONSYMMETRIC = ['LM', 'SM', 'LR', 'SR', 'LI', 'SI'] as const;
export type WhichNonSymmetric = (typeof WHICH_NONSYMMETRIC)[number];
