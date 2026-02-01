/**
 * ARPACK WASM Module Type Definitions
 *
 * ARPACK (ARnoldi PACKage) is a collection of Fortran77 subroutines designed
 * to solve large-scale sparse eigenvalue problems using the Implicitly
 * Restarted Arnoldi Method (IRAM) / Lanczos Method.
 *
 * Key features:
 * - Computes a few eigenvalues/eigenvectors of large sparse matrices
 * - Memory efficient: only requires matrix-vector products, not the full matrix
 * - Supports symmetric and non-symmetric problems
 * - Reverse communication interface for maximum flexibility
 *
 * This module supports double precision real eigenvalue problems:
 * - dsaupd/dseupd: Symmetric eigenvalues (Lanczos method)
 * - dnaupd/dneupd: Non-symmetric eigenvalues (Arnoldi method)
 *
 * Typical usage pattern (reverse communication):
 * 1. Initialize parameters and call dsaupd/dnaupd with ido=0
 * 2. If ido != 99, perform requested operation:
 *    - ido=1 or -1: compute y = Op*x (matrix-vector product)
 *    - ido=2: compute y = B*x (for generalized problems)
 * 3. Call dsaupd/dnaupd again with same parameters
 * 4. Repeat steps 2-3 until ido=99 (convergence)
 * 5. Call dseupd/dneupd to extract eigenvalues/eigenvectors
 *
 * @see http://www.caam.rice.edu/software/ARPACK/
 */

/**
 * ARPACK WebAssembly Module Interface
 *
 * All functions use Fortran calling conventions:
 * - Trailing underscore in function names
 * - All parameters passed by pointer
 * - Column-major array storage
 */
export interface ARPACKModule {
  // ============================================================
  // DOUBLE PRECISION SYMMETRIC EIGENVALUE SOLVER (LANCZOS)
  // ============================================================
  // For symmetric matrices: A = A^T
  // Uses the Lanczos method which is more efficient than general Arnoldi
  // Eigenvalues are guaranteed to be real
  // ============================================================

  /**
   * DSAUPD - Symmetric Arnoldi Update Iteration
   *
   * Reverse communication interface for computing eigenvalues and eigenvectors
   * of a real symmetric operator using the Implicitly Restarted Lanczos Method.
   *
   * This function implements one step of the Lanczos iteration with implicit
   * restarts. It must be called repeatedly in a loop until convergence (ido=99).
   *
   * Problem types:
   * - Standard:    A*x = lambda*x       (bmat='I')
   * - Generalized: A*x = lambda*B*x     (bmat='G', B symmetric positive definite)
   *
   * Reverse communication actions (based on ido output):
   * - ido = -1: Compute y = Op*x where Op = inv(B)*A (first call)
   * - ido =  1: Compute y = Op*x, but B*x is available in workd(ipntr(3))
   * - ido =  2: Compute y = B*x
   * - ido = 99: Converged, call dseupd to extract results
   *
   * @param idoPtr - Reverse communication flag (input/output)
   *                 Input: 0 for initial call, previous value for subsequent calls
   *                 Output: action code (-1, 1, 2, or 99)
   * @param bmatPtr - Problem type (input): 'I'=standard, 'G'=generalized
   * @param nPtr - Dimension of the eigenproblem (input)
   * @param whichPtr - Which eigenvalues to compute (input):
   *                   'LA'=largest algebraic, 'SA'=smallest algebraic,
   *                   'LM'=largest magnitude, 'SM'=smallest magnitude,
   *                   'BE'=both ends (half from each end)
   * @param nevPtr - Number of eigenvalues requested (input, 0 < nev < n)
   * @param tolPtr - Convergence tolerance (input, 0=machine precision)
   * @param residPtr - Initial/final residual vector, size n (input/output)
   * @param ncvPtr - Number of Lanczos vectors (input, nev < ncv <= n)
   *                 Recommended: ncv >= 2*nev
   * @param vPtr - Lanczos basis vectors, n×ncv column-major (output)
   * @param ldvPtr - Leading dimension of v, >= n (input)
   * @param iparamPtr - Integer parameters array of size 11 (input/output):
   *                    [0]=ishfts: shift strategy (1=exact shifts)
   *                    [2]=mxiter: max iterations
   *                    [3]=nb: block size (must be 1)
   *                    [4]=nconv: number of converged eigenvalues (output)
   *                    [6]=mode: 1=standard, 2=shift-invert, 3=Buckling, 4=Cayley
   * @param ipntrPtr - Pointer array of size 11 (output)
   *                   ipntr(1): location of x in workd for matrix-vector product
   *                   ipntr(2): location of y in workd to store result
   *                   ipntr(3): location of B*x when ido=1
   * @param workdPtr - Work array of size 3*n (workspace)
   * @param worklPtr - Work array of size ncv*(ncv+8) (workspace)
   * @param lworklPtr - Length of workl (input)
   * @param infoPtr - Status/error code (input/output)
   *                  Input: 0=random start, 1=user-supplied residual
   *                  Output: 0=success, negative=error
   */
  _dsaupd_(
    idoPtr: number,
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

  /**
   * DSEUPD - Symmetric Eigenvalue Post-Processing
   *
   * Extracts eigenvalues and (optionally) eigenvectors after dsaupd has
   * converged. This function returns the Ritz values (eigenvalue approximations)
   * and Ritz vectors (eigenvector approximations) that satisfy the convergence
   * criteria specified to dsaupd.
   *
   * Must be called after dsaupd returns with ido=99.
   *
   * For shift-invert mode (iparam[6]=2), the eigenvalues of the original
   * problem are recovered from the Ritz values using the shift sigma.
   *
   * @param rvecPtr - Compute eigenvectors? (input): 0=no, nonzero=yes
   * @param howmnyPtr - Eigenvector selection (input):
   *                    'A'=all nev eigenvectors,
   *                    'P'=compute Schur vectors,
   *                    'S'=use select array
   * @param selectPtr - Selection array of size ncv (input for howmny='S')
   *                    select[i] nonzero means include i-th Ritz vector
   * @param dPtr - Eigenvalues output array of size nev (output)
   *               Contains the nev converged eigenvalues
   * @param zPtr - Eigenvectors output, n×nev column-major (output)
   *               z(:,i) is the eigenvector for d(i)
   * @param ldzPtr - Leading dimension of z, >= n (input)
   * @param sigmaPtr - Shift value used in shift-invert mode (input)
   *                   Eigenvalues computed as: lambda = 1/theta + sigma
   * @param bmatPtr - Same as in dsaupd call
   * @param nPtr - Same as in dsaupd call
   * @param whichPtr - Same as in dsaupd call
   * @param nevPtr - Same as in dsaupd call
   * @param tolPtr - Same as in dsaupd call
   * @param residPtr - Same as in dsaupd call
   * @param ncvPtr - Same as in dsaupd call
   * @param vPtr - Same as in dsaupd call (preserved from dsaupd)
   * @param ldvPtr - Same as in dsaupd call
   * @param iparamPtr - Same as in dsaupd call
   * @param ipntrPtr - Same as in dsaupd call
   * @param workdPtr - Same as in dsaupd call
   * @param worklPtr - Same as in dsaupd call
   * @param lworklPtr - Same as in dsaupd call
   * @param infoPtr - Status code (output): 0=success, negative=error
   */
  _dseupd_(
    rvecPtr: number,
    howmnyPtr: number,
    selectPtr: number,
    dPtr: number,
    zPtr: number,
    ldzPtr: number,
    sigmaPtr: number,
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
  // DOUBLE PRECISION NON-SYMMETRIC EIGENVALUE SOLVER (ARNOLDI)
  // ============================================================
  // For general (non-symmetric) matrices
  // Uses the Arnoldi method
  // Eigenvalues may be complex (returned as real/imaginary pairs)
  // ============================================================

  /**
   * DNAUPD - Non-Symmetric Arnoldi Update Iteration
   *
   * Reverse communication interface for computing eigenvalues and eigenvectors
   * of a real non-symmetric operator using the Implicitly Restarted Arnoldi Method.
   *
   * Similar to dsaupd but for non-symmetric matrices. Eigenvalues may be complex
   * and are returned as (real, imaginary) pairs. Complex eigenvalues always
   * appear in conjugate pairs.
   *
   * Reverse communication actions (based on ido output):
   * - ido = -1: Compute y = Op*x (first call initialization)
   * - ido =  1: Compute y = Op*x, B*x available at ipntr(3) for generalized
   * - ido =  2: Compute y = B*x (for generalized problems)
   * - ido =  3: Compute shifts (advanced usage)
   * - ido = 99: Converged, call dneupd to extract results
   *
   * @param idoPtr - Reverse communication flag (input/output)
   * @param bmatPtr - Problem type: 'I'=standard, 'G'=generalized
   * @param nPtr - Dimension of the eigenproblem
   * @param whichPtr - Which eigenvalues to compute:
   *                   'LM'=largest magnitude, 'SM'=smallest magnitude,
   *                   'LR'=largest real part, 'SR'=smallest real part,
   *                   'LI'=largest imaginary part, 'SI'=smallest imaginary part
   * @param nevPtr - Number of eigenvalues requested (0 < nev < n-1)
   * @param tolPtr - Convergence tolerance (0=machine precision)
   * @param residPtr - Initial/final residual vector, size n
   * @param ncvPtr - Number of Arnoldi vectors (nev+2 <= ncv <= n)
   *                 Recommended: ncv >= 2*nev+1
   * @param vPtr - Arnoldi basis vectors, n×ncv column-major
   * @param ldvPtr - Leading dimension of v, >= n
   * @param iparamPtr - Integer parameters array of size 11
   * @param ipntrPtr - Pointer array of size 14 (larger than symmetric case)
   * @param workdPtr - Work array of size 3*n
   * @param worklPtr - Work array of size 3*ncv² + 6*ncv
   * @param lworklPtr - Length of workl
   * @param infoPtr - Status/error code
   */
  _dnaupd_(
    idoPtr: number,
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

  /**
   * DNEUPD - Non-Symmetric Eigenvalue Post-Processing
   *
   * Extracts eigenvalues and (optionally) eigenvectors after dnaupd has
   * converged. For non-symmetric problems, eigenvalues may be complex
   * and are returned as separate real and imaginary parts.
   *
   * Complex conjugate eigenvalue pairs share storage for eigenvectors:
   * - If eigenvalue j is real: z(:,j) is the real eigenvector
   * - If eigenvalue j and j+1 form a complex conjugate pair:
   *   z(:,j) is the real part, z(:,j+1) is the imaginary part
   *   Eigenvector for lambda_j = z(:,j) + i*z(:,j+1)
   *   Eigenvector for lambda_{j+1} = z(:,j) - i*z(:,j+1)
   *
   * Must be called after dnaupd returns with ido=99.
   *
   * @param rvecPtr - Compute eigenvectors? 0=no, nonzero=yes
   * @param howmnyPtr - Eigenvector selection: 'A'=all, 'P'=Schur vectors, 'S'=selected
   * @param selectPtr - Selection array of size ncv
   * @param drPtr - Real parts of eigenvalues, size nev+1 (output)
   *                Extra element for possible complex conjugate pair
   * @param diPtr - Imaginary parts of eigenvalues, size nev+1 (output)
   *                di=0 means real eigenvalue
   * @param zPtr - Eigenvectors output, n×(nev+1) column-major
   * @param ldzPtr - Leading dimension of z, >= n
   * @param sigmarPtr - Real part of shift (for shift-invert mode)
   * @param sigmaiPtr - Imaginary part of shift
   * @param workevPtr - Work array of size 3*ncv
   * @param bmatPtr - Same as in dnaupd call
   * @param nPtr - Same as in dnaupd call
   * @param whichPtr - Same as in dnaupd call
   * @param nevPtr - Same as in dnaupd call
   * @param tolPtr - Same as in dnaupd call
   * @param residPtr - Same as in dnaupd call
   * @param ncvPtr - Same as in dnaupd call
   * @param vPtr - Same as in dnaupd call
   * @param ldvPtr - Same as in dnaupd call
   * @param iparamPtr - Same as in dnaupd call
   * @param ipntrPtr - Same as in dnaupd call
   * @param workdPtr - Same as in dnaupd call
   * @param worklPtr - Same as in dnaupd call
   * @param lworklPtr - Same as in dnaupd call
   * @param infoPtr - Status code: 0=success, negative=error
   */
  _dneupd_(
    rvecPtr: number,
    howmnyPtr: number,
    selectPtr: number,
    drPtr: number,
    diPtr: number,
    zPtr: number,
    ldzPtr: number,
    sigmarPtr: number,
    sigmaiPtr: number,
    workevPtr: number,
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
  // COMPLEX DOUBLE PRECISION EIGENVALUE SOLVER (ARNOLDI)
  // ============================================================
  // For complex matrices (general or Hermitian)
  // Uses the Arnoldi method with complex arithmetic
  // All eigenvalues and eigenvectors are complex
  // ============================================================

  /**
   * ZNAUPD - Complex Non-Symmetric Arnoldi Update Iteration
   *
   * Reverse communication interface for computing eigenvalues and eigenvectors
   * of a complex (general or Hermitian) operator using the Implicitly Restarted
   * Arnoldi Method.
   *
   * Complex arrays are stored in interleaved format: [re0, im0, re1, im1, ...]
   * Each complex number takes 2 doubles (16 bytes).
   *
   * Problem types:
   * - Standard:    A*x = lambda*x       (bmat='I')
   * - Generalized: A*x = lambda*B*x     (bmat='G', B Hermitian positive definite)
   *
   * Reverse communication actions (based on ido output):
   * - ido = -1: Compute y = Op*x (first call initialization)
   * - ido =  1: Compute y = Op*x, B*x available at ipntr(3) for generalized
   * - ido =  2: Compute y = B*x (for generalized problems)
   * - ido =  3: Compute shifts (advanced usage)
   * - ido = 99: Converged, call zneupd to extract results
   *
   * @param idoPtr - Reverse communication flag (input/output)
   * @param bmatPtr - Problem type: 'I'=standard, 'G'=generalized
   * @param nPtr - Dimension of the eigenproblem
   * @param whichPtr - Which eigenvalues to compute:
   *                   'LM'=largest magnitude, 'SM'=smallest magnitude,
   *                   'LR'=largest real part, 'SR'=smallest real part,
   *                   'LI'=largest imaginary part, 'SI'=smallest imaginary part
   * @param nevPtr - Number of eigenvalues requested (0 < nev < n-1)
   * @param tolPtr - Convergence tolerance (0=machine precision)
   * @param residPtr - Initial/final residual vector, size 2*n (complex, interleaved)
   * @param ncvPtr - Number of Arnoldi vectors (nev+2 <= ncv <= n)
   * @param vPtr - Arnoldi basis vectors, 2*n*ncv (complex, interleaved, column-major)
   * @param ldvPtr - Leading dimension of v, >= n
   * @param iparamPtr - Integer parameters array of size 11
   * @param ipntrPtr - Pointer array of size 14
   * @param workdPtr - Complex work array of size 3*n (6*n doubles, interleaved)
   * @param worklPtr - Complex work array of size 3*ncv² + 5*ncv
   * @param lworklPtr - Length of workl (in complex numbers)
   * @param rworkPtr - Real work array of size ncv
   * @param infoPtr - Status/error code
   */
  _znaupd_(
    idoPtr: number,
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
    rworkPtr: number,
    infoPtr: number
  ): void;

  /**
   * ZNEUPD - Complex Eigenvalue Post-Processing
   *
   * Extracts eigenvalues and (optionally) eigenvectors after znaupd has
   * converged. All eigenvalues and eigenvectors are complex.
   *
   * Must be called after znaupd returns with ido=99.
   *
   * @param rvecPtr - Compute eigenvectors? 0=no, nonzero=yes
   * @param howmnyPtr - Eigenvector selection: 'A'=all, 'P'=Schur vectors, 'S'=selected
   * @param selectPtr - Selection array of size ncv (integer)
   * @param dPtr - Complex eigenvalues output, size 2*(nev+1) (interleaved)
   * @param zPtr - Complex eigenvectors output, 2*n*(nev+1) (interleaved, column-major)
   * @param ldzPtr - Leading dimension of z, >= n
   * @param sigmaPtr - Complex shift value, size 2 (real, imag)
   * @param workevPtr - Complex work array of size 2*ncv (4*ncv doubles)
   * @param bmatPtr - Same as in znaupd call
   * @param nPtr - Same as in znaupd call
   * @param whichPtr - Same as in znaupd call
   * @param nevPtr - Same as in znaupd call
   * @param tolPtr - Same as in znaupd call
   * @param residPtr - Same as in znaupd call
   * @param ncvPtr - Same as in znaupd call
   * @param vPtr - Same as in znaupd call
   * @param ldvPtr - Same as in znaupd call
   * @param iparamPtr - Same as in znaupd call
   * @param ipntrPtr - Same as in znaupd call
   * @param workdPtr - Same as in znaupd call
   * @param worklPtr - Same as in znaupd call
   * @param lworklPtr - Same as in znaupd call
   * @param rworkPtr - Same as in znaupd call
   * @param infoPtr - Status code: 0=success, negative=error
   */
  _zneupd_(
    rvecPtr: number,
    howmnyPtr: number,
    selectPtr: number,
    dPtr: number,
    zPtr: number,
    ldzPtr: number,
    sigmaPtr: number,
    workevPtr: number,
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
    rworkPtr: number,
    infoPtr: number
  ): void;

  // ============================================================
  // MEMORY MANAGEMENT
  // ============================================================

  /**
   * Allocate memory in the WASM heap.
   * @param size - Number of bytes to allocate
   * @returns Pointer to allocated memory, or 0 on failure
   */
  _malloc(size: number): number;

  /**
   * Free previously allocated memory.
   * @param ptr - Pointer returned by _malloc
   */
  _free(ptr: number): void;

  // ============================================================
  // HEAP VIEWS - Direct Access to WASM Memory
  // ============================================================

  /** Float64 view of WASM heap. Index = byteOffset / 8 for doubles. */
  HEAPF64: Float64Array;

  /** Float32 view of WASM heap. Index = byteOffset / 4 for floats. */
  HEAPF32: Float32Array;

  /** Int32 view of WASM heap. Index = byteOffset / 4 for 32-bit integers. */
  HEAP32: Int32Array;

  /** Int8 view of WASM heap. Index = byteOffset directly. */
  HEAP8: Int8Array;

  /** Uint8 view of WASM heap. Index = byteOffset directly. */
  HEAPU8: Uint8Array;

  // ============================================================
  // EMSCRIPTEN RUNTIME UTILITIES
  // ============================================================

  /**
   * Read a value from WASM memory at the specified pointer.
   * @param ptr - Memory address (byte offset)
   * @param type - Value type: 'i8', 'i16', 'i32', 'i64', 'float', 'double'
   * @returns The value at that address
   */
  getValue(ptr: number, type: string): number;

  /**
   * Write a value to WASM memory at the specified pointer.
   * @param ptr - Memory address (byte offset)
   * @param value - Value to write
   * @param type - Value type: 'i8', 'i16', 'i32', 'i64', 'float', 'double'
   */
  setValue(ptr: number, value: number, type: string): void;
}

/**
 * Options for configuring the ARPACK WASM module loader.
 */
export interface ARPACKModuleOptions {
  /**
   * Custom function to locate WASM files.
   * @param path - Filename being requested (e.g., 'arpack.wasm')
   * @param scriptDirectory - Directory of the JS loader script
   * @returns Full URL or path to the file
   */
  locateFile?: (path: string, scriptDirectory: string) => string;
}

/**
 * Factory function to create an ARPACK WASM module instance.
 *
 * @param options - Optional configuration
 * @returns Promise resolving to initialized ARPACK module
 *
 * @example
 * import createARPACKModule from 'arwasm';
 * const arpack = await createARPACKModule();
 * // Use reverse communication to find eigenvalues
 */
export type ARPACKModuleFactory = (
  options?: ARPACKModuleOptions
) => Promise<ARPACKModule>;

// ============================================================
// ARPACK ERROR MESSAGES
// ============================================================
// Error codes are returned in the 'info' parameter.
// Positive values indicate warnings, negative values indicate errors.
// ============================================================

/**
 * Error messages for DSAUPD (symmetric Arnoldi update).
 * Maps info return codes to human-readable descriptions.
 */
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

/**
 * Error messages for DNAUPD (non-symmetric Arnoldi update).
 * Maps info return codes to human-readable descriptions.
 */
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

/**
 * Error messages for DSEUPD (symmetric eigenvalue extraction).
 * Maps info return codes to human-readable descriptions.
 */
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

/**
 * Error messages for DNEUPD (non-symmetric eigenvalue extraction).
 * Maps info return codes to human-readable descriptions.
 */
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

/**
 * Valid 'which' parameter values for symmetric eigensolvers (dsaupd/dseupd).
 *
 * - 'LA': Largest Algebraic eigenvalues (most positive)
 * - 'SA': Smallest Algebraic eigenvalues (most negative)
 * - 'LM': Largest Magnitude eigenvalues (largest |λ|)
 * - 'SM': Smallest Magnitude eigenvalues (smallest |λ|)
 * - 'BE': Both Ends - half from each end of the spectrum
 */
export const WHICH_SYMMETRIC = ['LA', 'SA', 'LM', 'SM', 'BE'] as const;

/** Type for symmetric 'which' parameter values. */
export type WhichSymmetric = (typeof WHICH_SYMMETRIC)[number];

/**
 * Valid 'which' parameter values for non-symmetric eigensolvers (dnaupd/dneupd).
 *
 * - 'LM': Largest Magnitude eigenvalues (largest |λ|)
 * - 'SM': Smallest Magnitude eigenvalues (smallest |λ|)
 * - 'LR': Largest Real part eigenvalues
 * - 'SR': Smallest Real part eigenvalues
 * - 'LI': Largest Imaginary part eigenvalues
 * - 'SI': Smallest Imaginary part eigenvalues
 */
export const WHICH_NONSYMMETRIC = ['LM', 'SM', 'LR', 'SR', 'LI', 'SI'] as const;

/** Type for non-symmetric 'which' parameter values. */
export type WhichNonSymmetric = (typeof WHICH_NONSYMMETRIC)[number];

/**
 * Valid 'which' parameter values for complex eigensolvers (znaupd/zneupd).
 * Same as non-symmetric since complex eigenvalues can have any position in the complex plane.
 */
export const WHICH_COMPLEX = ['LM', 'SM', 'LR', 'SR', 'LI', 'SI'] as const;

/** Type for complex 'which' parameter values. */
export type WhichComplex = (typeof WHICH_COMPLEX)[number];

/**
 * Error messages for ZNAUPD (complex Arnoldi update).
 * Maps info return codes to human-readable descriptions.
 */
export const ZNAUPD_ERRORS: Record<number, string> = {
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
  [-10]: 'IPARAM(7) must be 1, 2, or 3.',
  [-11]: 'IPARAM(7) = 1 and BMAT = "G" are incompatible.',
  [-12]: 'IPARAM(1) must be 0 or 1.',
  [-9999]: 'Could not build an Arnoldi factorization.',
};

/**
 * Error messages for ZNEUPD (complex eigenvalue extraction).
 * Maps info return codes to human-readable descriptions.
 */
export const ZNEUPD_ERRORS: Record<number, string> = {
  0: 'Normal exit.',
  1: 'The Schur form computed by LAPACK routine could not be reordered.',
  [-1]: 'N must be positive.',
  [-2]: 'NEV must be positive.',
  [-3]: 'NCV must satisfy NEV+2 <= NCV <= N.',
  [-5]: 'WHICH must be one of "LM", "SM", "LR", "SR", "LI", "SI".',
  [-6]: 'BMAT must be "I" or "G".',
  [-7]: 'Length of WORKL not sufficient.',
  [-8]: 'Error return from LAPACK eigenvalue calculation.',
  [-9]: 'Error return from computation of eigenvectors.',
  [-10]: 'IPARAM(7) must be 1, 2, or 3.',
  [-11]: 'IPARAM(7) = 1 and BMAT = "G" are incompatible.',
  [-12]: 'HOWMNY = "S" not yet implemented.',
  [-13]: 'HOWMNY must be "A" or "P" if RVEC = true.',
  [-14]: 'ZNAUPD did not find any eigenvalues to sufficient accuracy.',
};
