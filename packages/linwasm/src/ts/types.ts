/**
 * LINPACK WASM Module Type Definitions
 *
 * LINPACK is a classic Fortran library for numerical linear algebra, developed
 * at Argonne National Laboratory in the 1970s-80s. It provides routines for:
 *
 * - Solving linear systems of equations (Ax = b)
 * - Matrix factorizations (LU, Cholesky, QR)
 * - Computing matrix inverses and determinants
 * - Estimating condition numbers
 * - Singular value decomposition (SVD)
 *
 * Matrix types supported:
 * - General (dense) matrices
 * - Banded matrices (efficient storage for matrices with limited bandwidth)
 * - Symmetric/Hermitian positive definite matrices
 * - Symmetric indefinite matrices
 * - Triangular matrices
 * - Tridiagonal matrices
 *
 * Precision variants:
 * - s: single precision real (float, 4 bytes)
 * - d: double precision real (double, 8 bytes)
 * - c: single precision complex (2 floats, 8 bytes)
 * - z: double precision complex (2 doubles, 16 bytes)
 *
 * Naming convention (e.g., DGEFA):
 * - D = double precision
 * - GE = general matrix
 * - FA = factor (LU factorization)
 *
 * All routines use Fortran calling conventions:
 * - Trailing underscore in function names
 * - All parameters passed by pointer (use WASM memory addresses)
 * - Column-major array storage
 *
 * @see http://www.netlib.org/linpack/
 */

/**
 * Base Emscripten module interface with memory management utilities.
 *
 * These methods are provided by Emscripten for interacting with the
 * WebAssembly memory heap and calling C functions.
 */
export interface EmscriptenModule {
  /** WebAssembly memory as Float64Array (for double precision) */
  HEAPF64: Float64Array;
  /** WebAssembly memory as Float32Array (for single precision) */
  HEAPF32: Float32Array;
  /** WebAssembly memory as Int32Array */
  HEAP32: Int32Array;
  /** WebAssembly memory as Uint8Array */
  HEAPU8: Uint8Array;

  /**
   * Allocate memory on the WebAssembly heap
   * @param size Number of bytes to allocate
   * @returns Pointer to the allocated memory
   */
  _malloc(size: number): number;

  /**
   * Free memory on the WebAssembly heap
   * @param ptr Pointer to memory to free
   */
  _free(ptr: number): void;

  /**
   * Get a value from the WebAssembly heap
   * @param ptr Pointer to read from
   * @param type Type of value ('i8', 'i16', 'i32', 'float', 'double')
   */
  getValue(ptr: number, type: string): number;

  /**
   * Set a value on the WebAssembly heap
   * @param ptr Pointer to write to
   * @param value Value to write
   * @param type Type of value ('i8', 'i16', 'i32', 'float', 'double')
   */
  setValue(ptr: number, value: number, type: string): void;

  /**
   * Call a C function
   * @param name Function name
   * @param returnType Return type
   * @param argTypes Argument types
   * @param args Arguments
   */
  ccall(
    name: string,
    returnType: string | null,
    argTypes: string[],
    args: (number | string)[]
  ): number | void;

  /**
   * Wrap a C function as a JavaScript function
   * @param name Function name
   * @param returnType Return type
   * @param argTypes Argument types
   */
  cwrap(
    name: string,
    returnType: string | null,
    argTypes: string[]
  ): (...args: number[]) => number | void;

  /**
   * Add a function to the WebAssembly table
   * @param func JavaScript function
   * @param sig Function signature
   */
  addFunction(func: Function, sig: string): number;

  /**
   * Remove a function from the WebAssembly table
   * @param ptr Function pointer
   */
  removeFunction(ptr: number): void;
}

/**
 * LINPACK WASM Module Interface
 *
 * All exported LINPACK routines using f2c calling conventions.
 * Parameters are WASM memory pointers - use module._malloc() and HEAP views.
 */
export interface LINPACKModule extends EmscriptenModule {
  // ============================================================================
  // DOUBLE PRECISION GENERAL MATRIX ROUTINES (DGE*)
  // ============================================================================
  // For general (dense) n×n matrices stored in column-major order.
  // Uses LU factorization with partial pivoting: P*A = L*U
  // ============================================================================

  /**
   * DGEFA - Factor a double precision general matrix.
   *
   * Computes the LU factorization of an n×n matrix A using Gaussian
   * elimination with partial pivoting: P*A = L*U
   *
   * The factors L and U overwrite A:
   * - L is unit lower triangular (diagonal = 1, stored below diagonal)
   * - U is upper triangular (stored on and above diagonal)
   *
   * @param a - Matrix A, n×n in column-major order (input/output)
   *            On output, contains L and U factors
   * @param lda - Leading dimension of A (≥ n)
   * @param n - Order of the matrix
   * @param ipvt - Pivot indices array of size n (output)
   *               ipvt[k] = row interchanged with row k
   * @param info - Status (output): 0=success, k>0=U(k,k)=0 (singular)
   */
  _dgefa_(a: number, lda: number, n: number, ipvt: number, info: number): void;

  /**
   * DGESL - Solve a double precision general system.
   *
   * Solves A*x = b or A'*x = b using the LU factorization from DGEFA.
   *
   * @param a - LU factors from dgefa (input)
   * @param lda - Leading dimension of A
   * @param n - Order of the matrix
   * @param ipvt - Pivot indices from dgefa (input)
   * @param b - Right-hand side vector (input), solution x (output)
   * @param job - 0=solve A*x=b, nonzero=solve A'*x=b (transpose)
   */
  _dgesl_(
    a: number,
    lda: number,
    n: number,
    ipvt: number,
    b: number,
    job: number
  ): void;

  /**
   * DGECO - Factor and estimate condition number.
   *
   * Computes LU factorization and estimates the reciprocal condition
   * number rcond = 1/||A|| * ||A^{-1}||. Small rcond indicates near-singularity.
   *
   * @param a - Matrix A (input), LU factors (output)
   * @param lda - Leading dimension of A
   * @param n - Order of the matrix
   * @param ipvt - Pivot indices (output)
   * @param rcond - Reciprocal condition number (output), 0 < rcond < 1
   *                rcond ≈ 0 means singular, rcond ≈ 1 means well-conditioned
   * @param z - Work array of size n
   */
  _dgeco_(
    a: number,
    lda: number,
    n: number,
    ipvt: number,
    rcond: number,
    z: number
  ): void;

  /**
   * DGEDI - Compute determinant and/or inverse.
   *
   * Uses the LU factorization from DGEFA or DGECO to compute
   * the determinant and/or inverse of the matrix.
   *
   * @param a - LU factors (input), inverse A^{-1} (output if job includes 1)
   * @param lda - Leading dimension of A
   * @param n - Order of the matrix
   * @param ipvt - Pivot indices from dgefa/dgeco
   * @param det - Determinant as det[0]*10^det[1] (output if job includes 10)
   * @param work - Work array of size n
   * @param job - 11=both, 01=inverse only, 10=determinant only
   */
  _dgedi_(
    a: number,
    lda: number,
    n: number,
    ipvt: number,
    det: number,
    work: number,
    job: number
  ): void;

  // ============================================================================
  // DOUBLE PRECISION BANDED MATRIX ROUTINES (DGB*)
  // ============================================================================
  // For banded matrices with ml subdiagonals and mu superdiagonals.
  // Storage: abd[i,j] = A[i-j+ml+mu, j] in band format
  // Leading dimension = 2*ml + mu + 1
  // ============================================================================

  /**
   * DGBFA - Factor a double precision general band matrix.
   *
   * Computes LU factorization of an n×n band matrix with ml subdiagonals
   * and mu superdiagonals.
   *
   * @param abd - Band matrix in band storage (input/output)
   *              Storage: abd(i,j) = A(i-j+ml+mu, j)
   * @param lda - Leading dimension of abd (≥ 2*ml + mu + 1)
   * @param n - Order of the matrix
   * @param ml - Number of diagonals below the main diagonal
   * @param mu - Number of diagonals above the main diagonal
   * @param ipvt - Pivot indices (output)
   * @param info - Status (output): 0=success, k>0=singular at k
   */
  _dgbfa_(
    abd: number,
    lda: number,
    n: number,
    ml: number,
    mu: number,
    ipvt: number,
    info: number
  ): void;

  /**
   * DGBSL - Solve a double precision general band system.
   *
   * Solves A*x = b or A'*x = b using LU factorization from DGBFA.
   *
   * @param abd - LU factors in band storage from dgbfa
   * @param lda - Leading dimension of abd
   * @param n - Order of the matrix
   * @param ml - Number of subdiagonals
   * @param mu - Number of superdiagonals
   * @param ipvt - Pivot indices from dgbfa
   * @param b - Right-hand side (input), solution (output)
   * @param job - 0=solve A*x=b, nonzero=solve A'*x=b
   */
  _dgbsl_(
    abd: number,
    lda: number,
    n: number,
    ml: number,
    mu: number,
    ipvt: number,
    b: number,
    job: number
  ): void;

  /**
   * DGBCO - Factor and estimate condition of band matrix.
   *
   * @param abd - Band matrix (input), LU factors (output)
   * @param lda - Leading dimension of abd
   * @param n - Order of the matrix
   * @param ml - Number of subdiagonals
   * @param mu - Number of superdiagonals
   * @param ipvt - Pivot indices (output)
   * @param rcond - Reciprocal condition number (output)
   * @param z - Work array of size n
   */
  _dgbco_(
    abd: number,
    lda: number,
    n: number,
    ml: number,
    mu: number,
    ipvt: number,
    rcond: number,
    z: number
  ): void;

  /**
   * DGBDI - Compute determinant of factored band matrix.
   *
   * @param abd - LU factors from dgbfa
   * @param lda - Leading dimension
   * @param n - Order of the matrix
   * @param ml - Number of subdiagonals
   * @param mu - Number of superdiagonals
   * @param ipvt - Pivot indices
   * @param det - Determinant as det[0]*10^det[1] (output)
   */
  _dgbdi_(abd: number, lda: number, n: number, ml: number, mu: number, ipvt: number, det: number): void;

  // ============================================================================
  // DOUBLE PRECISION POSITIVE DEFINITE MATRIX ROUTINES (DPO*)
  // ============================================================================
  // For symmetric positive definite matrices (A = A^T, x^T*A*x > 0 for all x≠0).
  // Uses Cholesky factorization: A = R^T * R (R upper triangular)
  // Only the upper triangle of A is accessed.
  // ============================================================================

  /**
   * DPOFA - Cholesky factorization of positive definite matrix.
   *
   * Computes the Cholesky factorization A = R^T * R where R is upper triangular.
   * Only the diagonal and upper triangle of A are used.
   *
   * @param a - Symmetric positive definite matrix (input), R factor (output)
   * @param lda - Leading dimension of A
   * @param n - Order of the matrix
   * @param info - Status (output): 0=success, k>0=not positive definite at k
   */
  _dpofa_(a: number, lda: number, n: number, info: number): void;

  /**
   * DPOSL - Solve positive definite system.
   *
   * Solves A*x = b using the Cholesky factorization from DPOFA.
   *
   * @param a - Cholesky factor R from dpofa
   * @param lda - Leading dimension of A
   * @param n - Order of the matrix
   * @param b - Right-hand side (input), solution (output)
   */
  _dposl_(a: number, lda: number, n: number, b: number): void;

  /**
   * DPOCO - Cholesky factorization with condition estimation.
   *
   * @param a - Positive definite matrix (input), R factor (output)
   * @param lda - Leading dimension
   * @param n - Order of the matrix
   * @param rcond - Reciprocal condition number (output)
   * @param z - Work array of size n
   * @param info - Status (output): 0=success, k>0=not positive definite
   */
  _dpoco_(a: number, lda: number, n: number, rcond: number, z: number, info: number): void;

  /**
   * DPODI - Determinant and/or inverse of positive definite matrix.
   *
   * @param a - Cholesky factor (input), inverse (output if job includes 1)
   * @param lda - Leading dimension
   * @param n - Order of the matrix
   * @param det - Determinant as det[0]*10^det[1] (output if job includes 10)
   * @param job - 11=both, 01=inverse only, 10=determinant only
   */
  _dpodi_(a: number, lda: number, n: number, det: number, job: number): void;

  // ============================================================================
  // DOUBLE PRECISION PACKED POSITIVE DEFINITE MATRIX ROUTINES (DPP*)
  // ============================================================================
  // For symmetric positive definite matrices in packed storage.
  // Packed storage: upper triangle stored by columns in 1D array of size n(n+1)/2
  // ap[k] = A[i,j] where k = i + j*(j-1)/2 for i ≤ j
  // ============================================================================

  /**
   * DPPFA - Cholesky factorization of packed positive definite matrix.
   *
   * @param ap - Packed upper triangle (input), R factor (output)
   * @param n - Order of the matrix
   * @param info - Status: 0=success, k>0=not positive definite
   */
  _dppfa_(ap: number, n: number, info: number): void;

  /**
   * DPPSL - Solve packed positive definite system.
   *
   * @param ap - Packed Cholesky factor from dppfa
   * @param n - Order of the matrix
   * @param b - Right-hand side (input), solution (output)
   */
  _dppsl_(ap: number, n: number, b: number): void;

  /**
   * DPPCO - Cholesky factorization with condition estimation (packed).
   *
   * @param ap - Packed upper triangle (input), R factor (output)
   * @param n - Order of the matrix
   * @param rcond - Reciprocal condition number (output)
   * @param z - Work array of size n
   * @param info - Status: 0=success, k>0=not positive definite
   */
  _dppco_(ap: number, n: number, rcond: number, z: number, info: number): void;

  /**
   * DPPDI - Determinant and/or inverse of packed positive definite matrix.
   *
   * @param ap - Packed Cholesky factor (input), inverse (output if job includes 1)
   * @param n - Order of the matrix
   * @param det - Determinant as det[0]*10^det[1] (output if job includes 10)
   * @param job - 11=both, 01=inverse only, 10=determinant only
   */
  _dppdi_(ap: number, n: number, det: number, job: number): void;

  // ============================================================================
  // DOUBLE PRECISION BANDED POSITIVE DEFINITE MATRIX ROUTINES (DPB*)
  // ============================================================================
  // For symmetric positive definite banded matrices.
  // Only the main diagonal and m superdiagonals are stored.
  // ============================================================================

  /**
   * DPBFA - Cholesky factorization of band positive definite matrix.
   *
   * @param abd - Band matrix storage (input), R factor (output)
   * @param lda - Leading dimension of abd (≥ m + 1)
   * @param n - Order of the matrix
   * @param m - Number of superdiagonals (bandwidth - 1)
   * @param info - Status: 0=success, k>0=not positive definite
   */
  _dpbfa_(abd: number, lda: number, n: number, m: number, info: number): void;

  /**
   * DPBSL - Solve band positive definite system.
   *
   * @param abd - Cholesky factor from dpbfa
   * @param lda - Leading dimension
   * @param n - Order of the matrix
   * @param m - Number of superdiagonals
   * @param b - Right-hand side (input), solution (output)
   */
  _dpbsl_(abd: number, lda: number, n: number, m: number, b: number): void;

  /**
   * DPBCO - Cholesky factorization with condition estimation (band).
   *
   * @param abd - Band matrix (input), R factor (output)
   * @param lda - Leading dimension
   * @param n - Order of the matrix
   * @param m - Number of superdiagonals
   * @param rcond - Reciprocal condition number (output)
   * @param z - Work array of size n
   * @param info - Status
   */
  _dpbco_(abd: number, lda: number, n: number, m: number, rcond: number, z: number, info: number): void;

  /**
   * DPBDI - Determinant of band positive definite matrix.
   *
   * @param abd - Cholesky factor from dpbfa
   * @param lda - Leading dimension
   * @param n - Order of the matrix
   * @param m - Number of superdiagonals
   * @param det - Determinant as det[0]*10^det[1] (output)
   */
  _dpbdi_(abd: number, lda: number, n: number, m: number, det: number): void;

  // ============================================================================
  // DOUBLE PRECISION SYMMETRIC INDEFINITE MATRIX ROUTINES (DSI*)
  // ============================================================================
  // For symmetric matrices that may be indefinite (not necessarily positive definite).
  // Uses Bunch-Kaufman diagonal pivoting factorization.
  // ============================================================================

  /**
   * DSIFA - Factor symmetric indefinite matrix.
   *
   * Computes the Bunch-Kaufman factorization: A = U*D*U' where
   * U is upper triangular and D is block diagonal (1×1 and 2×2 blocks).
   *
   * @param a - Symmetric matrix (input), factors U and D (output)
   * @param lda - Leading dimension
   * @param n - Order of the matrix
   * @param kpvt - Pivot information (output)
   * @param info - Status: 0=success, k>0=D(k,k)=0
   */
  _dsifa_(a: number, lda: number, n: number, kpvt: number, info: number): void;

  /**
   * DSISL - Solve symmetric indefinite system.
   *
   * @param a - Factors from dsifa
   * @param lda - Leading dimension
   * @param n - Order of the matrix
   * @param kpvt - Pivot information from dsifa
   * @param b - Right-hand side (input), solution (output)
   */
  _dsisl_(a: number, lda: number, n: number, kpvt: number, b: number): void;

  /**
   * DSICO - Factor with condition estimation (symmetric indefinite).
   *
   * @param a - Symmetric matrix (input), factors (output)
   * @param lda - Leading dimension
   * @param n - Order of the matrix
   * @param kpvt - Pivot information (output)
   * @param rcond - Reciprocal condition number (output)
   * @param z - Work array of size n
   */
  _dsico_(a: number, lda: number, n: number, kpvt: number, rcond: number, z: number): void;

  /**
   * DSIDI - Determinant, inertia, and/or inverse of symmetric indefinite matrix.
   *
   * @param a - Factors (input), inverse (output if job includes 1)
   * @param lda - Leading dimension
   * @param n - Order of the matrix
   * @param kpvt - Pivot information
   * @param det - Determinant as det[0]*10^det[1] (if job includes 10)
   * @param inert - Inertia: [neg, zero, pos] eigenvalue counts (if job includes 100)
   * @param work - Work array of size n
   * @param job - 111=all, 011=inverse+det, etc.
   */
  _dsidi_(a: number, lda: number, n: number, kpvt: number, det: number, inert: number, work: number, job: number): void;

  // ============================================================================
  // DOUBLE PRECISION PACKED SYMMETRIC INDEFINITE MATRIX ROUTINES (DSP*)
  // ============================================================================
  // Symmetric indefinite matrices in packed storage (upper triangle by columns).
  // ============================================================================

  /**
   * DSPFA - Factor packed symmetric indefinite matrix.
   * @param ap - Packed upper triangle (input), factors (output)
   * @param n - Order of the matrix
   * @param kpvt - Pivot information (output)
   * @param info - Status: 0=success, k>0=D(k,k)=0
   */
  _dspfa_(ap: number, n: number, kpvt: number, info: number): void;

  /**
   * DSPSL - Solve packed symmetric indefinite system.
   * @param ap - Factors from dspfa
   * @param n - Order of the matrix
   * @param kpvt - Pivot information
   * @param b - Right-hand side (input), solution (output)
   */
  _dspsl_(ap: number, n: number, kpvt: number, b: number): void;

  /**
   * DSPCO - Factor with condition estimation (packed symmetric indefinite).
   * @param ap - Packed matrix (input), factors (output)
   * @param n - Order of the matrix
   * @param kpvt - Pivot information (output)
   * @param rcond - Reciprocal condition number (output)
   * @param z - Work array of size n
   */
  _dspco_(ap: number, n: number, kpvt: number, rcond: number, z: number): void;

  /**
   * DSPDI - Determinant, inertia, and/or inverse (packed symmetric indefinite).
   * @param ap - Factors (input), inverse (output if requested)
   * @param n - Order of the matrix
   * @param kpvt - Pivot information
   * @param det - Determinant output
   * @param inert - Inertia output
   * @param work - Work array of size n
   * @param job - Computation selector
   */
  _dspdi_(ap: number, n: number, kpvt: number, det: number, inert: number, work: number, job: number): void;

  // ============================================================================
  // DOUBLE PRECISION TRIANGULAR MATRIX ROUTINES (DTR*)
  // ============================================================================
  // For triangular matrices (upper or lower).
  // ============================================================================

  /**
   * DTRCO - Estimate condition number of triangular matrix.
   * @param t - Triangular matrix
   * @param ldt - Leading dimension
   * @param n - Order of the matrix
   * @param rcond - Reciprocal condition number (output)
   * @param z - Work array of size n
   * @param job - 0=lower triangular, nonzero=upper triangular
   */
  _dtrco_(t: number, ldt: number, n: number, rcond: number, z: number, job: number): void;

  /**
   * DTRDI - Determinant and/or inverse of triangular matrix.
   * @param t - Triangular matrix (input), inverse (output if requested)
   * @param ldt - Leading dimension
   * @param n - Order of the matrix
   * @param det - Determinant output
   * @param job - (10*triangle + inverse): triangle 0=lower/1=upper, inverse 0=no/1=yes
   * @param info - Status: 0=success, k>0=singular at k
   */
  _dtrdi_(t: number, ldt: number, n: number, det: number, job: number, info: number): void;

  /**
   * DTRSL - Solve triangular system.
   * @param t - Triangular matrix
   * @param ldt - Leading dimension
   * @param n - Order of the matrix
   * @param b - Right-hand side (input), solution (output)
   * @param job - 00=T*x=b lower, 01=T'*x=b lower, 10=T*x=b upper, 11=T'*x=b upper
   * @param info - Status: 0=success, k>0=singular at k
   */
  _dtrsl_(t: number, ldt: number, n: number, b: number, job: number, info: number): void;

  // ============================================================================
  // DOUBLE PRECISION TRIDIAGONAL MATRIX ROUTINES (DGT*, DPT*)
  // ============================================================================
  // For tridiagonal matrices (efficient O(n) solvers).
  // Storage: c=subdiagonal, d=diagonal, e=superdiagonal
  // ============================================================================

  /**
   * DGTSL - Solve general tridiagonal system.
   *
   * Solves the system A*x = b where A is tridiagonal.
   * Uses Gaussian elimination with partial pivoting.
   *
   * @param n - Order of the matrix
   * @param c - Subdiagonal elements (n-1), modified on output
   * @param d - Diagonal elements (n), modified on output
   * @param e - Superdiagonal elements (n-1), modified on output
   * @param b - Right-hand side (input), solution (output)
   * @param info - Status: 0=success, k>0=singular at k
   */
  _dgtsl_(n: number, c: number, d: number, e: number, b: number, info: number): void;

  /**
   * DPTSL - Solve positive definite tridiagonal system.
   *
   * Solves A*x = b where A is symmetric positive definite tridiagonal.
   * Uses Cholesky factorization.
   *
   * @param n - Order of the matrix
   * @param d - Diagonal elements (n), modified on output
   * @param e - Off-diagonal elements (n-1), modified on output
   * @param b - Right-hand side (input), solution (output)
   */
  _dptsl_(n: number, d: number, e: number, b: number): void;

  // ============================================================================
  // DOUBLE PRECISION QR DECOMPOSITION ROUTINES (DQR*)
  // ============================================================================
  // Computes A = Q*R where Q is orthogonal and R is upper triangular.
  // Supports column pivoting: A*P = Q*R for better numerical stability.
  // ============================================================================

  /**
   * DQRDC - Compute QR decomposition with optional column pivoting.
   *
   * Computes the QR factorization of an n×p matrix X.
   * Q is represented implicitly via Householder vectors.
   *
   * @param x - Matrix X (input), R and Householder vectors (output)
   * @param ldx - Leading dimension of X
   * @param n - Number of rows
   * @param p - Number of columns
   * @param qraux - Auxiliary output for Q reconstruction (size p)
   * @param jpvt - Pivot indices (input/output), 0=free, >0=first, <0=last
   * @param work - Work array of size p
   * @param job - 0=no pivoting, nonzero=pivoting
   */
  _dqrdc_(x: number, ldx: number, n: number, p: number, qraux: number, jpvt: number, work: number, job: number): void;

  /**
   * DQRSL - Apply QR factorization to compute various quantities.
   *
   * Uses the QR factorization from DQRDC to compute Q*y, Q'*y,
   * the least squares solution, residuals, etc.
   *
   * @param x - QR factors from dqrdc
   * @param ldx - Leading dimension
   * @param n - Number of rows of original X
   * @param k - Number of columns to use (k ≤ min(n,p))
   * @param qraux - Auxiliary information from dqrdc
   * @param y - Vector to apply Q to (size n)
   * @param qy - Q*y output (if job includes 10000)
   * @param qty - Q'*y output (if job includes 1000)
   * @param b - Least squares solution (if job includes 100)
   * @param rsd - Residuals y - X*b (if job includes 10)
   * @param xb - Fitted values X*b (if job includes 1)
   * @param job - Bit flags selecting which outputs to compute
   * @param info - Status: 0=success, k>0=rank deficient at k
   */
  _dqrsl_(x: number, ldx: number, n: number, k: number, qraux: number, y: number, qy: number, qty: number, b: number, rsd: number, xb: number, job: number, info: number): void;

  // ============================================================================
  // DOUBLE PRECISION SINGULAR VALUE DECOMPOSITION (DSVDC)
  // ============================================================================
  // Computes A = U * S * V' where U and V are orthogonal, S is diagonal.
  // ============================================================================

  /**
   * DSVDC - Compute singular value decomposition.
   *
   * Computes the SVD of an n×p matrix X: X = U * S * V'
   * where U is n×min(n,p), S is diagonal (singular values), V is p×p.
   *
   * @param x - Matrix X (input), destroyed on output
   * @param ldx - Leading dimension of X
   * @param n - Number of rows
   * @param p - Number of columns
   * @param s - Singular values in descending order (output, size min(n+1,p))
   * @param e - Superdiagonal of bidiagonal form (output, for debugging)
   * @param u - Left singular vectors (output, n×n or n×min(n,p))
   * @param ldu - Leading dimension of U
   * @param v - Right singular vectors (output, p×p)
   * @param ldv - Leading dimension of V
   * @param work - Work array of size n
   * @param job - Controls what is computed:
   *              job = abc where:
   *              a=0: don't compute U, a=1: compute min(n,p) columns of U
   *              a≥2: compute all n columns of U
   *              b=0: don't compute V, b≠0: compute V
   *              c=0: no permutation, c≠0: permute columns
   * @param info - Status: 0=success, k>0=k singular values not computed
   */
  _dsvdc_(x: number, ldx: number, n: number, p: number, s: number, e: number, u: number, ldu: number, v: number, ldv: number, work: number, job: number, info: number): void;

  // ============================================================================
  // DOUBLE PRECISION CHOLESKY UPDATE ROUTINES (DCH*)
  // ============================================================================
  // For updating Cholesky factorizations when the matrix changes.
  // ============================================================================

  /**
   * DCHDC - Cholesky decomposition with column pivoting.
   * @param r - Matrix (input/output)
   * @param ldr - Leading dimension
   * @param p - Order of the matrix
   * @param x - Pivot permutation
   * @param z - Optional auxiliary matrix
   * @param ldz - Leading dimension of z
   * @param nz - Number of columns in z
   * @param c - Cosines of rotations (output)
   * @param s - Sines of rotations (output)
   * @param job - Control flag
   */
  _dchdc_(r: number, ldr: number, p: number, x: number, z: number, ldz: number, nz: number, c: number, s: number, job: number): void;

  /**
   * DCHDD - Cholesky rank-1 downdate.
   * Updates R to the Cholesky factor of A - x*x' where R'*R = A.
   * @param r - Cholesky factor (input/output)
   * @param ldr - Leading dimension
   * @param p - Order of the matrix
   * @param x - Vector to downdate by
   * @param z - Optional auxiliary matrix
   * @param ldz - Leading dimension of z
   * @param nz - Number of columns in z
   * @param y - Additional vector for y updates
   * @param rho - Scalar for rho updates
   * @param c - Cosines (output)
   * @param s - Sines (output)
   * @param info - Status: 0=success, -1=downdate would destroy positive definiteness
   */
  _dchdd_(r: number, ldr: number, p: number, x: number, z: number, ldz: number, nz: number, y: number, rho: number, c: number, s: number, info: number): void;

  /**
   * DCHEX - Cholesky exchange (reorder columns).
   * @param r - Cholesky factor (input/output)
   * @param ldr - Leading dimension
   * @param p - Order of the matrix
   * @param k - First column to exchange
   * @param l - Last column to exchange
   * @param z - Optional auxiliary matrix
   * @param ldz - Leading dimension of z
   * @param nz - Number of columns in z
   * @param c - Cosines (output)
   * @param s - Sines (output)
   * @param job - 1=right circular shift, 2=left circular shift
   */
  _dchex_(r: number, ldr: number, p: number, k: number, l: number, z: number, ldz: number, nz: number, c: number, s: number, job: number): void;

  /**
   * DCHUD - Cholesky rank-1 update.
   * Updates R to the Cholesky factor of A + x*x' where R'*R = A.
   * @param r - Cholesky factor (input/output)
   * @param ldr - Leading dimension
   * @param p - Order of the matrix
   * @param x - Vector to update by
   * @param z - Optional auxiliary matrix
   * @param ldz - Leading dimension of z
   * @param nz - Number of columns in z
   * @param y - Additional vector
   * @param rho - Scalar
   * @param c - Cosines (output)
   * @param s - Sines (output)
   */
  _dchud_(r: number, ldr: number, p: number, x: number, z: number, ldz: number, nz: number, y: number, rho: number, c: number, s: number): void;

  // ============================================================================
  // SINGLE PRECISION ROUTINES (SGE*, SGB*, SPO*, SQR*, SSVDC)
  // ============================================================================
  // Same algorithms as double precision but using float (4 bytes).
  // Use HEAPF32 for array access. ~7 decimal digits precision.
  // ============================================================================

  /** SGEFA - Factor single precision general matrix */
  _sgefa_(a: number, lda: number, n: number, ipvt: number, info: number): void;
  /** SGESL - Solve single precision general system */
  _sgesl_(a: number, lda: number, n: number, ipvt: number, b: number, job: number): void;
  /** SGECO - Factor with condition estimation (single precision) */
  _sgeco_(a: number, lda: number, n: number, ipvt: number, rcond: number, z: number): void;
  /** SGEDI - Determinant/inverse of general matrix (single precision) */
  _sgedi_(a: number, lda: number, n: number, ipvt: number, det: number, work: number, job: number): void;

  /** SGBFA - Factor band matrix (single precision) */
  _sgbfa_(abd: number, lda: number, n: number, ml: number, mu: number, ipvt: number, info: number): void;
  /** SGBSL - Solve band system (single precision) */
  _sgbsl_(abd: number, lda: number, n: number, ml: number, mu: number, ipvt: number, b: number, job: number): void;
  /** SGBCO - Band matrix condition estimation (single precision) */
  _sgbco_(abd: number, lda: number, n: number, ml: number, mu: number, ipvt: number, rcond: number, z: number): void;
  /** SGBDI - Band matrix determinant (single precision) */
  _sgbdi_(abd: number, lda: number, n: number, ml: number, mu: number, ipvt: number, det: number): void;

  /** SPOFA - Cholesky factorization (single precision) */
  _spofa_(a: number, lda: number, n: number, info: number): void;
  /** SPOSL - Solve positive definite system (single precision) */
  _sposl_(a: number, lda: number, n: number, b: number): void;
  /** SPOCO - Cholesky with condition estimation (single precision) */
  _spoco_(a: number, lda: number, n: number, rcond: number, z: number, info: number): void;
  /** SPODI - Determinant/inverse of positive definite (single precision) */
  _spodi_(a: number, lda: number, n: number, det: number, job: number): void;

  /** SQRDC - QR decomposition (single precision) */
  _sqrdc_(x: number, ldx: number, n: number, p: number, qraux: number, jpvt: number, work: number, job: number): void;
  /** SQRSL - Apply QR factorization (single precision) */
  _sqrsl_(x: number, ldx: number, n: number, k: number, qraux: number, y: number, qy: number, qty: number, b: number, rsd: number, xb: number, job: number, info: number): void;

  /** SSVDC - Singular value decomposition (single precision) */
  _ssvdc_(x: number, ldx: number, n: number, p: number, s: number, e: number, u: number, ldu: number, v: number, ldv: number, work: number, job: number, info: number): void;

  // ============================================================================
  // COMPLEX SINGLE PRECISION ROUTINES (CGE*, CPO*, CQR*, CSVDC)
  // ============================================================================
  // Complex single: each element is 2 floats (real, imag) = 8 bytes.
  // Use HEAPF32 with pairs of values for complex numbers.
  // ============================================================================

  /** CGEFA - Factor complex single general matrix */
  _cgefa_(a: number, lda: number, n: number, ipvt: number, info: number): void;
  /** CGESL - Solve complex single general system */
  _cgesl_(a: number, lda: number, n: number, ipvt: number, b: number, job: number): void;
  /** CGECO - Factor with condition estimation (complex single) */
  _cgeco_(a: number, lda: number, n: number, ipvt: number, rcond: number, z: number): void;
  /** CGEDI - Determinant/inverse (complex single) */
  _cgedi_(a: number, lda: number, n: number, ipvt: number, det: number, work: number, job: number): void;

  /** CPOFA - Cholesky factorization (complex single Hermitian) */
  _cpofa_(a: number, lda: number, n: number, info: number): void;
  /** CPOSL - Solve Hermitian positive definite system (complex single) */
  _cposl_(a: number, lda: number, n: number, b: number): void;

  /** CQRDC - QR decomposition (complex single) */
  _cqrdc_(x: number, ldx: number, n: number, p: number, qraux: number, jpvt: number, work: number, job: number): void;
  /** CSVDC - Singular value decomposition (complex single) */
  _csvdc_(x: number, ldx: number, n: number, p: number, s: number, e: number, u: number, ldu: number, v: number, ldv: number, work: number, job: number, info: number): void;

  // ============================================================================
  // COMPLEX DOUBLE PRECISION ROUTINES (ZGE*, ZPO*, ZQR*, ZSVDC)
  // ============================================================================
  // Complex double: each element is 2 doubles (real, imag) = 16 bytes.
  // Use HEAPF64 with pairs of values for complex numbers.
  // ============================================================================

  /** ZGEFA - Factor complex double general matrix */
  _zgefa_(a: number, lda: number, n: number, ipvt: number, info: number): void;
  /** ZGESL - Solve complex double general system */
  _zgesl_(a: number, lda: number, n: number, ipvt: number, b: number, job: number): void;
  /** ZGECO - Factor with condition estimation (complex double) */
  _zgeco_(a: number, lda: number, n: number, ipvt: number, rcond: number, z: number): void;
  /** ZGEDI - Determinant/inverse (complex double) */
  _zgedi_(a: number, lda: number, n: number, ipvt: number, det: number, work: number, job: number): void;

  /** ZPOFA - Cholesky factorization (complex double Hermitian) */
  _zpofa_(a: number, lda: number, n: number, info: number): void;
  /** ZPOSL - Solve Hermitian positive definite system (complex double) */
  _zposl_(a: number, lda: number, n: number, b: number): void;

  /** ZQRDC - QR decomposition (complex double) */
  _zqrdc_(x: number, ldx: number, n: number, p: number, qraux: number, jpvt: number, work: number, job: number): void;
  /** ZSVDC - Singular value decomposition (complex double) */
  _zsvdc_(x: number, ldx: number, n: number, p: number, s: number, e: number, u: number, ldu: number, v: number, ldv: number, work: number, job: number, info: number): void;

  // ============================================================================
  // BLAS LEVEL 1 ROUTINES - VECTOR OPERATIONS
  // ============================================================================
  // Basic vector operations: copy, scale, add, dot product, norms, rotations.
  // incx/incy = stride between elements (1 for contiguous, negative for reverse)
  // ============================================================================

  // Double precision BLAS Level 1
  /** DAXPY: y := alpha*x + y */
  _daxpy_(n: number, da: number, dx: number, incx: number, dy: number, incy: number): void;
  /** DCOPY: y := x */
  _dcopy_(n: number, dx: number, incx: number, dy: number, incy: number): void;
  /** DSCAL: x := alpha*x */
  _dscal_(n: number, da: number, dx: number, incx: number): void;
  /** DSWAP: swap x and y */
  _dswap_(n: number, dx: number, incx: number, dy: number, incy: number): void;
  /** DDOT: returns x'*y (dot product) */
  _ddot_(n: number, dx: number, incx: number, dy: number, incy: number): number;
  /** DASUM: returns sum(|x_i|) (1-norm of x) */
  _dasum_(n: number, dx: number, incx: number): number;
  /** DNRM2: returns sqrt(x'*x) (Euclidean norm) */
  _dnrm2_(n: number, dx: number, incx: number): number;
  /** DROT: apply Givens rotation [c s; -s c] to (x,y) */
  _drot_(n: number, dx: number, incx: number, dy: number, incy: number, c: number, s: number): void;
  /** DROTG: construct Givens rotation parameters from (a,b) */
  _drotg_(da: number, db: number, c: number, s: number): void;
  /** IDAMAX: returns index of max |x_i| (1-based) */
  _idamax_(n: number, dx: number, incx: number): number;

  // Single precision BLAS Level 1
  /** SAXPY: y := alpha*x + y (single) */
  _saxpy_(n: number, sa: number, sx: number, incx: number, sy: number, incy: number): void;
  /** SCOPY: y := x (single) */
  _scopy_(n: number, sx: number, incx: number, sy: number, incy: number): void;
  /** SSCAL: x := alpha*x (single) */
  _sscal_(n: number, sa: number, sx: number, incx: number): void;
  /** SSWAP: swap x and y (single) */
  _sswap_(n: number, sx: number, incx: number, sy: number, incy: number): void;
  /** SDOT: returns x'*y (single) */
  _sdot_(n: number, sx: number, incx: number, sy: number, incy: number): number;
  /** SASUM: returns sum(|x_i|) (single) */
  _sasum_(n: number, sx: number, incx: number): number;
  /** SNRM2: returns sqrt(x'*x) (single) */
  _snrm2_(n: number, sx: number, incx: number): number;
  /** SROT: apply Givens rotation (single) */
  _srot_(n: number, sx: number, incx: number, sy: number, incy: number, c: number, s: number): void;
  /** SROTG: construct Givens rotation (single) */
  _srotg_(sa: number, sb: number, c: number, s: number): void;
  /** ISAMAX: returns index of max |x_i| (single) */
  _isamax_(n: number, sx: number, incx: number): number;

  // Complex single BLAS Level 1
  /** CAXPY: y := alpha*x + y (complex single) */
  _caxpy_(n: number, ca: number, cx: number, incx: number, cy: number, incy: number): void;
  /** CCOPY: y := x (complex single) */
  _ccopy_(n: number, cx: number, incx: number, cy: number, incy: number): void;
  /** CSCAL: x := alpha*x (complex single) */
  _cscal_(n: number, ca: number, cx: number, incx: number): void;
  /** CSWAP: swap x and y (complex single) */
  _cswap_(n: number, cx: number, incx: number, cy: number, incy: number): void;
  /** CDOTC: returns conj(x)'*y (conjugated dot product) */
  _cdotc_(n: number, cx: number, incx: number, cy: number, incy: number): number;
  /** CDOTU: returns x'*y (unconjugated dot product) */
  _cdotu_(n: number, cx: number, incx: number, cy: number, incy: number): number;
  /** SCNRM2: returns sqrt(conj(x)'*x) (complex single Euclidean norm) */
  _scnrm2_(n: number, cx: number, incx: number): number;
  /** ICAMAX: returns index of max |x_i| (complex single) */
  _icamax_(n: number, cx: number, incx: number): number;

  // Complex double BLAS Level 1
  /** ZAXPY: y := alpha*x + y (complex double) */
  _zaxpy_(n: number, za: number, zx: number, incx: number, zy: number, incy: number): void;
  /** ZCOPY: y := x (complex double) */
  _zcopy_(n: number, zx: number, incx: number, zy: number, incy: number): void;
  /** ZSCAL: x := alpha*x (complex double) */
  _zscal_(n: number, za: number, zx: number, incx: number): void;
  /** ZSWAP: swap x and y (complex double) */
  _zswap_(n: number, zx: number, incx: number, zy: number, incy: number): void;
  /** ZDOTC: returns conj(x)'*y (complex double conjugated) */
  _zdotc_(n: number, zx: number, incx: number, zy: number, incy: number): number;
  /** ZDOTU: returns x'*y (complex double unconjugated) */
  _zdotu_(n: number, zx: number, incx: number, zy: number, incy: number): number;
  /** DZNRM2: returns sqrt(conj(x)'*x) (complex double Euclidean norm) */
  _dznrm2_(n: number, zx: number, incx: number): number;
  /** IZAMAX: returns index of max |x_i| (complex double) */
  _izamax_(n: number, zx: number, incx: number): number;
}

/**
 * Factory function that creates a LINPACK WASM module instance.
 *
 * @param config - Optional Emscripten configuration options
 * @returns Promise resolving to initialized LINPACK module
 *
 * @example
 * import createLINPACKModule from 'linwasm';
 *
 * const linpack = await createLINPACKModule();
 *
 * // Solve a 3x3 system Ax = b
 * const n = 3;
 * const aPtr = linpack._malloc(n * n * 8);
 * const bPtr = linpack._malloc(n * 8);
 * const ipvtPtr = linpack._malloc(n * 4);
 * const infoPtr = linpack._malloc(4);
 *
 * // ... populate A and b in column-major order ...
 *
 * // Factor A
 * linpack._dgefa_(aPtr, n, n, ipvtPtr, infoPtr);
 * // Solve
 * linpack._dgesl_(aPtr, n, n, ipvtPtr, bPtr, 0);
 *
 * // ... read solution from bPtr ...
 *
 * linpack._free(aPtr);
 * // ... free other memory ...
 */
/**
 * Configuration options for the LINPACK module factory.
 */
export interface LINPACKModuleConfig extends Partial<EmscriptenModule> {
  /**
   * Custom function to locate WASM assets
   * @param path The path to the file being requested
   * @returns The resolved URL or path to the file
   */
  locateFile?: (path: string) => string;
}

export interface LINPACKModuleFactory {
  (config?: LINPACKModuleConfig): Promise<LINPACKModule>;
}
