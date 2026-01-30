/**
 * LAPACK WebAssembly Module Types
 *
 * LAPACK (Linear Algebra PACKage) provides routines for:
 * - Solving systems of linear equations
 * - Eigenvalue problems
 * - Singular value decomposition
 * - Matrix factorizations (LU, Cholesky, QR)
 *
 * This module bundles BLAS internally, so no separate BLAS module is needed.
 *
 * All functions use the Fortran calling convention:
 * - Parameters are passed by pointer (number = pointer to WASM memory)
 * - Function names have trailing underscore (e.g., dgetrf_)
 * - Matrices are stored in column-major order
 * - Character arguments are passed as integers (e.g., 'N'=78, 'T'=84)
 */

export interface LAPACKModule {
  // ============================================================
  // Double Precision Linear Systems
  // ============================================================

  /**
   * DGETRF - LU factorization of a general M-by-N matrix A
   * A = P * L * U
   *
   * @param m - Number of rows (pointer to integer)
   * @param n - Number of columns (pointer to integer)
   * @param a - Matrix A, overwritten with L and U (pointer to double array)
   * @param lda - Leading dimension of A (pointer to integer)
   * @param ipiv - Pivot indices (pointer to integer array)
   * @param info - 0 = success, >0 = singular matrix (pointer to integer)
   */
  _dgetrf_(m: number, n: number, a: number, lda: number, ipiv: number, info: number): void;

  /**
   * DGETRS - Solve A*X = B using LU factorization from DGETRF
   *
   * @param trans - 'N' for A*X=B, 'T' for A^T*X=B (pointer to char)
   * @param n - Order of matrix A (pointer to integer)
   * @param nrhs - Number of right-hand sides (pointer to integer)
   * @param a - LU factorization from DGETRF (pointer to double array)
   * @param lda - Leading dimension of A (pointer to integer)
   * @param ipiv - Pivot indices from DGETRF (pointer to integer array)
   * @param b - Right-hand side matrix B, overwritten with solution X (pointer)
   * @param ldb - Leading dimension of B (pointer to integer)
   * @param info - 0 = success (pointer to integer)
   */
  _dgetrs_(
    trans: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    ipiv: number,
    b: number,
    ldb: number,
    info: number
  ): void;

  /**
   * DGETRI - Compute the inverse of a matrix using LU factorization from DGETRF
   *
   * Computes the inverse of a matrix A using the LU factorization computed by DGETRF.
   * The matrix A must have been previously factored by DGETRF.
   *
   * @param n - Order of the matrix A (pointer to integer)
   * @param a - On entry, L and U factors from DGETRF. On exit, the inverse of A (pointer to double array)
   * @param lda - Leading dimension of A (pointer to integer)
   * @param ipiv - Pivot indices from DGETRF (pointer to integer array)
   * @param work - Workspace array of dimension LWORK (pointer to double array)
   * @param lwork - Size of workspace. For optimal performance, LWORK >= N*NB where NB is optimal block size.
   *                If LWORK = -1, a workspace query is performed (pointer to integer)
   * @param info - 0 = success, <0 = illegal argument, >0 = singular matrix (pointer to integer)
   */
  _dgetri_(
    n: number,
    a: number,
    lda: number,
    ipiv: number,
    work: number,
    lwork: number,
    info: number
  ): void;

  /**
   * DGESV - Solve A*X = B for general matrices (driver routine)
   *
   * Computes the solution to a real system of linear equations A*X = B,
   * where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
   * Uses LU decomposition with partial pivoting and row interchanges.
   * This is a driver routine that combines DGETRF and DGETRS.
   *
   * @param n - Number of linear equations (order of A) (pointer to integer)
   * @param nrhs - Number of right-hand sides (columns of B) (pointer to integer)
   * @param a - On entry, the N-by-N coefficient matrix A. On exit, the L and U factors (pointer to double array)
   * @param lda - Leading dimension of A (pointer to integer)
   * @param ipiv - Pivot indices defining permutation matrix P (pointer to integer array)
   * @param b - On entry, N-by-NRHS right-hand side matrix B. On exit, the solution matrix X (pointer to double array)
   * @param ldb - Leading dimension of B (pointer to integer)
   * @param info - 0 = success, <0 = illegal argument, >0 = U(i,i) is zero, factorization complete but singular (pointer to integer)
   */
  _dgesv_(
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    ipiv: number,
    b: number,
    ldb: number,
    info: number
  ): void;

  /**
   * DGECON - Estimate the reciprocal condition number of a general matrix
   *
   * Estimates the reciprocal of the condition number of a general real matrix A,
   * in either the 1-norm or the infinity-norm, using the LU factorization computed by DGETRF.
   * The condition number is defined as κ(A) = ||A|| * ||A^(-1)||.
   *
   * @param norm - '1' or 'O' for 1-norm, 'I' for infinity-norm (pointer to char as integer)
   * @param n - Order of the matrix A (pointer to integer)
   * @param a - LU factorization from DGETRF (pointer to double array)
   * @param lda - Leading dimension of A (pointer to integer)
   * @param anorm - The norm of the original matrix A (1-norm or inf-norm matching norm parameter) (pointer to double)
   * @param rcond - Output: reciprocal condition number. If rcond is small, the matrix is nearly singular (pointer to double)
   * @param work - Workspace of dimension 4*N (pointer to double array)
   * @param iwork - Integer workspace of dimension N (pointer to integer array)
   * @param info - 0 = success, <0 = illegal argument (pointer to integer)
   */
  _dgecon_(
    norm: number,
    n: number,
    a: number,
    lda: number,
    anorm: number,
    rcond: number,
    work: number,
    iwork: number,
    info: number
  ): void;

  // ============================================================
  // Cholesky Factorization (Symmetric Positive Definite Matrices)
  // ============================================================

  /**
   * DPOTRF - Cholesky factorization of a symmetric positive definite matrix
   *
   * Computes the Cholesky factorization of a real symmetric positive definite matrix A.
   * The factorization has the form: A = U^T * U (if UPLO = 'U') or A = L * L^T (if UPLO = 'L')
   * where U is upper triangular and L is lower triangular.
   *
   * @param uplo - 'U' = upper triangle stored, 'L' = lower triangle stored (pointer to char as integer, 'U'=85, 'L'=76)
   * @param n - Order of the matrix A (pointer to integer)
   * @param a - On entry, symmetric matrix. On exit, the Cholesky factor U or L (pointer to double array)
   * @param lda - Leading dimension of A (pointer to integer)
   * @param info - 0 = success, <0 = illegal argument, >0 = matrix not positive definite (pointer to integer)
   */
  _dpotrf_(uplo: number, n: number, a: number, lda: number, info: number): void;

  /**
   * DPOTRS - Solve A*X = B using Cholesky factorization from DPOTRF
   *
   * Solves a system of linear equations A*X = B with a symmetric positive definite matrix A
   * using the Cholesky factorization A = U^T*U or A = L*L^T computed by DPOTRF.
   *
   * @param uplo - 'U' if upper triangle factor, 'L' if lower triangle factor (pointer to char as integer)
   * @param n - Order of the matrix A (pointer to integer)
   * @param nrhs - Number of right-hand sides (columns of B) (pointer to integer)
   * @param a - Cholesky factorization from DPOTRF (pointer to double array)
   * @param lda - Leading dimension of A (pointer to integer)
   * @param b - On entry, right-hand side matrix B. On exit, the solution matrix X (pointer to double array)
   * @param ldb - Leading dimension of B (pointer to integer)
   * @param info - 0 = success, <0 = illegal argument (pointer to integer)
   */
  _dpotrs_(
    uplo: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    info: number
  ): void;

  /**
   * DPOTRI - Compute inverse using Cholesky factorization from DPOTRF
   *
   * Computes the inverse of a real symmetric positive definite matrix A
   * using the Cholesky factorization A = U^T*U or A = L*L^T computed by DPOTRF.
   *
   * @param uplo - 'U' if upper triangle factor, 'L' if lower triangle factor (pointer to char as integer)
   * @param n - Order of the matrix A (pointer to integer)
   * @param a - On entry, Cholesky factor. On exit, the inverse of A (pointer to double array)
   * @param lda - Leading dimension of A (pointer to integer)
   * @param info - 0 = success, <0 = illegal argument, >0 = element (i,i) of factor is zero (pointer to integer)
   */
  _dpotri_(uplo: number, n: number, a: number, lda: number, info: number): void;

  /**
   * DPOSV - Solve A*X = B for positive definite matrices (driver routine)
   *
   * Computes the solution to a real system of linear equations A*X = B,
   * where A is an N-by-N symmetric positive definite matrix and X and B are N-by-NRHS matrices.
   * Uses Cholesky factorization. This is a driver routine that combines DPOTRF and DPOTRS.
   *
   * @param uplo - 'U' = upper triangle stored, 'L' = lower triangle stored (pointer to char as integer)
   * @param n - Number of linear equations (order of A) (pointer to integer)
   * @param nrhs - Number of right-hand sides (columns of B) (pointer to integer)
   * @param a - On entry, symmetric positive definite matrix. On exit, the Cholesky factor (pointer to double array)
   * @param lda - Leading dimension of A (pointer to integer)
   * @param b - On entry, right-hand side matrix B. On exit, the solution matrix X (pointer to double array)
   * @param ldb - Leading dimension of B (pointer to integer)
   * @param info - 0 = success, <0 = illegal argument, >0 = not positive definite (pointer to integer)
   */
  _dposv_(
    uplo: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    info: number
  ): void;

  /**
   * DPOCON - Estimate reciprocal condition number of a positive definite matrix
   *
   * Estimates the reciprocal of the condition number of a symmetric positive definite matrix
   * using the Cholesky factorization A = U^T*U or A = L*L^T computed by DPOTRF.
   *
   * @param uplo - 'U' if upper triangle factor, 'L' if lower triangle factor (pointer to char as integer)
   * @param n - Order of the matrix A (pointer to integer)
   * @param a - Cholesky factorization from DPOTRF (pointer to double array)
   * @param lda - Leading dimension of A (pointer to integer)
   * @param anorm - The 1-norm or infinity-norm of the original matrix A (pointer to double)
   * @param rcond - Output: reciprocal condition number (pointer to double)
   * @param work - Workspace of dimension 3*N (pointer to double array)
   * @param iwork - Integer workspace of dimension N (pointer to integer array)
   * @param info - 0 = success, <0 = illegal argument (pointer to integer)
   */
  _dpocon_(
    uplo: number,
    n: number,
    a: number,
    lda: number,
    anorm: number,
    rcond: number,
    work: number,
    iwork: number,
    info: number
  ): void;

  // ============================================================
  // Double Precision Eigenvalues
  // ============================================================

  /**
   * DGEEV - Compute eigenvalues and optionally eigenvectors of a general matrix
   *
   * @param jobvl - 'N' = don't compute left eigenvectors, 'V' = compute (pointer to char)
   * @param jobvr - 'N' = don't compute right eigenvectors, 'V' = compute (pointer to char)
   * @param n - Order of the matrix (pointer to integer)
   * @param a - Input matrix, destroyed on output (pointer to double array)
   * @param lda - Leading dimension of A (pointer to integer)
   * @param wr - Real parts of eigenvalues (pointer to double array)
   * @param wi - Imaginary parts of eigenvalues (pointer to double array)
   * @param vl - Left eigenvectors if jobvl='V' (pointer to double array)
   * @param ldvl - Leading dimension of VL (pointer to integer)
   * @param vr - Right eigenvectors if jobvr='V' (pointer to double array)
   * @param ldvr - Leading dimension of VR (pointer to integer)
   * @param work - Workspace (pointer to double array)
   * @param lwork - Size of workspace, -1 for query (pointer to integer)
   * @param info - 0 = success (pointer to integer)
   */
  _dgeev_(
    jobvl: number,
    jobvr: number,
    n: number,
    a: number,
    lda: number,
    wr: number,
    wi: number,
    vl: number,
    ldvl: number,
    vr: number,
    ldvr: number,
    work: number,
    lwork: number,
    info: number
  ): void;

  /**
   * DSYEV - Compute eigenvalues and optionally eigenvectors of a symmetric matrix
   *
   * @param jobz - 'N' = eigenvalues only, 'V' = eigenvalues and eigenvectors
   * @param uplo - 'U' = upper triangle, 'L' = lower triangle
   */
  _dsyev_(
    jobz: number,
    uplo: number,
    n: number,
    a: number,
    lda: number,
    w: number,
    work: number,
    lwork: number,
    info: number
  ): void;

  /**
   * DSYEVD - Compute eigenvalues of symmetric matrix using divide and conquer
   *
   * Computes all eigenvalues and, optionally, eigenvectors of a real symmetric matrix A
   * using a divide and conquer algorithm. This is generally faster than DSYEV for large matrices.
   *
   * @param jobz - 'N' = eigenvalues only, 'V' = eigenvalues and eigenvectors (pointer to char as integer)
   * @param uplo - 'U' = upper triangle stored, 'L' = lower triangle stored (pointer to char as integer)
   * @param n - Order of the matrix A (pointer to integer)
   * @param a - On entry, symmetric matrix. On exit, if JOBZ='V', the orthonormal eigenvectors (pointer to double array)
   * @param lda - Leading dimension of A (pointer to integer)
   * @param w - Output: eigenvalues in ascending order (pointer to double array of size N)
   * @param work - Workspace. On exit, WORK(1) contains optimal LWORK (pointer to double array)
   * @param lwork - Size of WORK. If LWORK=-1, workspace query is performed (pointer to integer)
   * @param iwork - Integer workspace (pointer to integer array)
   * @param liwork - Size of IWORK. If LIWORK=-1, workspace query is performed (pointer to integer)
   * @param info - 0 = success, <0 = illegal argument, >0 = algorithm failed to converge (pointer to integer)
   */
  _dsyevd_(
    jobz: number,
    uplo: number,
    n: number,
    a: number,
    lda: number,
    w: number,
    work: number,
    lwork: number,
    iwork: number,
    liwork: number,
    info: number
  ): void;

  /**
   * DSYEVR - Compute selected eigenvalues and eigenvectors using RRR algorithm (fastest)
   *
   * Computes selected eigenvalues and, optionally, eigenvectors of a real symmetric matrix A
   * using the Relatively Robust Representations (RRR) algorithm. This is the fastest method
   * and can compute a subset of eigenvalues/eigenvectors.
   *
   * @param jobz - 'N' = eigenvalues only, 'V' = eigenvalues and eigenvectors (pointer to char as integer)
   * @param range - 'A' = all, 'V' = in interval (VL,VU], 'I' = indices IL to IU (pointer to char as integer)
   * @param uplo - 'U' = upper triangle stored, 'L' = lower triangle stored (pointer to char as integer)
   * @param n - Order of the matrix A (pointer to integer)
   * @param a - On entry, symmetric matrix. On exit, destroyed (pointer to double array)
   * @param lda - Leading dimension of A (pointer to integer)
   * @param vl - Lower bound of interval if RANGE='V' (pointer to double)
   * @param vu - Upper bound of interval if RANGE='V' (pointer to double)
   * @param il - Index of smallest eigenvalue if RANGE='I' (pointer to integer)
   * @param iu - Index of largest eigenvalue if RANGE='I' (pointer to integer)
   * @param abstol - Absolute error tolerance for eigenvalues. Use DLAMCH('S') for best accuracy (pointer to double)
   * @param m - Output: total number of eigenvalues found (pointer to integer)
   * @param w - Output: eigenvalues in ascending order (pointer to double array)
   * @param z - Output: eigenvectors if JOBZ='V' (pointer to double array)
   * @param ldz - Leading dimension of Z (pointer to integer)
   * @param isuppz - Support of eigenvectors in Z (pointer to integer array of size 2*max(1,M))
   * @param work - Workspace (pointer to double array)
   * @param lwork - Size of WORK. If LWORK=-1, workspace query (pointer to integer)
   * @param iwork - Integer workspace (pointer to integer array)
   * @param liwork - Size of IWORK. If LIWORK=-1, workspace query (pointer to integer)
   * @param info - 0 = success (pointer to integer)
   */
  _dsyevr_(
    jobz: number,
    range: number,
    uplo: number,
    n: number,
    a: number,
    lda: number,
    vl: number,
    vu: number,
    il: number,
    iu: number,
    abstol: number,
    m: number,
    w: number,
    z: number,
    ldz: number,
    isuppz: number,
    work: number,
    lwork: number,
    iwork: number,
    liwork: number,
    info: number
  ): void;

  /**
   * DGEEVX - Expert driver for eigenvalues with balancing and condition estimation
   *
   * Computes eigenvalues and optionally eigenvectors of a general real matrix A.
   * Additionally performs balancing to improve accuracy and computes reciprocal
   * condition numbers for eigenvalues and eigenvectors.
   *
   * @param balanc - 'N' = no balancing, 'P' = permute, 'S' = scale, 'B' = both (pointer to char as integer)
   * @param jobvl - 'N' = don't compute left eigenvectors, 'V' = compute (pointer to char as integer)
   * @param jobvr - 'N' = don't compute right eigenvectors, 'V' = compute (pointer to char as integer)
   * @param sense - 'N' = no condition numbers, 'E' = for eigenvalues, 'V' = for eigenvectors, 'B' = both (pointer to char as integer)
   * @param n - Order of the matrix A (pointer to integer)
   * @param a - Input matrix, destroyed on output (pointer to double array)
   * @param lda - Leading dimension of A (pointer to integer)
   * @param wr - Real parts of eigenvalues (pointer to double array)
   * @param wi - Imaginary parts of eigenvalues (pointer to double array)
   * @param vl - Left eigenvectors if JOBVL='V' (pointer to double array)
   * @param ldvl - Leading dimension of VL (pointer to integer)
   * @param vr - Right eigenvectors if JOBVR='V' (pointer to double array)
   * @param ldvr - Leading dimension of VR (pointer to integer)
   * @param ilo - Output: index where balancing began (pointer to integer)
   * @param ihi - Output: index where balancing ended (pointer to integer)
   * @param scale - Output: permutations and scaling factors from balancing (pointer to double array)
   * @param abnrm - Output: 1-norm of balanced matrix (pointer to double)
   * @param rconde - Output: reciprocal condition numbers for eigenvalues (pointer to double array)
   * @param rcondv - Output: reciprocal condition numbers for eigenvectors (pointer to double array)
   * @param work - Workspace (pointer to double array)
   * @param lwork - Size of workspace (pointer to integer)
   * @param iwork - Integer workspace (pointer to integer array)
   * @param info - 0 = success (pointer to integer)
   */
  _dgeevx_(
    balanc: number,
    jobvl: number,
    jobvr: number,
    sense: number,
    n: number,
    a: number,
    lda: number,
    wr: number,
    wi: number,
    vl: number,
    ldvl: number,
    vr: number,
    ldvr: number,
    ilo: number,
    ihi: number,
    scale: number,
    abnrm: number,
    rconde: number,
    rcondv: number,
    work: number,
    lwork: number,
    iwork: number,
    info: number
  ): void;

  // ============================================================
  // Double Precision SVD
  // ============================================================

  /**
   * DGESVD - Compute singular value decomposition A = U * S * V^T
   *
   * @param jobu - 'A' = all U, 'S' = first min(m,n) columns, 'O' = overwrite A, 'N' = none
   * @param jobvt - 'A' = all V^T, 'S' = first min(m,n) rows, 'O' = overwrite A, 'N' = none
   * @param m - Number of rows of A (pointer to integer)
   * @param n - Number of columns of A (pointer to integer)
   * @param a - Input matrix, destroyed on output (pointer to double array)
   * @param lda - Leading dimension of A (pointer to integer)
   * @param s - Singular values in descending order (pointer to double array)
   * @param u - Left singular vectors (pointer to double array)
   * @param ldu - Leading dimension of U (pointer to integer)
   * @param vt - Right singular vectors (V^T) (pointer to double array)
   * @param ldvt - Leading dimension of VT (pointer to integer)
   * @param work - Workspace (pointer to double array)
   * @param lwork - Size of workspace, -1 for query (pointer to integer)
   * @param info - 0 = success (pointer to integer)
   */
  _dgesvd_(
    jobu: number,
    jobvt: number,
    m: number,
    n: number,
    a: number,
    lda: number,
    s: number,
    u: number,
    ldu: number,
    vt: number,
    ldvt: number,
    work: number,
    lwork: number,
    info: number
  ): void;

  /**
   * DGESDD - SVD using divide and conquer algorithm (faster for large matrices)
   *
   * Computes the singular value decomposition (SVD) of a real M-by-N matrix A
   * using a divide and conquer algorithm. This is generally faster than DGESVD
   * for large matrices but uses more workspace.
   *
   * A = U * SIGMA * V^T
   *
   * @param jobz - Specifies options for computing U and V^T:
   *               'A' = all M columns of U and all N rows of V^T
   *               'S' = first min(M,N) columns of U and rows of V^T (the singular vectors)
   *               'O' = if M >= N, first N columns of U overwritten to A and all rows of V^T computed;
   *                     if M < N, all columns of U computed and first M rows of V^T overwritten to A
   *               'N' = no columns of U or rows of V^T computed (pointer to char as integer)
   * @param m - Number of rows of A (pointer to integer)
   * @param n - Number of columns of A (pointer to integer)
   * @param a - Input matrix, destroyed on output (or overwritten with U/V^T if JOBZ='O') (pointer to double array)
   * @param lda - Leading dimension of A (pointer to integer)
   * @param s - Output: singular values in descending order (pointer to double array of size min(M,N))
   * @param u - Output: left singular vectors (pointer to double array)
   * @param ldu - Leading dimension of U (pointer to integer)
   * @param vt - Output: right singular vectors (V^T) (pointer to double array)
   * @param ldvt - Leading dimension of VT (pointer to integer)
   * @param work - Workspace (pointer to double array)
   * @param lwork - Size of workspace, -1 for query (pointer to integer)
   * @param iwork - Integer workspace of dimension 8*min(M,N) (pointer to integer array)
   * @param info - 0 = success, <0 = illegal argument, >0 = DBDSDC did not converge (pointer to integer)
   */
  _dgesdd_(
    jobz: number,
    m: number,
    n: number,
    a: number,
    lda: number,
    s: number,
    u: number,
    ldu: number,
    vt: number,
    ldvt: number,
    work: number,
    lwork: number,
    iwork: number,
    info: number
  ): void;

  // ============================================================
  // Double Precision QR Factorization
  // ============================================================

  /**
   * DGEQRF - QR factorization of an M-by-N matrix
   *
   * Computes a QR factorization of a real M-by-N matrix A: A = Q * R
   * where Q is an M-by-M orthogonal matrix and R is an M-by-N upper triangular matrix.
   * The matrix Q is represented as a product of elementary Householder reflectors
   * stored in the lower triangle of A and in the TAU array.
   *
   * @param m - Number of rows of A (pointer to integer)
   * @param n - Number of columns of A (pointer to integer)
   * @param a - On entry, M-by-N matrix. On exit, R in upper triangle; Householder vectors below diagonal (pointer to double array)
   * @param lda - Leading dimension of A (pointer to integer)
   * @param tau - Output: scalar factors of Householder reflectors (pointer to double array of size min(M,N))
   * @param work - Workspace (pointer to double array)
   * @param lwork - Size of workspace, -1 for query. Optimal is N*NB (pointer to integer)
   * @param info - 0 = success, <0 = illegal argument (pointer to integer)
   */
  _dgeqrf_(
    m: number,
    n: number,
    a: number,
    lda: number,
    tau: number,
    work: number,
    lwork: number,
    info: number
  ): void;

  /**
   * DGEQP3 - QR factorization with column pivoting
   *
   * Computes a QR factorization with column pivoting of an M-by-N matrix A: A * P = Q * R
   * where P is a permutation matrix, Q is orthogonal, and R is upper triangular.
   * Column pivoting produces a rank-revealing factorization.
   *
   * @param m - Number of rows of A (pointer to integer)
   * @param n - Number of columns of A (pointer to integer)
   * @param a - On entry, M-by-N matrix. On exit, R in upper triangle; Householder vectors below diagonal (pointer to double array)
   * @param lda - Leading dimension of A (pointer to integer)
   * @param jpvt - On entry, if JPVT(j)!=0, column j is permuted to front. On exit, permutation order (pointer to integer array)
   * @param tau - Output: scalar factors of Householder reflectors (pointer to double array)
   * @param work - Workspace (pointer to double array)
   * @param lwork - Size of workspace, -1 for query (pointer to integer)
   * @param info - 0 = success, <0 = illegal argument (pointer to integer)
   */
  _dgeqp3_(
    m: number,
    n: number,
    a: number,
    lda: number,
    jpvt: number,
    tau: number,
    work: number,
    lwork: number,
    info: number
  ): void;

  /**
   * DORGQR - Generate orthogonal matrix Q from QR factorization
   *
   * Generates an M-by-N real matrix Q with orthonormal columns from the
   * implicit representation stored by DGEQRF. The first K columns of Q
   * are the result of the K Householder reflectors.
   *
   * @param m - Number of rows of Q to generate (pointer to integer)
   * @param n - Number of columns of Q to generate, N <= M (pointer to integer)
   * @param k - Number of Householder reflectors, K <= N (pointer to integer)
   * @param a - On entry, Householder vectors from DGEQRF. On exit, the M-by-N matrix Q (pointer to double array)
   * @param lda - Leading dimension of A (pointer to integer)
   * @param tau - Scalar factors of Householder reflectors from DGEQRF (pointer to double array)
   * @param work - Workspace (pointer to double array)
   * @param lwork - Size of workspace, -1 for query (pointer to integer)
   * @param info - 0 = success, <0 = illegal argument (pointer to integer)
   */
  _dorgqr_(
    m: number,
    n: number,
    k: number,
    a: number,
    lda: number,
    tau: number,
    work: number,
    lwork: number,
    info: number
  ): void;

  /**
   * DORMQR - Multiply matrix by orthogonal Q from QR factorization
   *
   * Overwrites the general real M-by-N matrix C with:
   * - SIDE='L': Q*C or Q^T*C
   * - SIDE='R': C*Q or C*Q^T
   * where Q is the orthogonal matrix from DGEQRF. This is more efficient
   * than explicitly forming Q and then multiplying.
   *
   * @param side - 'L' = apply Q or Q^T from left, 'R' = from right (pointer to char as integer)
   * @param trans - 'N' = apply Q, 'T' = apply Q^T (pointer to char as integer)
   * @param m - Number of rows of C (pointer to integer)
   * @param n - Number of columns of C (pointer to integer)
   * @param k - Number of Householder reflectors (pointer to integer)
   * @param a - Householder vectors from DGEQRF (pointer to double array)
   * @param lda - Leading dimension of A (pointer to integer)
   * @param tau - Scalar factors of Householder reflectors from DGEQRF (pointer to double array)
   * @param c - On entry, M-by-N matrix C. On exit, overwritten with result (pointer to double array)
   * @param ldc - Leading dimension of C (pointer to integer)
   * @param work - Workspace (pointer to double array)
   * @param lwork - Size of workspace, -1 for query (pointer to integer)
   * @param info - 0 = success, <0 = illegal argument (pointer to integer)
   */
  _dormqr_(
    side: number,
    trans: number,
    m: number,
    n: number,
    k: number,
    a: number,
    lda: number,
    tau: number,
    c: number,
    ldc: number,
    work: number,
    lwork: number,
    info: number
  ): void;

  // ============================================================
  // Double Precision Least Squares
  // ============================================================

  /**
   * DGELS - Solve overdetermined/underdetermined linear systems using QR/LQ factorization
   *
   * Solves overdetermined or underdetermined real linear systems involving an M-by-N matrix A:
   * - If M >= N and TRANS='N': minimize ||B - A*X||_2 (least squares)
   * - If M < N and TRANS='N': find minimum norm solution to underdetermined system A*X = B
   * - TRANS='T' solves the transposed problems
   * Assumes A has full rank. Uses QR factorization for M >= N, LQ for M < N.
   *
   * @param trans - 'N' = no transpose, 'T' = transpose (pointer to char as integer)
   * @param m - Number of rows of A (pointer to integer)
   * @param n - Number of columns of A (pointer to integer)
   * @param nrhs - Number of right-hand sides (columns of B) (pointer to integer)
   * @param a - On entry, M-by-N matrix A. On exit, QR or LQ factorization (pointer to double array)
   * @param lda - Leading dimension of A (pointer to integer)
   * @param b - On entry, M-by-NRHS (or N-by-NRHS) matrix B. On exit, solution X (pointer to double array)
   * @param ldb - Leading dimension of B (pointer to integer)
   * @param work - Workspace (pointer to double array)
   * @param lwork - Size of workspace, -1 for query (pointer to integer)
   * @param info - 0 = success, <0 = illegal argument (pointer to integer)
   */
  _dgels_(
    trans: number,
    m: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    work: number,
    lwork: number,
    info: number
  ): void;

  /**
   * DGELSD - Solve least squares using SVD with divide and conquer
   *
   * Computes the minimum-norm solution to a real linear least squares problem:
   * minimize ||B - A*X||_2
   * using the singular value decomposition (SVD) of A with divide-and-conquer algorithm.
   * A is an M-by-N matrix which may be rank-deficient.
   *
   * @param m - Number of rows of A (pointer to integer)
   * @param n - Number of columns of A (pointer to integer)
   * @param nrhs - Number of right-hand sides (pointer to integer)
   * @param a - M-by-N matrix A, destroyed on output (pointer to double array)
   * @param lda - Leading dimension of A (pointer to integer)
   * @param b - On entry, M-by-NRHS matrix B. On exit, N-by-NRHS solution X (pointer to double array)
   * @param ldb - Leading dimension of B (pointer to integer)
   * @param s - Output: singular values of A in descending order (pointer to double array)
   * @param rcond - Used to determine effective rank. Singular values S(i) <= RCOND*S(1) treated as zero.
   *                If RCOND < 0, machine precision is used (pointer to double)
   * @param rank - Output: effective rank of A (number of singular values > RCOND*S(1)) (pointer to integer)
   * @param work - Workspace (pointer to double array)
   * @param lwork - Size of workspace, -1 for query (pointer to integer)
   * @param iwork - Integer workspace (pointer to integer array)
   * @param info - 0 = success, <0 = illegal argument, >0 = SVD failed to converge (pointer to integer)
   */
  _dgelsd_(
    m: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    s: number,
    rcond: number,
    rank: number,
    work: number,
    lwork: number,
    iwork: number,
    info: number
  ): void;

  /**
   * DGELSS - Solve least squares using SVD
   *
   * Computes the minimum-norm solution to a real linear least squares problem:
   * minimize ||B - A*X||_2
   * using the singular value decomposition (SVD) of A. This is the standard
   * SVD-based solver (DGELSD is usually faster for large problems).
   *
   * @param m - Number of rows of A (pointer to integer)
   * @param n - Number of columns of A (pointer to integer)
   * @param nrhs - Number of right-hand sides (pointer to integer)
   * @param a - M-by-N matrix A, destroyed on output (pointer to double array)
   * @param lda - Leading dimension of A (pointer to integer)
   * @param b - On entry, M-by-NRHS matrix B. On exit, N-by-NRHS solution X (pointer to double array)
   * @param ldb - Leading dimension of B (pointer to integer)
   * @param s - Output: singular values of A in descending order (pointer to double array)
   * @param rcond - Used to determine effective rank. If RCOND < 0, machine precision is used (pointer to double)
   * @param rank - Output: effective rank of A (pointer to integer)
   * @param work - Workspace (pointer to double array)
   * @param lwork - Size of workspace, -1 for query (pointer to integer)
   * @param info - 0 = success, <0 = illegal argument, >0 = SVD failed to converge (pointer to integer)
   */
  _dgelss_(
    m: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    s: number,
    rcond: number,
    rank: number,
    work: number,
    lwork: number,
    info: number
  ): void;

  /**
   * DGELSY - Solve least squares using complete orthogonal factorization
   *
   * Computes the minimum-norm solution to a real linear least squares problem:
   * minimize ||B - A*X||_2
   * using a complete orthogonal factorization of A with column pivoting.
   * This method is generally faster than SVD-based methods for rank-deficient problems.
   *
   * @param m - Number of rows of A (pointer to integer)
   * @param n - Number of columns of A (pointer to integer)
   * @param nrhs - Number of right-hand sides (pointer to integer)
   * @param a - M-by-N matrix A, destroyed on output (pointer to double array)
   * @param lda - Leading dimension of A (pointer to integer)
   * @param b - On entry, M-by-NRHS matrix B. On exit, N-by-NRHS solution X (pointer to double array)
   * @param ldb - Leading dimension of B (pointer to integer)
   * @param jpvt - On entry, if JPVT(i)!=0, column i is initial column of permutation.
   *               On exit, if JPVT(i)=k, column k was moved to position i (pointer to integer array)
   * @param rcond - Used to determine effective rank. If RCOND < 0, machine precision is used (pointer to double)
   * @param rank - Output: effective rank of A (pointer to integer)
   * @param work - Workspace (pointer to double array)
   * @param lwork - Size of workspace, -1 for query (pointer to integer)
   * @param info - 0 = success, <0 = illegal argument (pointer to integer)
   */
  _dgelsy_(
    m: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    jpvt: number,
    rcond: number,
    rank: number,
    work: number,
    lwork: number,
    info: number
  ): void;

  // ============================================================
  // Double Precision Triangular
  // ============================================================

  /**
   * DTRTRS - Solve triangular system A*X = B or A^T*X = B
   *
   * Solves a triangular system of the form A*X = B, A^T*X = B, or A^H*X = B,
   * where A is a triangular matrix of order N, and B is an N-by-NRHS matrix.
   * A check is made to verify that A is nonsingular.
   *
   * @param uplo - 'U' = A is upper triangular, 'L' = A is lower triangular (pointer to char as integer)
   * @param trans - 'N' = solve A*X=B, 'T' = solve A^T*X=B, 'C' = solve A^H*X=B (pointer to char as integer)
   * @param diag - 'N' = A is not unit triangular, 'U' = A is unit triangular (pointer to char as integer)
   * @param n - Order of the matrix A (pointer to integer)
   * @param nrhs - Number of right-hand sides (columns of B) (pointer to integer)
   * @param a - The triangular matrix A (pointer to double array)
   * @param lda - Leading dimension of A (pointer to integer)
   * @param b - On entry, right-hand side matrix B. On exit, solution matrix X (pointer to double array)
   * @param ldb - Leading dimension of B (pointer to integer)
   * @param info - 0 = success, <0 = illegal argument, >0 = A is singular (pointer to integer)
   */
  _dtrtrs_(
    uplo: number,
    trans: number,
    diag: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    info: number
  ): void;

  /**
   * DTRTRI - Compute inverse of a triangular matrix
   *
   * Computes the inverse of a real upper or lower triangular matrix A.
   *
   * @param uplo - 'U' = A is upper triangular, 'L' = A is lower triangular (pointer to char as integer)
   * @param diag - 'N' = A is not unit triangular, 'U' = A is unit triangular (pointer to char as integer)
   * @param n - Order of the matrix A (pointer to integer)
   * @param a - On entry, triangular matrix A. On exit, the inverse of A (pointer to double array)
   * @param lda - Leading dimension of A (pointer to integer)
   * @param info - 0 = success, <0 = illegal argument, >0 = A is singular (pointer to integer)
   */
  _dtrtri_(uplo: number, diag: number, n: number, a: number, lda: number, info: number): void;

  // ============================================================
  // Double Precision Utilities
  // ============================================================

  /**
   * DLAMCH - Determine double precision machine parameters
   *
   * Returns machine-dependent parameters for double precision floating-point arithmetic.
   * Essential for writing portable numerical software.
   *
   * @param cmach - Character specifying the parameter to return (pointer to char as integer):
   *                'E' or 'e' = relative machine precision (epsilon), eps such that 1+eps > 1
   *                'S' or 's' = safe minimum, smallest number such that 1/sfmin does not overflow
   *                'B' or 'b' = base of the machine (typically 2)
   *                'P' or 'p' = eps * base
   *                'N' or 'n' = number of base digits in the mantissa
   *                'R' or 'r' = 1 if rounding occurs in addition, 0 otherwise
   *                'M' or 'm' = minimum exponent before underflow
   *                'U' or 'u' = underflow threshold = base^(emin-1)
   *                'L' or 'l' = maximum exponent before overflow
   *                'O' or 'o' = overflow threshold = base^emax * (1-eps)
   * @returns The requested machine parameter value
   */
  _dlamch_(cmach: number): number;

  /**
   * DLANGE - Compute the norm of a general matrix
   *
   * Returns the value of the specified norm of a general real M-by-N matrix A.
   *
   * @param norm - Specifies which norm to compute (pointer to char as integer):
   *               'M' or 'm' = max(abs(A(i,j))) - max absolute element
   *               '1', 'O' or 'o' = 1-norm (maximum column sum)
   *               'I' or 'i' = infinity-norm (maximum row sum)
   *               'F', 'f', 'E' or 'e' = Frobenius norm = sqrt(sum(A(i,j)^2))
   * @param m - Number of rows (pointer to integer)
   * @param n - Number of columns (pointer to integer)
   * @param a - The M-by-N matrix (pointer to double array)
   * @param lda - Leading dimension of A (pointer to integer)
   * @param work - Workspace of dimension M (only needed for infinity-norm) (pointer to double array)
   * @returns The computed norm value
   */
  _dlange_(norm: number, m: number, n: number, a: number, lda: number, work: number): number;

  /**
   * DLANSY - Compute the norm of a symmetric matrix
   *
   * Returns the value of the specified norm of a real symmetric matrix A.
   *
   * @param norm - Specifies which norm to compute (pointer to char as integer):
   *               'M' = max absolute element, '1'/'O' = 1-norm = infinity-norm (symmetric),
   *               'I' = infinity-norm = 1-norm, 'F'/'E' = Frobenius norm
   * @param uplo - 'U' = upper triangle stored, 'L' = lower triangle stored (pointer to char as integer)
   * @param n - Order of the matrix (pointer to integer)
   * @param a - The symmetric matrix (pointer to double array)
   * @param lda - Leading dimension of A (pointer to integer)
   * @param work - Workspace of dimension N (only needed for 1-norm/infinity-norm) (pointer to double array)
   * @returns The computed norm value
   */
  _dlansy_(norm: number, uplo: number, n: number, a: number, lda: number, work: number): number;

  /**
   * DLASWP - Apply row permutations to a matrix
   *
   * Performs a series of row interchanges on a general rectangular matrix.
   * Typically used to apply pivoting from LU factorization.
   *
   * @param n - Number of columns of the matrix (pointer to integer)
   * @param a - The matrix to permute (pointer to double array)
   * @param lda - Leading dimension of A (pointer to integer)
   * @param k1 - First element of IPIV to use (pointer to integer)
   * @param k2 - Last element of IPIV to use (pointer to integer)
   * @param ipiv - Pivot indices. Row i is interchanged with row IPIV(i) (pointer to integer array)
   * @param incx - Increment between successive values of IPIV.
   *               If INCX > 0, pivots applied in order K1..K2.
   *               If INCX < 0, pivots applied in reverse order K2..K1 (pointer to integer)
   */
  _dlaswp_(n: number, a: number, lda: number, k1: number, k2: number, ipiv: number, incx: number): void;

  /**
   * DLASSQ - Update a sum of squares represented in scaled form
   *
   * Returns values scale and sumsq such that:
   *   (scale^2) * sumsq = x(1)^2 + ... + x(n)^2 + (scale_in^2) * sumsq_in
   *
   * This representation prevents overflow when computing the sum of squares of large values.
   * The final 2-norm can be computed as scale * sqrt(sumsq).
   *
   * @param n - Number of elements in vector X (pointer to integer)
   * @param x - Input vector (pointer to double array)
   * @param incx - Increment between elements of X (pointer to integer)
   * @param scale - On entry, initial scale. On exit, updated scale (pointer to double)
   * @param sumsq - On entry, initial sum of squares. On exit, updated sum of squares (pointer to double)
   */
  _dlassq_(n: number, x: number, incx: number, scale: number, sumsq: number): void;

  // ============================================================
  // Single Precision (S prefix)
  // ============================================================
  // Single precision versions have identical signatures to double precision (D prefix).
  // Use HEAPF32 instead of HEAPF64 for array access.
  // Trade-off: Half the memory usage, faster on some hardware, but reduced precision (~7 decimal digits vs ~15).

  /** SGETRF - Single precision LU factorization. See DGETRF for details. */
  _sgetrf_(m: number, n: number, a: number, lda: number, ipiv: number, info: number): void;
  /** SGETRS - Single precision triangular solve using LU factors. See DGETRS for details. */
  _sgetrs_(
    trans: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    ipiv: number,
    b: number,
    ldb: number,
    info: number
  ): void;
  /** SGETRI - Single precision matrix inverse using LU factors. See DGETRI for details. */
  _sgetri_(
    n: number,
    a: number,
    lda: number,
    ipiv: number,
    work: number,
    lwork: number,
    info: number
  ): void;
  /** SGESV - Single precision driver for A*X=B. See DGESV for details. */
  _sgesv_(
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    ipiv: number,
    b: number,
    ldb: number,
    info: number
  ): void;
  /** SGECON - Single precision condition number estimate. See DGECON for details. */
  _sgecon_(
    norm: number,
    n: number,
    a: number,
    lda: number,
    anorm: number,
    rcond: number,
    work: number,
    iwork: number,
    info: number
  ): void;

  /** SPOTRF - Single precision Cholesky factorization. See DPOTRF for details. */
  _spotrf_(uplo: number, n: number, a: number, lda: number, info: number): void;
  /** SPOTRS - Single precision solve using Cholesky factors. See DPOTRS for details. */
  _spotrs_(
    uplo: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    info: number
  ): void;
  /** SPOTRI - Single precision inverse using Cholesky factors. See DPOTRI for details. */
  _spotri_(uplo: number, n: number, a: number, lda: number, info: number): void;
  /** SPOSV - Single precision driver for positive definite A*X=B. See DPOSV for details. */
  _sposv_(
    uplo: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    info: number
  ): void;
  /** SPOCON - Single precision condition number for positive definite matrix. See DPOCON for details. */
  _spocon_(
    uplo: number,
    n: number,
    a: number,
    lda: number,
    anorm: number,
    rcond: number,
    work: number,
    iwork: number,
    info: number
  ): void;

  /** SGEEV - Single precision eigenvalues/eigenvectors of general matrix. See DGEEV for details. */
  _sgeev_(
    jobvl: number,
    jobvr: number,
    n: number,
    a: number,
    lda: number,
    wr: number,
    wi: number,
    vl: number,
    ldvl: number,
    vr: number,
    ldvr: number,
    work: number,
    lwork: number,
    info: number
  ): void;
  /** SSYEV - Single precision symmetric eigenvalues. See DSYEV for details. */
  _ssyev_(
    jobz: number,
    uplo: number,
    n: number,
    a: number,
    lda: number,
    w: number,
    work: number,
    lwork: number,
    info: number
  ): void;
  /** SSYEVD - Single precision symmetric eigenvalues (divide and conquer). See DSYEVD for details. */
  _ssyevd_(
    jobz: number,
    uplo: number,
    n: number,
    a: number,
    lda: number,
    w: number,
    work: number,
    lwork: number,
    iwork: number,
    liwork: number,
    info: number
  ): void;
  /** SSYEVR - Single precision selected symmetric eigenvalues (RRR). See DSYEVR for details. */
  _ssyevr_(
    jobz: number,
    range: number,
    uplo: number,
    n: number,
    a: number,
    lda: number,
    vl: number,
    vu: number,
    il: number,
    iu: number,
    abstol: number,
    m: number,
    w: number,
    z: number,
    ldz: number,
    isuppz: number,
    work: number,
    lwork: number,
    iwork: number,
    liwork: number,
    info: number
  ): void;
  /** SGEEVX - Single precision expert eigenvalue driver. See DGEEVX for details. */
  _sgeevx_(
    balanc: number,
    jobvl: number,
    jobvr: number,
    sense: number,
    n: number,
    a: number,
    lda: number,
    wr: number,
    wi: number,
    vl: number,
    ldvl: number,
    vr: number,
    ldvr: number,
    ilo: number,
    ihi: number,
    scale: number,
    abnrm: number,
    rconde: number,
    rcondv: number,
    work: number,
    lwork: number,
    iwork: number,
    info: number
  ): void;

  /** SGESVD - Single precision SVD. See DGESVD for details. */
  _sgesvd_(
    jobu: number,
    jobvt: number,
    m: number,
    n: number,
    a: number,
    lda: number,
    s: number,
    u: number,
    ldu: number,
    vt: number,
    ldvt: number,
    work: number,
    lwork: number,
    info: number
  ): void;
  /** SGESDD - Single precision SVD (divide and conquer). See DGESDD for details. */
  _sgesdd_(
    jobz: number,
    m: number,
    n: number,
    a: number,
    lda: number,
    s: number,
    u: number,
    ldu: number,
    vt: number,
    ldvt: number,
    work: number,
    lwork: number,
    iwork: number,
    info: number
  ): void;

  /** SGEQRF - Single precision QR factorization. See DGEQRF for details. */
  _sgeqrf_(
    m: number,
    n: number,
    a: number,
    lda: number,
    tau: number,
    work: number,
    lwork: number,
    info: number
  ): void;
  /** SGEQP3 - Single precision QR with column pivoting. See DGEQP3 for details. */
  _sgeqp3_(
    m: number,
    n: number,
    a: number,
    lda: number,
    jpvt: number,
    tau: number,
    work: number,
    lwork: number,
    info: number
  ): void;
  /** SORGQR - Single precision generate Q from QR. See DORGQR for details. */
  _sorgqr_(
    m: number,
    n: number,
    k: number,
    a: number,
    lda: number,
    tau: number,
    work: number,
    lwork: number,
    info: number
  ): void;
  /** SORMQR - Single precision multiply by Q from QR. See DORMQR for details. */
  _sormqr_(
    side: number,
    trans: number,
    m: number,
    n: number,
    k: number,
    a: number,
    lda: number,
    tau: number,
    c: number,
    ldc: number,
    work: number,
    lwork: number,
    info: number
  ): void;

  /** SGELS - Single precision least squares via QR/LQ. See DGELS for details. */
  _sgels_(
    trans: number,
    m: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    work: number,
    lwork: number,
    info: number
  ): void;
  /** SGELSD - Single precision least squares via SVD (divide and conquer). See DGELSD for details. */
  _sgelsd_(
    m: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    s: number,
    rcond: number,
    rank: number,
    work: number,
    lwork: number,
    iwork: number,
    info: number
  ): void;
  /** SGELSS - Single precision least squares via SVD. See DGELSS for details. */
  _sgelss_(
    m: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    s: number,
    rcond: number,
    rank: number,
    work: number,
    lwork: number,
    info: number
  ): void;
  /** SGELSY - Single precision least squares via complete orthogonal factorization. See DGELSY for details. */
  _sgelsy_(
    m: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    jpvt: number,
    rcond: number,
    rank: number,
    work: number,
    lwork: number,
    info: number
  ): void;

  /** STRTRS - Single precision triangular solve. See DTRTRS for details. */
  _strtrs_(
    uplo: number,
    trans: number,
    diag: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    info: number
  ): void;
  /** STRTRI - Single precision triangular matrix inverse. See DTRTRI for details. */
  _strtri_(uplo: number, diag: number, n: number, a: number, lda: number, info: number): void;

  /** SLAMCH - Single precision machine parameters. See DLAMCH for details. */
  _slamch_(cmach: number): number;
  /** SLANGE - Single precision matrix norm. See DLANGE for details. */
  _slange_(norm: number, m: number, n: number, a: number, lda: number, work: number): number;
  /** SLANSY - Single precision symmetric matrix norm. See DLANSY for details. */
  _slansy_(norm: number, uplo: number, n: number, a: number, lda: number, work: number): number;
  /** SLASWP - Single precision row permutations. See DLASWP for details. */
  _slaswp_(n: number, a: number, lda: number, k1: number, k2: number, ipiv: number, incx: number): void;
  /** SLASSQ - Single precision scaled sum of squares. See DLASSQ for details. */
  _slassq_(n: number, x: number, incx: number, scale: number, sumsq: number): void;

  // ============================================================
  // Complex Double Precision (Z prefix)
  // ============================================================
  // Complex numbers are stored as interleaved real/imaginary pairs in HEAPF64.
  // Each complex number occupies 16 bytes (two 8-byte doubles).
  // Arrays use column-major order with complex elements.
  // RWORK arrays contain real workspace for operations that need real intermediate results.

  /** ZGETRF - Complex double precision LU factorization. See DGETRF for details. */
  _zgetrf_(m: number, n: number, a: number, lda: number, ipiv: number, info: number): void;
  /** ZGETRS - Complex double precision triangular solve using LU factors. See DGETRS for details. */
  _zgetrs_(
    trans: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    ipiv: number,
    b: number,
    ldb: number,
    info: number
  ): void;
  /** ZGETRI - Complex double precision matrix inverse using LU factors. See DGETRI for details. */
  _zgetri_(
    n: number,
    a: number,
    lda: number,
    ipiv: number,
    work: number,
    lwork: number,
    info: number
  ): void;
  /** ZGESV - Complex double precision driver for A*X=B. See DGESV for details. */
  _zgesv_(
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    ipiv: number,
    b: number,
    ldb: number,
    info: number
  ): void;
  /**
   * ZGECON - Complex double precision condition number estimate
   * @param rwork - Real workspace of dimension 2*N (NOT complex)
   */
  _zgecon_(
    norm: number,
    n: number,
    a: number,
    lda: number,
    anorm: number,
    rcond: number,
    work: number,
    rwork: number,
    info: number
  ): void;

  /** ZPOTRF - Complex Hermitian positive definite Cholesky factorization. See DPOTRF for details. */
  _zpotrf_(uplo: number, n: number, a: number, lda: number, info: number): void;
  /** ZPOTRS - Complex solve using Cholesky factors. See DPOTRS for details. */
  _zpotrs_(
    uplo: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    info: number
  ): void;
  /** ZPOTRI - Complex inverse using Cholesky factors. See DPOTRI for details. */
  _zpotri_(uplo: number, n: number, a: number, lda: number, info: number): void;
  /** ZPOSV - Complex driver for Hermitian positive definite A*X=B. See DPOSV for details. */
  _zposv_(
    uplo: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    info: number
  ): void;
  /** ZPOCON - Complex condition number for Hermitian positive definite matrix. See DPOCON for details. */
  _zpocon_(
    uplo: number,
    n: number,
    a: number,
    lda: number,
    anorm: number,
    rcond: number,
    work: number,
    rwork: number,
    info: number
  ): void;

  /**
   * ZGEEV - Complex double precision eigenvalues/eigenvectors of general matrix
   *
   * Unlike DGEEV, eigenvalues are returned in a single complex array W (not separate WR/WI).
   * @param w - Output: complex eigenvalues (pointer to complex double array of size N)
   * @param rwork - Real workspace of dimension 2*N
   */
  _zgeev_(
    jobvl: number,
    jobvr: number,
    n: number,
    a: number,
    lda: number,
    w: number,
    vl: number,
    ldvl: number,
    vr: number,
    ldvr: number,
    work: number,
    lwork: number,
    rwork: number,
    info: number
  ): void;
  /**
   * ZHEEV - Complex Hermitian matrix eigenvalues (real eigenvalues)
   *
   * Hermitian matrices have real eigenvalues, so W is a real array.
   * @param w - Output: real eigenvalues in ascending order (pointer to REAL double array)
   * @param rwork - Real workspace of dimension max(1, 3*N-2)
   */
  _zheev_(
    jobz: number,
    uplo: number,
    n: number,
    a: number,
    lda: number,
    w: number,
    work: number,
    lwork: number,
    rwork: number,
    info: number
  ): void;
  /** ZHEEVD - Complex Hermitian eigenvalues (divide and conquer). See DSYEVD for details. */
  _zheevd_(
    jobz: number,
    uplo: number,
    n: number,
    a: number,
    lda: number,
    w: number,
    work: number,
    lwork: number,
    rwork: number,
    lrwork: number,
    iwork: number,
    liwork: number,
    info: number
  ): void;
  /** ZHEEVR - Complex Hermitian selected eigenvalues (RRR). See DSYEVR for details. */
  _zheevr_(
    jobz: number,
    range: number,
    uplo: number,
    n: number,
    a: number,
    lda: number,
    vl: number,
    vu: number,
    il: number,
    iu: number,
    abstol: number,
    m: number,
    w: number,
    z: number,
    ldz: number,
    isuppz: number,
    work: number,
    lwork: number,
    rwork: number,
    lrwork: number,
    iwork: number,
    liwork: number,
    info: number
  ): void;
  /** ZGEEVX - Complex expert eigenvalue driver with balancing and conditioning. See DGEEVX for details. */
  _zgeevx_(
    balanc: number,
    jobvl: number,
    jobvr: number,
    sense: number,
    n: number,
    a: number,
    lda: number,
    w: number,
    vl: number,
    ldvl: number,
    vr: number,
    ldvr: number,
    ilo: number,
    ihi: number,
    scale: number,
    abnrm: number,
    rconde: number,
    rcondv: number,
    work: number,
    lwork: number,
    rwork: number,
    info: number
  ): void;

  /**
   * ZGESVD - Complex double precision SVD
   *
   * Singular values are always real and returned in a real array S.
   * @param s - Output: real singular values in descending order (pointer to REAL double array)
   * @param rwork - Real workspace
   */
  _zgesvd_(
    jobu: number,
    jobvt: number,
    m: number,
    n: number,
    a: number,
    lda: number,
    s: number,
    u: number,
    ldu: number,
    vt: number,
    ldvt: number,
    work: number,
    lwork: number,
    rwork: number,
    info: number
  ): void;
  /** ZGESDD - Complex SVD (divide and conquer). See DGESDD for details. */
  _zgesdd_(
    jobz: number,
    m: number,
    n: number,
    a: number,
    lda: number,
    s: number,
    u: number,
    ldu: number,
    vt: number,
    ldvt: number,
    work: number,
    lwork: number,
    rwork: number,
    iwork: number,
    info: number
  ): void;

  /** ZGEQRF - Complex QR factorization. See DGEQRF for details. */
  _zgeqrf_(
    m: number,
    n: number,
    a: number,
    lda: number,
    tau: number,
    work: number,
    lwork: number,
    info: number
  ): void;
  /** ZGEQP3 - Complex QR with column pivoting. See DGEQP3 for details. */
  _zgeqp3_(
    m: number,
    n: number,
    a: number,
    lda: number,
    jpvt: number,
    tau: number,
    work: number,
    lwork: number,
    rwork: number,
    info: number
  ): void;
  /** ZUNGQR - Complex generate unitary Q from QR. See DORGQR for details. */
  _zungqr_(
    m: number,
    n: number,
    k: number,
    a: number,
    lda: number,
    tau: number,
    work: number,
    lwork: number,
    info: number
  ): void;
  /**
   * ZUNMQR - Complex multiply by unitary Q from QR
   *
   * For complex matrices, TRANS='C' applies Q^H (conjugate transpose) instead of Q^T.
   */
  _zunmqr_(
    side: number,
    trans: number,
    m: number,
    n: number,
    k: number,
    a: number,
    lda: number,
    tau: number,
    c: number,
    ldc: number,
    work: number,
    lwork: number,
    info: number
  ): void;

  /** ZGELS - Complex least squares via QR/LQ. See DGELS for details. */
  _zgels_(
    trans: number,
    m: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    work: number,
    lwork: number,
    info: number
  ): void;
  /** ZGELSD - Complex least squares via SVD (divide and conquer). See DGELSD for details. */
  _zgelsd_(
    m: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    s: number,
    rcond: number,
    rank: number,
    work: number,
    lwork: number,
    rwork: number,
    iwork: number,
    info: number
  ): void;
  /** ZGELSS - Complex least squares via SVD. See DGELSS for details. */
  _zgelss_(
    m: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    s: number,
    rcond: number,
    rank: number,
    work: number,
    lwork: number,
    rwork: number,
    info: number
  ): void;
  /** ZGELSY - Complex least squares via complete orthogonal factorization. See DGELSY for details. */
  _zgelsy_(
    m: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    jpvt: number,
    rcond: number,
    rank: number,
    work: number,
    lwork: number,
    rwork: number,
    info: number
  ): void;

  /** ZTRTRS - Complex triangular solve. See DTRTRS for details. */
  _ztrtrs_(
    uplo: number,
    trans: number,
    diag: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    info: number
  ): void;
  /** ZTRTRI - Complex triangular matrix inverse. See DTRTRI for details. */
  _ztrtri_(uplo: number, diag: number, n: number, a: number, lda: number, info: number): void;

  /** ZLANGE - Complex matrix norm. See DLANGE for details. */
  _zlange_(norm: number, m: number, n: number, a: number, lda: number, work: number): number;
  /** ZLANHE - Complex Hermitian matrix norm. See DLANSY for details. */
  _zlanhe_(norm: number, uplo: number, n: number, a: number, lda: number, work: number): number;
  /** ZLASWP - Complex row permutations. See DLASWP for details. */
  _zlaswp_(n: number, a: number, lda: number, k1: number, k2: number, ipiv: number, incx: number): void;
  /** ZLASSQ - Complex scaled sum of squares. See DLASSQ for details. */
  _zlassq_(n: number, x: number, incx: number, scale: number, sumsq: number): void;

  // ============================================================
  // Complex Single Precision (C prefix)
  // ============================================================
  // Complex single precision uses HEAPF32 with interleaved real/imaginary pairs.
  // Each complex number occupies 8 bytes (two 4-byte floats).
  // Same API as Z prefix functions but with single precision.

  /** CGETRF - Complex single precision LU factorization. See ZGETRF/DGETRF for details. */
  _cgetrf_(m: number, n: number, a: number, lda: number, ipiv: number, info: number): void;
  /** CGETRS - Complex single precision triangular solve. See ZGETRS for details. */
  _cgetrs_(
    trans: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    ipiv: number,
    b: number,
    ldb: number,
    info: number
  ): void;
  /** CGETRI - Complex single precision matrix inverse. See ZGETRI for details. */
  _cgetri_(
    n: number,
    a: number,
    lda: number,
    ipiv: number,
    work: number,
    lwork: number,
    info: number
  ): void;
  /** CGESV - Complex single precision driver for A*X=B. See ZGESV for details. */
  _cgesv_(
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    ipiv: number,
    b: number,
    ldb: number,
    info: number
  ): void;
  /** CGECON - Complex single precision condition number estimate. See ZGECON for details. */
  _cgecon_(
    norm: number,
    n: number,
    a: number,
    lda: number,
    anorm: number,
    rcond: number,
    work: number,
    rwork: number,
    info: number
  ): void;

  /** CPOTRF - Complex single Hermitian positive definite Cholesky. See ZPOTRF for details. */
  _cpotrf_(uplo: number, n: number, a: number, lda: number, info: number): void;
  /** CPOTRS - Complex single solve using Cholesky. See ZPOTRS for details. */
  _cpotrs_(
    uplo: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    info: number
  ): void;
  /** CPOTRI - Complex single inverse using Cholesky. See ZPOTRI for details. */
  _cpotri_(uplo: number, n: number, a: number, lda: number, info: number): void;
  /** CPOSV - Complex single driver for Hermitian positive definite A*X=B. See ZPOSV for details. */
  _cposv_(
    uplo: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    info: number
  ): void;
  /** CPOCON - Complex single condition number for Hermitian positive definite. See ZPOCON for details. */
  _cpocon_(
    uplo: number,
    n: number,
    a: number,
    lda: number,
    anorm: number,
    rcond: number,
    work: number,
    rwork: number,
    info: number
  ): void;

  /** CGEEV - Complex single eigenvalues/eigenvectors of general matrix. See ZGEEV for details. */
  _cgeev_(
    jobvl: number,
    jobvr: number,
    n: number,
    a: number,
    lda: number,
    w: number,
    vl: number,
    ldvl: number,
    vr: number,
    ldvr: number,
    work: number,
    lwork: number,
    rwork: number,
    info: number
  ): void;
  /** CHEEV - Complex single Hermitian eigenvalues. See ZHEEV for details. */
  _cheev_(
    jobz: number,
    uplo: number,
    n: number,
    a: number,
    lda: number,
    w: number,
    work: number,
    lwork: number,
    rwork: number,
    info: number
  ): void;
  /** CHEEVD - Complex single Hermitian eigenvalues (divide and conquer). See ZHEEVD for details. */
  _cheevd_(
    jobz: number,
    uplo: number,
    n: number,
    a: number,
    lda: number,
    w: number,
    work: number,
    lwork: number,
    rwork: number,
    lrwork: number,
    iwork: number,
    liwork: number,
    info: number
  ): void;
  /** CHEEVR - Complex single selected Hermitian eigenvalues (RRR). See ZHEEVR for details. */
  _cheevr_(
    jobz: number,
    range: number,
    uplo: number,
    n: number,
    a: number,
    lda: number,
    vl: number,
    vu: number,
    il: number,
    iu: number,
    abstol: number,
    m: number,
    w: number,
    z: number,
    ldz: number,
    isuppz: number,
    work: number,
    lwork: number,
    rwork: number,
    lrwork: number,
    iwork: number,
    liwork: number,
    info: number
  ): void;
  /** CGEEVX - Complex single expert eigenvalue driver. See ZGEEVX for details. */
  _cgeevx_(
    balanc: number,
    jobvl: number,
    jobvr: number,
    sense: number,
    n: number,
    a: number,
    lda: number,
    w: number,
    vl: number,
    ldvl: number,
    vr: number,
    ldvr: number,
    ilo: number,
    ihi: number,
    scale: number,
    abnrm: number,
    rconde: number,
    rcondv: number,
    work: number,
    lwork: number,
    rwork: number,
    info: number
  ): void;

  /** CGESVD - Complex single SVD. See ZGESVD for details. */
  _cgesvd_(
    jobu: number,
    jobvt: number,
    m: number,
    n: number,
    a: number,
    lda: number,
    s: number,
    u: number,
    ldu: number,
    vt: number,
    ldvt: number,
    work: number,
    lwork: number,
    rwork: number,
    info: number
  ): void;
  /** CGESDD - Complex single SVD (divide and conquer). See ZGESDD for details. */
  _cgesdd_(
    jobz: number,
    m: number,
    n: number,
    a: number,
    lda: number,
    s: number,
    u: number,
    ldu: number,
    vt: number,
    ldvt: number,
    work: number,
    lwork: number,
    rwork: number,
    iwork: number,
    info: number
  ): void;

  /** CGEQRF - Complex single QR factorization. See ZGEQRF for details. */
  _cgeqrf_(
    m: number,
    n: number,
    a: number,
    lda: number,
    tau: number,
    work: number,
    lwork: number,
    info: number
  ): void;
  /** CGEQP3 - Complex single QR with column pivoting. See ZGEQP3 for details. */
  _cgeqp3_(
    m: number,
    n: number,
    a: number,
    lda: number,
    jpvt: number,
    tau: number,
    work: number,
    lwork: number,
    rwork: number,
    info: number
  ): void;
  /** CUNGQR - Complex single generate unitary Q from QR. See ZUNGQR for details. */
  _cungqr_(
    m: number,
    n: number,
    k: number,
    a: number,
    lda: number,
    tau: number,
    work: number,
    lwork: number,
    info: number
  ): void;
  /** CUNMQR - Complex single multiply by unitary Q from QR. See ZUNMQR for details. */
  _cunmqr_(
    side: number,
    trans: number,
    m: number,
    n: number,
    k: number,
    a: number,
    lda: number,
    tau: number,
    c: number,
    ldc: number,
    work: number,
    lwork: number,
    info: number
  ): void;

  /** CGELS - Complex single least squares via QR/LQ. See ZGELS for details. */
  _cgels_(
    trans: number,
    m: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    work: number,
    lwork: number,
    info: number
  ): void;
  /** CGELSD - Complex single least squares via SVD (divide and conquer). See ZGELSD for details. */
  _cgelsd_(
    m: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    s: number,
    rcond: number,
    rank: number,
    work: number,
    lwork: number,
    rwork: number,
    iwork: number,
    info: number
  ): void;
  /** CGELSS - Complex single least squares via SVD. See ZGELSS for details. */
  _cgelss_(
    m: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    s: number,
    rcond: number,
    rank: number,
    work: number,
    lwork: number,
    rwork: number,
    info: number
  ): void;
  /** CGELSY - Complex single least squares via complete orthogonal factorization. See ZGELSY for details. */
  _cgelsy_(
    m: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    jpvt: number,
    rcond: number,
    rank: number,
    work: number,
    lwork: number,
    rwork: number,
    info: number
  ): void;

  /** CTRTRS - Complex single triangular solve. See ZTRTRS for details. */
  _ctrtrs_(
    uplo: number,
    trans: number,
    diag: number,
    n: number,
    nrhs: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    info: number
  ): void;
  /** CTRTRI - Complex single triangular matrix inverse. See ZTRTRI for details. */
  _ctrtri_(uplo: number, diag: number, n: number, a: number, lda: number, info: number): void;

  /** CLANGE - Complex single matrix norm. See ZLANGE for details. */
  _clange_(norm: number, m: number, n: number, a: number, lda: number, work: number): number;
  /** CLANHE - Complex single Hermitian matrix norm. See ZLANHE for details. */
  _clanhe_(norm: number, uplo: number, n: number, a: number, lda: number, work: number): number;
  /** CLASWP - Complex single row permutations. See ZLASWP for details. */
  _claswp_(n: number, a: number, lda: number, k1: number, k2: number, ipiv: number, incx: number): void;
  /** CLASSQ - Complex single scaled sum of squares. See ZLASSQ for details. */
  _classq_(n: number, x: number, incx: number, scale: number, sumsq: number): void;

  // ============================================================
  // BLAS Level 1 (Vector operations)
  // ============================================================
  // Level 1 BLAS perform O(n) operations on vectors.
  // INCX/INCY are strides (use 1 for contiguous arrays).
  // All parameters are pointers except where noted.

  /**
   * DAXPY - Double precision y = alpha*x + y
   *
   * Computes a constant times a vector plus a vector.
   * @param n - Number of elements (pointer to integer)
   * @param alpha - Scalar multiplier (pointer to double)
   * @param x - Input vector X (pointer to double array)
   * @param incx - Stride for X (pointer to integer)
   * @param y - Input/output vector Y, overwritten with result (pointer to double array)
   * @param incy - Stride for Y (pointer to integer)
   */
  _daxpy_(n: number, alpha: number, x: number, incx: number, y: number, incy: number): void;
  /** DCOPY - Copy vector: y = x */
  _dcopy_(n: number, x: number, incx: number, y: number, incy: number): void;
  /** DSCAL - Scale vector: x = alpha*x */
  _dscal_(n: number, alpha: number, x: number, incx: number): void;
  /** DSWAP - Swap vectors: interchange x and y */
  _dswap_(n: number, x: number, incx: number, y: number, incy: number): void;
  /**
   * DDOT - Double precision dot product
   * @returns x^T * y (sum of x[i]*y[i])
   */
  _ddot_(n: number, x: number, incx: number, y: number, incy: number): number;
  /**
   * DNRM2 - Double precision Euclidean norm
   * @returns ||x||_2 = sqrt(sum of x[i]^2)
   */
  _dnrm2_(n: number, x: number, incx: number): number;
  /**
   * DASUM - Sum of absolute values
   * @returns sum of |x[i]|
   */
  _dasum_(n: number, x: number, incx: number): number;
  /**
   * IDAMAX - Index of maximum absolute value (1-indexed)
   * @returns Index k such that |x[k]| = max|x[i]| (Fortran 1-based indexing)
   */
  _idamax_(n: number, x: number, incx: number): number;

  // Single precision versions - same signatures as double precision
  /** SAXPY - Single precision y = alpha*x + y */
  _saxpy_(n: number, alpha: number, x: number, incx: number, y: number, incy: number): void;
  /** SCOPY - Single precision copy: y = x */
  _scopy_(n: number, x: number, incx: number, y: number, incy: number): void;
  /** SSCAL - Single precision scale: x = alpha*x */
  _sscal_(n: number, alpha: number, x: number, incx: number): void;
  /** SSWAP - Single precision swap vectors */
  _sswap_(n: number, x: number, incx: number, y: number, incy: number): void;
  /** SDOT - Single precision dot product */
  _sdot_(n: number, x: number, incx: number, y: number, incy: number): number;
  /** SNRM2 - Single precision Euclidean norm */
  _snrm2_(n: number, x: number, incx: number): number;
  /** SASUM - Single precision sum of absolute values */
  _sasum_(n: number, x: number, incx: number): number;
  /** ISAMAX - Single precision index of max absolute value (1-indexed) */
  _isamax_(n: number, x: number, incx: number): number;

  // Complex double precision - note different calling conventions for dot products
  /** ZAXPY - Complex double y = alpha*x + y. Alpha is pointer to complex (2 doubles). */
  _zaxpy_(n: number, alpha: number, x: number, incx: number, y: number, incy: number): void;
  /** ZCOPY - Complex double copy: y = x */
  _zcopy_(n: number, x: number, incx: number, y: number, incy: number): void;
  /** ZSCAL - Complex double scale: x = alpha*x. Alpha is pointer to complex. */
  _zscal_(n: number, alpha: number, x: number, incx: number): void;
  /** ZSWAP - Complex double swap vectors */
  _zswap_(n: number, x: number, incx: number, y: number, incy: number): void;
  /**
   * ZDOTC - Complex double conjugated dot product: x^H * y
   *
   * Returns result via pointer (Fortran complex return convention).
   * @param ret - Pointer where result will be stored (2 doubles: real, imag)
   */
  _zdotc_(ret: number, n: number, x: number, incx: number, y: number, incy: number): void;
  /**
   * ZDOTU - Complex double unconjugated dot product: x^T * y
   * @param ret - Pointer where result will be stored (2 doubles)
   */
  _zdotu_(ret: number, n: number, x: number, incx: number, y: number, incy: number): void;
  /** DZNRM2 - Complex double Euclidean norm (returns real double) */
  _dznrm2_(n: number, x: number, incx: number): number;
  /** DZASUM - Complex double sum of |Re(x[i])| + |Im(x[i])| (returns real double) */
  _dzasum_(n: number, x: number, incx: number): number;
  /** IZAMAX - Complex double index of max |Re(x[i])| + |Im(x[i])| (1-indexed) */
  _izamax_(n: number, x: number, incx: number): number;

  // Complex single precision
  /** CAXPY - Complex single y = alpha*x + y */
  _caxpy_(n: number, alpha: number, x: number, incx: number, y: number, incy: number): void;
  /** CCOPY - Complex single copy: y = x */
  _ccopy_(n: number, x: number, incx: number, y: number, incy: number): void;
  /** CSCAL - Complex single scale: x = alpha*x */
  _cscal_(n: number, alpha: number, x: number, incx: number): void;
  /** CSWAP - Complex single swap vectors */
  _cswap_(n: number, x: number, incx: number, y: number, incy: number): void;
  /** CDOTC - Complex single conjugated dot product: x^H * y. Result via pointer. */
  _cdotc_(ret: number, n: number, x: number, incx: number, y: number, incy: number): void;
  /** CDOTU - Complex single unconjugated dot product: x^T * y. Result via pointer. */
  _cdotu_(ret: number, n: number, x: number, incx: number, y: number, incy: number): void;
  /** SCNRM2 - Complex single Euclidean norm (returns real float) */
  _scnrm2_(n: number, x: number, incx: number): number;
  /** SCASUM - Complex single sum of |Re| + |Im| (returns real float) */
  _scasum_(n: number, x: number, incx: number): number;
  /** ICAMAX - Complex single index of max |Re| + |Im| (1-indexed) */
  _icamax_(n: number, x: number, incx: number): number;

  // ============================================================
  // BLAS Level 2 (Matrix-Vector operations)
  // ============================================================
  // Level 2 BLAS perform O(n^2) operations involving a matrix and vector(s).
  // Matrices are stored in column-major order.

  /**
   * DGEMV - Double precision general matrix-vector multiply
   *
   * Computes y = alpha*op(A)*x + beta*y where op(A) = A or A^T.
   *
   * @param trans - 'N' for A*x, 'T' for A^T*x, 'C' for A^H*x (pointer to char as integer)
   * @param m - Number of rows of A (pointer to integer)
   * @param n - Number of columns of A (pointer to integer)
   * @param alpha - Scalar multiplier for A*x (pointer to double)
   * @param a - M-by-N matrix A in column-major order (pointer to double array)
   * @param lda - Leading dimension of A (pointer to integer)
   * @param x - Input vector (pointer to double array)
   * @param incx - Stride for x (pointer to integer)
   * @param beta - Scalar multiplier for y (pointer to double)
   * @param y - Input/output vector, overwritten with result (pointer to double array)
   * @param incy - Stride for y (pointer to integer)
   */
  _dgemv_(
    trans: number,
    m: number,
    n: number,
    alpha: number,
    a: number,
    lda: number,
    x: number,
    incx: number,
    beta: number,
    y: number,
    incy: number
  ): void;
  /**
   * DTRSV - Solve triangular system A*x = b
   *
   * Solves A*x = b where A is triangular. x overwrites b.
   * @param uplo - 'U' = upper triangular, 'L' = lower triangular
   * @param trans - 'N' = A*x=b, 'T' = A^T*x=b
   * @param diag - 'N' = non-unit diagonal, 'U' = unit diagonal
   */
  _dtrsv_(
    uplo: number,
    trans: number,
    diag: number,
    n: number,
    a: number,
    lda: number,
    x: number,
    incx: number
  ): void;
  /**
   * DTRMV - Triangular matrix-vector multiply: x = A*x or x = A^T*x
   */
  _dtrmv_(
    uplo: number,
    trans: number,
    diag: number,
    n: number,
    a: number,
    lda: number,
    x: number,
    incx: number
  ): void;
  /**
   * DGER - Rank-1 update: A = alpha*x*y^T + A
   */
  _dger_(
    m: number,
    n: number,
    alpha: number,
    x: number,
    incx: number,
    y: number,
    incy: number,
    a: number,
    lda: number
  ): void;
  /**
   * DSYR - Symmetric rank-1 update: A = alpha*x*x^T + A
   *
   * Only the upper or lower triangle is updated.
   */
  _dsyr_(uplo: number, n: number, alpha: number, x: number, incx: number, a: number, lda: number): void;
  /**
   * DSYR2 - Symmetric rank-2 update: A = alpha*x*y^T + alpha*y*x^T + A
   */
  _dsyr2_(
    uplo: number,
    n: number,
    alpha: number,
    x: number,
    incx: number,
    y: number,
    incy: number,
    a: number,
    lda: number
  ): void;
  /**
   * DSYMV - Symmetric matrix-vector multiply: y = alpha*A*x + beta*y
   *
   * A is symmetric, only upper or lower triangle is referenced.
   */
  _dsymv_(
    uplo: number,
    n: number,
    alpha: number,
    a: number,
    lda: number,
    x: number,
    incx: number,
    beta: number,
    y: number,
    incy: number
  ): void;

  // Single precision Level 2 BLAS
  /** SGEMV - Single precision general matrix-vector multiply */
  _sgemv_(
    trans: number,
    m: number,
    n: number,
    alpha: number,
    a: number,
    lda: number,
    x: number,
    incx: number,
    beta: number,
    y: number,
    incy: number
  ): void;
  /** STRSV - Single precision triangular solve */
  _strsv_(
    uplo: number,
    trans: number,
    diag: number,
    n: number,
    a: number,
    lda: number,
    x: number,
    incx: number
  ): void;
  /** STRMV - Single precision triangular matrix-vector multiply */
  _strmv_(
    uplo: number,
    trans: number,
    diag: number,
    n: number,
    a: number,
    lda: number,
    x: number,
    incx: number
  ): void;
  /** SGER - Single precision rank-1 update */
  _sger_(
    m: number,
    n: number,
    alpha: number,
    x: number,
    incx: number,
    y: number,
    incy: number,
    a: number,
    lda: number
  ): void;
  /** SSYR - Single precision symmetric rank-1 update */
  _ssyr_(uplo: number, n: number, alpha: number, x: number, incx: number, a: number, lda: number): void;
  /** SSYR2 - Single precision symmetric rank-2 update */
  _ssyr2_(
    uplo: number,
    n: number,
    alpha: number,
    x: number,
    incx: number,
    y: number,
    incy: number,
    a: number,
    lda: number
  ): void;
  /** SSYMV - Single precision symmetric matrix-vector multiply */
  _ssymv_(
    uplo: number,
    n: number,
    alpha: number,
    a: number,
    lda: number,
    x: number,
    incx: number,
    beta: number,
    y: number,
    incy: number
  ): void;

  // Complex double precision Level 2 BLAS
  /** ZGEMV - Complex double general matrix-vector multiply */
  _zgemv_(
    trans: number,
    m: number,
    n: number,
    alpha: number,
    a: number,
    lda: number,
    x: number,
    incx: number,
    beta: number,
    y: number,
    incy: number
  ): void;
  /** ZTRSV - Complex double triangular solve */
  _ztrsv_(
    uplo: number,
    trans: number,
    diag: number,
    n: number,
    a: number,
    lda: number,
    x: number,
    incx: number
  ): void;
  /** ZTRMV - Complex double triangular matrix-vector multiply */
  _ztrmv_(
    uplo: number,
    trans: number,
    diag: number,
    n: number,
    a: number,
    lda: number,
    x: number,
    incx: number
  ): void;
  /**
   * ZGERC - Complex double conjugate rank-1 update: A = alpha*x*y^H + A
   */
  _zgerc_(
    m: number,
    n: number,
    alpha: number,
    x: number,
    incx: number,
    y: number,
    incy: number,
    a: number,
    lda: number
  ): void;
  /**
   * ZGERU - Complex double unconjugated rank-1 update: A = alpha*x*y^T + A
   */
  _zgeru_(
    m: number,
    n: number,
    alpha: number,
    x: number,
    incx: number,
    y: number,
    incy: number,
    a: number,
    lda: number
  ): void;
  /**
   * ZHER - Hermitian rank-1 update: A = alpha*x*x^H + A
   *
   * Alpha must be real (pointer to double, not complex).
   */
  _zher_(uplo: number, n: number, alpha: number, x: number, incx: number, a: number, lda: number): void;
  /**
   * ZHER2 - Hermitian rank-2 update: A = alpha*x*y^H + conj(alpha)*y*x^H + A
   */
  _zher2_(
    uplo: number,
    n: number,
    alpha: number,
    x: number,
    incx: number,
    y: number,
    incy: number,
    a: number,
    lda: number
  ): void;
  /**
   * ZHEMV - Hermitian matrix-vector multiply: y = alpha*A*x + beta*y
   */
  _zhemv_(
    uplo: number,
    n: number,
    alpha: number,
    a: number,
    lda: number,
    x: number,
    incx: number,
    beta: number,
    y: number,
    incy: number
  ): void;

  // Complex single precision Level 2 BLAS
  /** CGEMV - Complex single general matrix-vector multiply */
  _cgemv_(
    trans: number,
    m: number,
    n: number,
    alpha: number,
    a: number,
    lda: number,
    x: number,
    incx: number,
    beta: number,
    y: number,
    incy: number
  ): void;
  /** CTRSV - Complex single triangular solve */
  _ctrsv_(
    uplo: number,
    trans: number,
    diag: number,
    n: number,
    a: number,
    lda: number,
    x: number,
    incx: number
  ): void;
  /** CTRMV - Complex single triangular matrix-vector multiply */
  _ctrmv_(
    uplo: number,
    trans: number,
    diag: number,
    n: number,
    a: number,
    lda: number,
    x: number,
    incx: number
  ): void;
  /** CGERC - Complex single conjugate rank-1 update */
  _cgerc_(
    m: number,
    n: number,
    alpha: number,
    x: number,
    incx: number,
    y: number,
    incy: number,
    a: number,
    lda: number
  ): void;
  /** CGERU - Complex single unconjugated rank-1 update */
  _cgeru_(
    m: number,
    n: number,
    alpha: number,
    x: number,
    incx: number,
    y: number,
    incy: number,
    a: number,
    lda: number
  ): void;
  /** CHER - Complex single Hermitian rank-1 update (alpha is real) */
  _cher_(uplo: number, n: number, alpha: number, x: number, incx: number, a: number, lda: number): void;
  /** CHER2 - Complex single Hermitian rank-2 update */
  _cher2_(
    uplo: number,
    n: number,
    alpha: number,
    x: number,
    incx: number,
    y: number,
    incy: number,
    a: number,
    lda: number
  ): void;
  /** CHEMV - Complex single Hermitian matrix-vector multiply */
  _chemv_(
    uplo: number,
    n: number,
    alpha: number,
    a: number,
    lda: number,
    x: number,
    incx: number,
    beta: number,
    y: number,
    incy: number
  ): void;

  // ============================================================
  // BLAS Level 3 (Matrix-Matrix operations)
  // ============================================================
  // Level 3 BLAS perform O(n^3) operations on matrices.
  // These are the most computationally intensive operations.

  /**
   * DGEMM - Double precision general matrix-matrix multiply
   *
   * Computes C = alpha*op(A)*op(B) + beta*C where op(X) = X, X^T, or X^H.
   * This is the most important BLAS routine - the foundation of most linear algebra.
   *
   * @param transa - 'N' = A, 'T' = A^T, 'C' = A^H (pointer to char as integer)
   * @param transb - 'N' = B, 'T' = B^T, 'C' = B^H (pointer to char as integer)
   * @param m - Number of rows of op(A) and C (pointer to integer)
   * @param n - Number of columns of op(B) and C (pointer to integer)
   * @param k - Number of columns of op(A) and rows of op(B) (pointer to integer)
   * @param alpha - Scalar multiplier for A*B (pointer to double)
   * @param a - Matrix A (pointer to double array)
   * @param lda - Leading dimension of A (pointer to integer)
   * @param b - Matrix B (pointer to double array)
   * @param ldb - Leading dimension of B (pointer to integer)
   * @param beta - Scalar multiplier for C (pointer to double)
   * @param c - Input/output matrix C (pointer to double array)
   * @param ldc - Leading dimension of C (pointer to integer)
   */
  _dgemm_(
    transa: number,
    transb: number,
    m: number,
    n: number,
    k: number,
    alpha: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    beta: number,
    c: number,
    ldc: number
  ): void;
  /**
   * DTRSM - Solve triangular matrix equation
   *
   * Solves op(A)*X = alpha*B or X*op(A) = alpha*B where A is triangular.
   * X overwrites B.
   *
   * @param side - 'L' = left (A*X=B), 'R' = right (X*A=B) (pointer to char as integer)
   * @param uplo - 'U' = upper triangular, 'L' = lower triangular
   * @param transa - 'N' = A, 'T' = A^T, 'C' = A^H
   * @param diag - 'N' = non-unit diagonal, 'U' = unit diagonal
   */
  _dtrsm_(
    side: number,
    uplo: number,
    transa: number,
    diag: number,
    m: number,
    n: number,
    alpha: number,
    a: number,
    lda: number,
    b: number,
    ldb: number
  ): void;
  /**
   * DTRMM - Triangular matrix-matrix multiply
   *
   * Computes B = alpha*op(A)*B or B = alpha*B*op(A) where A is triangular.
   */
  _dtrmm_(
    side: number,
    uplo: number,
    transa: number,
    diag: number,
    m: number,
    n: number,
    alpha: number,
    a: number,
    lda: number,
    b: number,
    ldb: number
  ): void;
  /**
   * DSYRK - Symmetric rank-k update
   *
   * Computes C = alpha*A*A^T + beta*C or C = alpha*A^T*A + beta*C
   * where C is symmetric. Only upper or lower triangle is updated.
   */
  _dsyrk_(
    uplo: number,
    trans: number,
    n: number,
    k: number,
    alpha: number,
    a: number,
    lda: number,
    beta: number,
    c: number,
    ldc: number
  ): void;
  /**
   * DSYR2K - Symmetric rank-2k update
   *
   * Computes C = alpha*A*B^T + alpha*B*A^T + beta*C (or transposed version)
   */
  _dsyr2k_(
    uplo: number,
    trans: number,
    n: number,
    k: number,
    alpha: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    beta: number,
    c: number,
    ldc: number
  ): void;

  // Single precision Level 3 BLAS
  /** SGEMM - Single precision general matrix-matrix multiply. See DGEMM for details. */
  _sgemm_(
    transa: number,
    transb: number,
    m: number,
    n: number,
    k: number,
    alpha: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    beta: number,
    c: number,
    ldc: number
  ): void;
  /** STRSM - Single precision triangular solve. See DTRSM for details. */
  _strsm_(
    side: number,
    uplo: number,
    transa: number,
    diag: number,
    m: number,
    n: number,
    alpha: number,
    a: number,
    lda: number,
    b: number,
    ldb: number
  ): void;
  /** STRMM - Single precision triangular matrix-matrix multiply. See DTRMM for details. */
  _strmm_(
    side: number,
    uplo: number,
    transa: number,
    diag: number,
    m: number,
    n: number,
    alpha: number,
    a: number,
    lda: number,
    b: number,
    ldb: number
  ): void;
  /** SSYRK - Single precision symmetric rank-k update. See DSYRK for details. */
  _ssyrk_(
    uplo: number,
    trans: number,
    n: number,
    k: number,
    alpha: number,
    a: number,
    lda: number,
    beta: number,
    c: number,
    ldc: number
  ): void;
  /** SSYR2K - Single precision symmetric rank-2k update. See DSYR2K for details. */
  _ssyr2k_(
    uplo: number,
    trans: number,
    n: number,
    k: number,
    alpha: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    beta: number,
    c: number,
    ldc: number
  ): void;

  // Complex double precision Level 3 BLAS
  /** ZGEMM - Complex double general matrix-matrix multiply. See DGEMM for details. */
  _zgemm_(
    transa: number,
    transb: number,
    m: number,
    n: number,
    k: number,
    alpha: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    beta: number,
    c: number,
    ldc: number
  ): void;
  /** ZTRSM - Complex double triangular solve. See DTRSM for details. */
  _ztrsm_(
    side: number,
    uplo: number,
    transa: number,
    diag: number,
    m: number,
    n: number,
    alpha: number,
    a: number,
    lda: number,
    b: number,
    ldb: number
  ): void;
  /** ZTRMM - Complex double triangular matrix-matrix multiply. See DTRMM for details. */
  _ztrmm_(
    side: number,
    uplo: number,
    transa: number,
    diag: number,
    m: number,
    n: number,
    alpha: number,
    a: number,
    lda: number,
    b: number,
    ldb: number
  ): void;
  /**
   * ZHERK - Hermitian rank-k update
   *
   * Computes C = alpha*A*A^H + beta*C or C = alpha*A^H*A + beta*C
   * where C is Hermitian. Alpha and beta are real.
   */
  _zherk_(
    uplo: number,
    trans: number,
    n: number,
    k: number,
    alpha: number,
    a: number,
    lda: number,
    beta: number,
    c: number,
    ldc: number
  ): void;
  /**
   * ZHER2K - Hermitian rank-2k update
   *
   * Computes C = alpha*A*B^H + conj(alpha)*B*A^H + beta*C
   * Beta is real.
   */
  _zher2k_(
    uplo: number,
    trans: number,
    n: number,
    k: number,
    alpha: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    beta: number,
    c: number,
    ldc: number
  ): void;

  // Complex single precision Level 3 BLAS
  /** CGEMM - Complex single general matrix-matrix multiply. See DGEMM for details. */
  _cgemm_(
    transa: number,
    transb: number,
    m: number,
    n: number,
    k: number,
    alpha: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    beta: number,
    c: number,
    ldc: number
  ): void;
  /** CTRSM - Complex single triangular solve. See DTRSM for details. */
  _ctrsm_(
    side: number,
    uplo: number,
    transa: number,
    diag: number,
    m: number,
    n: number,
    alpha: number,
    a: number,
    lda: number,
    b: number,
    ldb: number
  ): void;
  /** CTRMM - Complex single triangular matrix-matrix multiply. See DTRMM for details. */
  _ctrmm_(
    side: number,
    uplo: number,
    transa: number,
    diag: number,
    m: number,
    n: number,
    alpha: number,
    a: number,
    lda: number,
    b: number,
    ldb: number
  ): void;
  /** CHERK - Complex single Hermitian rank-k update. See ZHERK for details. */
  _cherk_(
    uplo: number,
    trans: number,
    n: number,
    k: number,
    alpha: number,
    a: number,
    lda: number,
    beta: number,
    c: number,
    ldc: number
  ): void;
  /** CHER2K - Complex single Hermitian rank-2k update. See ZHER2K for details. */
  _cher2k_(
    uplo: number,
    trans: number,
    n: number,
    k: number,
    alpha: number,
    a: number,
    lda: number,
    b: number,
    ldb: number,
    beta: number,
    c: number,
    ldc: number
  ): void;

  // ============================================================
  // Utility
  // ============================================================

  /**
   * LSAME - Test if two characters are the same regardless of case
   *
   * @param ca - First character (pointer to char as integer)
   * @param cb - Second character (pointer to char as integer)
   * @returns 1 if same (ignoring case), 0 otherwise
   */
  _lsame_(ca: number, cb: number): number;

  // ============================================================
  // Memory Management
  // ============================================================

  /**
   * Allocate memory in the WASM heap
   *
   * @param size - Number of bytes to allocate
   * @returns Pointer to allocated memory (byte offset in WASM linear memory)
   *
   * @example
   * // Allocate space for 10 doubles
   * const ptr = module._malloc(10 * 8);  // 8 bytes per double
   * // ... use the memory ...
   * module._free(ptr);
   */
  _malloc(size: number): number;

  /**
   * Free previously allocated WASM heap memory
   *
   * @param ptr - Pointer returned by _malloc
   */
  _free(ptr: number): void;

  // ============================================================
  // Emscripten Runtime
  // ============================================================

  /**
   * Float64Array view of WASM linear memory
   *
   * Used to read/write double precision values. Index = ptr / 8.
   * @example
   * module.HEAPF64[ptr / 8] = 3.14;  // Write double at ptr
   * const value = module.HEAPF64[ptr / 8];  // Read double at ptr
   */
  HEAPF64: Float64Array;

  /**
   * Float32Array view of WASM linear memory
   *
   * Used to read/write single precision values. Index = ptr / 4.
   */
  HEAPF32: Float32Array;

  /**
   * Int32Array view of WASM linear memory
   *
   * Used to read/write 32-bit integers. Index = ptr / 4.
   */
  HEAP32: Int32Array;

  /**
   * Uint8Array view of WASM linear memory
   *
   * Used for byte-level access. Index = ptr directly.
   */
  HEAPU8: Uint8Array;

  /**
   * Call a C function by name with automatic type marshalling
   *
   * @param name - C function name (without underscore prefix)
   * @param returnType - Return type: 'number', 'string', 'array', or null for void
   * @param argTypes - Array of argument types: 'number', 'string', 'array'
   * @param args - Array of argument values
   * @returns Function return value
   *
   * @example
   * const result = module.ccall('dnrm2', 'number', ['number', 'number', 'number'], [nPtr, xPtr, incxPtr]);
   */
  ccall(
    name: string,
    returnType: string | null,
    argTypes: string[],
    args: (number | string)[]
  ): number | void;

  /**
   * Create a JavaScript wrapper function for a C function
   *
   * @param name - C function name (without underscore prefix)
   * @param returnType - Return type: 'number', 'string', or null for void
   * @param argTypes - Array of argument types
   * @returns JavaScript function that calls the C function
   *
   * @example
   * const dnrm2 = module.cwrap('dnrm2', 'number', ['number', 'number', 'number']);
   * const norm = dnrm2(nPtr, xPtr, incxPtr);
   */
  cwrap(
    name: string,
    returnType: string | null,
    argTypes: string[]
  ): (...args: (number | string)[]) => number | void;

  /**
   * Read a value from WASM memory
   *
   * @param ptr - Pointer to memory location
   * @param type - Type string: 'i8', 'i16', 'i32', 'i64', 'float', 'double', '*'
   * @returns The value at that memory location
   *
   * @example
   * const intValue = module.getValue(ptr, 'i32');
   * const doubleValue = module.getValue(ptr, 'double');
   */
  getValue(ptr: number, type: string): number;

  /**
   * Write a value to WASM memory
   *
   * @param ptr - Pointer to memory location
   * @param value - Value to write
   * @param type - Type string: 'i8', 'i16', 'i32', 'i64', 'float', 'double', '*'
   *
   * @example
   * module.setValue(ptr, 42, 'i32');
   * module.setValue(ptr, 3.14159, 'double');
   */
  setValue(ptr: number, value: number, type: string): void;
}

/**
 * Options for loading the LAPACK WASM module
 */
export interface LAPACKModuleOptions {
  /**
   * Custom function to locate WASM files
   *
   * Override this to load the .wasm file from a different location (e.g., CDN).
   *
   * @param path - The filename being requested (e.g., 'lapack.wasm')
   * @param scriptDirectory - The directory where the JS loader is located
   * @returns The full URL or path to the file
   *
   * @example
   * const module = await createLAPACKModule({
   *   locateFile: (path) => `https://cdn.example.com/wasm/${path}`
   * });
   */
  locateFile?: (path: string, scriptDirectory: string) => string;
}

/**
 * Factory function type for creating a LAPACK module instance
 *
 * @param options - Optional configuration for module loading
 * @returns Promise that resolves to the initialized LAPACK module
 *
 * @example
 * import createLAPACKModule from 'lawasm';
 *
 * const lapack = await createLAPACKModule();
 *
 * // Allocate memory for a 3x3 matrix and perform LU factorization
 * const n = 3;
 * const aPtr = lapack._malloc(n * n * 8);
 * const ipivPtr = lapack._malloc(n * 4);
 * const infoPtr = lapack._malloc(4);
 * // ... set up matrix data ...
 * lapack._dgetrf_(nPtr, nPtr, aPtr, ldaPtr, ipivPtr, infoPtr);
 */
export type LAPACKModuleFactory = (options?: LAPACKModuleOptions) => Promise<LAPACKModule>;
