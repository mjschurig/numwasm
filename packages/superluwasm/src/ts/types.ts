/**
 * SuperLU WASM Module Type Definitions
 *
 * SuperLU is a library for the direct solution of large, sparse,
 * nonsymmetric systems of linear equations using LU factorization
 * with partial pivoting and triangular solves.
 *
 * Key features:
 * - Sparse matrix storage (CSC/CSR formats)
 * - Supernodal factorization for cache efficiency
 * - Column/row permutations for fill-in reduction
 * - Support for iterative refinement
 * - ILU preconditioning support
 *
 * This module supports all four precision variants:
 * - s: single precision real (float, 4 bytes)
 * - d: double precision real (double, 8 bytes)
 * - c: single precision complex (2 floats, 8 bytes)
 * - z: double precision complex (2 doubles, 16 bytes)
 *
 * Sparse Matrix Formats:
 * - CSC (Compressed Sparse Column): colptr[n+1], rowind[nnz], nzval[nnz]
 * - CSR (Compressed Sparse Row): rowptr[m+1], colind[nnz], nzval[nnz]
 */

// ============================================================
// SuperLU Options Enumerations
// ============================================================

/**
 * Factorization options - controls what factorization work to perform
 */
export enum SuperLUFactOption {
  /** Perform factorization from scratch (default for new matrix) */
  DOFACT = 0,
  /** Re-use column permutation from previous factorization (matrix has same sparsity pattern) */
  SamePattern = 1,
  /** Re-use both column and row permutations (same pattern, similar values) */
  SamePattern_SameRowPerm = 2,
  /** Matrix is already factored, just do triangular solve */
  FACTORED = 3,
}

/**
 * Row permutation strategy - controls how rows are reordered
 */
export enum SuperLURowPerm {
  /** No row permutation (natural ordering) */
  NOROWPERM = 0,
  /** Use MC64 algorithm to maximize diagonal entries */
  LargeDiag_MC64 = 1,
  /** Use AWPM (approximate weight perfect matching) algorithm */
  LargeDiag_AWPM = 2,
  /** Use user-supplied row permutation */
  MY_PERMR = 3,
}

/**
 * Column permutation strategy - controls how columns are reordered for fill-in reduction
 */
export enum SuperLUColPerm {
  /** Natural ordering (no permutation) - fastest but may have high fill-in */
  NATURAL = 0,
  /** Minimum degree ordering on A^T * A - good for unsymmetric matrices */
  MMD_ATA = 1,
  /** Minimum degree ordering on A^T + A - good for structurally symmetric matrices */
  MMD_AT_PLUS_A = 2,
  /** COLAMD approximate minimum degree (recommended default) */
  COLAMD = 3,
  /** METIS nested dissection on A^T + A */
  METIS_AT_PLUS_A = 4,
  /** ParMETIS parallel ordering (not available in WASM) */
  PARMETIS = 5,
  /** Zoltan hypergraph partitioning (not available in WASM) */
  ZOLTAN = 6,
  /** Use user-supplied column permutation */
  MY_PERMC = 7,
}

/**
 * Transpose options for triangular solve
 */
export enum SuperLUTrans {
  /** Solve A * X = B (no transpose) */
  NOTRANS = 0,
  /** Solve A^T * X = B (transpose) */
  TRANS = 1,
  /** Solve A^H * X = B (conjugate transpose, for complex matrices) */
  CONJ = 2,
}

/**
 * ILU drop rule options - controls which elements are dropped in incomplete factorization
 */
export enum SuperLUILUDrop {
  /** Basic threshold dropping */
  DROP_BASIC = 0,
  /** Drop by rows */
  DROP_PROWS = 1,
  /** Drop by columns */
  DROP_COLUMN = 2,
  /** Drop by area */
  DROP_AREA = 3,
  /** Dynamic dropping strategy */
  DROP_DYNAMIC = 4,
  /** Interpolation-based dropping */
  DROP_INTERP = 5,
}

/**
 * Modified ILU options
 */
export enum SuperLUMiluT {
  /** Standard ILU (no modification) */
  SILU = 0,
  /** MILU variant 1: add dropped elements to diagonal */
  SMILU_1 = 1,
  /** MILU variant 2 */
  SMILU_2 = 2,
  /** MILU variant 3 */
  SMILU_3 = 3,
}

/**
 * Matrix norm types for condition number estimation
 */
export enum SuperLUNorm {
  /** 1-norm (maximum column sum of absolute values) */
  ONE_NORM = 0,
  /** 2-norm (largest singular value - expensive to compute) */
  TWO_NORM = 1,
  /** Infinity-norm (maximum row sum of absolute values) */
  INF_NORM = 2,
}

// ============================================================
// SuperMatrix Type Enumerations
// ============================================================

/**
 * SuperMatrix storage types - specifies how matrix data is organized in memory
 */
export enum SuperLUStype {
  /** Compressed Sparse Column (CSC) format - standard input format */
  SLU_NC = 0,
  /** CSC with column permutation applied */
  SLU_NCP = 1,
  /** Compressed Sparse Row (CSR) format */
  SLU_NR = 2,
  /** Supernode Column format (output L factor) */
  SLU_SC = 3,
  /** Supernode Column with permutation */
  SLU_SCP = 4,
  /** Supernode Row format */
  SLU_SR = 5,
  /** Dense matrix in column-major (Fortran) order */
  SLU_DN = 6,
  /** Distributed CSR format (for parallel SuperLU) */
  SLU_NR_loc = 7,
}

/**
 * Data type of matrix elements
 */
export enum SuperLUDtype {
  /** Single precision real (float, 4 bytes) */
  SLU_S = 0,
  /** Double precision real (double, 8 bytes) */
  SLU_D = 1,
  /** Single precision complex (2 floats, 8 bytes total) */
  SLU_C = 2,
  /** Double precision complex (2 doubles, 16 bytes total) */
  SLU_Z = 3,
}

/**
 * Mathematical type of matrix - affects how matrix is interpreted
 */
export enum SuperLUMtype {
  /** General (no special structure) */
  SLU_GE = 0,
  /** Lower triangular with unit diagonal */
  SLU_TRLU = 1,
  /** Lower triangular with non-unit diagonal */
  SLU_TRLL = 2,
  /** Upper triangular with unit diagonal */
  SLU_TRU = 3,
  /** Upper triangular with non-unit diagonal */
  SLU_TRL = 4,
  /** Symmetric, lower triangle stored */
  SLU_SYL = 5,
  /** Symmetric, upper triangle stored */
  SLU_SYU = 6,
  /** Hermitian, lower triangle stored */
  SLU_HEL = 7,
  /** Hermitian, upper triangle stored */
  SLU_HEU = 8,
}

/**
 * SuperLU WebAssembly Module Interface
 *
 * All pointer parameters (ending in Ptr) are memory addresses in the WASM heap.
 * Use module._malloc() to allocate memory and module.HEAP* views to read/write data.
 */
export interface SuperLUModule {
  // ============================================================
  // Double Precision (d) Routines - Most Commonly Used
  // ============================================================

  /**
   * DGSSV - Simple driver for solving A*X = B (double precision)
   *
   * This is the simplest way to solve a sparse linear system. It performs:
   * 1. Column permutation for fill reduction
   * 2. Row permutation for numerical stability
   * 3. LU factorization: Pr * A * Pc = L * U
   * 4. Triangular solves to get X
   *
   * @param optionsPtr - Pointer to superlu_options_t structure (use set_default_options to initialize)
   * @param APtr - Pointer to SuperMatrix A in CSC format (input, may be modified)
   * @param permCPtr - Pointer to column permutation array of size n (output)
   * @param permRPtr - Pointer to row permutation array of size m (output)
   * @param LPtr - Pointer to SuperMatrix L factor (output, must be destroyed after use)
   * @param UPtr - Pointer to SuperMatrix U factor (output, must be destroyed after use)
   * @param BPtr - Pointer to SuperMatrix B (input: RHS, output: solution X)
   * @param statPtr - Pointer to SuperLUStat_t for performance statistics
   * @param infoPtr - Pointer to int: 0=success, <0=illegal argument, >0=singular matrix
   */
  _dgssv(
    optionsPtr: number,
    APtr: number,
    permCPtr: number,
    permRPtr: number,
    LPtr: number,
    UPtr: number,
    BPtr: number,
    statPtr: number,
    infoPtr: number
  ): void;

  /**
   * DGSTRF - LU factorization of a sparse matrix (double precision)
   *
   * Computes the LU factorization of a general sparse matrix A using
   * supernodal techniques for cache efficiency. This is the computational
   * core of SuperLU.
   *
   * Factorization: Pr * A * Pc = L * U
   *
   * @param optionsPtr - Pointer to superlu_options_t structure
   * @param APtr - Pointer to SuperMatrix A in CSC format (permuted columns)
   * @param relax - Supernode relaxation parameter (use sp_ienv(2) for default)
   * @param panelSize - Panel size for blocking (use sp_ienv(1) for default)
   * @param etreePtr - Pointer to elimination tree array of size n
   * @param workPtr - Pointer to workspace, or 0 for automatic allocation
   * @param lwork - Size of workspace in bytes, or 0 for automatic
   * @param permCPtr - Pointer to column permutation (input)
   * @param permRPtr - Pointer to row permutation (output)
   * @param LPtr - Pointer to SuperMatrix L (output)
   * @param UPtr - Pointer to SuperMatrix U (output)
   * @param GluPtr - Pointer to GlobalLU_t structure for L/U storage info
   * @param statPtr - Pointer to SuperLUStat_t
   * @param infoPtr - Pointer to int: 0=success, >0=number of zero pivots encountered
   */
  _dgstrf(
    optionsPtr: number,
    APtr: number,
    relax: number,
    panelSize: number,
    etreePtr: number,
    workPtr: number,
    lwork: number,
    permCPtr: number,
    permRPtr: number,
    LPtr: number,
    UPtr: number,
    GluPtr: number,
    statPtr: number,
    infoPtr: number
  ): void;

  /**
   * DGSTRS - Triangular solve using LU factors (double precision)
   *
   * Solves one of the systems:
   * - A * X = B (trans = NOTRANS)
   * - A^T * X = B (trans = TRANS)
   * - A^H * X = B (trans = CONJ)
   *
   * using the L and U factors from dgstrf.
   *
   * @param trans - Transpose option: NOTRANS=0, TRANS=1, CONJ=2
   * @param LPtr - Pointer to SuperMatrix L from dgstrf
   * @param UPtr - Pointer to SuperMatrix U from dgstrf
   * @param permCPtr - Pointer to column permutation from dgstrf
   * @param permRPtr - Pointer to row permutation from dgstrf
   * @param BPtr - Pointer to SuperMatrix B (input: RHS, output: solution)
   * @param statPtr - Pointer to SuperLUStat_t
   * @param infoPtr - Pointer to int: 0=success
   */
  _dgstrs(
    trans: number,
    LPtr: number,
    UPtr: number,
    permCPtr: number,
    permRPtr: number,
    BPtr: number,
    statPtr: number,
    infoPtr: number
  ): void;

  /**
   * DGSSVX - Expert driver for solving A*X = B (double precision)
   *
   * Full-featured driver that includes:
   * - Optional equilibration (row/column scaling)
   * - Condition number estimation
   * - Iterative refinement for improved accuracy
   * - Forward/backward error bounds
   * - Option to reuse factorization for multiple solves
   *
   * @param optionsPtr - Pointer to superlu_options_t (controls all behavior)
   * @param APtr - Pointer to SuperMatrix A
   * @param permCPtr - Pointer to column permutation
   * @param permRPtr - Pointer to row permutation
   * @param etreePtr - Pointer to elimination tree
   * @param equedPtr - Pointer to char: equilibration type ('N', 'R', 'C', 'B')
   * @param RPtr - Pointer to row scale factors (size m)
   * @param CPtr - Pointer to column scale factors (size n)
   * @param LPtr - Pointer to SuperMatrix L
   * @param UPtr - Pointer to SuperMatrix U
   * @param workPtr - Pointer to workspace
   * @param lwork - Workspace size
   * @param BPtr - Pointer to SuperMatrix B (right-hand side)
   * @param XPtr - Pointer to SuperMatrix X (solution)
   * @param rpgPtr - Pointer to double: reciprocal pivot growth factor
   * @param rcondPtr - Pointer to double: reciprocal condition number estimate
   * @param ferrPtr - Pointer to forward error bound array (size nrhs)
   * @param berrPtr - Pointer to backward error bound array (size nrhs)
   * @param GluPtr - Pointer to GlobalLU_t
   * @param memUsagePtr - Pointer to mem_usage_t structure
   * @param statPtr - Pointer to SuperLUStat_t
   * @param infoPtr - Pointer to int: 0=success
   */
  _dgssvx(
    optionsPtr: number,
    APtr: number,
    permCPtr: number,
    permRPtr: number,
    etreePtr: number,
    equedPtr: number,
    RPtr: number,
    CPtr: number,
    LPtr: number,
    UPtr: number,
    workPtr: number,
    lwork: number,
    BPtr: number,
    XPtr: number,
    rpgPtr: number,
    rcondPtr: number,
    ferrPtr: number,
    berrPtr: number,
    GluPtr: number,
    memUsagePtr: number,
    statPtr: number,
    infoPtr: number
  ): void;

  /**
   * DGSISX - ILU-based driver for solving A*X = B (double precision)
   *
   * Uses incomplete LU factorization for iterative solver preconditioning.
   * The ILU factors can be used as preconditioners in Krylov methods.
   *
   * @param optionsPtr - Pointer to superlu_options_t (set ILU options)
   * @param APtr - Pointer to SuperMatrix A
   * @param permCPtr - Pointer to column permutation
   * @param permRPtr - Pointer to row permutation
   * @param etreePtr - Pointer to elimination tree
   * @param equedPtr - Pointer to equilibration type
   * @param RPtr - Pointer to row scale factors
   * @param CPtr - Pointer to column scale factors
   * @param LPtr - Pointer to SuperMatrix L (ILU factor)
   * @param UPtr - Pointer to SuperMatrix U (ILU factor)
   * @param workPtr - Pointer to workspace
   * @param lwork - Workspace size
   * @param BPtr - Pointer to SuperMatrix B
   * @param XPtr - Pointer to SuperMatrix X
   * @param rpgPtr - Pointer to reciprocal pivot growth
   * @param rcondPtr - Pointer to condition number estimate
   * @param GluPtr - Pointer to GlobalLU_t
   * @param memUsagePtr - Pointer to mem_usage_t
   * @param statPtr - Pointer to SuperLUStat_t
   * @param infoPtr - Pointer to int
   */
  _dgsisx(
    optionsPtr: number,
    APtr: number,
    permCPtr: number,
    permRPtr: number,
    etreePtr: number,
    equedPtr: number,
    RPtr: number,
    CPtr: number,
    LPtr: number,
    UPtr: number,
    workPtr: number,
    lwork: number,
    BPtr: number,
    XPtr: number,
    rpgPtr: number,
    rcondPtr: number,
    GluPtr: number,
    memUsagePtr: number,
    statPtr: number,
    infoPtr: number
  ): void;

  /**
   * DGSITRF - Incomplete LU factorization (double precision)
   *
   * Computes an incomplete LU factorization with dropping based on
   * fill-in level or drop tolerance. Useful for preconditioning.
   *
   * @param optionsPtr - Pointer to options (ILU parameters: drop tolerance, fill ratio)
   * @param APtr - Pointer to SuperMatrix A
   * @param relax - Relaxation parameter
   * @param panelSize - Panel size
   * @param etreePtr - Pointer to elimination tree
   * @param workPtr - Pointer to workspace
   * @param lwork - Workspace size
   * @param permCPtr - Pointer to column permutation
   * @param permRPtr - Pointer to row permutation
   * @param LPtr - Pointer to incomplete L factor
   * @param UPtr - Pointer to incomplete U factor
   * @param GluPtr - Pointer to GlobalLU_t
   * @param statPtr - Pointer to SuperLUStat_t
   * @param infoPtr - Pointer to int
   */
  _dgsitrf(
    optionsPtr: number,
    APtr: number,
    relax: number,
    panelSize: number,
    etreePtr: number,
    workPtr: number,
    lwork: number,
    permCPtr: number,
    permRPtr: number,
    LPtr: number,
    UPtr: number,
    GluPtr: number,
    statPtr: number,
    infoPtr: number
  ): void;

  /**
   * DGSCON - Estimate reciprocal condition number (double precision)
   *
   * Estimates 1/κ(A) using the LU factorization. A small value indicates
   * the matrix is nearly singular.
   *
   * @param normPtr - Pointer to norm type (ONE_NORM or INF_NORM)
   * @param LPtr - Pointer to L factor from dgstrf
   * @param UPtr - Pointer to U factor from dgstrf
   * @param anorm - The norm of original matrix A (compute before factorization)
   * @param rcondPtr - Pointer to output reciprocal condition number
   * @param statPtr - Pointer to SuperLUStat_t
   * @param infoPtr - Pointer to int
   */
  _dgscon(
    normPtr: number,
    LPtr: number,
    UPtr: number,
    anorm: number,
    rcondPtr: number,
    statPtr: number,
    infoPtr: number
  ): void;

  /**
   * DGSEQU - Compute row and column scaling factors (double precision)
   *
   * Computes scaling factors R and C such that diag(R)*A*diag(C) has
   * entries with absolute value close to 1. Improves numerical stability.
   *
   * @param APtr - Pointer to SuperMatrix A
   * @param RPtr - Pointer to output row scale factors (size m)
   * @param CPtr - Pointer to output column scale factors (size n)
   * @param rowcndPtr - Pointer to output row condition: min(R)/max(R)
   * @param colcndPtr - Pointer to output column condition: min(C)/max(C)
   * @param amaxPtr - Pointer to output max absolute element
   * @param infoPtr - Pointer to int
   */
  _dgsequ(
    APtr: number,
    RPtr: number,
    CPtr: number,
    rowcndPtr: number,
    colcndPtr: number,
    amaxPtr: number,
    infoPtr: number
  ): void;

  /**
   * DGSRFS - Iterative refinement of solution (double precision)
   *
   * Improves the computed solution to A*X = B and provides error bounds.
   * Uses the original matrix A and the LU factors.
   *
   * @param trans - Transpose option
   * @param APtr - Pointer to original SuperMatrix A
   * @param LPtr - Pointer to L factor
   * @param UPtr - Pointer to U factor
   * @param permCPtr - Pointer to column permutation
   * @param permRPtr - Pointer to row permutation
   * @param equedPtr - Pointer to equilibration type used
   * @param RPtr - Pointer to row scales
   * @param CPtr - Pointer to column scales
   * @param BPtr - Pointer to original RHS B
   * @param XPtr - Pointer to solution X (input/output: refined)
   * @param ferrPtr - Pointer to output forward error bounds
   * @param berrPtr - Pointer to output backward error bounds
   * @param statPtr - Pointer to SuperLUStat_t
   * @param infoPtr - Pointer to int
   */
  _dgsrfs(
    trans: number,
    APtr: number,
    LPtr: number,
    UPtr: number,
    permCPtr: number,
    permRPtr: number,
    equedPtr: number,
    RPtr: number,
    CPtr: number,
    BPtr: number,
    XPtr: number,
    ferrPtr: number,
    berrPtr: number,
    statPtr: number,
    infoPtr: number
  ): void;

  /**
   * DLAQGS - Apply equilibration scaling to matrix (double precision)
   *
   * Scales the matrix A by diag(R)*A*diag(C) using scale factors
   * computed by dgsequ.
   *
   * @param APtr - Pointer to SuperMatrix A (modified in place)
   * @param RPtr - Pointer to row scale factors
   * @param CPtr - Pointer to column scale factors
   * @param rowcnd - Row condition number from dgsequ
   * @param colcnd - Column condition number from dgsequ
   * @param amax - Max absolute element from dgsequ
   * @param equedPtr - Pointer to output: scaling applied ('N', 'R', 'C', or 'B')
   */
  _dlaqgs(
    APtr: number,
    RPtr: number,
    CPtr: number,
    rowcnd: number,
    colcnd: number,
    amax: number,
    equedPtr: number
  ): void;

  /**
   * dCreate_CompCol_Matrix - Create a sparse matrix in CSC format (double precision)
   *
   * Creates a SuperMatrix wrapper around existing CSC data. The data arrays
   * are NOT copied - they must remain valid for the lifetime of the matrix.
   *
   * CSC format: colptr[j] to colptr[j+1]-1 gives indices into rowind/nzval for column j
   *
   * @param APtr - Pointer to SuperMatrix structure to initialize (output)
   * @param m - Number of rows
   * @param n - Number of columns
   * @param nnz - Number of nonzeros
   * @param nzvalPtr - Pointer to double array of nonzero values (size nnz)
   * @param rowindPtr - Pointer to int array of row indices (size nnz)
   * @param colptrPtr - Pointer to int array of column pointers (size n+1)
   * @param stype - Storage type (use SLU_NC for standard CSC)
   * @param dtype - Data type (use SLU_D for double)
   * @param mtype - Matrix type (use SLU_GE for general)
   */
  _dCreate_CompCol_Matrix(
    APtr: number,
    m: number,
    n: number,
    nnz: number,
    nzvalPtr: number,
    rowindPtr: number,
    colptrPtr: number,
    stype: number,
    dtype: number,
    mtype: number
  ): void;

  /**
   * dCreate_CompRow_Matrix - Create a sparse matrix in CSR format (double precision)
   *
   * CSR format: rowptr[i] to rowptr[i+1]-1 gives indices into colind/nzval for row i
   *
   * @param APtr - Pointer to SuperMatrix structure (output)
   * @param m - Number of rows
   * @param n - Number of columns
   * @param nnz - Number of nonzeros
   * @param nzvalPtr - Pointer to nonzero values (size nnz)
   * @param colindPtr - Pointer to column indices (size nnz)
   * @param rowptrPtr - Pointer to row pointers (size m+1)
   * @param stype - Storage type (use SLU_NR for CSR)
   * @param dtype - Data type (use SLU_D)
   * @param mtype - Matrix type (use SLU_GE)
   */
  _dCreate_CompRow_Matrix(
    APtr: number,
    m: number,
    n: number,
    nnz: number,
    nzvalPtr: number,
    colindPtr: number,
    rowptrPtr: number,
    stype: number,
    dtype: number,
    mtype: number
  ): void;

  /**
   * dCreate_Dense_Matrix - Create a dense matrix (double precision)
   *
   * Creates a SuperMatrix wrapper around column-major dense data.
   *
   * @param XPtr - Pointer to SuperMatrix structure (output)
   * @param m - Number of rows
   * @param n - Number of columns
   * @param xPtr - Pointer to double array of values (size m*n in column-major order)
   * @param ldx - Leading dimension (usually m)
   * @param stype - Storage type (use SLU_DN for dense)
   * @param dtype - Data type (use SLU_D)
   * @param mtype - Matrix type (use SLU_GE)
   */
  _dCreate_Dense_Matrix(
    XPtr: number,
    m: number,
    n: number,
    xPtr: number,
    ldx: number,
    stype: number,
    dtype: number,
    mtype: number
  ): void;

  /** Copy a compressed column matrix. B = A (deep copy) */
  _dCopy_CompCol_Matrix(APtr: number, BPtr: number): void;

  /**
   * Copy a dense matrix
   * @param M - Number of rows to copy
   * @param N - Number of columns to copy
   * @param XPtr - Pointer to source data
   * @param ldx - Leading dimension of source
   * @param YPtr - Pointer to destination data
   * @param ldy - Leading dimension of destination
   */
  _dCopy_Dense_Matrix(M: number, N: number, XPtr: number, ldx: number, YPtr: number, ldy: number): void;

  /**
   * dCompRow_to_CompCol - Convert CSR to CSC format (double precision)
   *
   * Transposes the sparse matrix storage format. Output arrays are allocated.
   *
   * @param m - Number of rows
   * @param n - Number of columns
   * @param nnz - Number of nonzeros
   * @param aPtr - Pointer to CSR values
   * @param colindPtr - Pointer to CSR column indices
   * @param rowptrPtr - Pointer to CSR row pointers
   * @param atPtrPtr - Pointer to pointer: output CSC values (allocated)
   * @param rowindPtrPtr - Pointer to pointer: output CSC row indices (allocated)
   * @param colptrPtrPtr - Pointer to pointer: output CSC column pointers (allocated)
   */
  _dCompRow_to_CompCol(
    m: number,
    n: number,
    nnz: number,
    aPtr: number,
    colindPtr: number,
    rowptrPtr: number,
    atPtrPtr: number,
    rowindPtrPtr: number,
    colptrPtrPtr: number
  ): void;

  /** Allocate storage for sparse matrix A (n columns, nnz nonzeros) */
  _dallocateA(n: number, nnz: number, aPtr: number, asub: number, xa: number): void;

  /** Print compressed column matrix to stdout (for debugging) */
  _dPrint_CompCol_Matrix(whatPtr: number, APtr: number): void;

  /** Print dense matrix to stdout (for debugging) */
  _dPrint_Dense_Matrix(whatPtr: number, APtr: number): void;

  /** Generate a random test solution vector */
  _dGenXtrue(n: number, nrhs: number, xPtr: number, ldx: number): void;

  /** Compute B = A * X to create RHS for testing */
  _dFillRHS(trans: number, nrhs: number, xPtr: number, ldx: number, APtr: number, BPtr: number): void;

  /**
   * Compute infinity-norm relative error: max|X - Xtrue| / max|Xtrue|
   * @returns The relative error value
   */
  _dinf_norm_error(nrhs: number, XPtr: number, ldx: number, xtruePtr: number, ldxtrue: number): number;

  /** Query memory usage of L and U factors */
  _dQuerySpace(LPtr: number, UPtr: number, memUsagePtr: number): number;

  /** Calculate memory usage estimate */
  _dmemory_usage(nzlmax: number, nzumax: number, nzlumax: number, n: number): number;

  // ============================================================
  // Single Precision (s) Routines
  // ============================================================
  // Single precision versions have identical signatures to double precision.
  // Use HEAPF32 instead of HEAPF64 for array access.
  // Trade-off: Half the memory, potentially faster, but ~7 decimal digits precision.

  /** SGSSV - Simple driver for A*X=B (single precision). See dgssv for details. */
  _sgssv(optionsPtr: number, APtr: number, permCPtr: number, permRPtr: number, LPtr: number, UPtr: number, BPtr: number, statPtr: number, infoPtr: number): void;
  /** SGSTRF - LU factorization (single precision). See dgstrf for details. */
  _sgstrf(optionsPtr: number, APtr: number, relax: number, panelSize: number, etreePtr: number, workPtr: number, lwork: number, permCPtr: number, permRPtr: number, LPtr: number, UPtr: number, GluPtr: number, statPtr: number, infoPtr: number): void;
  /** SGSTRS - Triangular solve (single precision). See dgstrs for details. */
  _sgstrs(trans: number, LPtr: number, UPtr: number, permCPtr: number, permRPtr: number, BPtr: number, statPtr: number, infoPtr: number): void;
  /** SGSSVX - Expert driver (single precision). See dgssvx for details. */
  _sgssvx(optionsPtr: number, APtr: number, permCPtr: number, permRPtr: number, etreePtr: number, equedPtr: number, RPtr: number, CPtr: number, LPtr: number, UPtr: number, workPtr: number, lwork: number, BPtr: number, XPtr: number, rpgPtr: number, rcondPtr: number, ferrPtr: number, berrPtr: number, GluPtr: number, memUsagePtr: number, statPtr: number, infoPtr: number): void;
  /** SGSISX - ILU driver (single precision). See dgsisx for details. */
  _sgsisx(optionsPtr: number, APtr: number, permCPtr: number, permRPtr: number, etreePtr: number, equedPtr: number, RPtr: number, CPtr: number, LPtr: number, UPtr: number, workPtr: number, lwork: number, BPtr: number, XPtr: number, rpgPtr: number, rcondPtr: number, GluPtr: number, memUsagePtr: number, statPtr: number, infoPtr: number): void;
  /** SGSITRF - ILU factorization (single precision). See dgsitrf for details. */
  _sgsitrf(optionsPtr: number, APtr: number, relax: number, panelSize: number, etreePtr: number, workPtr: number, lwork: number, permCPtr: number, permRPtr: number, LPtr: number, UPtr: number, GluPtr: number, statPtr: number, infoPtr: number): void;
  /** SGSCON - Condition number (single precision). See dgscon for details. */
  _sgscon(normPtr: number, LPtr: number, UPtr: number, anorm: number, rcondPtr: number, statPtr: number, infoPtr: number): void;
  /** SGSEQU - Equilibration (single precision). See dgsequ for details. */
  _sgsequ(APtr: number, RPtr: number, CPtr: number, rowcndPtr: number, colcndPtr: number, amaxPtr: number, infoPtr: number): void;
  /** SGSRFS - Iterative refinement (single precision). See dgsrfs for details. */
  _sgsrfs(trans: number, APtr: number, LPtr: number, UPtr: number, permCPtr: number, permRPtr: number, equedPtr: number, RPtr: number, CPtr: number, BPtr: number, XPtr: number, ferrPtr: number, berrPtr: number, statPtr: number, infoPtr: number): void;
  /** SLAQGS - Apply equilibration (single precision). See dlaqgs for details. */
  _slaqgs(APtr: number, RPtr: number, CPtr: number, rowcnd: number, colcnd: number, amax: number, equedPtr: number): void;
  /** Create CSC matrix (single precision). Use dtype=SLU_S. */
  _sCreate_CompCol_Matrix(APtr: number, m: number, n: number, nnz: number, nzvalPtr: number, rowindPtr: number, colptrPtr: number, stype: number, dtype: number, mtype: number): void;
  /** Create CSR matrix (single precision). Use dtype=SLU_S. */
  _sCreate_CompRow_Matrix(APtr: number, m: number, n: number, nnz: number, nzvalPtr: number, colindPtr: number, rowptrPtr: number, stype: number, dtype: number, mtype: number): void;
  /** Create dense matrix (single precision). Use dtype=SLU_S. */
  _sCreate_Dense_Matrix(XPtr: number, m: number, n: number, xPtr: number, ldx: number, stype: number, dtype: number, mtype: number): void;
  /** Copy CSC matrix (single precision). */
  _sCopy_CompCol_Matrix(APtr: number, BPtr: number): void;
  /** Copy dense matrix (single precision). */
  _sCopy_Dense_Matrix(M: number, N: number, XPtr: number, ldx: number, YPtr: number, ldy: number): void;
  /** Convert CSR to CSC (single precision). */
  _sCompRow_to_CompCol(m: number, n: number, nnz: number, aPtr: number, colindPtr: number, rowptrPtr: number, atPtrPtr: number, rowindPtrPtr: number, colptrPtrPtr: number): void;
  /** Allocate sparse storage (single precision). */
  _sallocateA(n: number, nnz: number, aPtr: number, asub: number, xa: number): void;
  /** Print CSC matrix (single precision). */
  _sPrint_CompCol_Matrix(whatPtr: number, APtr: number): void;
  /** Print dense matrix (single precision). */
  _sPrint_Dense_Matrix(whatPtr: number, APtr: number): void;
  /** Generate test solution (single precision). */
  _sGenXtrue(n: number, nrhs: number, xPtr: number, ldx: number): void;
  /** Compute RHS for testing (single precision). */
  _sFillRHS(trans: number, nrhs: number, xPtr: number, ldx: number, APtr: number, BPtr: number): void;
  /** Infinity norm error (single precision). */
  _sinf_norm_error(nrhs: number, XPtr: number, ldx: number, xtruePtr: number, ldxtrue: number): number;
  /** Query memory usage (single precision). */
  _sQuerySpace(LPtr: number, UPtr: number, memUsagePtr: number): number;
  /** Calculate memory usage (single precision). */
  _smemory_usage(nzlmax: number, nzumax: number, nzlumax: number, n: number): number;

  // ============================================================
  // Single Precision Complex (c) Routines
  // ============================================================
  // Complex single precision: each element is 2 floats (real, imag) = 8 bytes.
  // Use HEAPF32 with pairs of values for complex numbers.
  // Use dtype=SLU_C when creating matrices.

  /** CGSSV - Simple driver for A*X=B (complex single). See dgssv for details. */
  _cgssv(optionsPtr: number, APtr: number, permCPtr: number, permRPtr: number, LPtr: number, UPtr: number, BPtr: number, statPtr: number, infoPtr: number): void;
  /** CGSTRF - LU factorization (complex single). See dgstrf for details. */
  _cgstrf(optionsPtr: number, APtr: number, relax: number, panelSize: number, etreePtr: number, workPtr: number, lwork: number, permCPtr: number, permRPtr: number, LPtr: number, UPtr: number, GluPtr: number, statPtr: number, infoPtr: number): void;
  /** CGSTRS - Triangular solve (complex single). trans=CONJ applies A^H. */
  _cgstrs(trans: number, LPtr: number, UPtr: number, permCPtr: number, permRPtr: number, BPtr: number, statPtr: number, infoPtr: number): void;
  /** CGSSVX - Expert driver (complex single). See dgssvx for details. */
  _cgssvx(optionsPtr: number, APtr: number, permCPtr: number, permRPtr: number, etreePtr: number, equedPtr: number, RPtr: number, CPtr: number, LPtr: number, UPtr: number, workPtr: number, lwork: number, BPtr: number, XPtr: number, rpgPtr: number, rcondPtr: number, ferrPtr: number, berrPtr: number, GluPtr: number, memUsagePtr: number, statPtr: number, infoPtr: number): void;
  /** CGSISX - ILU driver (complex single). See dgsisx for details. */
  _cgsisx(optionsPtr: number, APtr: number, permCPtr: number, permRPtr: number, etreePtr: number, equedPtr: number, RPtr: number, CPtr: number, LPtr: number, UPtr: number, workPtr: number, lwork: number, BPtr: number, XPtr: number, rpgPtr: number, rcondPtr: number, GluPtr: number, memUsagePtr: number, statPtr: number, infoPtr: number): void;
  /** CGSITRF - ILU factorization (complex single). See dgsitrf for details. */
  _cgsitrf(optionsPtr: number, APtr: number, relax: number, panelSize: number, etreePtr: number, workPtr: number, lwork: number, permCPtr: number, permRPtr: number, LPtr: number, UPtr: number, GluPtr: number, statPtr: number, infoPtr: number): void;
  /** CGSCON - Condition number (complex single). See dgscon for details. */
  _cgscon(normPtr: number, LPtr: number, UPtr: number, anorm: number, rcondPtr: number, statPtr: number, infoPtr: number): void;
  /** CGSEQU - Equilibration (complex single). Scale factors are real. */
  _cgsequ(APtr: number, RPtr: number, CPtr: number, rowcndPtr: number, colcndPtr: number, amaxPtr: number, infoPtr: number): void;
  /** CGSRFS - Iterative refinement (complex single). See dgsrfs for details. */
  _cgsrfs(trans: number, APtr: number, LPtr: number, UPtr: number, permCPtr: number, permRPtr: number, equedPtr: number, RPtr: number, CPtr: number, BPtr: number, XPtr: number, ferrPtr: number, berrPtr: number, statPtr: number, infoPtr: number): void;
  /** CLAQGS - Apply equilibration (complex single). See dlaqgs for details. */
  _claqgs(APtr: number, RPtr: number, CPtr: number, rowcnd: number, colcnd: number, amax: number, equedPtr: number): void;
  /** Create CSC matrix (complex single). Use dtype=SLU_C. */
  _cCreate_CompCol_Matrix(APtr: number, m: number, n: number, nnz: number, nzvalPtr: number, rowindPtr: number, colptrPtr: number, stype: number, dtype: number, mtype: number): void;
  /** Create CSR matrix (complex single). Use dtype=SLU_C. */
  _cCreate_CompRow_Matrix(APtr: number, m: number, n: number, nnz: number, nzvalPtr: number, colindPtr: number, rowptrPtr: number, stype: number, dtype: number, mtype: number): void;
  /** Create dense matrix (complex single). Use dtype=SLU_C. */
  _cCreate_Dense_Matrix(XPtr: number, m: number, n: number, xPtr: number, ldx: number, stype: number, dtype: number, mtype: number): void;
  /** Copy CSC matrix (complex single). */
  _cCopy_CompCol_Matrix(APtr: number, BPtr: number): void;
  /** Copy dense matrix (complex single). */
  _cCopy_Dense_Matrix(M: number, N: number, XPtr: number, ldx: number, YPtr: number, ldy: number): void;
  /** Convert CSR to CSC (complex single). */
  _cCompRow_to_CompCol(m: number, n: number, nnz: number, aPtr: number, colindPtr: number, rowptrPtr: number, atPtrPtr: number, rowindPtrPtr: number, colptrPtrPtr: number): void;
  /** Allocate sparse storage (complex single). */
  _callocateA(n: number, nnz: number, aPtr: number, asub: number, xa: number): void;
  /** Print CSC matrix (complex single). */
  _cPrint_CompCol_Matrix(whatPtr: number, APtr: number): void;
  /** Print dense matrix (complex single). */
  _cPrint_Dense_Matrix(whatPtr: number, APtr: number): void;
  /** Generate test solution (complex single). */
  _cGenXtrue(n: number, nrhs: number, xPtr: number, ldx: number): void;
  /** Compute RHS for testing (complex single). */
  _cFillRHS(trans: number, nrhs: number, xPtr: number, ldx: number, APtr: number, BPtr: number): void;
  /** Infinity norm error (complex single). */
  _cinf_norm_error(nrhs: number, XPtr: number, ldx: number, xtruePtr: number, ldxtrue: number): number;
  /** Query memory usage (complex single). */
  _cQuerySpace(LPtr: number, UPtr: number, memUsagePtr: number): number;
  /** Calculate memory usage (complex single). */
  _cmemory_usage(nzlmax: number, nzumax: number, nzlumax: number, n: number): number;

  // ============================================================
  // Double Precision Complex (z) Routines
  // ============================================================
  // Complex double precision: each element is 2 doubles (real, imag) = 16 bytes.
  // Use HEAPF64 with pairs of values for complex numbers.
  // Use dtype=SLU_Z when creating matrices.

  /** ZGSSV - Simple driver for A*X=B (complex double). See dgssv for details. */
  _zgssv(optionsPtr: number, APtr: number, permCPtr: number, permRPtr: number, LPtr: number, UPtr: number, BPtr: number, statPtr: number, infoPtr: number): void;
  /** ZGSTRF - LU factorization (complex double). See dgstrf for details. */
  _zgstrf(optionsPtr: number, APtr: number, relax: number, panelSize: number, etreePtr: number, workPtr: number, lwork: number, permCPtr: number, permRPtr: number, LPtr: number, UPtr: number, GluPtr: number, statPtr: number, infoPtr: number): void;
  /** ZGSTRS - Triangular solve (complex double). trans=CONJ applies A^H. */
  _zgstrs(trans: number, LPtr: number, UPtr: number, permCPtr: number, permRPtr: number, BPtr: number, statPtr: number, infoPtr: number): void;
  /** ZGSSVX - Expert driver (complex double). See dgssvx for details. */
  _zgssvx(optionsPtr: number, APtr: number, permCPtr: number, permRPtr: number, etreePtr: number, equedPtr: number, RPtr: number, CPtr: number, LPtr: number, UPtr: number, workPtr: number, lwork: number, BPtr: number, XPtr: number, rpgPtr: number, rcondPtr: number, ferrPtr: number, berrPtr: number, GluPtr: number, memUsagePtr: number, statPtr: number, infoPtr: number): void;
  /** ZGSISX - ILU driver (complex double). See dgsisx for details. */
  _zgsisx(optionsPtr: number, APtr: number, permCPtr: number, permRPtr: number, etreePtr: number, equedPtr: number, RPtr: number, CPtr: number, LPtr: number, UPtr: number, workPtr: number, lwork: number, BPtr: number, XPtr: number, rpgPtr: number, rcondPtr: number, GluPtr: number, memUsagePtr: number, statPtr: number, infoPtr: number): void;
  /** ZGSITRF - ILU factorization (complex double). See dgsitrf for details. */
  _zgsitrf(optionsPtr: number, APtr: number, relax: number, panelSize: number, etreePtr: number, workPtr: number, lwork: number, permCPtr: number, permRPtr: number, LPtr: number, UPtr: number, GluPtr: number, statPtr: number, infoPtr: number): void;
  /** ZGSCON - Condition number (complex double). See dgscon for details. */
  _zgscon(normPtr: number, LPtr: number, UPtr: number, anorm: number, rcondPtr: number, statPtr: number, infoPtr: number): void;
  /** ZGSEQU - Equilibration (complex double). Scale factors are real. */
  _zgsequ(APtr: number, RPtr: number, CPtr: number, rowcndPtr: number, colcndPtr: number, amaxPtr: number, infoPtr: number): void;
  /** ZGSRFS - Iterative refinement (complex double). See dgsrfs for details. */
  _zgsrfs(trans: number, APtr: number, LPtr: number, UPtr: number, permCPtr: number, permRPtr: number, equedPtr: number, RPtr: number, CPtr: number, BPtr: number, XPtr: number, ferrPtr: number, berrPtr: number, statPtr: number, infoPtr: number): void;
  /** ZLAQGS - Apply equilibration (complex double). See dlaqgs for details. */
  _zlaqgs(APtr: number, RPtr: number, CPtr: number, rowcnd: number, colcnd: number, amax: number, equedPtr: number): void;
  /** Create CSC matrix (complex double). Use dtype=SLU_Z. */
  _zCreate_CompCol_Matrix(APtr: number, m: number, n: number, nnz: number, nzvalPtr: number, rowindPtr: number, colptrPtr: number, stype: number, dtype: number, mtype: number): void;
  /** Create CSR matrix (complex double). Use dtype=SLU_Z. */
  _zCreate_CompRow_Matrix(APtr: number, m: number, n: number, nnz: number, nzvalPtr: number, colindPtr: number, rowptrPtr: number, stype: number, dtype: number, mtype: number): void;
  /** Create dense matrix (complex double). Use dtype=SLU_Z. */
  _zCreate_Dense_Matrix(XPtr: number, m: number, n: number, xPtr: number, ldx: number, stype: number, dtype: number, mtype: number): void;
  /** Copy CSC matrix (complex double). */
  _zCopy_CompCol_Matrix(APtr: number, BPtr: number): void;
  /** Copy dense matrix (complex double). */
  _zCopy_Dense_Matrix(M: number, N: number, XPtr: number, ldx: number, YPtr: number, ldy: number): void;
  /** Convert CSR to CSC (complex double). */
  _zCompRow_to_CompCol(m: number, n: number, nnz: number, aPtr: number, colindPtr: number, rowptrPtr: number, atPtrPtr: number, rowindPtrPtr: number, colptrPtrPtr: number): void;
  /** Allocate sparse storage (complex double). */
  _zallocateA(n: number, nnz: number, aPtr: number, asub: number, xa: number): void;
  /** Print CSC matrix (complex double). */
  _zPrint_CompCol_Matrix(whatPtr: number, APtr: number): void;
  /** Print dense matrix (complex double). */
  _zPrint_Dense_Matrix(whatPtr: number, APtr: number): void;
  /** Generate test solution (complex double). */
  _zGenXtrue(n: number, nrhs: number, xPtr: number, ldx: number): void;
  /** Compute RHS for testing (complex double). */
  _zFillRHS(trans: number, nrhs: number, xPtr: number, ldx: number, APtr: number, BPtr: number): void;
  /** Infinity norm error (complex double). */
  _zinf_norm_error(nrhs: number, XPtr: number, ldx: number, xtruePtr: number, ldxtrue: number): number;
  /** Query memory usage (complex double). */
  _zQuerySpace(LPtr: number, UPtr: number, memUsagePtr: number): number;
  /** Calculate memory usage (complex double). */
  _zmemory_usage(nzlmax: number, nzumax: number, nzlumax: number, n: number): number;

  // ============================================================
  // Common Utility Functions
  // ============================================================

  /**
   * Set default SuperLU options
   *
   * Initializes a superlu_options_t structure with recommended defaults:
   * - Fact = DOFACT (perform factorization)
   * - Equil = YES (equilibrate if needed)
   * - ColPerm = COLAMD (column approximate minimum degree)
   * - Trans = NOTRANS (no transpose)
   * - IterRefine = NOREFINE (no iterative refinement)
   * - PrintStat = YES (print statistics)
   * - etc.
   *
   * @param optionsPtr - Pointer to superlu_options_t structure to initialize
   */
  _set_default_options(optionsPtr: number): void;

  /**
   * Initialize statistics structure
   *
   * Must be called before using any SuperLU solve/factor routine.
   * The structure tracks timing and operation counts.
   *
   * @param statPtr - Pointer to SuperLUStat_t structure
   */
  _StatInit(statPtr: number): void;

  /**
   * Free statistics structure
   *
   * Releases memory allocated by StatInit.
   *
   * @param statPtr - Pointer to SuperLUStat_t structure
   */
  _StatFree(statPtr: number): void;

  /**
   * Destroy compressed column matrix
   *
   * Frees the storage arrays (nzval, rowind, colptr) and the Store structure.
   * Call this for CSC matrices after use.
   *
   * @param APtr - Pointer to SuperMatrix in CSC format
   */
  _Destroy_CompCol_Matrix(APtr: number): void;

  /**
   * Destroy compressed row matrix
   *
   * Frees the storage arrays (nzval, colind, rowptr) and the Store structure.
   * Call this for CSR matrices after use.
   *
   * @param APtr - Pointer to SuperMatrix in CSR format
   */
  _Destroy_CompRow_Matrix(APtr: number): void;

  /**
   * Destroy dense matrix
   *
   * Frees the dense storage array and the Store structure.
   * Call this for dense matrices (typically RHS/solution) after use.
   *
   * @param APtr - Pointer to SuperMatrix in dense format
   */
  _Destroy_Dense_Matrix(APtr: number): void;

  /**
   * Destroy supernode matrix
   *
   * Frees the supernode storage used by L and U factors.
   * Call this for L/U matrices after use.
   *
   * @param APtr - Pointer to SuperMatrix in supernode format
   */
  _Destroy_SuperNode_Matrix(APtr: number): void;

  /**
   * Destroy only the Store component of a SuperMatrix
   *
   * Frees the Store structure but not the data arrays.
   * Use when data arrays were allocated separately.
   *
   * @param APtr - Pointer to SuperMatrix
   */
  _Destroy_SuperMatrix_Store(APtr: number): void;

  /**
   * Compute column permutation for fill reduction
   *
   * Computes a column permutation to reduce fill-in during factorization.
   *
   * @param ispec - Permutation algorithm:
   *                0 = NATURAL (no permutation)
   *                1 = MMD_ATA (minimum degree on A^T*A)
   *                2 = MMD_AT_PLUS_A (minimum degree on A^T+A)
   *                3 = COLAMD (column approximate minimum degree)
   * @param APtr - Pointer to SuperMatrix A
   * @param permCPtr - Pointer to output permutation array (size n)
   */
  _get_perm_c(ispec: number, APtr: number, permCPtr: number): void;

  /**
   * Symbolic factorization / preordering
   *
   * Performs symbolic factorization to determine the structure of L and U.
   * Computes the elimination tree and applies column permutation.
   *
   * @param optionsPtr - Pointer to options structure
   * @param APtr - Pointer to original SuperMatrix A
   * @param permCPtr - Pointer to column permutation
   * @param etreePtr - Pointer to output elimination tree (size n)
   * @param ACPtr - Pointer to output permuted matrix A*Pc
   */
  _sp_preorder(optionsPtr: number, APtr: number, permCPtr: number, etreePtr: number, ACPtr: number): void;

  /**
   * Get SuperLU environment parameter
   *
   * Returns internal tuning parameters.
   *
   * @param ispec - Parameter to query:
   *                1 = panel size (default: platform-dependent)
   *                2 = relaxation parameter for supernode amalgamation
   *                3 = maximum supernode size
   * @returns The parameter value
   */
  _sp_ienv(ispec: number): number;

  /**
   * Report input error
   *
   * Prints an error message for illegal input arguments.
   *
   * @param srname - Pointer to routine name string
   * @param info - Error code (usually negative of bad argument position)
   */
  _input_error(srname: number, info: number): void;

  /**
   * Abort with error message
   *
   * Prints error message and terminates (throws in WASM context).
   *
   * @param msgPtr - Pointer to error message string
   */
  _superlu_abort_and_exit(msgPtr: number): void;

  // ============================================================
  // Memory Allocation Utilities
  // ============================================================
  // SuperLU provides typed allocation functions for convenience.
  // These handle element size automatically.

  /**
   * Allocate array of n integers (uninitialized)
   * @param n - Number of integers
   * @returns Pointer to allocated array
   */
  _intMalloc(n: number): number;

  /**
   * Allocate array of n integers (zero-initialized)
   * @param n - Number of integers
   * @returns Pointer to allocated array
   */
  _intCalloc(n: number): number;

  /**
   * Allocate array of n floats (uninitialized)
   * @param n - Number of floats
   * @returns Pointer to allocated array
   */
  _floatMalloc(n: number): number;

  /**
   * Allocate array of n floats (zero-initialized)
   * @param n - Number of floats
   * @returns Pointer to allocated array
   */
  _floatCalloc(n: number): number;

  /**
   * Allocate array of n doubles (uninitialized)
   * @param n - Number of doubles
   * @returns Pointer to allocated array
   */
  _doubleMalloc(n: number): number;

  /**
   * Allocate array of n doubles (zero-initialized)
   * @param n - Number of doubles
   * @returns Pointer to allocated array
   */
  _doubleCalloc(n: number): number;

  /**
   * Allocate array of n single-precision complex numbers (uninitialized)
   * Each complex number is 2 floats (8 bytes total).
   * @param n - Number of complex numbers
   * @returns Pointer to allocated array
   */
  _singlecomplexMalloc(n: number): number;

  /**
   * Allocate array of n single-precision complex numbers (zero-initialized)
   * @param n - Number of complex numbers
   * @returns Pointer to allocated array
   */
  _singlecomplexCalloc(n: number): number;

  /**
   * Allocate array of n double-precision complex numbers (uninitialized)
   * Each complex number is 2 doubles (16 bytes total).
   * @param n - Number of complex numbers
   * @returns Pointer to allocated array
   */
  _doublecomplexMalloc(n: number): number;

  /**
   * Allocate array of n double-precision complex numbers (zero-initialized)
   * @param n - Number of complex numbers
   * @returns Pointer to allocated array
   */
  _doublecomplexCalloc(n: number): number;

  /**
   * SuperLU's internal allocation (with memory tracking)
   * @param size - Number of bytes to allocate
   * @returns Pointer to allocated memory
   */
  _superlu_malloc(size: number): number;

  /**
   * Free memory allocated by superlu_malloc
   * @param ptr - Pointer to free
   */
  _superlu_free(ptr: number): void;

  // ============================================================
  // Standard C Memory Functions
  // ============================================================

  /**
   * Allocate memory in the WASM heap
   * @param size - Number of bytes to allocate
   * @returns Pointer (byte offset) to allocated memory, or 0 on failure
   */
  _malloc(size: number): number;

  /**
   * Free memory allocated by _malloc
   * @param ptr - Pointer to free
   */
  _free(ptr: number): void;

  // ============================================================
  // Emscripten Runtime Methods
  // ============================================================

  /**
   * Read a value from WASM memory
   * @param ptr - Memory address
   * @param type - Type: 'i8', 'i16', 'i32', 'i64', 'float', 'double'
   * @returns The value at that address
   */
  getValue(ptr: number, type: string): number;

  /**
   * Write a value to WASM memory
   * @param ptr - Memory address
   * @param value - Value to write
   * @param type - Type: 'i8', 'i16', 'i32', 'i64', 'float', 'double'
   */
  setValue(ptr: number, value: number, type: string): void;

  // ============================================================
  // Heap Array Views
  // ============================================================

  /** Float64Array view of WASM memory. Index = ptr/8 for doubles. */
  HEAPF64: Float64Array;

  /** Float32Array view of WASM memory. Index = ptr/4 for floats. */
  HEAPF32: Float32Array;

  /** Int32Array view of WASM memory. Index = ptr/4 for 32-bit ints. */
  HEAP32: Int32Array;

  /** Int8Array view of WASM memory. Index = ptr for bytes. */
  HEAP8: Int8Array;

  /** Uint8Array view of WASM memory. Index = ptr for unsigned bytes. */
  HEAPU8: Uint8Array;

  /** Uint32Array view of WASM memory. Index = ptr/4 for 32-bit unsigned ints. */
  HEAPU32: Uint32Array;

  // ============================================================
  // Function Table Management
  // ============================================================

  /**
   * Add a JavaScript function to the WASM function table
   *
   * Used for callbacks (e.g., custom comparison functions).
   *
   * @param func - JavaScript function to add
   * @param signature - Signature string: 'v'=void, 'i'=i32, 'j'=i64, 'f'=f32, 'd'=f64
   *                    First char is return type, rest are parameter types.
   *                    Example: 'iiii' = returns int, takes 3 ints
   * @returns Function pointer (index in table)
   */
  addFunction(func: Function, signature: string): number;

  /**
   * Remove a function from the WASM function table
   * @param ptr - Function pointer returned by addFunction
   */
  removeFunction(ptr: number): void;
}

/**
 * Options for loading the SuperLU WASM module
 */
export interface SuperLUModuleOptions {
  /**
   * Custom function to locate WASM files
   *
   * @param path - Filename being requested (e.g., 'superlu.wasm')
   * @param scriptDirectory - Directory where the JS loader is located
   * @returns Full URL or path to the file
   *
   * @example
   * const module = await createSuperLUModule({
   *   locateFile: (path) => `https://cdn.example.com/wasm/${path}`
   * });
   */
  locateFile?: (path: string, scriptDirectory: string) => string;
}

/**
 * Factory function to create a SuperLU WASM module instance
 *
 * @param options - Optional configuration
 * @returns Promise resolving to initialized SuperLU module
 *
 * @example
 * import createSuperLUModule from 'superluwasm';
 *
 * const slu = await createSuperLUModule();
 *
 * // Solve a sparse system Ax = b
 * // 1. Allocate and populate CSC matrix A
 * // 2. Create dense matrix B for RHS
 * // 3. Call slu._dgssv(...) to solve
 * // 4. Read solution from B
 * // 5. Destroy matrices with Destroy_*_Matrix
 */
export type SuperLUModuleFactory = (
  options?: SuperLUModuleOptions
) => Promise<SuperLUModule>;
