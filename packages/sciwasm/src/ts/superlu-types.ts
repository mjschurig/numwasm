/**
 * SuperLU WASM module type definitions
 *
 * SuperLU is a library for the direct solution of large, sparse,
 * nonsymmetric systems of linear equations.
 *
 * This module supports all four precision variants:
 * - s: single precision real (float)
 * - d: double precision real (double)
 * - c: single precision complex
 * - z: double precision complex
 */

/**
 * SuperLU options structure constants
 */
export const enum SuperLUFactOption {
  DOFACT = 0,
  SamePattern = 1,
  SamePattern_SameRowPerm = 2,
  FACTORED = 3,
}

export const enum SuperLURowPerm {
  NOROWPERM = 0,
  LargeDiag_MC64 = 1,
  LargeDiag_AWPM = 2,
  MY_PERMR = 3,
}

export const enum SuperLUColPerm {
  NATURAL = 0,
  MMD_ATA = 1,
  MMD_AT_PLUS_A = 2,
  COLAMD = 3,
  METIS_AT_PLUS_A = 4,
  PARMETIS = 5,
  ZOLTAN = 6,
  MY_PERMC = 7,
}

export const enum SuperLUTrans {
  NOTRANS = 0,
  TRANS = 1,
  CONJ = 2,
}

export const enum SuperLUILUDrop {
  DROP_BASIC = 0,
  DROP_PROWS = 1,
  DROP_COLUMN = 2,
  DROP_AREA = 3,
  DROP_DYNAMIC = 4,
  DROP_INTERP = 5,
}

export const enum SuperLUMiluT {
  SILU = 0,
  SMILU_1 = 1,
  SMILU_2 = 2,
  SMILU_3 = 3,
}

export const enum SuperLUNorm {
  ONE_NORM = 0,
  TWO_NORM = 1,
  INF_NORM = 2,
}

/**
 * SuperMatrix storage types
 */
export const enum SuperLUStype {
  SLU_NC = 0,    // Column-wise, no supernode
  SLU_NCP = 1,   // Column-wise, column permutation
  SLU_NR = 2,    // Row-wise, no supernode
  SLU_SC = 3,    // Column-wise, supernode
  SLU_SCP = 4,   // Supernode, column permutation
  SLU_SR = 5,    // Row-wise, supernode
  SLU_DN = 6,    // Dense (Fortran style)
  SLU_NR_loc = 7, // Distributed CSR
}

export const enum SuperLUDtype {
  SLU_S = 0,  // Single precision real
  SLU_D = 1,  // Double precision real
  SLU_C = 2,  // Single precision complex
  SLU_Z = 3,  // Double precision complex
}

export const enum SuperLUMtype {
  SLU_GE = 0,   // General
  SLU_TRLU = 1, // Lower triangular, unit diagonal
  SLU_TRLL = 2, // Lower triangular
  SLU_TRU = 3,  // Upper triangular, unit diagonal
  SLU_TRL = 4,  // Upper triangular
  SLU_SYL = 5,  // Symmetric, lower triangle stored
  SLU_SYU = 6,  // Symmetric, upper triangle stored
  SLU_HEL = 7,  // Hermitian, lower triangle stored
  SLU_HEU = 8,  // Hermitian, upper triangle stored
}

export interface SuperLUModule {
  // ============================================================
  // Double precision (d) routines - most commonly used
  // ============================================================

  /** Simple driver for solving Ax=b with double precision */
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

  /** LU factorization for double precision */
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

  /** Triangular solve for double precision */
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

  /** Expert driver for double precision */
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

  /** ILU simple driver */
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

  /** ILU factorization */
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

  /** Condition number estimation */
  _dgscon(
    normPtr: number,
    LPtr: number,
    UPtr: number,
    anorm: number,
    rcondPtr: number,
    statPtr: number,
    infoPtr: number
  ): void;

  /** Row/column equilibration */
  _dgsequ(
    APtr: number,
    RPtr: number,
    CPtr: number,
    rowcndPtr: number,
    colcndPtr: number,
    amaxPtr: number,
    infoPtr: number
  ): void;

  /** Iterative refinement */
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

  /** Apply equilibration scaling */
  _dlaqgs(
    APtr: number,
    RPtr: number,
    CPtr: number,
    rowcnd: number,
    colcnd: number,
    amax: number,
    equedPtr: number
  ): void;

  /** Create compressed column matrix */
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

  /** Create compressed row matrix */
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

  /** Create dense matrix */
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

  /** Copy compressed column matrix */
  _dCopy_CompCol_Matrix(APtr: number, BPtr: number): void;

  /** Copy dense matrix */
  _dCopy_Dense_Matrix(M: number, N: number, XPtr: number, ldx: number, YPtr: number, ldy: number): void;

  /** Convert CSR to CSC */
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

  /** Allocate matrix storage */
  _dallocateA(n: number, nnz: number, aPtr: number, asub: number, xa: number): void;

  /** Print compressed column matrix */
  _dPrint_CompCol_Matrix(whatPtr: number, APtr: number): void;

  /** Print dense matrix */
  _dPrint_Dense_Matrix(whatPtr: number, APtr: number): void;

  /** Generate test RHS */
  _dGenXtrue(n: number, nrhs: number, xPtr: number, ldx: number): void;

  /** Fill RHS for testing */
  _dFillRHS(trans: number, nrhs: number, xPtr: number, ldx: number, APtr: number, BPtr: number): void;

  /** Compute infinity norm error */
  _dinf_norm_error(nrhs: number, XPtr: number, ldx: number, xtruePtr: number, ldxtrue: number): number;

  /** Query memory usage */
  _dQuerySpace(LPtr: number, UPtr: number, memUsagePtr: number): number;

  /** Calculate memory usage */
  _dmemory_usage(nzlmax: number, nzumax: number, nzlumax: number, n: number): number;

  // ============================================================
  // Single precision (s) routines
  // ============================================================

  _sgssv(optionsPtr: number, APtr: number, permCPtr: number, permRPtr: number, LPtr: number, UPtr: number, BPtr: number, statPtr: number, infoPtr: number): void;
  _sgstrf(optionsPtr: number, APtr: number, relax: number, panelSize: number, etreePtr: number, workPtr: number, lwork: number, permCPtr: number, permRPtr: number, LPtr: number, UPtr: number, GluPtr: number, statPtr: number, infoPtr: number): void;
  _sgstrs(trans: number, LPtr: number, UPtr: number, permCPtr: number, permRPtr: number, BPtr: number, statPtr: number, infoPtr: number): void;
  _sgssvx(optionsPtr: number, APtr: number, permCPtr: number, permRPtr: number, etreePtr: number, equedPtr: number, RPtr: number, CPtr: number, LPtr: number, UPtr: number, workPtr: number, lwork: number, BPtr: number, XPtr: number, rpgPtr: number, rcondPtr: number, ferrPtr: number, berrPtr: number, GluPtr: number, memUsagePtr: number, statPtr: number, infoPtr: number): void;
  _sgsisx(optionsPtr: number, APtr: number, permCPtr: number, permRPtr: number, etreePtr: number, equedPtr: number, RPtr: number, CPtr: number, LPtr: number, UPtr: number, workPtr: number, lwork: number, BPtr: number, XPtr: number, rpgPtr: number, rcondPtr: number, GluPtr: number, memUsagePtr: number, statPtr: number, infoPtr: number): void;
  _sgsitrf(optionsPtr: number, APtr: number, relax: number, panelSize: number, etreePtr: number, workPtr: number, lwork: number, permCPtr: number, permRPtr: number, LPtr: number, UPtr: number, GluPtr: number, statPtr: number, infoPtr: number): void;
  _sgscon(normPtr: number, LPtr: number, UPtr: number, anorm: number, rcondPtr: number, statPtr: number, infoPtr: number): void;
  _sgsequ(APtr: number, RPtr: number, CPtr: number, rowcndPtr: number, colcndPtr: number, amaxPtr: number, infoPtr: number): void;
  _sgsrfs(trans: number, APtr: number, LPtr: number, UPtr: number, permCPtr: number, permRPtr: number, equedPtr: number, RPtr: number, CPtr: number, BPtr: number, XPtr: number, ferrPtr: number, berrPtr: number, statPtr: number, infoPtr: number): void;
  _slaqgs(APtr: number, RPtr: number, CPtr: number, rowcnd: number, colcnd: number, amax: number, equedPtr: number): void;
  _sCreate_CompCol_Matrix(APtr: number, m: number, n: number, nnz: number, nzvalPtr: number, rowindPtr: number, colptrPtr: number, stype: number, dtype: number, mtype: number): void;
  _sCreate_CompRow_Matrix(APtr: number, m: number, n: number, nnz: number, nzvalPtr: number, colindPtr: number, rowptrPtr: number, stype: number, dtype: number, mtype: number): void;
  _sCreate_Dense_Matrix(XPtr: number, m: number, n: number, xPtr: number, ldx: number, stype: number, dtype: number, mtype: number): void;
  _sCopy_CompCol_Matrix(APtr: number, BPtr: number): void;
  _sCopy_Dense_Matrix(M: number, N: number, XPtr: number, ldx: number, YPtr: number, ldy: number): void;
  _sCompRow_to_CompCol(m: number, n: number, nnz: number, aPtr: number, colindPtr: number, rowptrPtr: number, atPtrPtr: number, rowindPtrPtr: number, colptrPtrPtr: number): void;
  _sallocateA(n: number, nnz: number, aPtr: number, asub: number, xa: number): void;
  _sPrint_CompCol_Matrix(whatPtr: number, APtr: number): void;
  _sPrint_Dense_Matrix(whatPtr: number, APtr: number): void;
  _sGenXtrue(n: number, nrhs: number, xPtr: number, ldx: number): void;
  _sFillRHS(trans: number, nrhs: number, xPtr: number, ldx: number, APtr: number, BPtr: number): void;
  _sinf_norm_error(nrhs: number, XPtr: number, ldx: number, xtruePtr: number, ldxtrue: number): number;
  _sQuerySpace(LPtr: number, UPtr: number, memUsagePtr: number): number;
  _smemory_usage(nzlmax: number, nzumax: number, nzlumax: number, n: number): number;

  // ============================================================
  // Single precision complex (c) routines
  // ============================================================

  _cgssv(optionsPtr: number, APtr: number, permCPtr: number, permRPtr: number, LPtr: number, UPtr: number, BPtr: number, statPtr: number, infoPtr: number): void;
  _cgstrf(optionsPtr: number, APtr: number, relax: number, panelSize: number, etreePtr: number, workPtr: number, lwork: number, permCPtr: number, permRPtr: number, LPtr: number, UPtr: number, GluPtr: number, statPtr: number, infoPtr: number): void;
  _cgstrs(trans: number, LPtr: number, UPtr: number, permCPtr: number, permRPtr: number, BPtr: number, statPtr: number, infoPtr: number): void;
  _cgssvx(optionsPtr: number, APtr: number, permCPtr: number, permRPtr: number, etreePtr: number, equedPtr: number, RPtr: number, CPtr: number, LPtr: number, UPtr: number, workPtr: number, lwork: number, BPtr: number, XPtr: number, rpgPtr: number, rcondPtr: number, ferrPtr: number, berrPtr: number, GluPtr: number, memUsagePtr: number, statPtr: number, infoPtr: number): void;
  _cgsisx(optionsPtr: number, APtr: number, permCPtr: number, permRPtr: number, etreePtr: number, equedPtr: number, RPtr: number, CPtr: number, LPtr: number, UPtr: number, workPtr: number, lwork: number, BPtr: number, XPtr: number, rpgPtr: number, rcondPtr: number, GluPtr: number, memUsagePtr: number, statPtr: number, infoPtr: number): void;
  _cgsitrf(optionsPtr: number, APtr: number, relax: number, panelSize: number, etreePtr: number, workPtr: number, lwork: number, permCPtr: number, permRPtr: number, LPtr: number, UPtr: number, GluPtr: number, statPtr: number, infoPtr: number): void;
  _cgscon(normPtr: number, LPtr: number, UPtr: number, anorm: number, rcondPtr: number, statPtr: number, infoPtr: number): void;
  _cgsequ(APtr: number, RPtr: number, CPtr: number, rowcndPtr: number, colcndPtr: number, amaxPtr: number, infoPtr: number): void;
  _cgsrfs(trans: number, APtr: number, LPtr: number, UPtr: number, permCPtr: number, permRPtr: number, equedPtr: number, RPtr: number, CPtr: number, BPtr: number, XPtr: number, ferrPtr: number, berrPtr: number, statPtr: number, infoPtr: number): void;
  _claqgs(APtr: number, RPtr: number, CPtr: number, rowcnd: number, colcnd: number, amax: number, equedPtr: number): void;
  _cCreate_CompCol_Matrix(APtr: number, m: number, n: number, nnz: number, nzvalPtr: number, rowindPtr: number, colptrPtr: number, stype: number, dtype: number, mtype: number): void;
  _cCreate_CompRow_Matrix(APtr: number, m: number, n: number, nnz: number, nzvalPtr: number, colindPtr: number, rowptrPtr: number, stype: number, dtype: number, mtype: number): void;
  _cCreate_Dense_Matrix(XPtr: number, m: number, n: number, xPtr: number, ldx: number, stype: number, dtype: number, mtype: number): void;
  _cCopy_CompCol_Matrix(APtr: number, BPtr: number): void;
  _cCopy_Dense_Matrix(M: number, N: number, XPtr: number, ldx: number, YPtr: number, ldy: number): void;
  _cCompRow_to_CompCol(m: number, n: number, nnz: number, aPtr: number, colindPtr: number, rowptrPtr: number, atPtrPtr: number, rowindPtrPtr: number, colptrPtrPtr: number): void;
  _callocateA(n: number, nnz: number, aPtr: number, asub: number, xa: number): void;
  _cPrint_CompCol_Matrix(whatPtr: number, APtr: number): void;
  _cPrint_Dense_Matrix(whatPtr: number, APtr: number): void;
  _cGenXtrue(n: number, nrhs: number, xPtr: number, ldx: number): void;
  _cFillRHS(trans: number, nrhs: number, xPtr: number, ldx: number, APtr: number, BPtr: number): void;
  _cinf_norm_error(nrhs: number, XPtr: number, ldx: number, xtruePtr: number, ldxtrue: number): number;
  _cQuerySpace(LPtr: number, UPtr: number, memUsagePtr: number): number;
  _cmemory_usage(nzlmax: number, nzumax: number, nzlumax: number, n: number): number;

  // ============================================================
  // Double precision complex (z) routines
  // ============================================================

  _zgssv(optionsPtr: number, APtr: number, permCPtr: number, permRPtr: number, LPtr: number, UPtr: number, BPtr: number, statPtr: number, infoPtr: number): void;
  _zgstrf(optionsPtr: number, APtr: number, relax: number, panelSize: number, etreePtr: number, workPtr: number, lwork: number, permCPtr: number, permRPtr: number, LPtr: number, UPtr: number, GluPtr: number, statPtr: number, infoPtr: number): void;
  _zgstrs(trans: number, LPtr: number, UPtr: number, permCPtr: number, permRPtr: number, BPtr: number, statPtr: number, infoPtr: number): void;
  _zgssvx(optionsPtr: number, APtr: number, permCPtr: number, permRPtr: number, etreePtr: number, equedPtr: number, RPtr: number, CPtr: number, LPtr: number, UPtr: number, workPtr: number, lwork: number, BPtr: number, XPtr: number, rpgPtr: number, rcondPtr: number, ferrPtr: number, berrPtr: number, GluPtr: number, memUsagePtr: number, statPtr: number, infoPtr: number): void;
  _zgsisx(optionsPtr: number, APtr: number, permCPtr: number, permRPtr: number, etreePtr: number, equedPtr: number, RPtr: number, CPtr: number, LPtr: number, UPtr: number, workPtr: number, lwork: number, BPtr: number, XPtr: number, rpgPtr: number, rcondPtr: number, GluPtr: number, memUsagePtr: number, statPtr: number, infoPtr: number): void;
  _zgsitrf(optionsPtr: number, APtr: number, relax: number, panelSize: number, etreePtr: number, workPtr: number, lwork: number, permCPtr: number, permRPtr: number, LPtr: number, UPtr: number, GluPtr: number, statPtr: number, infoPtr: number): void;
  _zgscon(normPtr: number, LPtr: number, UPtr: number, anorm: number, rcondPtr: number, statPtr: number, infoPtr: number): void;
  _zgsequ(APtr: number, RPtr: number, CPtr: number, rowcndPtr: number, colcndPtr: number, amaxPtr: number, infoPtr: number): void;
  _zgsrfs(trans: number, APtr: number, LPtr: number, UPtr: number, permCPtr: number, permRPtr: number, equedPtr: number, RPtr: number, CPtr: number, BPtr: number, XPtr: number, ferrPtr: number, berrPtr: number, statPtr: number, infoPtr: number): void;
  _zlaqgs(APtr: number, RPtr: number, CPtr: number, rowcnd: number, colcnd: number, amax: number, equedPtr: number): void;
  _zCreate_CompCol_Matrix(APtr: number, m: number, n: number, nnz: number, nzvalPtr: number, rowindPtr: number, colptrPtr: number, stype: number, dtype: number, mtype: number): void;
  _zCreate_CompRow_Matrix(APtr: number, m: number, n: number, nnz: number, nzvalPtr: number, colindPtr: number, rowptrPtr: number, stype: number, dtype: number, mtype: number): void;
  _zCreate_Dense_Matrix(XPtr: number, m: number, n: number, xPtr: number, ldx: number, stype: number, dtype: number, mtype: number): void;
  _zCopy_CompCol_Matrix(APtr: number, BPtr: number): void;
  _zCopy_Dense_Matrix(M: number, N: number, XPtr: number, ldx: number, YPtr: number, ldy: number): void;
  _zCompRow_to_CompCol(m: number, n: number, nnz: number, aPtr: number, colindPtr: number, rowptrPtr: number, atPtrPtr: number, rowindPtrPtr: number, colptrPtrPtr: number): void;
  _zallocateA(n: number, nnz: number, aPtr: number, asub: number, xa: number): void;
  _zPrint_CompCol_Matrix(whatPtr: number, APtr: number): void;
  _zPrint_Dense_Matrix(whatPtr: number, APtr: number): void;
  _zGenXtrue(n: number, nrhs: number, xPtr: number, ldx: number): void;
  _zFillRHS(trans: number, nrhs: number, xPtr: number, ldx: number, APtr: number, BPtr: number): void;
  _zinf_norm_error(nrhs: number, XPtr: number, ldx: number, xtruePtr: number, ldxtrue: number): number;
  _zQuerySpace(LPtr: number, UPtr: number, memUsagePtr: number): number;
  _zmemory_usage(nzlmax: number, nzumax: number, nzlumax: number, n: number): number;

  // ============================================================
  // Common utility functions
  // ============================================================

  /** Set default SuperLU options */
  _set_default_options(optionsPtr: number): void;

  /** Initialize statistics structure */
  _StatInit(statPtr: number): void;

  /** Free statistics structure */
  _StatFree(statPtr: number): void;

  /** Destroy compressed column matrix */
  _Destroy_CompCol_Matrix(APtr: number): void;

  /** Destroy compressed row matrix */
  _Destroy_CompRow_Matrix(APtr: number): void;

  /** Destroy dense matrix */
  _Destroy_Dense_Matrix(APtr: number): void;

  /** Destroy supernode matrix */
  _Destroy_SuperNode_Matrix(APtr: number): void;

  /** Destroy matrix store only */
  _Destroy_SuperMatrix_Store(APtr: number): void;

  /** Get column permutation */
  _get_perm_c(ispec: number, APtr: number, permCPtr: number): void;

  /** Symbolic factorization */
  _sp_preorder(optionsPtr: number, APtr: number, permCPtr: number, etreePtr: number, ACPtr: number): void;

  /** Get environment parameter */
  _sp_ienv(ispec: number): number;

  /** Report input error */
  _input_error(srname: number, info: number): void;

  /** Abort with error message */
  _superlu_abort_and_exit(msgPtr: number): void;

  // ============================================================
  // Memory allocation utilities
  // ============================================================

  _intMalloc(n: number): number;
  _intCalloc(n: number): number;
  _floatMalloc(n: number): number;
  _floatCalloc(n: number): number;
  _doubleMalloc(n: number): number;
  _doubleCalloc(n: number): number;
  _singlecomplexMalloc(n: number): number;
  _singlecomplexCalloc(n: number): number;
  _doublecomplexMalloc(n: number): number;
  _doublecomplexCalloc(n: number): number;
  _superlu_malloc(size: number): number;
  _superlu_free(ptr: number): void;

  // Standard C memory
  _malloc(size: number): number;
  _free(ptr: number): void;

  // Emscripten runtime methods
  getValue(ptr: number, type: string): number;
  setValue(ptr: number, value: number, type: string): void;

  // Heap views
  HEAPF64: Float64Array;
  HEAPF32: Float32Array;
  HEAP32: Int32Array;
  HEAP8: Int8Array;
  HEAPU8: Uint8Array;
  HEAPU32: Uint32Array;

  // Function table management
  addFunction(func: Function, signature: string): number;
  removeFunction(ptr: number): void;
}

export interface SuperLUModuleOptions {
  locateFile?: (path: string, scriptDirectory: string) => string;
}

export type SuperLUModuleFactory = (
  options?: SuperLUModuleOptions
) => Promise<SuperLUModule>;
