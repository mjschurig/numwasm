/**
 * WASM module type definitions for sciwasm
 */

export interface SciWasmModule {
  // Nelder-Mead optimization
  _nelder_mead_minimize(
    n: number,
    x0Ptr: number,
    xOutPtr: number,
    simOutPtr: number,
    fsimOutPtr: number,
    resultPtr: number,
    xatol: number,
    fatol: number,
    maxiter: number,
    maxfev: number,
    adaptive: number,
    hasBounds: number,
    lowerPtr: number,
    upperPtr: number,
    initialSimplexPtr: number,
    funcPtr: number
  ): number;

  // BFGS optimization
  _bfgs_minimize(
    n: number,
    x0Ptr: number,
    xOutPtr: number,
    hessInvOutPtr: number,
    jacOutPtr: number,
    resultPtr: number,
    gtol: number,
    maxiter: number,
    c1: number,
    c2: number,
    hasJac: number,
    eps: number,
    funcPtr: number,
    gradPtr: number
  ): number;

  // L-BFGS-B optimization (reverse communication via setulb)
  _setulb(
    n: number, m: number,
    xPtr: number, lPtr: number, uPtr: number, nbdPtr: number,
    fPtr: number, gPtr: number,
    factr: number, pgtol: number,
    waPtr: number, iwaPtr: number,
    taskPtr: number, lsavePtr: number,
    isavePtr: number, dsavePtr: number,
    maxls: number, lnTaskPtr: number
  ): void;

  // QUADPACK adaptive quadrature
  _wasm_dqagse(
    fcn: number, a: number, b: number,
    epsabs: number, epsrel: number, limit: number,
    result: number, abserr: number, neval: number, ier: number,
    alist: number, blist: number, rlist: number, elist: number,
    iord: number, last: number,
  ): void;

  _wasm_dqagie(
    fcn: number, bound: number, inf: number,
    epsabs: number, epsrel: number, limit: number,
    result: number, abserr: number, neval: number, ier: number,
    alist: number, blist: number, rlist: number, elist: number,
    iord: number, last: number,
  ): void;

  // Special functions (gamma)
  _wasm_gamma(x: number): number;
  _wasm_gammaln(x: number): number;
  _wasm_rgamma(x: number): number;

  // Special functions (combinatorial)
  _wasm_binom(n: number, k: number): number;
  _wasm_binom_exact(n: number, k: number): number;
  _wasm_poch(x: number, m: number): number;
  _wasm_perm_exact(n: number, k: number): number;

  // Sparse matrix operations (sparsetools)
  _sp_csr_matvec_f64(n_row: number, n_col: number, Ap: number, Aj: number, Ax: number, Xx: number, Yx: number): void;
  _sp_csr_matvecs_f64(n_row: number, n_col: number, n_vecs: number, Ap: number, Aj: number, Ax: number, Xx: number, Yx: number): void;
  _sp_csr_tocsc_f64(n_row: number, n_col: number, Ap: number, Aj: number, Ax: number, Bp: number, Bi: number, Bx: number): void;
  _sp_csr_todense_f64(n_row: number, n_col: number, Ap: number, Aj: number, Ax: number, Bx: number): void;
  _sp_csr_diagonal_f64(k: number, n_row: number, n_col: number, Ap: number, Aj: number, Ax: number, Yx: number): void;
  _sp_csr_sort_indices_f64(n_row: number, Ap: number, Aj: number, Ax: number): void;
  _sp_csr_has_sorted_indices(n_row: number, Ap: number, Aj: number): number;
  _sp_csr_has_canonical_format(n_row: number, Ap: number, Aj: number): number;
  _sp_csr_sum_duplicates_f64(n_row: number, n_col: number, Ap: number, Aj: number, Ax: number): void;
  _sp_csr_eliminate_zeros_f64(n_row: number, n_col: number, Ap: number, Aj: number, Ax: number): void;
  _sp_csr_matmat_maxnnz(n_row: number, n_col: number, Ap: number, Aj: number, Bp: number, Bj: number): number;
  _sp_csr_matmat_f64(n_row: number, n_col: number, Ap: number, Aj: number, Ax: number, Bp: number, Bj: number, Bx: number, Cp: number, Cj: number, Cx: number): void;
  _sp_csr_plus_csr_f64(n_row: number, n_col: number, Ap: number, Aj: number, Ax: number, Bp: number, Bj: number, Bx: number, Cp: number, Cj: number, Cx: number): void;
  _sp_csr_minus_csr_f64(n_row: number, n_col: number, Ap: number, Aj: number, Ax: number, Bp: number, Bj: number, Bx: number, Cp: number, Cj: number, Cx: number): void;
  _sp_csr_elmul_csr_f64(n_row: number, n_col: number, Ap: number, Aj: number, Ax: number, Bp: number, Bj: number, Bx: number, Cp: number, Cj: number, Cx: number): void;
  _sp_csr_eldiv_csr_f64(n_row: number, n_col: number, Ap: number, Aj: number, Ax: number, Bp: number, Bj: number, Bx: number, Cp: number, Cj: number, Cx: number): void;
  _sp_csr_scale_rows_f64(n_row: number, n_col: number, Ap: number, Aj: number, Ax: number, Xx: number): void;
  _sp_csr_scale_columns_f64(n_row: number, n_col: number, Ap: number, Aj: number, Ax: number, Xx: number): void;
  _sp_expandptr(n_row: number, Ap: number, Bi: number): void;
  _sp_csr_row_index_f64(n_row_idx: number, rows: number, Ap: number, Aj: number, Ax: number, Bj: number, Bx: number): void;
  _sp_csr_row_slice_f64(start: number, stop: number, step: number, Ap: number, Aj: number, Ax: number, Bj: number, Bx: number): void;
  _sp_csr_column_index1(n_idx: number, col_idxs: number, n_row: number, n_col: number, Ap: number, Aj: number, col_offsets: number, Bp: number): void;
  _sp_csr_column_index2_f64(col_order: number, col_offsets: number, nnz: number, Aj: number, Ax: number, Bj: number, Bx: number): void;
  _sp_csr_sample_offsets(n_row: number, n_col: number, Ap: number, Aj: number, n_samples: number, Bi: number, Bj: number, Bp: number): number;
  _sp_csr_sample_values_f64(n_row: number, n_col: number, Ap: number, Aj: number, Ax: number, n_samples: number, Bi: number, Bj: number, Bx: number): void;
  _sp_get_csr_submatrix_nnz(n_row: number, n_col: number, Ap: number, Aj: number, ir0: number, ir1: number, ic0: number, ic1: number): number;
  _sp_get_csr_submatrix_f64(n_row: number, n_col: number, Ap: number, Aj: number, Ax: number, ir0: number, ir1: number, ic0: number, ic1: number, Bp: number, Bj: number, Bx: number): void;
  _sp_csc_matvec_f64(n_row: number, n_col: number, Ap: number, Ai: number, Ax: number, Xx: number, Yx: number): void;
  _sp_csc_matvecs_f64(n_row: number, n_col: number, n_vecs: number, Ap: number, Ai: number, Ax: number, Xx: number, Yx: number): void;
  _sp_coo_tocsr_f64(n_row: number, n_col: number, nnz: number, Ai: number, Aj: number, Ax: number, Bp: number, Bj: number, Bx: number): void;
  _sp_coo_todense_f64(n_row: number, n_col: number, nnz: number, Ai: number, Aj: number, Ax: number, Bx: number, fortran: number): void;
  _sp_coo_matvec_f64(nnz: number, Ai: number, Aj: number, Ax: number, Xx: number, Yx: number): void;

  // DIA sparse matrix operations
  _sp_dia_matvec_f64(n_row: number, n_col: number, n_diags: number, L: number,
    offsets: number, diags: number, Xx: number, Yx: number): void;
  _sp_dia_tocsr_f64(n_rows: number, n_cols: number, n_diags: number, L: number,
    offsets: number, data: number, order: number,
    csr_data: number, indices: number, indptr: number): number;

  // BSR (Block Sparse Row) matrix operations
  _sp_bsr_matvec_f64(n_brow: number, n_bcol: number, R: number, C: number,
    Ap: number, Aj: number, Ax: number, Xx: number, Yx: number): void;
  _sp_bsr_tocsr_f64(n_brow: number, n_bcol: number, R: number, C: number,
    Ap: number, Aj: number, Ax: number, Bp: number, Bj: number, Bx: number): void;
  _sp_bsr_transpose_f64(n_brow: number, n_bcol: number, R: number, C: number,
    Ap: number, Aj: number, Ax: number, Bp: number, Bi: number, Bx: number): void;
  _sp_bsr_diagonal_f64(k: number, n_brow: number, n_bcol: number, R: number, C: number,
    Ap: number, Aj: number, Ax: number, Yx: number): void;

  // Sparse matrix construction functions
  _sp_coo_tril_f64_entry(nnz: number, k: number,
    row_in: number, col_in: number, data_in: number,
    row_out: number, col_out: number, data_out: number): number;
  _sp_coo_triu_f64_entry(nnz: number, k: number,
    row_in: number, col_in: number, data_in: number,
    row_out: number, col_out: number, data_out: number): number;
  _sp_csr_hstack_f64_entry(n_blocks: number, n_row: number, n_col: number,
    Ap_cat: number, Aj_cat: number, Ax_cat: number,
    Bp: number, Bj: number, Bx: number): void;
  _sp_coo_vstack_f64_entry(n_blocks: number, n_row: number, nnz_per_block: number,
    row_cat: number, col_cat: number, data_cat: number,
    row_out: number, col_out: number, data_out: number): void;
  _sp_coo_hstack_f64_entry(n_blocks: number, n_col: number, nnz_per_block: number,
    row_cat: number, col_cat: number, data_cat: number,
    row_out: number, col_out: number, data_out: number): void;
  _sp_coo_block_diag_f64_entry(n_blocks: number, n_row: number, n_col: number,
    nnz_per_block: number, row_cat: number, col_cat: number, data_cat: number,
    row_out: number, col_out: number, data_out: number): void;
  _sp_coo_kron_f64_entry(nnz_A: number, nnz_B: number, B_nrow: number, B_ncol: number,
    A_row: number, A_col: number, A_data: number,
    B_row: number, B_col: number, B_data: number,
    out_row: number, out_col: number, out_data: number): void;
  _sp_coo_random_f64_entry(nnz: number, n_row: number, n_col: number,
    flat_indices: number, random_values: number,
    row_out: number, col_out: number, data_out: number): void;

  _sp_malloc(size: number): number;
  _sp_free(ptr: number): void;

  // KDTree spatial data structure
  _kdtree_build(
    dataPtr: number,
    n: number,
    m: number,
    leafsize: number,
    balanced: number,
    compact: number
  ): number;
  _kdtree_query_knn(
    treePtr: number,
    xPtr: number,
    nQueries: number,
    k: number,
    p: number,
    eps: number,
    distanceUpperBound: number,
    distancesOutPtr: number,
    indicesOutPtr: number
  ): number;
  _kdtree_query_ball_point(
    treePtr: number,
    xPtr: number,
    nQueries: number,
    rPtr: number,
    p: number,
    eps: number,
    returnLength: number,
    sortOutput: number,
    resultsOutPtr: number,
    countsOutPtr: number
  ): number;
  _kdtree_free(treePtr: number): void;

  // Memory management
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

  // Function table management (for JS→C callbacks)
  addFunction(func: Function, signature: string): number;
  removeFunction(ptr: number): void;
}

export interface WasmModuleOptions {
  locateFile?: (path: string, scriptDirectory: string) => string;
}

export type WasmModuleFactory = (
  options?: WasmModuleOptions
) => Promise<SciWasmModule>;
