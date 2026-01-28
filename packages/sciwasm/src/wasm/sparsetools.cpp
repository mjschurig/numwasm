/**
 * WASM entry points for scipy sparse matrix operations.
 *
 * Instantiates C++ template functions from sparsetools for
 * I=int32_t, T=double (the primary use case for browser-based sparse math).
 *
 * Each function is exported via EMSCRIPTEN_KEEPALIVE with C linkage.
 */

#include <cstdint>
#include <cstring>
#include <vector>
#include <emscripten.h>

#include "npy_compat.h"
#include "bool_ops.h"
#include "complex_ops.h"
#include "csr.h"
#include "csc.h"
#include "coo.h"

extern "C" {

/* ===== CSR kernels ===== */

EMSCRIPTEN_KEEPALIVE
void sp_csr_matvec_f64(int32_t n_row, int32_t n_col,
    const int32_t* Ap, const int32_t* Aj, const double* Ax,
    const double* Xx, double* Yx)
{
    csr_matvec(n_row, n_col, Ap, Aj, Ax, Xx, Yx);
}

EMSCRIPTEN_KEEPALIVE
void sp_csr_matvecs_f64(int32_t n_row, int32_t n_col, int32_t n_vecs,
    const int32_t* Ap, const int32_t* Aj, const double* Ax,
    const double* Xx, double* Yx)
{
    csr_matvecs(n_row, n_col, n_vecs, Ap, Aj, Ax, Xx, Yx);
}

EMSCRIPTEN_KEEPALIVE
void sp_csr_tocsc_f64(int32_t n_row, int32_t n_col,
    const int32_t* Ap, const int32_t* Aj, const double* Ax,
    int32_t* Bp, int32_t* Bi, double* Bx)
{
    csr_tocsc(n_row, n_col, Ap, Aj, Ax, Bp, Bi, Bx);
}

EMSCRIPTEN_KEEPALIVE
void sp_csr_todense_f64(int32_t n_row, int32_t n_col,
    const int32_t* Ap, const int32_t* Aj, const double* Ax,
    double* Bx)
{
    csr_todense(n_row, n_col, Ap, Aj, Ax, Bx);
}

EMSCRIPTEN_KEEPALIVE
void sp_csr_diagonal_f64(int32_t k, int32_t n_row, int32_t n_col,
    const int32_t* Ap, const int32_t* Aj, const double* Ax,
    double* Yx)
{
    csr_diagonal(k, n_row, n_col, Ap, Aj, Ax, Yx);
}

EMSCRIPTEN_KEEPALIVE
void sp_csr_sort_indices_f64(int32_t n_row,
    const int32_t* Ap, int32_t* Aj, double* Ax)
{
    csr_sort_indices(n_row, Ap, Aj, Ax);
}

EMSCRIPTEN_KEEPALIVE
int32_t sp_csr_has_sorted_indices(int32_t n_row,
    const int32_t* Ap, const int32_t* Aj)
{
    return csr_has_sorted_indices(n_row, Ap, Aj) ? 1 : 0;
}

EMSCRIPTEN_KEEPALIVE
int32_t sp_csr_has_canonical_format(int32_t n_row,
    const int32_t* Ap, const int32_t* Aj)
{
    return csr_has_canonical_format(n_row, Ap, Aj) ? 1 : 0;
}

EMSCRIPTEN_KEEPALIVE
void sp_csr_sum_duplicates_f64(int32_t n_row, int32_t n_col,
    int32_t* Ap, int32_t* Aj, double* Ax)
{
    csr_sum_duplicates(n_row, n_col, Ap, Aj, Ax);
}

EMSCRIPTEN_KEEPALIVE
void sp_csr_eliminate_zeros_f64(int32_t n_row, int32_t n_col,
    int32_t* Ap, int32_t* Aj, double* Ax)
{
    csr_eliminate_zeros(n_row, n_col, Ap, Aj, Ax);
}

EMSCRIPTEN_KEEPALIVE
int32_t sp_csr_matmat_maxnnz(int32_t n_row, int32_t n_col,
    const int32_t* Ap, const int32_t* Aj,
    const int32_t* Bp, const int32_t* Bj)
{
    return (int32_t)csr_matmat_maxnnz(n_row, n_col, Ap, Aj, Bp, Bj);
}

EMSCRIPTEN_KEEPALIVE
void sp_csr_matmat_f64(int32_t n_row, int32_t n_col,
    const int32_t* Ap, const int32_t* Aj, const double* Ax,
    const int32_t* Bp, const int32_t* Bj, const double* Bx,
    int32_t* Cp, int32_t* Cj, double* Cx)
{
    csr_matmat(n_row, n_col, Ap, Aj, Ax, Bp, Bj, Bx, Cp, Cj, Cx);
}

EMSCRIPTEN_KEEPALIVE
void sp_csr_plus_csr_f64(int32_t n_row, int32_t n_col,
    const int32_t* Ap, const int32_t* Aj, const double* Ax,
    const int32_t* Bp, const int32_t* Bj, const double* Bx,
    int32_t* Cp, int32_t* Cj, double* Cx)
{
    csr_plus_csr(n_row, n_col, Ap, Aj, Ax, Bp, Bj, Bx, Cp, Cj, Cx);
}

EMSCRIPTEN_KEEPALIVE
void sp_csr_minus_csr_f64(int32_t n_row, int32_t n_col,
    const int32_t* Ap, const int32_t* Aj, const double* Ax,
    const int32_t* Bp, const int32_t* Bj, const double* Bx,
    int32_t* Cp, int32_t* Cj, double* Cx)
{
    csr_minus_csr(n_row, n_col, Ap, Aj, Ax, Bp, Bj, Bx, Cp, Cj, Cx);
}

EMSCRIPTEN_KEEPALIVE
void sp_csr_elmul_csr_f64(int32_t n_row, int32_t n_col,
    const int32_t* Ap, const int32_t* Aj, const double* Ax,
    const int32_t* Bp, const int32_t* Bj, const double* Bx,
    int32_t* Cp, int32_t* Cj, double* Cx)
{
    csr_elmul_csr(n_row, n_col, Ap, Aj, Ax, Bp, Bj, Bx, Cp, Cj, Cx);
}

EMSCRIPTEN_KEEPALIVE
void sp_csr_eldiv_csr_f64(int32_t n_row, int32_t n_col,
    const int32_t* Ap, const int32_t* Aj, const double* Ax,
    const int32_t* Bp, const int32_t* Bj, const double* Bx,
    int32_t* Cp, int32_t* Cj, double* Cx)
{
    csr_eldiv_csr(n_row, n_col, Ap, Aj, Ax, Bp, Bj, Bx, Cp, Cj, Cx);
}

EMSCRIPTEN_KEEPALIVE
void sp_csr_scale_rows_f64(int32_t n_row, int32_t n_col,
    const int32_t* Ap, const int32_t* Aj, double* Ax,
    const double* Xx)
{
    csr_scale_rows(n_row, n_col, Ap, Aj, Ax, Xx);
}

EMSCRIPTEN_KEEPALIVE
void sp_csr_scale_columns_f64(int32_t n_row, int32_t n_col,
    const int32_t* Ap, const int32_t* Aj, double* Ax,
    const double* Xx)
{
    csr_scale_columns(n_row, n_col, Ap, Aj, Ax, Xx);
}

EMSCRIPTEN_KEEPALIVE
void sp_expandptr(int32_t n_row, const int32_t* Ap, int32_t* Bi)
{
    expandptr(n_row, Ap, Bi);
}

EMSCRIPTEN_KEEPALIVE
void sp_csr_row_index_f64(int32_t n_row_idx,
    const int32_t* rows,
    const int32_t* Ap, const int32_t* Aj, const double* Ax,
    int32_t* Bj, double* Bx)
{
    csr_row_index(n_row_idx, rows, Ap, Aj, Ax, Bj, Bx);
}

EMSCRIPTEN_KEEPALIVE
void sp_csr_row_slice_f64(int32_t start, int32_t stop, int32_t step,
    const int32_t* Ap, const int32_t* Aj, const double* Ax,
    int32_t* Bj, double* Bx)
{
    csr_row_slice(start, stop, step, Ap, Aj, Ax, Bj, Bx);
}

EMSCRIPTEN_KEEPALIVE
void sp_csr_column_index1(int32_t n_idx,
    const int32_t* col_idxs,
    int32_t n_row, int32_t n_col,
    const int32_t* Ap, const int32_t* Aj,
    int32_t* col_offsets, int32_t* Bp)
{
    csr_column_index1(n_idx, col_idxs, n_row, n_col, Ap, Aj, col_offsets, Bp);
}

EMSCRIPTEN_KEEPALIVE
void sp_csr_column_index2_f64(const int32_t* col_order,
    const int32_t* col_offsets,
    int32_t nnz,
    const int32_t* Aj, const double* Ax,
    int32_t* Bj, double* Bx)
{
    csr_column_index2(col_order, col_offsets, nnz, Aj, Ax, Bj, Bx);
}

EMSCRIPTEN_KEEPALIVE
int32_t sp_csr_sample_offsets(int32_t n_row, int32_t n_col,
    const int32_t* Ap, const int32_t* Aj,
    int32_t n_samples,
    const int32_t* Bi, const int32_t* Bj,
    int32_t* Bp)
{
    return csr_sample_offsets(n_row, n_col, Ap, Aj, n_samples, Bi, Bj, Bp);
}

EMSCRIPTEN_KEEPALIVE
void sp_csr_sample_values_f64(int32_t n_row, int32_t n_col,
    const int32_t* Ap, const int32_t* Aj, const double* Ax,
    int32_t n_samples,
    const int32_t* Bi, const int32_t* Bj,
    double* Bx)
{
    csr_sample_values(n_row, n_col, Ap, Aj, Ax, n_samples, Bi, Bj, Bx);
}

/*
 * get_csr_submatrix wrapper — since the C++ version uses std::vector* outputs,
 * we provide a two-step interface:
 *   1. sp_get_csr_submatrix_nnz: compute output nnz so caller can allocate
 *   2. sp_get_csr_submatrix_f64: fill preallocated output arrays
 */
EMSCRIPTEN_KEEPALIVE
int32_t sp_get_csr_submatrix_nnz(int32_t n_row, int32_t n_col,
    const int32_t* Ap, const int32_t* Aj,
    int32_t ir0, int32_t ir1, int32_t ic0, int32_t ic1)
{
    int32_t new_n_row = ir1 - ir0;
    int32_t new_nnz = 0;
    for(int32_t i = 0; i < new_n_row; i++){
        int32_t row_start = Ap[ir0+i];
        int32_t row_end   = Ap[ir0+i+1];
        for(int32_t jj = row_start; jj < row_end; jj++){
            if ((Aj[jj] >= ic0) && (Aj[jj] < ic1)) {
                new_nnz++;
            }
        }
    }
    return new_nnz;
}

EMSCRIPTEN_KEEPALIVE
void sp_get_csr_submatrix_f64(int32_t n_row, int32_t n_col,
    const int32_t* Ap, const int32_t* Aj, const double* Ax,
    int32_t ir0, int32_t ir1, int32_t ic0, int32_t ic1,
    int32_t* Bp, int32_t* Bj, double* Bx)
{
    int32_t new_n_row = ir1 - ir0;
    int32_t kk = 0;

    Bp[0] = 0;
    for(int32_t i = 0; i < new_n_row; i++){
        int32_t row_start = Ap[ir0+i];
        int32_t row_end   = Ap[ir0+i+1];
        for(int32_t jj = row_start; jj < row_end; jj++){
            if ((Aj[jj] >= ic0) && (Aj[jj] < ic1)) {
                Bj[kk] = Aj[jj] - ic0;
                Bx[kk] = Ax[jj];
                kk++;
            }
        }
        Bp[i+1] = kk;
    }
}

/* ===== CSC kernels ===== */

EMSCRIPTEN_KEEPALIVE
void sp_csc_matvec_f64(int32_t n_row, int32_t n_col,
    const int32_t* Ap, const int32_t* Ai, const double* Ax,
    const double* Xx, double* Yx)
{
    csc_matvec(n_row, n_col, Ap, Ai, Ax, Xx, Yx);
}

EMSCRIPTEN_KEEPALIVE
void sp_csc_matvecs_f64(int32_t n_row, int32_t n_col, int32_t n_vecs,
    const int32_t* Ap, const int32_t* Ai, const double* Ax,
    const double* Xx, double* Yx)
{
    csc_matvecs(n_row, n_col, n_vecs, Ap, Ai, Ax, Xx, Yx);
}

/* ===== COO kernels ===== */

EMSCRIPTEN_KEEPALIVE
void sp_coo_tocsr_f64(int32_t n_row, int32_t n_col, int32_t nnz,
    const int32_t* Ai, const int32_t* Aj, const double* Ax,
    int32_t* Bp, int32_t* Bj, double* Bx)
{
    coo_tocsr(n_row, n_col, nnz, Ai, Aj, Ax, Bp, Bj, Bx);
}

EMSCRIPTEN_KEEPALIVE
void sp_coo_todense_f64(int32_t n_row, int32_t n_col, int64_t nnz,
    const int32_t* Ai, const int32_t* Aj, const double* Ax,
    double* Bx, int32_t fortran)
{
    coo_todense(n_row, n_col, nnz, Ai, Aj, Ax, Bx, fortran);
}

EMSCRIPTEN_KEEPALIVE
void sp_coo_matvec_f64(int64_t nnz,
    const int32_t* Ai, const int32_t* Aj, const double* Ax,
    const double* Xx, double* Yx)
{
    coo_matvec(nnz, Ai, Aj, Ax, Xx, Yx);
}

/* ===== Memory helpers ===== */

EMSCRIPTEN_KEEPALIVE
void* sp_malloc(int32_t size) {
    return malloc(size);
}

EMSCRIPTEN_KEEPALIVE
void sp_free(void* ptr) {
    free(ptr);
}

} // extern "C"
