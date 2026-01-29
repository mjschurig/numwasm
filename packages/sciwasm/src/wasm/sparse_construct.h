/**
 * sparse_construct.h - Sparse matrix construction utilities
 *
 * C++ implementations for sparse matrix construction operations:
 * - tril/triu: triangle extraction
 * - hstack/vstack: stacking operations
 * - block_diag: block diagonal construction
 * - kron: Kronecker product
 * - random: random sparse matrix generation
 */

#ifndef SPARSE_CONSTRUCT_H
#define SPARSE_CONSTRUCT_H

#include <cstdint>
#include <vector>
#include <algorithm>

/**
 * Extract lower triangle from COO matrix
 *
 * Keep elements where row + k >= col (on or below k-th diagonal)
 *
 * @param nnz Number of non-zeros in input
 * @param k Diagonal offset (0=main, >0=above, <0=below)
 * @param row_in Input row indices
 * @param col_in Input column indices
 * @param data_in Input data values
 * @param row_out Output row indices (pre-allocated, size nnz)
 * @param col_out Output column indices (pre-allocated, size nnz)
 * @param data_out Output data values (pre-allocated, size nnz)
 * @return Number of elements in output
 */
inline int32_t sp_coo_tril_f64(
    int32_t nnz,
    int32_t k,
    const int32_t* row_in,
    const int32_t* col_in,
    const double* data_in,
    int32_t* row_out,
    int32_t* col_out,
    double* data_out)
{
    int32_t out_idx = 0;
    for (int32_t i = 0; i < nnz; i++) {
        if (row_in[i] + k >= col_in[i]) {
            row_out[out_idx] = row_in[i];
            col_out[out_idx] = col_in[i];
            data_out[out_idx] = data_in[i];
            out_idx++;
        }
    }
    return out_idx;
}

/**
 * Extract upper triangle from COO matrix
 *
 * Keep elements where row + k <= col (on or above k-th diagonal)
 *
 * @param nnz Number of non-zeros in input
 * @param k Diagonal offset (0=main, >0=above, <0=below)
 * @param row_in Input row indices
 * @param col_in Input column indices
 * @param data_in Input data values
 * @param row_out Output row indices (pre-allocated, size nnz)
 * @param col_out Output column indices (pre-allocated, size nnz)
 * @param data_out Output data values (pre-allocated, size nnz)
 * @return Number of elements in output
 */
inline int32_t sp_coo_triu_f64(
    int32_t nnz,
    int32_t k,
    const int32_t* row_in,
    const int32_t* col_in,
    const double* data_in,
    int32_t* row_out,
    int32_t* col_out,
    double* data_out)
{
    int32_t out_idx = 0;
    for (int32_t i = 0; i < nnz; i++) {
        if (row_in[i] + k <= col_in[i]) {
            row_out[out_idx] = row_in[i];
            col_out[out_idx] = col_in[i];
            data_out[out_idx] = data_in[i];
            out_idx++;
        }
    }
    return out_idx;
}

/**
 * Horizontal stack of CSR matrices
 *
 * Stack multiple CSR matrices with the same number of rows side by side.
 * Based on scipy's csr_hstack (csr.h:1649-1692).
 *
 * @param n_blocks Number of blocks to stack
 * @param n_row Number of rows (same for all blocks)
 * @param n_col Array of column counts per block [n_blocks]
 * @param Ap_cat Concatenated indptr arrays [(n_row+1) * n_blocks]
 * @param Aj_cat Concatenated indices arrays
 * @param Ax_cat Concatenated data arrays
 * @param Bp Output indptr [n_row + 1]
 * @param Bj Output indices
 * @param Bx Output data
 */
inline void sp_csr_hstack_f64(
    int32_t n_blocks,
    int32_t n_row,
    const int32_t* n_col,
    const int32_t* Ap_cat,
    const int32_t* Aj_cat,
    const double* Ax_cat,
    int32_t* Bp,
    int32_t* Bj,
    double* Bx)
{
    // Compute column offsets and locate block boundaries
    std::vector<int32_t> col_offset(n_blocks);
    std::vector<const int32_t*> bAp(n_blocks);
    std::vector<const int32_t*> bAj(n_blocks);
    std::vector<const double*> bAx(n_blocks);

    col_offset[0] = 0;
    bAp[0] = Ap_cat;
    bAj[0] = Aj_cat;
    bAx[0] = Ax_cat;

    for (int32_t b = 1; b < n_blocks; b++) {
        col_offset[b] = col_offset[b - 1] + n_col[b - 1];
        bAp[b] = bAp[b - 1] + (n_row + 1);
        bAj[b] = bAj[b - 1] + bAp[b - 1][n_row];
        bAx[b] = bAx[b - 1] + bAp[b - 1][n_row];
    }

    // Build output matrix
    Bp[0] = 0;
    int32_t s = 0;
    for (int32_t i = 0; i < n_row; i++) {
        for (int32_t b = 0; b < n_blocks; b++) {
            int32_t jj_start = bAp[b][i];
            int32_t jj_end = bAp[b][i + 1];
            int32_t offset = col_offset[b];
            for (int32_t jj = jj_start; jj < jj_end; jj++) {
                Bj[s] = bAj[b][jj] + offset;
                Bx[s] = bAx[b][jj];
                s++;
            }
        }
        Bp[i + 1] = s;
    }
}

/**
 * Vertical stack of COO matrices
 *
 * Stack multiple COO matrices with the same number of columns on top of each other.
 *
 * @param n_blocks Number of blocks to stack
 * @param n_row Array of row counts per block [n_blocks]
 * @param nnz_per_block Array of nnz counts per block [n_blocks]
 * @param row_cat Concatenated row arrays
 * @param col_cat Concatenated column arrays
 * @param data_cat Concatenated data arrays
 * @param row_out Output row indices
 * @param col_out Output column indices
 * @param data_out Output data values
 */
inline void sp_coo_vstack_f64(
    int32_t n_blocks,
    const int32_t* n_row,
    const int32_t* nnz_per_block,
    const int32_t* row_cat,
    const int32_t* col_cat,
    const double* data_cat,
    int32_t* row_out,
    int32_t* col_out,
    double* data_out)
{
    int32_t row_offset = 0;
    int32_t in_idx = 0;
    int32_t out_idx = 0;

    for (int32_t b = 0; b < n_blocks; b++) {
        for (int32_t i = 0; i < nnz_per_block[b]; i++) {
            row_out[out_idx] = row_cat[in_idx] + row_offset;
            col_out[out_idx] = col_cat[in_idx];
            data_out[out_idx] = data_cat[in_idx];
            in_idx++;
            out_idx++;
        }
        row_offset += n_row[b];
    }
}

/**
 * Horizontal stack of COO matrices
 *
 * Stack multiple COO matrices with the same number of rows side by side.
 *
 * @param n_blocks Number of blocks to stack
 * @param n_col Array of column counts per block [n_blocks]
 * @param nnz_per_block Array of nnz counts per block [n_blocks]
 * @param row_cat Concatenated row arrays
 * @param col_cat Concatenated column arrays
 * @param data_cat Concatenated data arrays
 * @param row_out Output row indices
 * @param col_out Output column indices
 * @param data_out Output data values
 */
inline void sp_coo_hstack_f64(
    int32_t n_blocks,
    const int32_t* n_col,
    const int32_t* nnz_per_block,
    const int32_t* row_cat,
    const int32_t* col_cat,
    const double* data_cat,
    int32_t* row_out,
    int32_t* col_out,
    double* data_out)
{
    int32_t col_offset = 0;
    int32_t in_idx = 0;
    int32_t out_idx = 0;

    for (int32_t b = 0; b < n_blocks; b++) {
        for (int32_t i = 0; i < nnz_per_block[b]; i++) {
            row_out[out_idx] = row_cat[in_idx];
            col_out[out_idx] = col_cat[in_idx] + col_offset;
            data_out[out_idx] = data_cat[in_idx];
            in_idx++;
            out_idx++;
        }
        col_offset += n_col[b];
    }
}

/**
 * Block diagonal construction from COO matrices
 *
 * Place matrices along the main diagonal of a larger matrix.
 *
 * @param n_blocks Number of blocks
 * @param n_row Array of row counts per block [n_blocks]
 * @param n_col Array of column counts per block [n_blocks]
 * @param nnz_per_block Array of nnz counts per block [n_blocks]
 * @param row_cat Concatenated row arrays
 * @param col_cat Concatenated column arrays
 * @param data_cat Concatenated data arrays
 * @param row_out Output row indices
 * @param col_out Output column indices
 * @param data_out Output data values
 */
inline void sp_coo_block_diag_f64(
    int32_t n_blocks,
    const int32_t* n_row,
    const int32_t* n_col,
    const int32_t* nnz_per_block,
    const int32_t* row_cat,
    const int32_t* col_cat,
    const double* data_cat,
    int32_t* row_out,
    int32_t* col_out,
    double* data_out)
{
    int32_t row_offset = 0;
    int32_t col_offset = 0;
    int32_t in_idx = 0;
    int32_t out_idx = 0;

    for (int32_t b = 0; b < n_blocks; b++) {
        for (int32_t i = 0; i < nnz_per_block[b]; i++) {
            row_out[out_idx] = row_cat[in_idx] + row_offset;
            col_out[out_idx] = col_cat[in_idx] + col_offset;
            data_out[out_idx] = data_cat[in_idx];
            in_idx++;
            out_idx++;
        }
        row_offset += n_row[b];
        col_offset += n_col[b];
    }
}

/**
 * Kronecker product of two COO matrices
 *
 * Computes A ⊗ B where output[i*B_nrow + p, j*B_ncol + q] = A[i,j] * B[p,q]
 * Output nnz = nnz_A * nnz_B
 *
 * @param nnz_A Number of non-zeros in A
 * @param nnz_B Number of non-zeros in B
 * @param B_nrow Number of rows in B
 * @param B_ncol Number of columns in B
 * @param A_row Row indices of A
 * @param A_col Column indices of A
 * @param A_data Data values of A
 * @param B_row Row indices of B
 * @param B_col Column indices of B
 * @param B_data Data values of B
 * @param out_row Output row indices [nnz_A * nnz_B]
 * @param out_col Output column indices [nnz_A * nnz_B]
 * @param out_data Output data values [nnz_A * nnz_B]
 */
inline void sp_coo_kron_f64(
    int32_t nnz_A,
    int32_t nnz_B,
    int32_t B_nrow,
    int32_t B_ncol,
    const int32_t* A_row,
    const int32_t* A_col,
    const double* A_data,
    const int32_t* B_row,
    const int32_t* B_col,
    const double* B_data,
    int32_t* out_row,
    int32_t* out_col,
    double* out_data)
{
    int32_t out_idx = 0;
    for (int32_t a = 0; a < nnz_A; a++) {
        int32_t ai = A_row[a];
        int32_t aj = A_col[a];
        double av = A_data[a];
        for (int32_t b = 0; b < nnz_B; b++) {
            out_row[out_idx] = ai * B_nrow + B_row[b];
            out_col[out_idx] = aj * B_ncol + B_col[b];
            out_data[out_idx] = av * B_data[b];
            out_idx++;
        }
    }
}

/**
 * Generate COO arrays from flat indices
 *
 * Converts flat indices to (row, col) coordinates.
 * Random index generation is done in TypeScript, this just converts.
 *
 * @param nnz Number of non-zeros
 * @param n_row Number of rows
 * @param n_col Number of columns
 * @param flat_indices Flat indices in [0, n_row*n_col)
 * @param random_values Random values for data
 * @param row_out Output row indices
 * @param col_out Output column indices
 * @param data_out Output data values
 */
inline void sp_coo_random_f64(
    int32_t nnz,
    int32_t n_row,
    int32_t n_col,
    const int32_t* flat_indices,
    const double* random_values,
    int32_t* row_out,
    int32_t* col_out,
    double* data_out)
{
    for (int32_t i = 0; i < nnz; i++) {
        int32_t flat = flat_indices[i];
        row_out[i] = flat / n_col;
        col_out[i] = flat % n_col;
        data_out[i] = random_values[i];
    }
}

#endif // SPARSE_CONSTRUCT_H
