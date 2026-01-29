#ifndef __DIA_H__
#define __DIA_H__

/**
 * DIA (DIAgonal) sparse matrix operations.
 * From scipy/sparse/sparsetools/dia.h — adapted for Emscripten/WASM.
 */

#include <algorithm>
#include "npy_compat.h"

using std::min;
using std::max;

/*
 * Compute Y += A * X for DIA matrix A and dense vectors X, Y
 *
 * Input Arguments:
 *   I  n_row            - number of rows in A
 *   I  n_col            - number of columns in A
 *   I  n_diags          - number of diagonals
 *   I  L                - length of each diagonal
 *   I  offsets[n_diags] - diagonal offsets
 *   T  diags[n_diags,L] - nonzeros
 *   T  Xx[n_col]        - input vector
 *
 * Output Arguments:
 *   T  Yx[n_row]        - output vector
 *
 * Note:
 *   Output array Yx must be preallocated
 *   Negative offsets correspond to lower diagonals
 *   Positive offsets correspond to upper diagonals
 */
template <class I, class T>
void dia_matvec(const I n_row,
                const I n_col,
                const I n_diags,
                const I L,
                const I offsets[],
                const T diags[],
                const T Xx[],
                      T Yx[])
{
    for (I i = 0; i < n_diags; i++){
        const I k = offsets[i];  //diagonal offset

        const I i_start = max<I>(0, -k);
        const I j_start = max<I>(0, k);
        const I j_end   = min<I>(min<I>(n_row + k, n_col), L);

        const I N = j_end - j_start;  //number of elements to process

        const T * diag = diags + (npy_intp)i * L + j_start;
        const T * x = Xx + j_start;
              T * y = Yx + i_start;

        for (I n = 0; n < N; n++) {
            y[n] += diag[n] * x[n];
        }
    }
}


/*
 * Compute B = A for DIA matrix A, CSR matrix B
 *
 * Input Arguments:
 *   I  n_rows            - number of rows in A
 *   I  n_cols            - number of columns in A
 *   I  n_diags           - number of diagonals in A
 *   I  L                 - length of each diagonal in A
 *   I  offsets[n_diags]  - diagonal offsets in A
 *   T  data[n_diags,L]   - diagonals data of A (in C order)
 *   I  order[n_diags]    - indices for traversing offsets[] in ascending order
 *
 * Output Arguments:
 *   T  csr_data[max_nnz] - CSR format data array for B
 *   I  indices[max_nnz]  - CSR format index array for B
 *   I  indptr[n_rows+1]  - CSR format index pointer array for B
 *
 * Return Value:
 *   I  nnz               - number of (non-zero) values stored in B
 *
 * Note:
 *   Output arrays csr_data, indices and indptr must be preallocated and have
 *   sufficient size (with max_nnz >= nnz == A.count_nonzero()), then resulting
 *   arrays csr_data and indices must be truncated to the actual size returned
 *   as nnz
 *
 * Note:
 *   Output has canonical CSR format (sorted indices and no duplicates)
 */
template <class I, class T>
I dia_tocsr(const I n_rows,
            const I n_cols,
            const I n_diags,
            const I L,
            const I offsets[],
            const T data[],
            const I order[],
                  T csr_data[],
                  I indices[],
                  I indptr[])
{
    const I j_end = min(L, n_cols); // columns limit
    indptr[0] = 0;
    I nnz = 0;
    // loop over rows
    for (I i = 0; i < n_rows; ++i) {
        // loop over offsets in ascending order
        for (I k = 0; k < n_diags; ++k) {
            const I n = order[k], // index of diagonal
                    j = i + offsets[n]; // column
            if (j < 0 || j >= j_end)
                continue;
            const T x = data[npy_intp(L) * n + j];
            if (x != 0) {
                indices[nnz] = j;
                csr_data[nnz] = x;
                ++nnz;
            }
        }
        indptr[1 + i] = nnz;
    }
    return nnz; // actual output lengths
}


#endif
