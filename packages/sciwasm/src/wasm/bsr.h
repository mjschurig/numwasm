#ifndef __BSR_H__
#define __BSR_H__

/**
 * BSR (Block Sparse Row) matrix operations.
 * From scipy/sparse/sparsetools/bsr.h — adapted for Emscripten/WASM.
 */

#include <vector>
#include <algorithm>
#include <cassert>

#include "npy_compat.h"
#include "csr.h"
#include "dense.h"

static inline npy_intp diagonal_size(const npy_intp k,
                                     const npy_intp rows,
                                     const npy_intp cols)
{
    return std::min(rows + std::min(k, (npy_intp)0),
                    cols - std::max(k, (npy_intp)0));
}


template <class I, class T>
void bsr_diagonal(const I k,
                  const I n_brow,
                  const I n_bcol,
                  const I R,
                  const I C,
                  const I Ap[],
                  const I Aj[],
                  const T Ax[],
                        T Yx[])
{
    const npy_intp RC = R * C;
    const npy_intp D = diagonal_size(k, (npy_intp)n_brow * R,
                                        (npy_intp)n_bcol * C);
    const npy_intp first_row = (k >= 0) ? 0 : -k;
    /* First and next-to-last brows of the diagonal. */
    const npy_intp first_brow = first_row / R;
    const npy_intp last_brow = (first_row + D - 1) / R + 1;

    for (npy_intp brow = first_brow; brow < last_brow; ++brow) {
        /* First and next-to-last bcols of the diagonal in this brow. */
        const npy_intp first_bcol = (brow * R + k) / C;
        const npy_intp last_bcol = ((brow + 1) * R + k - 1) / C + 1;

        for (npy_intp jj = Ap[brow]; jj < Ap[brow + 1]; ++jj) {
            const npy_intp bcol = Aj[jj];

            if (first_bcol <= bcol && bcol < last_bcol) {
                /*
                 * Compute and extract diagonal of block corresponding to the
                 * k-th overall diagonal and add it to output in right place.
                 */
                const npy_intp block_k = brow * R + k - bcol * C;
                const npy_intp block_D = diagonal_size(block_k, R, C);
                const npy_intp block_first_row = (block_k >= 0) ? 0 : -block_k;
                const npy_intp Y_idx = brow * R + block_first_row - first_row;
                const npy_intp Ax_idx = RC * jj +
                                        ((block_k >= 0) ? block_k :
                                                          -C * block_k);

                for (npy_intp kk = 0; kk < block_D; ++kk) {
                    Yx[Y_idx + kk] += Ax[Ax_idx + kk * (C + 1)];
                }

            }
        }
    }
}


/*
 * Compute transpose(A) BSR matrix A
 *
 * Input Arguments:
 *   I  n_brow        - number of row blocks in A
 *   I  n_bcol        - number of column blocks in A
 *   I  R             - rows per block
 *   I  C             - columns per block
 *   I  Ap[n_brow+1]  - row pointer
 *   I  Aj[nblk(A)]   - column indices
 *   T  Ax[nnz(A)]    - nonzeros
 *
 * Output Arguments:
 *   I  Bp[n_col+1]   - row pointer
 *   I  Bj[nblk(A)]   - column indices
 *   T  Bx[nnz(A)]    - nonzeros
 *
 * Note:
 *   Output arrays Bp, Bj, Bx must be preallocated
 */
template <class I, class T>
void bsr_transpose(const I n_brow,
                   const I n_bcol,
                   const I R,
                   const I C,
                   const I Ap[],
                   const I Aj[],
                   const T Ax[],
                         I Bp[],
                         I Bj[],
                         T Bx[])
{
    const I nblks = Ap[n_brow];
    const npy_intp RC    = (npy_intp)R*C;

    //compute permutation of blocks using transpose(CSR)
    std::vector<I> perm_in (nblks);
    std::vector<I> perm_out(nblks);

    for(I i = 0; i < nblks; i++)
        perm_in[i] = i;

    csr_tocsc(n_brow, n_bcol, Ap, Aj, &perm_in[0], Bp, Bj, &perm_out[0]);

    for(I i = 0; i < nblks; i++){
        const T * Ax_blk = Ax + RC * perm_out[i];
              T * Bx_blk = Bx + RC * i;
        for(I r = 0; r < R; r++){
            for(I c = 0; c < C; c++){
                Bx_blk[(npy_intp)c * R + r] = Ax_blk[(npy_intp)r * C + c];
            }
        }
    }
}


/*
 * Compute B = A for BSR matrix A, CSR matrix B
 *
 * Input Arguments:
 *   I  n_brow          - number of block rows in A
 *   I  n_bcol          - number of block columns in A
 *   I  R               - row blocksize
 *   I  C               - column blocksize
 *   I  Ap[n_brow+1]    - block row pointer
 *   I  Aj[nnz(A)]      - block column indices
 *   T  Ax[nnz(A)]      - nonzero blocks
 *
 * Output Arguments:
 *   I  Bp[n_brow*R + 1]- row pointer
 *   I  Bj[nnz(B)]      - column indices
 *   T  Bx[nnz(B)]      - nonzero values
 *
 * Note:
 *   Output arrays must be preallocated
 */
template <class I, class T>
void bsr_tocsr(const I n_brow, const I n_bcol, const I R, const I C,
               const I Ap[], const I Aj[], const T Ax[],
                     I Bp[],       I Bj[],       T Bx[])
{
    // number of elements per block
    const I RC = R*C;
    // nnz
    const I nnz = Ap[n_brow] * RC;
    // last element in Bp is always nnz
    Bp[n_brow * R] = nnz;
    // loop for block row
    for(I brow = 0; brow < n_brow; brow++){
        // size of block rows
        const I brow_size = Ap[brow + 1] - Ap[brow];
        // size of row in csr
        const I row_size = C * brow_size;
        // loop of rows inside block
        for(I r = 0; r < R; r++){
            // csr row number
            const I row = R * brow + r;
            Bp[row] = RC * Ap[brow] + r * row_size;
            // loop for block column
            // block index inside row as loop variable
            for (I bjj = 0; bjj < brow_size; bjj++)
            {
                const I b_ind = Ap[brow] + bjj;
                // block column number
                const I bcol = Aj[b_ind];
                // loop for columns inside block
                for (I c = 0; c < C; c++)
                {
                    // bsr data index in Ax
                    // Ax is in C order
                    const I b_data_ind = RC * b_ind + C * r + c;
                    // csr column number
                    const I col = C * bcol + c;
                    // csr data anc col index in Bj and Bx
                    // start from Bp[row], offset by current bjj*C and c
                    const I data_ind = Bp[row] + bjj * C + c;
                    // assign col and data to Bj and Bx
                    Bj[data_ind] = col;
                    Bx[data_ind] = Ax[b_data_ind];
                }
            }
        }
    }
}


template <class I, class T>
void bsr_matvec(const I n_brow,
                const I n_bcol,
                const I R,
                const I C,
                const I Ap[],
                const I Aj[],
                const T Ax[],
                const T Xx[],
                      T Yx[])
{
    assert(R > 0 && C > 0);

    if( R == 1 && C == 1 ){
        //use CSR for 1x1 blocksize
        csr_matvec(n_brow, n_bcol, Ap, Aj, Ax, Xx, Yx);
        return;
    }

    const npy_intp RC = (npy_intp)R*C;
    for(I i = 0; i < n_brow; i++){
        T * y = Yx + (npy_intp)R * i;
        for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
            const I j = Aj[jj];
            const T * A = Ax + RC * jj;
            const T * x = Xx + (npy_intp)C * j;
            gemv(R, C, A, x, y); // y += A*x
        }
    }
}


#endif
