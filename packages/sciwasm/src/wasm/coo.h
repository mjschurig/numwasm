#ifndef __COO_H__
#define __COO_H__

/*
 * COO sparse matrix kernels.
 * From scipy/sparse/sparsetools/coo.h — adapted for Emscripten/WASM.
 */

#include <algorithm>
#include "npy_compat.h"

/*
 * Compute B = A for COO matrix A, CSR matrix B
 */
template <class I, class T>
void coo_tocsr(const I n_row,
               const I n_col,
               const I nnz,
               const I Ai[],
               const I Aj[],
               const T Ax[],
                     I Bp[],
                     I Bj[],
                     T Bx[])
{
    //compute number of non-zero entries per row of A
    std::fill(Bp, Bp + n_row, 0);

    for (I n = 0; n < nnz; n++){
        Bp[Ai[n]]++;
    }

    //cumsum the nnz per row to get Bp[]
    for(I i = 0, cumsum = 0; i < n_row; i++){
        I temp = Bp[i];
        Bp[i] = cumsum;
        cumsum += temp;
    }
    Bp[n_row] = nnz;

    //write Aj,Ax into Bj,Bx
    for(I n = 0; n < nnz; n++){
        I row  = Ai[n];
        I dest = Bp[row];

        Bj[dest] = Aj[n];
        Bx[dest] = Ax[n];

        Bp[row]++;
    }

    for(I i = 0, last = 0; i <= n_row; i++){
        I temp = Bp[i];
        Bp[i]  = last;
        last   = temp;
    }

    //now Bp,Bj,Bx form a CSR representation (with possible duplicates)
}

/*
 * Compute B += A for COO matrix A, dense matrix B
 */
template <class I, class T>
void coo_todense(const I n_row,
                 const I n_col,
                 const npy_int64 nnz,
                 const I Ai[],
                 const I Aj[],
                 const T Ax[],
                       T Bx[],
                 const int fortran)
{
    if (!fortran) {
        for(npy_int64 n = 0; n < nnz; n++){
            Bx[ (npy_intp)n_col * Ai[n] + Aj[n] ] += Ax[n];
        }
    }
    else {
        for(npy_int64 n = 0; n < nnz; n++){
            Bx[ (npy_intp)n_row * Aj[n] + Ai[n] ] += Ax[n];
        }
    }
}


/*
 * Compute Y += A*X for COO matrix A and dense vectors X,Y
 */
template <class I, class T>
void coo_matvec(const npy_int64 nnz,
                const I Ai[],
                const I Aj[],
                const T Ax[],
                const T Xx[],
                      T Yx[])
{
    for(npy_int64 n = 0; n < nnz; n++){
        Yx[Ai[n]] += Ax[n] * Xx[Aj[n]];
    }
}

#endif
