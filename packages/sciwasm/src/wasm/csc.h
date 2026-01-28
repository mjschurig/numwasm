#ifndef __CSC_H__
#define __CSC_H__

/*
 * CSC sparse matrix kernels.
 * From scipy/sparse/sparsetools/csc.h — adapted for Emscripten/WASM.
 */

#include "csr.h"

/*
 * Compute Y += A*X for CSC matrix A and dense vectors X,Y
 */
template <class I, class T>
void csc_matvec(const I n_row,
                const I n_col,
                const I Ap[],
                const I Ai[],
                const T Ax[],
                const T Xx[],
                      T Yx[])
{
    for(I j = 0; j < n_col; j++){
        I col_start = Ap[j];
        I col_end   = Ap[j+1];

        for(I ii = col_start; ii < col_end; ii++){
            I i    = Ai[ii];
            Yx[i] += Ax[ii] * Xx[j];
        }
    }
}


/*
 * Compute Y += A*X for CSC matrix A and dense block vectors X,Y
 */
template <class I, class T>
void csc_matvecs(const I n_row,
                 const I n_col,
                 const I n_vecs,
                 const I Ap[],
                 const I Ai[],
                 const T Ax[],
                 const T Xx[],
                       T Yx[])
{
    for(I j = 0; j < n_col; j++){
        for(I ii = Ap[j]; ii < Ap[j+1]; ii++){
            const I i = Ai[ii];
            axpy(n_vecs, Ax[ii], Xx + (npy_intp)n_vecs * j, Yx + (npy_intp)n_vecs * i);
        }
    }
}

#endif
