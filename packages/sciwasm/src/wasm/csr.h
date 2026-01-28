#ifndef __CSR_H__
#define __CSR_H__

/*
 * CSR sparse matrix kernels.
 * From scipy/sparse/sparsetools/csr.h — adapted for Emscripten/WASM.
 * Removed SPTOOLS_CSR_EXTERN_TEMPLATE macro (not needed for WASM).
 */

#include <set>
#include <vector>
#include <algorithm>
#include <functional>
#include <numeric>
#include <stdexcept>

#include "npy_compat.h"
#include "util.h"
#include "dense.h"

/*
 * Extract k-th diagonal of CSR matrix A
 */
template <class I, class T>
void csr_diagonal(const I k,
                  const I n_row,
                  const I n_col,
                  const I Ap[],
                  const I Aj[],
                  const T Ax[],
                        T Yx[])
{
    const I first_row = (k >= 0) ? 0 : -k;
    const I first_col = (k >= 0) ? k : 0;
    const I N = std::min(n_row - first_row, n_col - first_col);

    for (I i = 0; i < N; ++i) {
        const I row = first_row + i;
        const I col = first_col + i;
        const I row_begin = Ap[row];
        const I row_end = Ap[row + 1];

        T diag = 0;
        for (I j = row_begin; j < row_end; ++j) {
            if (Aj[j] == col) {
                diag += Ax[j];
            }
        }

        Yx[i] = diag;
    }
}


/*
 * Expand a compressed row pointer into a row array
 */
template <class I>
void expandptr(const I n_row,
               const I Ap[],
                     I Bi[])
{
    for(I i = 0; i < n_row; i++){
        for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
            Bi[jj] = i;
        }
    }
}


/*
 * Scale the rows of a CSR matrix *in place*
 *   A[i,:] *= X[i]
 */
template <class I, class T>
void csr_scale_rows(const I n_row,
                    const I n_col,
                    const I Ap[],
                    const I Aj[],
                          T Ax[],
                    const T Xx[])
{
    for(I i = 0; i < n_row; i++){
        for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
            Ax[jj] *= Xx[i];
        }
    }
}


/*
 * Scale the columns of a CSR matrix *in place*
 *   A[:,i] *= X[i]
 */
template <class I, class T>
void csr_scale_columns(const I n_row,
                       const I n_col,
                       const I Ap[],
                       const I Aj[],
                             T Ax[],
                       const T Xx[])
{
    const I nnz = Ap[n_row];
    for(I i = 0; i < nnz; i++){
        Ax[i] *= Xx[Aj[i]];
    }
}


/*
 * Compute B += A for CSR matrix A, C-contiguous dense matrix B
 */
template <class I, class T>
void csr_todense(const I n_row,
                 const I n_col,
                 const I Ap[],
                 const I Aj[],
                 const T Ax[],
                       T Bx[])
{
    T * Bx_row = Bx;
    for(I i = 0; i < n_row; i++){
        for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
            Bx_row[Aj[jj]] += Ax[jj];
        }
        Bx_row += (npy_intp)n_col;
    }
}


/*
 * Determine whether the CSR column indices are in sorted order.
 */
template <class I>
bool csr_has_sorted_indices(const I n_row,
                            const I Ap[],
                            const I Aj[])
{
  for(I i = 0; i < n_row; i++){
      for(I jj = Ap[i]; jj < Ap[i+1] - 1; jj++){
          if(Aj[jj] > Aj[jj+1]){
              return false;
          }
      }
  }
  return true;
}


/*
 * Determine whether the matrix structure is canonical CSR.
 */
template <class I>
bool csr_has_canonical_format(const I n_row,
                              const I Ap[],
                              const I Aj[])
{
    for(I i = 0; i < n_row; i++){
        if (Ap[i] > Ap[i+1])
            return false;
        for(I jj = Ap[i] + 1; jj < Ap[i+1]; jj++){
            if( !(Aj[jj-1] < Aj[jj]) ){
                return false;
            }
        }
    }
    return true;
}


template< class T1, class T2 >
bool kv_pair_less(const std::pair<T1,T2>& x, const std::pair<T1,T2>& y){
    return x.first < y.first;
}

/*
 * Sort CSR column indices inplace
 */
template<class I, class T>
void csr_sort_indices(const I n_row,
                      const I Ap[],
                            I Aj[],
                            T Ax[])
{
    std::vector< std::pair<I,T> > temp;

    for(I i = 0; i < n_row; i++){
        I row_start = Ap[i];
        I row_end   = Ap[i+1];

        temp.resize(row_end - row_start);
        for (I jj = row_start, n = 0; jj < row_end; jj++, n++){
            temp[n].first  = Aj[jj];
            temp[n].second = Ax[jj];
        }

        std::sort(temp.begin(),temp.end(),kv_pair_less<I,T>);

        for(I jj = row_start, n = 0; jj < row_end; jj++, n++){
            Aj[jj] = temp[n].first;
            Ax[jj] = temp[n].second;
        }
    }
}


/*
 * Compute B = A for CSR matrix A, CSC matrix B
 * Also used for CSR->CSC, CSC->CSR, transpose operations.
 */
template <class I, class T>
void csr_tocsc(const I n_row,
               const I n_col,
               const I Ap[],
               const I Aj[],
               const T Ax[],
                     I Bp[],
                     I Bi[],
                     T Bx[])
{
    const I nnz = Ap[n_row];

    //compute number of non-zero entries per column of A
    std::fill(Bp, Bp + n_col, 0);

    for (I n = 0; n < nnz; n++){
        Bp[Aj[n]]++;
    }

    //cumsum the nnz per column to get Bp[]
    for(I col = 0, cumsum = 0; col < n_col; col++){
        I temp  = Bp[col];
        Bp[col] = cumsum;
        cumsum += temp;
    }
    Bp[n_col] = nnz;

    for(I row = 0; row < n_row; row++){
        for(I jj = Ap[row]; jj < Ap[row+1]; jj++){
            I col  = Aj[jj];
            I dest = Bp[col];

            Bi[dest] = row;
            Bx[dest] = Ax[jj];

            Bp[col]++;
        }
    }

    for(I col = 0, last = 0; col <= n_col; col++){
        I temp  = Bp[col];
        Bp[col] = last;
        last    = temp;
    }
}


/*
 * Compute the number of non-zeroes (nnz) in the result of C = A * B.
 */
template <class I>
npy_intp csr_matmat_maxnnz(const I n_row,
                           const I n_col,
                           const I Ap[],
                           const I Aj[],
                           const I Bp[],
                           const I Bj[])
{
    std::vector<I> mask(n_col, -1);

    npy_intp nnz = 0;
    for(I i = 0; i < n_row; i++){
        npy_intp row_nnz = 0;

        for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
            I j = Aj[jj];
            for(I kk = Bp[j]; kk < Bp[j+1]; kk++){
                I k = Bj[kk];
                if(mask[k] != i){
                    mask[k] = i;
                    row_nnz++;
                }
            }
        }

        npy_intp next_nnz = nnz + row_nnz;

        if (row_nnz > NPY_MAX_INTP - nnz) {
            throw std::overflow_error("nnz of the result is too large");
        }

        nnz = next_nnz;
    }

    return nnz;
}

/*
 * Compute CSR entries for matrix C = A*B.
 */
template <class I, class T>
void csr_matmat(const I n_row,
                const I n_col,
                const I Ap[],
                const I Aj[],
                const T Ax[],
                const I Bp[],
                const I Bj[],
                const T Bx[],
                      I Cp[],
                      I Cj[],
                      T Cx[])
{
    std::vector<I> next(n_col,-1);
    std::vector<T> sums(n_col, 0);

    I nnz = 0;

    Cp[0] = 0;

    for(I i = 0; i < n_row; i++){
        I head   = -2;
        I length =  0;

        I jj_start = Ap[i];
        I jj_end   = Ap[i+1];
        for(I jj = jj_start; jj < jj_end; jj++){
            I j = Aj[jj];
            T v = Ax[jj];

            I kk_start = Bp[j];
            I kk_end   = Bp[j+1];
            for(I kk = kk_start; kk < kk_end; kk++){
                I k = Bj[kk];

                sums[k] += v*Bx[kk];

                if(next[k] == -1){
                    next[k] = head;
                    head  = k;
                    length++;
                }
            }
        }

        for(I jj = 0; jj < length; jj++){

            if(sums[head] != 0){
                Cj[nnz] = head;
                Cx[nnz] = sums[head];
                nnz++;
            }

            I temp = head;
            head = next[head];

            next[temp] = -1; //clear arrays
            sums[temp] =  0;
        }

        Cp[i+1] = nnz;
    }
}


/*
 * Compute C = A (binary_op) B for CSR matrices (general, non-canonical).
 */
template <class I, class T, class T2, class binary_op>
void csr_binop_csr_general(const I n_row, const I n_col,
                           const I Ap[], const I Aj[], const T Ax[],
                           const I Bp[], const I Bj[], const T Bx[],
                                 I Cp[],       I Cj[],       T2 Cx[],
                           const binary_op& op)
{
    std::vector<I>  next(n_col,-1);
    std::vector<T> A_row(n_col, 0);
    std::vector<T> B_row(n_col, 0);

    I nnz = 0;
    Cp[0] = 0;

    for(I i = 0; i < n_row; i++){
        I head   = -2;
        I length =  0;

        I i_start = Ap[i];
        I i_end   = Ap[i+1];
        for(I jj = i_start; jj < i_end; jj++){
            I j = Aj[jj];
            A_row[j] += Ax[jj];
            if(next[j] == -1){
                next[j] = head;
                head = j;
                length++;
            }
        }

        i_start = Bp[i];
        i_end   = Bp[i+1];
        for(I jj = i_start; jj < i_end; jj++){
            I j = Bj[jj];
            B_row[j] += Bx[jj];
            if(next[j] == -1){
                next[j] = head;
                head = j;
                length++;
            }
        }

        for(I jj = 0; jj < length; jj++){
            T result = op(A_row[head], B_row[head]);
            if(result != 0){
                Cj[nnz] = head;
                Cx[nnz] = result;
                nnz++;
            }
            I temp = head;
            head = next[head];
            next[temp]  = -1;
            A_row[temp] =  0;
            B_row[temp] =  0;
        }

        Cp[i + 1] = nnz;
    }
}


/*
 * Compute C = A (binary_op) B for CSR matrices (canonical).
 */
template <class I, class T, class T2, class binary_op>
void csr_binop_csr_canonical(const I n_row, const I n_col,
                             const I Ap[], const I Aj[], const T Ax[],
                             const I Bp[], const I Bj[], const T Bx[],
                                   I Cp[],       I Cj[],       T2 Cx[],
                             const binary_op& op)
{
    Cp[0] = 0;
    I nnz = 0;

    for(I i = 0; i < n_row; i++){
        I A_pos = Ap[i];
        I B_pos = Bp[i];
        I A_end = Ap[i+1];
        I B_end = Bp[i+1];

        while(A_pos < A_end && B_pos < B_end){
            I A_j = Aj[A_pos];
            I B_j = Bj[B_pos];

            if(A_j == B_j){
                T result = op(Ax[A_pos],Bx[B_pos]);
                if(result != 0){
                    Cj[nnz] = A_j;
                    Cx[nnz] = result;
                    nnz++;
                }
                A_pos++;
                B_pos++;
            } else if (A_j < B_j) {
                T result = op(Ax[A_pos],0);
                if (result != 0){
                    Cj[nnz] = A_j;
                    Cx[nnz] = result;
                    nnz++;
                }
                A_pos++;
            } else {
                T result = op(0,Bx[B_pos]);
                if (result != 0){
                    Cj[nnz] = B_j;
                    Cx[nnz] = result;
                    nnz++;
                }
                B_pos++;
            }
        }

        while(A_pos < A_end){
            T result = op(Ax[A_pos],0);
            if (result != 0){
                Cj[nnz] = Aj[A_pos];
                Cx[nnz] = result;
                nnz++;
            }
            A_pos++;
        }
        while(B_pos < B_end){
            T result = op(0,Bx[B_pos]);
            if (result != 0){
                Cj[nnz] = Bj[B_pos];
                Cx[nnz] = result;
                nnz++;
            }
            B_pos++;
        }

        Cp[i+1] = nnz;
    }
}


/*
 * Compute C = A (binary_op) B for CSR matrices A,B.
 */
template <class I, class T, class T2, class binary_op>
void csr_binop_csr(const I n_row,
                   const I n_col,
                   const I Ap[],
                   const I Aj[],
                   const T Ax[],
                   const I Bp[],
                   const I Bj[],
                   const T Bx[],
                         I Cp[],
                         I Cj[],
                         T2 Cx[],
                   const binary_op& op)
{
    if (csr_has_canonical_format(n_row,Ap,Aj) && csr_has_canonical_format(n_row,Bp,Bj))
        csr_binop_csr_canonical(n_row, n_col, Ap, Aj, Ax, Bp, Bj, Bx, Cp, Cj, Cx, op);
    else
        csr_binop_csr_general(n_row, n_col, Ap, Aj, Ax, Bp, Bj, Bx, Cp, Cj, Cx, op);
}

/* element-wise binary operations */
template <class I, class T>
void csr_elmul_csr(const I n_row, const I n_col,
                   const I Ap[], const I Aj[], const T Ax[],
                   const I Bp[], const I Bj[], const T Bx[],
                         I Cp[],       I Cj[],       T Cx[])
{
    csr_binop_csr(n_row,n_col,Ap,Aj,Ax,Bp,Bj,Bx,Cp,Cj,Cx,std::multiplies<T>());
}

template <class I, class T>
void csr_eldiv_csr(const I n_row, const I n_col,
                   const I Ap[], const I Aj[], const T Ax[],
                   const I Bp[], const I Bj[], const T Bx[],
                         I Cp[],       I Cj[],       T Cx[])
{
    csr_binop_csr(n_row,n_col,Ap,Aj,Ax,Bp,Bj,Bx,Cp,Cj,Cx,safe_divides<T>());
}

template <class I, class T>
void csr_plus_csr(const I n_row, const I n_col,
                  const I Ap[], const I Aj[], const T Ax[],
                  const I Bp[], const I Bj[], const T Bx[],
                        I Cp[],       I Cj[],       T Cx[])
{
    csr_binop_csr(n_row,n_col,Ap,Aj,Ax,Bp,Bj,Bx,Cp,Cj,Cx,std::plus<T>());
}

template <class I, class T>
void csr_minus_csr(const I n_row, const I n_col,
                   const I Ap[], const I Aj[], const T Ax[],
                   const I Bp[], const I Bj[], const T Bx[],
                         I Cp[],       I Cj[],       T Cx[])
{
    csr_binop_csr(n_row,n_col,Ap,Aj,Ax,Bp,Bj,Bx,Cp,Cj,Cx,std::minus<T>());
}


/*
 * Sum together duplicate column entries in each row of CSR matrix A
 */
template <class I, class T>
void csr_sum_duplicates(const I n_row,
                        const I n_col,
                              I Ap[],
                              I Aj[],
                              T Ax[])
{
    I nnz = 0;
    I row_end = 0;
    for(I i = 0; i < n_row; i++){
        I jj = row_end;
        row_end = Ap[i+1];
        while( jj < row_end ){
            I j = Aj[jj];
            T x = Ax[jj];
            jj++;
            while( jj < row_end && Aj[jj] == j ){
                x += Ax[jj];
                jj++;
            }
            Aj[nnz] = j;
            Ax[nnz] = x;
            nnz++;
        }
        Ap[i+1] = nnz;
    }
}

/*
 * Eliminate zero entries from CSR matrix A
 */
template <class I, class T>
void csr_eliminate_zeros(const I n_row,
                         const I n_col,
                               I Ap[],
                               I Aj[],
                               T Ax[])
{
    I nnz = 0;
    I row_end = 0;
    for(I i = 0; i < n_row; i++){
        I jj = row_end;
        row_end = Ap[i+1];
        while( jj < row_end ){
            I j = Aj[jj];
            T x = Ax[jj];
            if(x != 0){
                Aj[nnz] = j;
                Ax[nnz] = x;
                nnz++;
            }
            jj++;
        }
        Ap[i+1] = nnz;
    }
}


/*
 * Compute Y += A*X for CSR matrix A and dense vectors X,Y
 */
template <class I, class T>
void csr_matvec(const I n_row,
                const I n_col,
                const I Ap[],
                const I Aj[],
                const T Ax[],
                const T Xx[],
                      T Yx[])
{
    for(I i = 0; i < n_row; i++){
        T sum = Yx[i];
        for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
            sum += Ax[jj] * Xx[Aj[jj]];
        }
        Yx[i] = sum;
    }
}


/*
 * Compute Y += A*X for CSR matrix A and dense block vectors X,Y
 */
template <class I, class T>
void csr_matvecs(const I n_row,
                 const I n_col,
                 const I n_vecs,
                 const I Ap[],
                 const I Aj[],
                 const T Ax[],
                 const T Xx[],
                       T Yx[])
{
    for(I i = 0; i < n_row; i++){
        T * y = Yx + (npy_intp)n_vecs * i;
        for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
            const I j = Aj[jj];
            const T a = Ax[jj];
            const T * x = Xx + (npy_intp)n_vecs * j;
            axpy(n_vecs, a, x, y);
        }
    }
}


/*
 * Extract submatrix from CSR matrix A.
 * Output arrays are allocated internally via std::vector.
 */
template<class I, class T>
void get_csr_submatrix(const I n_row,
                       const I n_col,
                       const I Ap[],
                       const I Aj[],
                       const T Ax[],
                       const I ir0,
                       const I ir1,
                       const I ic0,
                       const I ic1,
                       std::vector<I>* Bp,
                       std::vector<I>* Bj,
                       std::vector<T>* Bx)
{
    I new_n_row = ir1 - ir0;
    I new_nnz = 0;
    I kk = 0;

    for(I i = 0; i < new_n_row; i++){
        I row_start = Ap[ir0+i];
        I row_end   = Ap[ir0+i+1];
        for(I jj = row_start; jj < row_end; jj++){
            if ((Aj[jj] >= ic0) && (Aj[jj] < ic1)) {
                new_nnz++;
            }
        }
    }

    Bp->resize(new_n_row+1);
    Bj->resize(new_nnz);
    Bx->resize(new_nnz);

    (*Bp)[0] = 0;
    for(I i = 0; i < new_n_row; i++){
        I row_start = Ap[ir0+i];
        I row_end   = Ap[ir0+i+1];
        for(I jj = row_start; jj < row_end; jj++){
            if ((Aj[jj] >= ic0) && (Aj[jj] < ic1)) {
                (*Bj)[kk] = Aj[jj] - ic0;
                (*Bx)[kk] = Ax[jj];
                kk++;
            }
        }
        (*Bp)[i+1] = kk;
    }
}


/*
 * Slice rows given as an array of indices.
 */
template<class I, class T>
void csr_row_index(const I n_row_idx,
                   const I rows[],
                   const I Ap[],
                   const I Aj[],
                   const T Ax[],
                   I Bj[],
                   T Bx[])
{
    for(I i = 0; i < n_row_idx; i++){
        const I row = rows[i];
        const I row_start = Ap[row];
        const I row_end   = Ap[row+1];
        Bj = std::copy(Aj + row_start, Aj + row_end, Bj);
        Bx = std::copy(Ax + row_start, Ax + row_end, Bx);
    }
}


/*
 * Slice rows given as a (start, stop, step) tuple.
 */
template<class I, class T>
void csr_row_slice(const I start,
                   const I stop,
                   const I step,
                   const I Ap[],
                   const I Aj[],
                   const T Ax[],
                   I Bj[],
                   T Bx[])
{
    if (step > 0) {
        for(I row = start; row < stop; row += step){
            const I row_start = Ap[row];
            const I row_end   = Ap[row+1];
            Bj = std::copy(Aj + row_start, Aj + row_end, Bj);
            Bx = std::copy(Ax + row_start, Ax + row_end, Bx);
        }
    } else {
        for(I row = start; row > stop; row += step){
            const I row_start = Ap[row];
            const I row_end   = Ap[row+1];
            Bj = std::copy(Aj + row_start, Aj + row_end, Bj);
            Bx = std::copy(Ax + row_start, Ax + row_end, Bx);
        }
    }
}


/*
 * Slice columns given as an array of indices (pass 1).
 */
template<class I>
void csr_column_index1(const I n_idx,
                       const I col_idxs[],
                       const I n_row,
                       const I n_col,
                       const I Ap[],
                       const I Aj[],
                       I col_offsets[],
                       I Bp[])
{
    for(I jj = 0; jj < n_idx; jj++){
        const I j = col_idxs[jj];
        col_offsets[j]++;
    }

    I new_nnz = 0;
    Bp[0] = 0;
    for(I i = 0; i < n_row; i++){
        for(I jj = Ap[i]; jj < Ap[i+1]; jj++){
            new_nnz += col_offsets[Aj[jj]];
        }
        Bp[i+1] = new_nnz;
    }

    for(I j = 1; j < n_col; j++){
        col_offsets[j] += col_offsets[j - 1];
    }
}


/*
 * Slice columns given as an array of indices (pass 2).
 */
template<class I, class T>
void csr_column_index2(const I col_order[],
                       const I col_offsets[],
                       const I nnz,
                       const I Aj[],
                       const T Ax[],
                       I Bj[],
                       T Bx[])
{
    I n = 0;
    for(I jj = 0; jj < nnz; jj++){
        const I j = Aj[jj];
        const I offset = col_offsets[j];
        const I prev_offset = j == 0 ? 0 : col_offsets[j-1];
        if (offset != prev_offset) {
            const T v = Ax[jj];
            for(I k = prev_offset; k < offset; k++){
                Bj[n] = col_order[k];
                Bx[n] = v;
                n++;
            }
        }
    }
}


/*
 * Sample the matrix at specific locations
 */
template <class I, class T>
void csr_sample_values(const I n_row,
                       const I n_col,
                       const I Ap[],
                       const I Aj[],
                       const T Ax[],
                       const I n_samples,
                       const I Bi[],
                       const I Bj[],
                             T Bx[])
{
    const I nnz = Ap[n_row];
    const I threshold = nnz / 10;

    if (n_samples > threshold && csr_has_canonical_format(n_row, Ap, Aj))
    {
        for(I n = 0; n < n_samples; n++)
        {
            const I i = Bi[n] < 0 ? Bi[n] + n_row : Bi[n];
            const I j = Bj[n] < 0 ? Bj[n] + n_col : Bj[n];

            const I row_start = Ap[i];
            const I row_end   = Ap[i+1];

            if (row_start < row_end)
            {
                const I offset = std::lower_bound(Aj + row_start, Aj + row_end, j) - Aj;
                if (offset < row_end && Aj[offset] == j)
                    Bx[n] = Ax[offset];
                else
                    Bx[n] = 0;
            }
            else
            {
                Bx[n] = 0;
            }
        }
    }
    else
    {
        for(I n = 0; n < n_samples; n++)
        {
            const I i = Bi[n] < 0 ? Bi[n] + n_row : Bi[n];
            const I j = Bj[n] < 0 ? Bj[n] + n_col : Bj[n];

            const I row_start = Ap[i];
            const I row_end   = Ap[i+1];

            T x = 0;
            for(I jj = row_start; jj < row_end; jj++)
            {
                if (Aj[jj] == j)
                    x += Ax[jj];
            }
            Bx[n] = x;
        }
    }
}

/*
 * Determine the data offset at specific locations
 */
template <class I>
int csr_sample_offsets(const I n_row,
                       const I n_col,
                       const I Ap[],
                       const I Aj[],
                       const I n_samples,
                       const I Bi[],
                       const I Bj[],
                             I Bp[])
{
    const I nnz = Ap[n_row];
    const I threshold = nnz / 10;

    if (n_samples > threshold && csr_has_canonical_format(n_row, Ap, Aj))
    {
        for(I n = 0; n < n_samples; n++)
        {
            const I i = Bi[n] < 0 ? Bi[n] + n_row : Bi[n];
            const I j = Bj[n] < 0 ? Bj[n] + n_col : Bj[n];

            const I row_start = Ap[i];
            const I row_end   = Ap[i+1];

            if (row_start < row_end)
            {
                const I offset = std::lower_bound(Aj + row_start, Aj + row_end, j) - Aj;
                if (offset < row_end && Aj[offset] == j)
                    Bp[n] = offset;
                else
                    Bp[n] = -1;
            }
            else
            {
                Bp[n] = -1;
            }
        }
    }
    else
    {
        for(I n = 0; n < n_samples; n++)
        {
            const I i = Bi[n] < 0 ? Bi[n] + n_row : Bi[n];
            const I j = Bj[n] < 0 ? Bj[n] + n_col : Bj[n];

            const I row_start = Ap[i];
            const I row_end   = Ap[i+1];

            I offset = -1;
            for(I jj = row_start; jj < row_end; jj++)
            {
                if (Aj[jj] == j) {
                    offset = jj;
                    for (jj++; jj < row_end; jj++) {
                        if (Aj[jj] == j) {
                            offset = -2;
                            return 1;
                        }
                    }
                }
            }
            Bp[n] = offset;
        }
    }
    return 0;
}

#endif
