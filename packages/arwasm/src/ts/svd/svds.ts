/**
 * SVDS - Partial Singular Value Decomposition
 *
 * Computes the k largest or smallest singular values and vectors
 * of a sparse matrix using ARPACK. This is the WASM equivalent of
 * scipy.sparse.linalg.svds.
 *
 * The SVD of an m×n matrix A is: A = U * diag(s) * Vt
 * where U (m×k), s (k), Vt (k×n) are the truncated SVD components.
 *
 * The algorithm works by computing eigenvalues of A^T*A (for n < m)
 * or A*A^T (for m <= n), since singular values are the square roots
 * of these eigenvalues.
 */

import { eigs } from '../core/eigs.js';
import type { MatVecFunction, SvdsOptions, SvdsResult } from '../high-level-types.js';

/**
 * Compute partial singular value decomposition of a sparse matrix.
 *
 * Finds the k largest or smallest singular values and optionally the
 * corresponding left and right singular vectors.
 *
 * @param matvec - Function computing y = A*x (m×n matrix applied to n-vector)
 * @param matvecT - Function computing y = A^T*x (n×m matrix applied to m-vector)
 * @param m - Number of rows in the matrix
 * @param n - Number of columns in the matrix
 * @param k - Number of singular values to compute
 * @param options - Solver options
 * @returns Partial SVD: U, s, Vt
 *
 * @example
 * ```ts
 * import { svds } from 'arwasm';
 *
 * // Define a sparse matrix via matvec
 * const m = 1000, n = 500;
 * const matvec = (x: Float64Array): Float64Array => {
 *   // Your matrix-vector product A*x
 *   const y = new Float64Array(m);
 *   // ... compute y = A*x
 *   return y;
 * };
 * const matvecT = (x: Float64Array): Float64Array => {
 *   // Your transpose matrix-vector product A^T*x
 *   const y = new Float64Array(n);
 *   // ... compute y = A^T*x
 *   return y;
 * };
 *
 * // Find 6 largest singular values
 * const { U, s, Vt } = await svds(matvec, matvecT, m, n, 6);
 * console.log('Singular values:', s);
 * ```
 */
export async function svds(
  matvec: MatVecFunction,
  matvecT: MatVecFunction,
  m: number,
  n: number,
  k: number,
  options?: SvdsOptions
): Promise<SvdsResult> {
  // Validate inputs
  if (m < 1 || n < 1) {
    throw new Error('Matrix dimensions must be positive');
  }
  if (k < 1) {
    throw new Error('Number of singular values k must be positive');
  }

  const minDim = Math.min(m, n);
  if (k >= minDim) {
    throw new Error(`k must be less than min(m, n) (got k=${k}, min(m,n)=${minDim})`);
  }

  const {
    which = 'LM',
    tol = 0,
    ncv: userNcv,
    maxiter = minDim * 10,
    v0,
    return_singular_vectors = true,
  } = options ?? {};

  // Determine whether to compute U, Vt, or both
  let computeU = false;
  let computeVt = false;
  if (return_singular_vectors === true || return_singular_vectors === 'both') {
    computeU = true;
    computeVt = true;
  } else if (return_singular_vectors === 'u') {
    computeU = true;
  } else if (return_singular_vectors === 'vh') {
    computeVt = true;
  }

  // Strategy: compute eigenvalues of A^T*A (if n <= m) or A*A^T (if m < n)
  // This gives us squared singular values. Singular vectors come from eigenvectors.
  //
  // If n <= m: compute eigenpairs of A^T*A (n×n), get right singular vectors V
  //            then U = A*V / s
  // If m < n:  compute eigenpairs of A*A^T (m×m), get left singular vectors U
  //            then V = A^T*U / s

  const useATA = n <= m; // true: work with A^T*A, false: work with A*A^T
  const eigDim = useATA ? n : m;

  // Helper to ensure Float64Array (matvec may return number[])
  const toFloat64 = (arr: Float64Array | number[]): Float64Array => {
    return arr instanceof Float64Array ? arr : new Float64Array(arr);
  };

  // Create the symmetric matvec for the eigenvalue problem
  // For A^T*A: y = A^T*(A*x)
  // For A*A^T: y = A*(A^T*x)
  const symmetricMatvec: MatVecFunction = useATA
    ? (x: Float64Array) => {
        const Ax = toFloat64(matvec(x)); // m-vector
        return matvecT(Ax); // n-vector
      }
    : (x: Float64Array) => {
        const ATx = toFloat64(matvecT(x)); // n-vector
        return matvec(ATx); // m-vector
      };

  // Map 'which' for singular values to eigenvalue selection
  // For largest singular values, we want largest eigenvalues of A^T*A
  // For smallest singular values, we want smallest positive eigenvalues
  const eigWhich = which === 'LM' ? 'LM' : 'SM';

  // Compute eigenvalues/eigenvectors
  const needVectors = computeU || computeVt;
  const result = await eigs(symmetricMatvec, eigDim, k, {
    which: eigWhich,
    tol,
    ncv: userNcv,
    maxiter,
    v0: v0 as Float64Array | undefined,
    return_eigenvectors: needVectors,
  });

  if (!result.success) {
    return {
      s: new Float64Array(0),
      nconv: 0,
      niter: result.niter,
      nops: result.nops,
      success: false,
      message: result.message,
    };
  }

  // Convert eigenvalues to singular values (take square root)
  // Handle potential numerical issues (small negative eigenvalues should be 0)
  const nconv = result.nconv;
  const singularValues = new Float64Array(nconv);
  for (let i = 0; i < nconv; i++) {
    const eigval = result.eigenvalues[i];
    singularValues[i] = eigval > 0 ? Math.sqrt(eigval) : 0;
  }

  // Sort singular values in descending order and get sort indices
  const indices = Array.from({ length: nconv }, (_, i) => i);
  indices.sort((a, b) => singularValues[b] - singularValues[a]);

  const sortedS = new Float64Array(nconv);
  for (let i = 0; i < nconv; i++) {
    sortedS[i] = singularValues[indices[i]];
  }

  // Compute singular vectors if requested
  let U: Float64Array[] | undefined;
  let Vt: Float64Array[] | undefined;

  if (needVectors && result.eigenvectors) {
    // Reorder eigenvectors according to sorted singular values
    const sortedEigvecs: Float64Array[] = [];
    for (let i = 0; i < nconv; i++) {
      sortedEigvecs.push(result.eigenvectors[indices[i]]);
    }

    if (useATA) {
      // Eigenvectors of A^T*A are right singular vectors V
      // V columns are the eigenvectors (each is n-dimensional)
      // Vt rows are the eigenvectors
      if (computeVt) {
        Vt = sortedEigvecs; // Each Vt[i] is an n-vector (i-th row of Vt)
      }

      // Compute U from V: u_i = A*v_i / s_i
      if (computeU) {
        U = [];
        for (let i = 0; i < nconv; i++) {
          const v = sortedEigvecs[i];
          const Av = matvec(v);
          const s = sortedS[i];
          if (s > 1e-14) {
            const u = new Float64Array(m);
            for (let j = 0; j < m; j++) {
              u[j] = Av[j] / s;
            }
            U.push(u);
          } else {
            // Zero singular value - U vector is arbitrary unit vector
            U.push(new Float64Array(m));
          }
        }
      }
    } else {
      // Eigenvectors of A*A^T are left singular vectors U
      if (computeU) {
        U = sortedEigvecs; // Each U[i] is an m-vector
      }

      // Compute Vt from U: v_i = A^T*u_i / s_i
      if (computeVt) {
        Vt = [];
        for (let i = 0; i < nconv; i++) {
          const u = sortedEigvecs[i];
          const ATu = matvecT(u);
          const s = sortedS[i];
          if (s > 1e-14) {
            const v = new Float64Array(n);
            for (let j = 0; j < n; j++) {
              v[j] = ATu[j] / s;
            }
            Vt.push(v);
          } else {
            // Zero singular value
            Vt.push(new Float64Array(n));
          }
        }
      }
    }
  }

  return {
    U,
    s: sortedS,
    Vt,
    nconv,
    niter: result.niter,
    nops: result.nops,
    success: true,
    message: 'Converged successfully',
  };
}
