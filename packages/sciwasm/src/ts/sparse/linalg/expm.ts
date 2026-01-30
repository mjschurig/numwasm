/**
 * Sparse matrix exponential
 *
 * Provides functions for computing matrix exponentials:
 * - expm: Full matrix exponential using Pade approximation
 * - expm_multiply: Action of matrix exponential on a vector using Krylov methods
 */

import type { SparseMatrix } from '../base.js';
import type {
  LinearOperatorLike,
  ExpmOptions,
  ExpmMultiplyOptions,
  ExpmMultiplyResult,
} from './types.js';
import { aslinearoperator } from './interface.js';
import { DimensionMismatchError } from '../../errors.js';

// ============================================================
// Pade Coefficients
// ============================================================

// Pade coefficients for order 13 (from scipy)
const PADE_COEFFS_13 = [
  64764752532480000,
  32382376266240000,
  7771770303897600,
  1187353796428800,
  129060195264000,
  10559470521600,
  670442572800,
  33522128640,
  1323241920,
  40840800,
  960960,
  16380,
  182,
  1,
];

// Theta values for scaling (from scipy)
const THETA_13 = 5.371920351148152;

// ============================================================
// Dense Matrix Operations (for internal use)
// ============================================================

function denseZeros(n: number): number[][] {
  const result: number[][] = [];
  for (let i = 0; i < n; i++) {
    result.push(new Array(n).fill(0));
  }
  return result;
}

function denseIdentity(n: number): number[][] {
  const result = denseZeros(n);
  for (let i = 0; i < n; i++) {
    result[i][i] = 1;
  }
  return result;
}

function denseAdd(A: number[][], B: number[][]): number[][] {
  const n = A.length;
  const result = denseZeros(n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      result[i][j] = A[i][j] + B[i][j];
    }
  }
  return result;
}

function denseSub(A: number[][], B: number[][]): number[][] {
  const n = A.length;
  const result = denseZeros(n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      result[i][j] = A[i][j] - B[i][j];
    }
  }
  return result;
}

function denseScale(s: number, A: number[][]): number[][] {
  const n = A.length;
  const result = denseZeros(n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < n; j++) {
      result[i][j] = s * A[i][j];
    }
  }
  return result;
}

function denseMatmul(A: number[][], B: number[][]): number[][] {
  const n = A.length;
  const m = B[0].length;
  const k = A[0].length;
  const result = denseZeros(n);
  for (let i = 0; i < n; i++) {
    for (let j = 0; j < m; j++) {
      let sum = 0;
      for (let l = 0; l < k; l++) {
        sum += A[i][l] * B[l][j];
      }
      result[i][j] = sum;
    }
  }
  return result;
}

function denseOneNorm(A: number[][]): number {
  const n = A.length;
  const m = A[0].length;
  let maxSum = 0;
  for (let j = 0; j < m; j++) {
    let colSum = 0;
    for (let i = 0; i < n; i++) {
      colSum += Math.abs(A[i][j]);
    }
    maxSum = Math.max(maxSum, colSum);
  }
  return maxSum;
}

/**
 * Solve dense linear system A @ X = B using LU decomposition with partial pivoting
 * Returns X such that A @ X = B
 */
function denseLUSolve(A: number[][], B: number[][]): number[][] {
  const n = A.length;
  const nrhs = B[0].length;

  // Create copies
  const LU: number[][] = A.map(row => [...row]);
  const X: number[][] = B.map(row => [...row]);
  const perm: number[] = Array.from({ length: n }, (_, i) => i);

  // LU decomposition with partial pivoting
  for (let k = 0; k < n - 1; k++) {
    // Find pivot
    let maxVal = Math.abs(LU[k][k]);
    let maxIdx = k;
    for (let i = k + 1; i < n; i++) {
      if (Math.abs(LU[i][k]) > maxVal) {
        maxVal = Math.abs(LU[i][k]);
        maxIdx = i;
      }
    }

    // Swap rows if necessary
    if (maxIdx !== k) {
      [LU[k], LU[maxIdx]] = [LU[maxIdx], LU[k]];
      [X[k], X[maxIdx]] = [X[maxIdx], X[k]];
      [perm[k], perm[maxIdx]] = [perm[maxIdx], perm[k]];
    }

    // Check for singularity
    if (Math.abs(LU[k][k]) < 1e-15) {
      throw new Error('Matrix is singular or nearly singular');
    }

    // Elimination
    for (let i = k + 1; i < n; i++) {
      const factor = LU[i][k] / LU[k][k];
      LU[i][k] = factor;  // Store L factor
      for (let j = k + 1; j < n; j++) {
        LU[i][j] -= factor * LU[k][j];
      }
      for (let j = 0; j < nrhs; j++) {
        X[i][j] -= factor * X[k][j];
      }
    }
  }

  // Back substitution
  for (let j = 0; j < nrhs; j++) {
    for (let i = n - 1; i >= 0; i--) {
      for (let k = i + 1; k < n; k++) {
        X[i][j] -= LU[i][k] * X[k][j];
      }
      X[i][j] /= LU[i][i];
    }
  }

  return X;
}

// ============================================================
// Pade Approximation
// ============================================================

/**
 * Compute Pade approximation of exp(A) where ||A|| <= theta
 */
function padeApproximation13(A: number[][]): number[][] {
  const n = A.length;
  const I = denseIdentity(n);

  // Compute powers of A
  const A2 = denseMatmul(A, A);
  const A4 = denseMatmul(A2, A2);
  const A6 = denseMatmul(A4, A2);

  // Compute U and V for Pade approximant
  // U = A @ (A6 @ (b13*A6 + b11*A4 + b9*A2) + b7*A6 + b5*A4 + b3*A2 + b1*I)
  // V = A6 @ (b12*A6 + b10*A4 + b8*A2) + b6*A6 + b4*A4 + b2*A2 + b0*I

  const b = PADE_COEFFS_13;

  // Inner terms for U
  let U_inner = denseAdd(
    denseScale(b[13], A6),
    denseAdd(denseScale(b[11], A4), denseScale(b[9], A2))
  );
  U_inner = denseMatmul(A6, U_inner);
  U_inner = denseAdd(U_inner, denseScale(b[7], A6));
  U_inner = denseAdd(U_inner, denseScale(b[5], A4));
  U_inner = denseAdd(U_inner, denseScale(b[3], A2));
  U_inner = denseAdd(U_inner, denseScale(b[1], I));
  const U = denseMatmul(A, U_inner);

  // Terms for V
  let V_inner = denseAdd(
    denseScale(b[12], A6),
    denseAdd(denseScale(b[10], A4), denseScale(b[8], A2))
  );
  V_inner = denseMatmul(A6, V_inner);
  let V = denseAdd(V_inner, denseScale(b[6], A6));
  V = denseAdd(V, denseScale(b[4], A4));
  V = denseAdd(V, denseScale(b[2], A2));
  V = denseAdd(V, denseScale(b[0], I));

  // Solve (V - U) @ X = (V + U) for X = exp(A)
  const VmU = denseSub(V, U);
  const VpU = denseAdd(V, U);
  return denseLUSolve(VmU, VpU);
}

/**
 * Square a matrix s times: A -> A^(2^s)
 */
function squareMatrix(A: number[][], s: number): number[][] {
  let result = A;
  for (let i = 0; i < s; i++) {
    result = denseMatmul(result, result);
  }
  return result;
}

// ============================================================
// expm - Full Matrix Exponential
// ============================================================

/**
 * Compute the matrix exponential of a sparse matrix.
 *
 * Uses the Pade approximation with scaling and squaring method.
 * Note: The result is always a dense 2D array since exp(A) is generally dense.
 *
 * @param A - Sparse square matrix
 * @param options - Computation options
 * @returns Dense matrix representing exp(A)
 *
 * @example
 * ```typescript
 * const A = csr_matrix([
 *   [0, -1],
 *   [1, 0]
 * ]);
 * const expA = expm(A);
 * // expA is [[cos(1), -sin(1)], [sin(1), cos(1)]]
 * ```
 */
export function expm(
  A: SparseMatrix | LinearOperatorLike,
  _options?: ExpmOptions
): number[][] {
  // Convert to dense for the computation
  let Adense: number[][];

  if ('toArray' in A && typeof A.toArray === 'function') {
    Adense = (A as SparseMatrix).toArray();
  } else {
    throw new Error('Input must be a sparse matrix with toArray method');
  }

  const [m, n] = [Adense.length, Adense[0]?.length ?? 0];

  if (m !== n) {
    throw new DimensionMismatchError(`Matrix must be square, got ${m}x${n}`);
  }

  if (m === 0) {
    return [];
  }

  // Compute 1-norm for scaling
  const normA = denseOneNorm(Adense);

  if (normA === 0) {
    // exp(0) = I
    return denseIdentity(n);
  }

  // Determine scaling factor
  const s = Math.max(0, Math.ceil(Math.log2(normA / THETA_13)));

  // Scale the matrix
  const scaleFactor = Math.pow(2, -s);
  const Ascaled = denseScale(scaleFactor, Adense);

  // Compute Pade approximation
  let result = padeApproximation13(Ascaled);

  // Square s times to undo scaling
  result = squareMatrix(result, s);

  return result;
}

// ============================================================
// Krylov/Arnoldi Methods for expm_multiply
// ============================================================

/**
 * Arnoldi iteration to build orthonormal Krylov basis
 * Returns V (n x m+1) and H ((m+1) x m) such that A @ V_m = V_{m+1} @ H
 */
function arnoldiIteration(
  matvec: (x: Float64Array) => Float64Array,
  v0: Float64Array,
  m: number
): { V: Float64Array[]; H: number[][] } {
  const n = v0.length;

  // Initialize
  const V: Float64Array[] = [];
  const H: number[][] = [];
  for (let i = 0; i <= m; i++) {
    H.push(new Array(m).fill(0));
  }

  // Normalize starting vector
  let beta = Math.sqrt(v0.reduce((s, x) => s + x * x, 0));
  if (beta < 1e-15) {
    throw new Error('Starting vector is zero');
  }

  const v = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    v[i] = v0[i] / beta;
  }
  V.push(v);

  for (let j = 0; j < m; j++) {
    // w = A @ v_j
    let w = matvec(V[j]);

    // Orthogonalize against previous vectors (modified Gram-Schmidt)
    for (let i = 0; i <= j; i++) {
      let h = 0;
      for (let k = 0; k < n; k++) {
        h += V[i][k] * w[k];
      }
      H[i][j] = h;
      for (let k = 0; k < n; k++) {
        w[k] -= h * V[i][k];
      }
    }

    // Compute norm of w
    const wnorm = Math.sqrt(w.reduce((s, x) => s + x * x, 0));
    H[j + 1][j] = wnorm;

    // Check for happy breakdown
    if (wnorm < 1e-12) {
      // Happy breakdown: we found an invariant subspace
      break;
    }

    // Normalize and store
    const vNew = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      vNew[i] = w[i] / wnorm;
    }
    V.push(vNew);
  }

  return { V, H };
}

/**
 * Compute exp(t * H) for a small Hessenberg matrix
 * Uses direct Pade approximation (small matrix, so no scaling needed usually)
 */
function smallExpm(H: number[][], t: number): number[][] {
  const m = H.length;
  if (m === 0) return [];

  // Scale H by t
  const Ht: number[][] = [];
  for (let i = 0; i < m; i++) {
    Ht.push([]);
    for (let j = 0; j < m; j++) {
      Ht[i].push(t * (H[i]?.[j] ?? 0));
    }
  }

  // For small matrices, use Pade if the norm is reasonable
  const normHt = denseOneNorm(Ht);
  if (normHt < 10) {
    // Use direct Pade
    const s = Math.max(0, Math.ceil(Math.log2(normHt / THETA_13)));
    const scaleFactor = Math.pow(2, -s);
    const Hscaled = denseScale(scaleFactor, Ht);
    let result = padeApproximation13(Hscaled);
    return squareMatrix(result, s);
  } else {
    // Use Taylor series for very large norms (rare)
    let result = denseIdentity(m);
    let term = denseIdentity(m);
    for (let k = 1; k <= 30; k++) {
      term = denseScale(1 / k, denseMatmul(term, Ht));
      result = denseAdd(result, term);
      if (denseOneNorm(term) < 1e-15 * denseOneNorm(result)) {
        break;
      }
    }
    return result;
  }
}

// ============================================================
// expm_multiply - Action of Matrix Exponential on a Vector
// ============================================================

/**
 * Compute the action of the matrix exponential on a vector: exp(t*A) * v
 *
 * Uses the Krylov subspace method, which is more efficient than computing
 * the full matrix exponential for large sparse matrices.
 *
 * @param A - Sparse matrix or LinearOperator
 * @param v - Vector to multiply
 * @param options - Computation options
 * @returns Result containing exp(t*A) * v and convergence info
 *
 * @example
 * ```typescript
 * const A = csr_matrix([
 *   [0, -1],
 *   [1, 0]
 * ]);
 * const v = new Float64Array([1, 0]);
 * const result = expm_multiply(A, v, { t: Math.PI / 2 });
 * // result.result is approximately [0, 1]
 * ```
 */
export function expm_multiply(
  A: SparseMatrix | LinearOperatorLike,
  v: Float64Array,
  options?: ExpmMultiplyOptions
): ExpmMultiplyResult {
  const Aop = aslinearoperator(A);
  const [m, n] = Aop.shape;

  if (m !== n) {
    throw new DimensionMismatchError(`Matrix must be square, got ${m}x${n}`);
  }

  if (v.length !== n) {
    throw new DimensionMismatchError(
      `Vector length ${v.length} doesn't match matrix size ${n}`
    );
  }

  const t = options?.t ?? 1.0;
  const krylovDim = options?.krylovDim ?? Math.min(30, n);

  // Handle zero time
  if (Math.abs(t) < 1e-15) {
    return {
      result: new Float64Array(v),
      iterations: 0,
      errorEstimate: 0,
    };
  }

  // Compute ||v||
  const vNorm = Math.sqrt(v.reduce((s, x) => s + x * x, 0));
  if (vNorm < 1e-15) {
    return {
      result: new Float64Array(n),
      iterations: 0,
      errorEstimate: 0,
    };
  }

  // Run Arnoldi iteration
  const { V, H } = arnoldiIteration(Aop.matvec.bind(Aop), v, krylovDim);
  const actualDim = V.length - 1;

  // Truncate H to actual dimension
  const Hm: number[][] = [];
  for (let i = 0; i <= actualDim; i++) {
    Hm.push(H[i].slice(0, actualDim));
  }

  // Compute exp(t * H_m) for the (m x m) upper part
  const HmSquare: number[][] = Hm.slice(0, actualDim);
  const expHm = smallExpm(HmSquare, t);

  // Result = ||v|| * V_m @ exp(t*H_m) @ e_1
  // where e_1 = [1, 0, 0, ...]
  // So we need the first column of exp(t*H_m)
  const expHe1 = new Float64Array(actualDim);
  for (let i = 0; i < actualDim; i++) {
    expHe1[i] = expHm[i][0];
  }

  // Compute result = V_m @ (||v|| * expHe1)
  const result = new Float64Array(n);
  for (let j = 0; j < actualDim; j++) {
    const coeff = vNorm * expHe1[j];
    for (let i = 0; i < n; i++) {
      result[i] += coeff * V[j][i];
    }
  }

  // Error estimate: ||v|| * |h_{m+1,m}| * |e_m^T @ exp(t*H_m) @ e_1|
  const h_mp1_m = actualDim < krylovDim ? H[actualDim][actualDim - 1] : 0;
  const exp_last = actualDim > 0 ? expHm[actualDim - 1][0] : 0;
  const errorEstimate = Math.abs(vNorm * h_mp1_m * exp_last);

  return {
    result,
    iterations: actualDim,
    errorEstimate,
  };
}
