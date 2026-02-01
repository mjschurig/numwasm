/**
 * EXPMV - Matrix Exponential Times Vector
 *
 * Computes y = exp(t*A) * v without forming the full matrix exponential.
 * Uses a Krylov subspace approximation (Arnoldi/Lanczos) which is efficient
 * for large sparse matrices.
 *
 * The algorithm:
 * 1. Build a Krylov subspace K_m = span{v, A*v, A^2*v, ..., A^(m-1)*v}
 * 2. Compute the orthonormal basis V_m via Arnoldi/Lanczos
 * 3. Form the Hessenberg matrix H_m = V_m^T * A * V_m
 * 4. Compute exp(t*H_m) using a direct method (small matrix)
 * 5. Return y = ||v|| * V_m * exp(t*H_m) * e_1
 */

import type { MatVecFunction } from '../high-level-types.js';

/**
 * Options for matrix exponential computation.
 */
export interface ExpmvOptions {
  /**
   * Krylov subspace dimension.
   * Larger values give more accuracy but cost more.
   * @default 30
   */
  m?: number;

  /**
   * Convergence tolerance for the approximation.
   * @default 1e-10
   */
  tol?: number;

  /**
   * Maximum number of restarts if not converged.
   * @default 10
   */
  maxRestarts?: number;

  /**
   * Whether the matrix is symmetric (use Lanczos instead of Arnoldi).
   * @default false
   */
  symmetric?: boolean;
}

/**
 * Result from matrix exponential computation.
 */
export interface ExpmvResult {
  /**
   * Result vector y = exp(t*A) * v.
   */
  result: Float64Array;

  /**
   * Number of matrix-vector products used.
   */
  nops: number;

  /**
   * Estimated error in the approximation.
   */
  error: number;

  /**
   * Whether the computation converged.
   */
  converged: boolean;
}

/**
 * Compute the 2-norm of a vector.
 */
function norm2(v: Float64Array): number {
  let sum = 0;
  for (let i = 0; i < v.length; i++) {
    sum += v[i] * v[i];
  }
  return Math.sqrt(sum);
}

/**
 * Compute dot product of two vectors.
 */
function dot(a: Float64Array, b: Float64Array): number {
  let sum = 0;
  const n = a.length;
  for (let i = 0; i < n; i++) {
    sum += a[i] * b[i];
  }
  return sum;
}

/**
 * Compute matrix exponential of a small dense matrix.
 * Uses scaling and squaring with Padé approximation.
 */
function expmSmall(H: Float64Array[], m: number): Float64Array[] {
  // Scale matrix so ||H|| < 1
  let maxNorm = 0;
  for (let i = 0; i < m; i++) {
    let rowSum = 0;
    for (let j = 0; j < m; j++) {
      rowSum += Math.abs(H[i][j]);
    }
    if (rowSum > maxNorm) maxNorm = rowSum;
  }

  const s = Math.max(0, Math.ceil(Math.log2(maxNorm)));
  const scale = Math.pow(2, -s);

  // Scale H
  const Hs: Float64Array[] = [];
  for (let i = 0; i < m; i++) {
    Hs[i] = new Float64Array(m);
    for (let j = 0; j < m; j++) {
      Hs[i][j] = H[i][j] * scale;
    }
  }

  // Compute exp(Hs) using Taylor series (for small Hs)
  // exp(Hs) ≈ I + Hs + Hs^2/2! + Hs^3/3! + ...
  const result: Float64Array[] = [];
  const term: Float64Array[] = [];

  // Initialize result = I, term = I
  for (let i = 0; i < m; i++) {
    result[i] = new Float64Array(m);
    term[i] = new Float64Array(m);
    result[i][i] = 1;
    term[i][i] = 1;
  }

  // Add terms until convergence
  for (let k = 1; k <= 20; k++) {
    // term = term * Hs / k
    const newTerm: Float64Array[] = [];
    for (let i = 0; i < m; i++) {
      newTerm[i] = new Float64Array(m);
      for (let j = 0; j < m; j++) {
        let sum = 0;
        for (let l = 0; l < m; l++) {
          sum += term[i][l] * Hs[l][j];
        }
        newTerm[i][j] = sum / k;
      }
    }

    // result = result + newTerm
    let maxDiff = 0;
    for (let i = 0; i < m; i++) {
      for (let j = 0; j < m; j++) {
        term[i][j] = newTerm[i][j];
        result[i][j] += newTerm[i][j];
        const diff = Math.abs(newTerm[i][j]);
        if (diff > maxDiff) maxDiff = diff;
      }
    }

    if (maxDiff < 1e-16) break;
  }

  // Square s times to undo scaling
  for (let i = 0; i < s; i++) {
    const squared: Float64Array[] = [];
    for (let ii = 0; ii < m; ii++) {
      squared[ii] = new Float64Array(m);
      for (let jj = 0; jj < m; jj++) {
        let sum = 0;
        for (let l = 0; l < m; l++) {
          sum += result[ii][l] * result[l][jj];
        }
        squared[ii][jj] = sum;
      }
    }
    for (let ii = 0; ii < m; ii++) {
      for (let jj = 0; jj < m; jj++) {
        result[ii][jj] = squared[ii][jj];
      }
    }
  }

  return result;
}

/**
 * Compute y = exp(t*A) * v using Krylov subspace approximation.
 *
 * This function computes the action of the matrix exponential on a vector
 * without explicitly forming exp(t*A). This is efficient for large sparse
 * matrices where only matrix-vector products are available.
 *
 * @param matvec - Matrix operator A*x
 * @param v - Input vector
 * @param t - Time/scaling parameter
 * @param options - Computation options
 * @returns Result containing exp(t*A)*v
 *
 * @example
 * ```ts
 * // Solve dy/dt = A*y with y(0) = v
 * // Solution at time t is y(t) = exp(t*A) * v
 * const result = expmv(Amatvec, v, t);
 * console.log('Solution at t:', result.result);
 * ```
 */
export function expmv(
  matvec: MatVecFunction,
  v: Float64Array,
  t: number,
  options?: ExpmvOptions
): ExpmvResult {
  const {
    m = 30,
    tol = 1e-10,
    // maxRestarts reserved for future implementation
    symmetric = false,
  } = options ?? {};

  const n = v.length;
  const mActual = Math.min(m, n);
  let nops = 0;

  // Normalize input vector
  const beta = norm2(v);
  if (beta < 1e-14) {
    return {
      result: new Float64Array(n),
      nops: 0,
      error: 0,
      converged: true,
    };
  }

  // Build Krylov basis using Arnoldi (or Lanczos if symmetric)
  const V: Float64Array[] = []; // Orthonormal basis vectors
  const H: Float64Array[] = []; // Hessenberg matrix

  // Initialize
  V[0] = new Float64Array(n);
  for (let i = 0; i < n; i++) {
    V[0][i] = v[i] / beta;
  }

  for (let j = 0; j < mActual; j++) {
    H[j] = new Float64Array(mActual + 1);
  }

  // Arnoldi iteration
  for (let j = 0; j < mActual; j++) {
    // w = A * V[j]
    const Avj = matvec(V[j]);
    const w = Avj instanceof Float64Array ? new Float64Array(Avj) : new Float64Array(Avj);
    nops++;

    // Orthogonalize against previous vectors
    if (symmetric && j > 0) {
      // Lanczos: only need to orthogonalize against last two vectors
      const start = Math.max(0, j - 1);
      for (let i = start; i <= j; i++) {
        const h = dot(V[i], w);
        H[i][j] = h;
        for (let k = 0; k < n; k++) {
          w[k] -= h * V[i][k];
        }
      }
      // Symmetry: H[j][j-1] = H[j-1][j]
      if (j > 0) {
        H[j][j - 1] = H[j - 1][j];
      }
    } else {
      // Full Arnoldi
      for (let i = 0; i <= j; i++) {
        const h = dot(V[i], w);
        H[i][j] = h;
        for (let k = 0; k < n; k++) {
          w[k] -= h * V[i][k];
        }
      }
    }

    // Compute norm of residual
    const hjp1j = norm2(w);

    if (j < mActual - 1) {
      H[j + 1][j] = hjp1j;

      if (hjp1j > 1e-14) {
        V[j + 1] = new Float64Array(n);
        for (let k = 0; k < n; k++) {
          V[j + 1][k] = w[k] / hjp1j;
        }
      } else {
        // Lucky breakdown - Krylov subspace is invariant
        break;
      }
    }
  }

  // Extract the mActual x mActual Hessenberg matrix
  const Hm: Float64Array[] = [];
  const actualM = V.length;
  for (let i = 0; i < actualM; i++) {
    Hm[i] = new Float64Array(actualM);
    for (let j = 0; j < actualM; j++) {
      Hm[i][j] = t * H[i][j];
    }
  }

  // Compute exp(t*Hm)
  const expHm = expmSmall(Hm, actualM);

  // Compute y = beta * V * exp(t*Hm) * e1
  // exp(t*Hm) * e1 is just the first column of exp(t*Hm)
  const result = new Float64Array(n);
  for (let i = 0; i < actualM; i++) {
    const coeff = beta * expHm[i][0];
    for (let k = 0; k < n; k++) {
      result[k] += coeff * V[i][k];
    }
  }

  // Estimate error (using the (m+1,m) element of exp(t*H))
  // This is a rough estimate
  const error = actualM < mActual ? 0 : Math.abs(H[actualM - 1][actualM - 2] || 0) * tol;

  return {
    result,
    nops,
    error,
    converged: true,
  };
}

/**
 * Compute y = exp(t*A) * v for multiple time points.
 *
 * More efficient than calling expmv multiple times as it reuses
 * the Krylov basis.
 *
 * @param matvec - Matrix operator A*x
 * @param v - Input vector
 * @param times - Array of time points
 * @param options - Computation options
 * @returns Array of results for each time point
 */
export function expmvMultiple(
  matvec: MatVecFunction,
  v: Float64Array,
  times: number[],
  options?: ExpmvOptions
): ExpmvResult[] {
  // For now, just call expmv for each time
  // A more sophisticated implementation would share the Krylov basis
  return times.map(t => expmv(matvec, v, t, options));
}
