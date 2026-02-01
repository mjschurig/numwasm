/**
 * SQRTMV - Matrix Square Root Times Vector
 *
 * Computes y = A^(1/2) * v for symmetric positive definite matrices.
 *
 * Uses eigendecomposition: if A = V * D * V^T where D is diagonal,
 * then A^(1/2) = V * D^(1/2) * V^T
 *
 * For large matrices, uses a partial eigendecomposition via ARPACK
 * to approximate the result.
 */

import { eigs } from '../core/eigs.js';
import type { MatVecFunction } from '../high-level-types.js';

/**
 * Options for matrix square root computation.
 */
export interface SqrtmvOptions {
  /**
   * Number of eigenvalues/eigenvectors to use for approximation.
   * More eigenvalues give better accuracy but cost more.
   * @default 20
   */
  k?: number;

  /**
   * Convergence tolerance for eigenvalue computation.
   * @default 0 (machine precision)
   */
  tol?: number;

  /**
   * Maximum iterations for eigenvalue computation.
   */
  maxiter?: number;
}

/**
 * Result from matrix square root computation.
 */
export interface SqrtmvResult {
  /**
   * Result vector y = A^(1/2) * v.
   */
  result: Float64Array;

  /**
   * Number of eigenvalues used in the approximation.
   */
  k: number;

  /**
   * Whether the computation was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
}

/**
 * Compute y = A^(1/2) * v for symmetric positive definite A.
 *
 * Uses a partial eigendecomposition to approximate A^(1/2).
 * The approximation is: A^(1/2) ≈ Σ sqrt(λ_i) * v_i * v_i^T
 *
 * @param matvec - Symmetric positive definite matrix operator A*x
 * @param v - Input vector
 * @param options - Computation options
 * @returns Result containing A^(1/2)*v
 *
 * @example
 * ```ts
 * // Compute matrix square root times vector
 * const result = await sqrtmv(Amatvec, v);
 * if (result.success) {
 *   console.log('A^(1/2) * v:', result.result);
 * }
 * ```
 */
export async function sqrtmv(
  matvec: MatVecFunction,
  v: Float64Array,
  options?: SqrtmvOptions
): Promise<SqrtmvResult> {
  const {
    k = 20,
    tol = 0,
    maxiter,
  } = options ?? {};

  const n = v.length;
  const actualK = Math.min(k, n - 2);

  if (actualK < 1) {
    return {
      result: new Float64Array(n),
      k: 0,
      success: false,
      message: 'Matrix dimension too small for partial eigendecomposition',
    };
  }

  // Compute largest eigenvalues (these dominate the square root)
  const eigResult = await eigs(matvec, n, actualK, {
    which: 'LM',
    tol,
    maxiter,
    return_eigenvectors: true,
  });

  if (!eigResult.success || !eigResult.eigenvectors) {
    return {
      result: new Float64Array(n),
      k: 0,
      success: false,
      message: 'Eigendecomposition failed: ' + eigResult.message,
    };
  }

  // Check for negative eigenvalues (matrix not positive definite)
  for (let i = 0; i < eigResult.nconv; i++) {
    if (eigResult.eigenvalues[i] < -1e-10) {
      return {
        result: new Float64Array(n),
        k: eigResult.nconv,
        success: false,
        message: `Matrix is not positive definite (eigenvalue ${i} = ${eigResult.eigenvalues[i]})`,
      };
    }
  }

  // Compute y = Σ sqrt(λ_i) * (v_i^T * v) * v_i
  const result = new Float64Array(n);

  for (let i = 0; i < eigResult.nconv; i++) {
    const lambda = eigResult.eigenvalues[i];
    const sqrtLambda = Math.sqrt(Math.max(0, lambda));
    const vi = eigResult.eigenvectors[i];

    // Compute v_i^T * v
    let viTv = 0;
    for (let j = 0; j < n; j++) {
      viTv += vi[j] * v[j];
    }

    // Add sqrt(λ_i) * (v_i^T * v) * v_i to result
    const coeff = sqrtLambda * viTv;
    for (let j = 0; j < n; j++) {
      result[j] += coeff * vi[j];
    }
  }

  return {
    result,
    k: eigResult.nconv,
    success: true,
    message: `Approximation using ${eigResult.nconv} eigenvalues`,
  };
}

/**
 * Compute y = A^(-1/2) * v for symmetric positive definite A.
 *
 * Uses a partial eigendecomposition to approximate A^(-1/2).
 *
 * @param matvec - Symmetric positive definite matrix operator A*x
 * @param v - Input vector
 * @param options - Computation options
 * @returns Result containing A^(-1/2)*v
 */
export async function invsqrtmv(
  matvec: MatVecFunction,
  v: Float64Array,
  options?: SqrtmvOptions
): Promise<SqrtmvResult> {
  const {
    k = 20,
    tol = 0,
    maxiter,
  } = options ?? {};

  const n = v.length;
  const actualK = Math.min(k, n - 2);

  if (actualK < 1) {
    return {
      result: new Float64Array(n),
      k: 0,
      success: false,
      message: 'Matrix dimension too small for partial eigendecomposition',
    };
  }

  // For A^(-1/2), we need smallest eigenvalues (they become largest in inverse)
  // But we also need largest for accuracy
  // Compute both ends
  const eigResultLM = await eigs(matvec, n, Math.ceil(actualK / 2), {
    which: 'LM',
    tol,
    maxiter,
    return_eigenvectors: true,
  });

  const eigResultSM = await eigs(matvec, n, Math.floor(actualK / 2), {
    which: 'SM',
    tol,
    maxiter,
    return_eigenvectors: true,
  });

  if (!eigResultLM.success || !eigResultSM.success ||
      !eigResultLM.eigenvectors || !eigResultSM.eigenvectors) {
    return {
      result: new Float64Array(n),
      k: 0,
      success: false,
      message: 'Eigendecomposition failed',
    };
  }

  // Combine eigenvalues and eigenvectors
  const eigenvalues: number[] = [];
  const eigenvectors: Float64Array[] = [];

  for (let i = 0; i < eigResultLM.nconv; i++) {
    eigenvalues.push(eigResultLM.eigenvalues[i]);
    eigenvectors.push(eigResultLM.eigenvectors[i]);
  }

  for (let i = 0; i < eigResultSM.nconv; i++) {
    const lambda = eigResultSM.eigenvalues[i];
    // Avoid duplicates (eigenvalues that appear in both)
    const isDuplicate = eigenvalues.some(l => Math.abs(l - lambda) < 1e-10);
    if (!isDuplicate) {
      eigenvalues.push(lambda);
      eigenvectors.push(eigResultSM.eigenvectors[i]);
    }
  }

  // Check for non-positive eigenvalues
  for (let i = 0; i < eigenvalues.length; i++) {
    if (eigenvalues[i] < 1e-10) {
      return {
        result: new Float64Array(n),
        k: eigenvalues.length,
        success: false,
        message: `Matrix is not positive definite (eigenvalue = ${eigenvalues[i]})`,
      };
    }
  }

  // Compute y = Σ (1/sqrt(λ_i)) * (v_i^T * v) * v_i
  const result = new Float64Array(n);

  for (let i = 0; i < eigenvalues.length; i++) {
    const lambda = eigenvalues[i];
    const invSqrtLambda = 1 / Math.sqrt(lambda);
    const vi = eigenvectors[i];

    let viTv = 0;
    for (let j = 0; j < n; j++) {
      viTv += vi[j] * v[j];
    }

    const coeff = invSqrtLambda * viTv;
    for (let j = 0; j < n; j++) {
      result[j] += coeff * vi[j];
    }
  }

  return {
    result,
    k: eigenvalues.length,
    success: true,
    message: `Approximation using ${eigenvalues.length} eigenvalues`,
  };
}

/**
 * Compute y = A^p * v for symmetric positive definite A and any real power p.
 *
 * Generalizes sqrtmv (p=0.5) and invsqrtmv (p=-0.5).
 *
 * @param matvec - Symmetric positive definite matrix operator A*x
 * @param v - Input vector
 * @param p - Power to raise the matrix to
 * @param options - Computation options
 * @returns Result containing A^p*v
 */
export async function matpowv(
  matvec: MatVecFunction,
  v: Float64Array,
  p: number,
  options?: SqrtmvOptions
): Promise<SqrtmvResult> {
  const {
    k = 20,
    tol = 0,
    maxiter,
  } = options ?? {};

  const n = v.length;
  const actualK = Math.min(k, n - 2);

  // Determine which eigenvalues to compute
  const which = p >= 0 ? 'LM' : 'SM';

  const eigResult = await eigs(matvec, n, actualK, {
    which,
    tol,
    maxiter,
    return_eigenvectors: true,
  });

  if (!eigResult.success || !eigResult.eigenvectors) {
    return {
      result: new Float64Array(n),
      k: 0,
      success: false,
      message: 'Eigendecomposition failed: ' + eigResult.message,
    };
  }

  // Check for non-positive eigenvalues when p is not an integer
  if (p !== Math.floor(p)) {
    for (let i = 0; i < eigResult.nconv; i++) {
      if (eigResult.eigenvalues[i] < 1e-10) {
        return {
          result: new Float64Array(n),
          k: eigResult.nconv,
          success: false,
          message: `Matrix must be positive definite for non-integer powers`,
        };
      }
    }
  }

  // Compute y = Σ λ_i^p * (v_i^T * v) * v_i
  const result = new Float64Array(n);

  for (let i = 0; i < eigResult.nconv; i++) {
    const lambda = eigResult.eigenvalues[i];
    const lambdaP = Math.pow(Math.max(0, lambda), p);
    const vi = eigResult.eigenvectors[i];

    let viTv = 0;
    for (let j = 0; j < n; j++) {
      viTv += vi[j] * v[j];
    }

    const coeff = lambdaP * viTv;
    for (let j = 0; j < n; j++) {
      result[j] += coeff * vi[j];
    }
  }

  return {
    result,
    k: eigResult.nconv,
    success: true,
    message: `Approximation using ${eigResult.nconv} eigenvalues`,
  };
}
