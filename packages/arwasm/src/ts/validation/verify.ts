/**
 * Validation and Verification Functions
 *
 * Utilities for verifying eigenvalue/eigenvector results and
 * checking matrix properties.
 */

import type { MatVecFunction } from '../high-level-types.js';

/**
 * Options for verification functions.
 */
export interface VerifyOptions {
  /**
   * Tolerance for considering values as valid.
   * @default 1e-10
   */
  tol?: number;
}

/**
 * Result from eigenvalue verification.
 */
export interface VerifyEigsResult {
  /**
   * Residual norms ||A*v - λ*v|| for each eigenpair.
   */
  residuals: Float64Array;

  /**
   * Maximum residual norm across all eigenpairs.
   */
  maxResidual: number;

  /**
   * Relative residuals ||A*v - λ*v|| / |λ| for each eigenpair.
   */
  relativeResiduals: Float64Array;

  /**
   * Maximum relative residual.
   */
  maxRelativeResidual: number;

  /**
   * Maximum off-diagonal entry in V^T * V (measures orthogonality).
   * 0 = perfectly orthogonal, larger = less orthogonal.
   */
  orthogonality: number;

  /**
   * Whether the eigenpairs pass validation (maxResidual < tol).
   */
  isValid: boolean;
}

/**
 * Result from SVD verification.
 */
export interface VerifySvdsResult {
  /**
   * Residual norms ||A*v - σ*u|| for each singular triplet.
   */
  residualsAV: Float64Array;

  /**
   * Residual norms ||A^T*u - σ*v|| for each singular triplet.
   */
  residualsATU: Float64Array;

  /**
   * Maximum residual across all checks.
   */
  maxResidual: number;

  /**
   * Orthogonality of U vectors.
   */
  orthogonalityU: number;

  /**
   * Orthogonality of V vectors.
   */
  orthogonalityV: number;

  /**
   * Whether the SVD passes validation.
   */
  isValid: boolean;
}

/**
 * Result from symmetry check.
 */
export interface SymmetryCheckResult {
  /**
   * Whether the matrix appears symmetric within tolerance.
   */
  isSymmetric: boolean;

  /**
   * Maximum asymmetry found: max |x^T*A*y - y^T*A*x|.
   */
  maxAsymmetry: number;

  /**
   * Number of probe pairs used.
   */
  nProbes: number;
}

/**
 * Result from positive definiteness check.
 */
export interface PositiveDefiniteCheckResult {
  /**
   * Whether the matrix appears positive definite.
   */
  isPositiveDefinite: boolean;

  /**
   * Smallest eigenvalue found (negative if not PD).
   */
  smallestEigenvalue: number;

  /**
   * Whether the check was successful.
   */
  success: boolean;

  /**
   * Status message.
   */
  message: string;
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
 * Compute the dot product of two vectors.
 */
function dot(a: Float64Array, b: Float64Array): number {
  let sum = 0;
  const n = Math.min(a.length, b.length);
  for (let i = 0; i < n; i++) {
    sum += a[i] * b[i];
  }
  return sum;
}

/**
 * Verify eigenvalue/eigenvector pairs.
 *
 * For each eigenpair (λ, v), computes the residual ||A*v - λ*v||
 * and checks orthogonality of eigenvectors.
 *
 * @param matvec - Matrix operator A
 * @param eigenvalues - Computed eigenvalues
 * @param eigenvectors - Computed eigenvectors
 * @param options - Verification options
 * @returns Verification results
 *
 * @example
 * ```ts
 * const result = await eigs(matvec, n, 6);
 * const verify = verifyEigs(matvec, result.eigenvalues, result.eigenvectors!);
 * if (!verify.isValid) {
 *   console.warn('Eigenvalues may not be accurate:', verify.maxResidual);
 * }
 * ```
 */
export function verifyEigs(
  matvec: MatVecFunction,
  eigenvalues: Float64Array,
  eigenvectors: Float64Array[],
  options?: VerifyOptions
): VerifyEigsResult {
  const { tol = 1e-10 } = options ?? {};
  const k = eigenvalues.length;

  if (eigenvectors.length !== k) {
    throw new Error('eigenvalues and eigenvectors must have same length');
  }

  const residuals = new Float64Array(k);
  const relativeResiduals = new Float64Array(k);

  // Compute residuals for each eigenpair
  for (let i = 0; i < k; i++) {
    const lambda = eigenvalues[i];
    const v = eigenvectors[i];
    const n = v.length;

    // Compute A*v
    const Av = matvec(v);
    const AvArr = Av instanceof Float64Array ? Av : new Float64Array(Av);

    // Compute residual r = A*v - λ*v
    const r = new Float64Array(n);
    for (let j = 0; j < n; j++) {
      r[j] = AvArr[j] - lambda * v[j];
    }

    residuals[i] = norm2(r);
    relativeResiduals[i] = Math.abs(lambda) > 1e-14
      ? residuals[i] / Math.abs(lambda)
      : residuals[i];
  }

  // Compute orthogonality: max |v_i · v_j| for i ≠ j
  let orthogonality = 0;
  for (let i = 0; i < k; i++) {
    for (let j = i + 1; j < k; j++) {
      const d = Math.abs(dot(eigenvectors[i], eigenvectors[j]));
      if (d > orthogonality) {
        orthogonality = d;
      }
    }
  }

  const maxResidual = Math.max(...residuals);
  const maxRelativeResidual = Math.max(...relativeResiduals);

  return {
    residuals,
    maxResidual,
    relativeResiduals,
    maxRelativeResidual,
    orthogonality,
    isValid: maxResidual < tol,
  };
}

/**
 * Verify non-symmetric eigenvalue/eigenvector pairs.
 *
 * For non-symmetric matrices, eigenvalues may be complex.
 * This verifies pairs where eigenvaluesImag[i] != 0 indicates
 * complex conjugate pairs.
 *
 * @param matvec - Matrix operator A
 * @param eigenvaluesReal - Real parts of eigenvalues
 * @param eigenvaluesImag - Imaginary parts of eigenvalues
 * @param eigenvectors - Eigenvectors (see EignResult for format)
 * @param options - Verification options
 * @returns Verification results
 */
export function verifyEign(
  matvec: MatVecFunction,
  eigenvaluesReal: Float64Array,
  eigenvaluesImag: Float64Array,
  eigenvectors: Float64Array[],
  options?: VerifyOptions
): VerifyEigsResult {
  const { tol = 1e-10 } = options ?? {};
  const k = eigenvaluesReal.length;

  const residuals = new Float64Array(k);
  const relativeResiduals = new Float64Array(k);

  let i = 0;
  while (i < k) {
    const lambdaR = eigenvaluesReal[i];
    const lambdaI = eigenvaluesImag[i];

    if (Math.abs(lambdaI) < 1e-14) {
      // Real eigenvalue
      const v = eigenvectors[i];
      const n = v.length;
      const Av = matvec(v);
      const AvArr = Av instanceof Float64Array ? Av : new Float64Array(Av);

      const r = new Float64Array(n);
      for (let j = 0; j < n; j++) {
        r[j] = AvArr[j] - lambdaR * v[j];
      }

      residuals[i] = norm2(r);
      const lambdaMag = Math.abs(lambdaR);
      relativeResiduals[i] = lambdaMag > 1e-14 ? residuals[i] / lambdaMag : residuals[i];
      i++;
    } else {
      // Complex conjugate pair
      // v_real = eigenvectors[i], v_imag = eigenvectors[i+1]
      // Full eigenvector: v = v_real + j*v_imag
      // A*v = lambda*v where lambda = lambdaR + j*lambdaI
      const vR = eigenvectors[i];
      const vI = eigenvectors[i + 1];
      const n = vR.length;

      const AvR = matvec(vR);
      const AvI = matvec(vI);
      const AvRArr = AvR instanceof Float64Array ? AvR : new Float64Array(AvR);
      const AvIArr = AvI instanceof Float64Array ? AvI : new Float64Array(AvI);

      // Residual for first eigenvalue (lambdaR + j*lambdaI):
      // A*(vR + j*vI) - (lambdaR + j*lambdaI)*(vR + j*vI)
      // Real part: A*vR - lambdaR*vR + lambdaI*vI
      // Imag part: A*vI - lambdaR*vI - lambdaI*vR
      const rR = new Float64Array(n);
      const rI = new Float64Array(n);
      for (let j = 0; j < n; j++) {
        rR[j] = AvRArr[j] - lambdaR * vR[j] + lambdaI * vI[j];
        rI[j] = AvIArr[j] - lambdaR * vI[j] - lambdaI * vR[j];
      }

      const res = Math.sqrt(norm2(rR) ** 2 + norm2(rI) ** 2);
      const lambdaMag = Math.sqrt(lambdaR * lambdaR + lambdaI * lambdaI);

      residuals[i] = res;
      residuals[i + 1] = res; // Same for conjugate pair
      relativeResiduals[i] = lambdaMag > 1e-14 ? res / lambdaMag : res;
      relativeResiduals[i + 1] = relativeResiduals[i];

      i += 2;
    }
  }

  // Orthogonality check (only for real eigenvectors)
  let orthogonality = 0;
  for (let ii = 0; ii < k; ii++) {
    for (let jj = ii + 1; jj < k; jj++) {
      const d = Math.abs(dot(eigenvectors[ii], eigenvectors[jj]));
      if (d > orthogonality) {
        orthogonality = d;
      }
    }
  }

  const maxResidual = Math.max(...residuals);
  const maxRelativeResidual = Math.max(...relativeResiduals);

  return {
    residuals,
    maxResidual,
    relativeResiduals,
    maxRelativeResidual,
    orthogonality,
    isValid: maxResidual < tol,
  };
}

/**
 * Verify singular value decomposition.
 *
 * Checks that:
 * - A*v_i ≈ σ_i * u_i
 * - A^T*u_i ≈ σ_i * v_i
 * - U vectors are orthonormal
 * - V vectors are orthonormal
 *
 * @param matvec - Matrix operator A*x
 * @param matvecT - Transpose operator A^T*x
 * @param U - Left singular vectors
 * @param s - Singular values
 * @param Vt - Right singular vectors (rows)
 * @param options - Verification options
 * @returns Verification results
 */
export function verifySvds(
  matvec: MatVecFunction,
  matvecT: MatVecFunction,
  U: Float64Array[],
  s: Float64Array,
  Vt: Float64Array[],
  options?: VerifyOptions
): VerifySvdsResult {
  const { tol = 1e-10 } = options ?? {};
  const k = s.length;

  if (U.length !== k || Vt.length !== k) {
    throw new Error('U, s, and Vt must have consistent lengths');
  }

  const residualsAV = new Float64Array(k);
  const residualsATU = new Float64Array(k);

  for (let i = 0; i < k; i++) {
    const sigma = s[i];
    const u = U[i];
    const v = Vt[i]; // v is i-th row of Vt, which is i-th right singular vector

    // Check A*v ≈ σ*u
    const Av = matvec(v);
    const AvArr = Av instanceof Float64Array ? Av : new Float64Array(Av);
    let resAV = 0;
    for (let j = 0; j < u.length; j++) {
      const diff = AvArr[j] - sigma * u[j];
      resAV += diff * diff;
    }
    residualsAV[i] = Math.sqrt(resAV);

    // Check A^T*u ≈ σ*v
    const ATu = matvecT(u);
    const ATuArr = ATu instanceof Float64Array ? ATu : new Float64Array(ATu);
    let resATU = 0;
    for (let j = 0; j < v.length; j++) {
      const diff = ATuArr[j] - sigma * v[j];
      resATU += diff * diff;
    }
    residualsATU[i] = Math.sqrt(resATU);
  }

  // Check orthogonality of U
  let orthogonalityU = 0;
  for (let i = 0; i < k; i++) {
    for (let j = i + 1; j < k; j++) {
      const d = Math.abs(dot(U[i], U[j]));
      if (d > orthogonalityU) {
        orthogonalityU = d;
      }
    }
  }

  // Check orthogonality of V
  let orthogonalityV = 0;
  for (let i = 0; i < k; i++) {
    for (let j = i + 1; j < k; j++) {
      const d = Math.abs(dot(Vt[i], Vt[j]));
      if (d > orthogonalityV) {
        orthogonalityV = d;
      }
    }
  }

  const maxResidual = Math.max(
    Math.max(...residualsAV),
    Math.max(...residualsATU)
  );

  return {
    residualsAV,
    residualsATU,
    maxResidual,
    orthogonalityU,
    orthogonalityV,
    isValid: maxResidual < tol,
  };
}

/**
 * Check if a matrix is symmetric via random probing.
 *
 * Tests whether x^T*A*y ≈ y^T*A*x for random vectors x, y.
 * This is a probabilistic test that can detect asymmetry but
 * cannot prove symmetry.
 *
 * @param matvec - Matrix operator A
 * @param n - Matrix dimension
 * @param nProbes - Number of random probe pairs (default 10)
 * @param tol - Tolerance for asymmetry (default 1e-10)
 * @returns Symmetry check result
 *
 * @example
 * ```ts
 * const result = checkSymmetry(matvec, 100);
 * if (!result.isSymmetric) {
 *   console.log('Matrix is not symmetric, max asymmetry:', result.maxAsymmetry);
 * }
 * ```
 */
export function checkSymmetry(
  matvec: MatVecFunction,
  n: number,
  nProbes: number = 10,
  tol: number = 1e-10
): SymmetryCheckResult {
  let maxAsymmetry = 0;

  for (let probe = 0; probe < nProbes; probe++) {
    // Generate random vectors
    const x = new Float64Array(n);
    const y = new Float64Array(n);
    for (let i = 0; i < n; i++) {
      x[i] = Math.random() - 0.5;
      y[i] = Math.random() - 0.5;
    }

    // Compute A*x and A*y
    const Ax = matvec(x);
    const Ay = matvec(y);
    const AxArr = Ax instanceof Float64Array ? Ax : new Float64Array(Ax);
    const AyArr = Ay instanceof Float64Array ? Ay : new Float64Array(Ay);

    // Compute x^T*A*y and y^T*A*x
    const xTAy = dot(x, AyArr);
    const yTAx = dot(y, AxArr);

    const asymmetry = Math.abs(xTAy - yTAx);
    if (asymmetry > maxAsymmetry) {
      maxAsymmetry = asymmetry;
    }
  }

  return {
    isSymmetric: maxAsymmetry < tol,
    maxAsymmetry,
    nProbes,
  };
}

/**
 * Check if a matrix is positive definite.
 *
 * Computes the smallest eigenvalue using eigs and checks if it's positive.
 * Note: This requires the eigs function to be available.
 *
 * @param matvec - Matrix operator A (must be symmetric)
 * @param n - Matrix dimension
 * @param tol - Tolerance for positive eigenvalue (default 0)
 * @returns Positive definiteness check result
 */
export async function checkPositiveDefinite(
  matvec: MatVecFunction,
  n: number,
  tol: number = 0
): Promise<PositiveDefiniteCheckResult> {
  // Dynamic import to avoid circular dependency
  const { eigs } = await import('../core/eigs.js');

  try {
    // Find smallest algebraic eigenvalue
    const result = await eigs(matvec, n, 1, {
      which: 'SA',
      return_eigenvectors: false,
    });

    if (!result.success || result.nconv === 0) {
      return {
        isPositiveDefinite: false,
        smallestEigenvalue: NaN,
        success: false,
        message: 'Failed to compute smallest eigenvalue: ' + result.message,
      };
    }

    const smallestEigenvalue = result.eigenvalues[0];

    return {
      isPositiveDefinite: smallestEigenvalue > tol,
      smallestEigenvalue,
      success: true,
      message: smallestEigenvalue > tol
        ? 'Matrix is positive definite'
        : `Matrix is not positive definite (smallest eigenvalue: ${smallestEigenvalue})`,
    };
  } catch (error) {
    return {
      isPositiveDefinite: false,
      smallestEigenvalue: NaN,
      success: false,
      message: `Error checking positive definiteness: ${error}`,
    };
  }
}

/**
 * Check normalization of eigenvectors.
 *
 * Verifies that each eigenvector has unit norm.
 *
 * @param eigenvectors - Array of eigenvectors
 * @param tol - Tolerance for unit norm (default 1e-10)
 * @returns Object with normalization info
 */
export function checkNormalization(
  eigenvectors: Float64Array[],
  tol: number = 1e-10
): { isNormalized: boolean; norms: Float64Array; maxDeviation: number } {
  const k = eigenvectors.length;
  const norms = new Float64Array(k);

  for (let i = 0; i < k; i++) {
    norms[i] = norm2(eigenvectors[i]);
  }

  let maxDeviation = 0;
  for (let i = 0; i < k; i++) {
    const dev = Math.abs(norms[i] - 1);
    if (dev > maxDeviation) {
      maxDeviation = dev;
    }
  }

  return {
    isNormalized: maxDeviation < tol,
    norms,
    maxDeviation,
  };
}
