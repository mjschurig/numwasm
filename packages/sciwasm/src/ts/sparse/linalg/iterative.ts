/**
 * Iterative solvers for sparse linear systems
 *
 * Based on scipy.sparse.linalg._isolve.iterative
 */

import {
  LinearOperator,
  aslinearoperator,
  IdentityOperator,
} from './interface.js';
import type {
  IterativeSolverResult,
  IterativeSolverOptions,
  GMRESOptions,
  LinearOperatorLike,
} from './types.js';

/**
 * Helper vector operations (pure TypeScript, synchronous)
 */
function vecAdd(a: Float64Array, b: Float64Array): Float64Array {
  const result = new Float64Array(a.length);
  for (let i = 0; i < a.length; i++) {
    result[i] = a[i] + b[i];
  }
  return result;
}

function vecSub(a: Float64Array, b: Float64Array): Float64Array {
  const result = new Float64Array(a.length);
  for (let i = 0; i < a.length; i++) {
    result[i] = a[i] - b[i];
  }
  return result;
}

function vecScale(scalar: number, a: Float64Array): Float64Array {
  const result = new Float64Array(a.length);
  for (let i = 0; i < a.length; i++) {
    result[i] = scalar * a[i];
  }
  return result;
}

function vecAxpy(alpha: number, x: Float64Array, y: Float64Array): Float64Array {
  // y + alpha * x
  const result = new Float64Array(y.length);
  for (let i = 0; i < y.length; i++) {
    result[i] = y[i] + alpha * x[i];
  }
  return result;
}

function vecDot(a: Float64Array, b: Float64Array): number {
  let sum = 0;
  for (let i = 0; i < a.length; i++) {
    sum += a[i] * b[i];
  }
  return sum;
}

function vecNorm(a: Float64Array): number {
  let sumSq = 0;
  for (let i = 0; i < a.length; i++) {
    sumSq += a[i] * a[i];
  }
  return Math.sqrt(sumSq);
}

function vecCopy(a: Float64Array): Float64Array {
  return new Float64Array(a);
}

/**
 * Conjugate Gradient method for solving Ax = b
 *
 * Solves the linear system Ax = b for x, where A is a symmetric
 * positive-definite matrix.
 *
 * @param A - Sparse matrix or LinearOperator (must be symmetric positive-definite)
 * @param b - Right-hand side vector
 * @param options - Solver options
 * @returns Solver result with solution x and convergence info
 *
 * @example
 * ```typescript
 * // Create a symmetric positive-definite matrix
 * const A = csr_matrix([
 *   [4, 1, 0],
 *   [1, 3, 1],
 *   [0, 1, 2]
 * ]);
 * const b = new Float64Array([1, 2, 3]);
 *
 * const result = cg(A, b);
 * console.log(result.x);    // Solution
 * console.log(result.info); // 0 if converged
 * ```
 */
export function cg(
  A: LinearOperatorLike,
  b: Float64Array,
  options: IterativeSolverOptions = {}
): IterativeSolverResult {
  const Aop = aslinearoperator(A);
  const [m, n] = Aop.shape;

  if (m !== n) {
    throw new Error(`Matrix must be square, got ${m}x${n}`);
  }

  if (b.length !== n) {
    throw new Error(`Dimension mismatch: A is ${n}x${n}, b has length ${b.length}`);
  }

  // Solver parameters
  const tol = options.tol ?? 1e-5;
  const atol = options.atol ?? 0;
  const maxiter = options.maxiter ?? n;
  const callback = options.callback;

  // Preconditioner (default: identity)
  let Mop: LinearOperator;
  if (options.M) {
    Mop = aslinearoperator(options.M);
  } else {
    Mop = new IdentityOperator(n);
  }

  // Initial guess
  let x: Float64Array;
  if (options.x0) {
    x = vecCopy(options.x0);
  } else {
    x = new Float64Array(n);
  }

  // Compute initial residual: r = b - A @ x
  let r = vecSub(b, Aop.matvec(x));

  // Convergence threshold
  const bnorm = vecNorm(b);
  const threshold = Math.max(tol * bnorm, atol);

  // Check for trivial solution
  let rnorm = vecNorm(r);
  if (rnorm <= threshold) {
    return {
      x,
      info: 0,
      iterations: 0,
      residualNorm: rnorm,
    };
  }

  // Apply preconditioner: z = M @ r
  let z = Mop.matvec(r);

  // p = z
  let p = vecCopy(z);

  // rho = r.dot(z)
  let rho = vecDot(r, z);

  let iterations = 0;

  for (let k = 0; k < maxiter; k++) {
    iterations = k + 1;

    // Ap = A @ p
    const Ap = Aop.matvec(p);

    // alpha = rho / (p.dot(Ap))
    const pAp = vecDot(p, Ap);
    if (Math.abs(pAp) < 1e-16) {
      // Breakdown
      return {
        x,
        info: -1,
        iterations,
        residualNorm: rnorm,
      };
    }
    const alpha = rho / pAp;

    // x = x + alpha * p
    x = vecAxpy(alpha, p, x);

    // r = r - alpha * Ap
    r = vecAxpy(-alpha, Ap, r);

    // Check convergence
    rnorm = vecNorm(r);
    if (rnorm <= threshold) {
      if (callback) callback(x);
      return {
        x,
        info: 0,
        iterations,
        residualNorm: rnorm,
      };
    }

    // z = M @ r
    z = Mop.matvec(r);

    // rho_new = r.dot(z)
    const rho_new = vecDot(r, z);

    // beta = rho_new / rho
    const beta = rho_new / rho;
    rho = rho_new;

    // p = z + beta * p
    p = vecAxpy(beta, p, z);

    if (callback) callback(x);
  }

  // Did not converge
  return {
    x,
    info: maxiter,
    iterations,
    residualNorm: rnorm,
  };
}

/**
 * BiConjugate Gradient Stabilized method for solving Ax = b
 *
 * Solves the linear system Ax = b for x, where A is a general
 * (possibly non-symmetric) matrix.
 *
 * @param A - Sparse matrix or LinearOperator
 * @param b - Right-hand side vector
 * @param options - Solver options
 * @returns Solver result with solution x and convergence info
 *
 * @example
 * ```typescript
 * const A = csr_matrix([
 *   [4, 1, 0],
 *   [1, 3, 1],
 *   [0, 1, 2]
 * ]);
 * const b = new Float64Array([1, 2, 3]);
 *
 * const result = bicgstab(A, b);
 * console.log(result.x);    // Solution
 * console.log(result.info); // 0 if converged
 * ```
 */
export function bicgstab(
  A: LinearOperatorLike,
  b: Float64Array,
  options: IterativeSolverOptions = {}
): IterativeSolverResult {
  const Aop = aslinearoperator(A);
  const [m, n] = Aop.shape;

  if (m !== n) {
    throw new Error(`Matrix must be square, got ${m}x${n}`);
  }

  if (b.length !== n) {
    throw new Error(`Dimension mismatch: A is ${n}x${n}, b has length ${b.length}`);
  }

  // Solver parameters
  const tol = options.tol ?? 1e-5;
  const atol = options.atol ?? 0;
  const maxiter = options.maxiter ?? n;
  const callback = options.callback;

  // Preconditioner (default: identity)
  let Mop: LinearOperator;
  if (options.M) {
    Mop = aslinearoperator(options.M);
  } else {
    Mop = new IdentityOperator(n);
  }

  // Initial guess
  let x: Float64Array;
  if (options.x0) {
    x = vecCopy(options.x0);
  } else {
    x = new Float64Array(n);
  }

  // Compute initial residual: r = b - A @ x
  let r = vecSub(b, Aop.matvec(x));

  // Convergence threshold
  const bnorm = vecNorm(b);
  const threshold = Math.max(tol * bnorm, atol);

  // Check for trivial solution
  let rnorm = vecNorm(r);
  if (rnorm <= threshold) {
    return {
      x,
      info: 0,
      iterations: 0,
      residualNorm: rnorm,
    };
  }

  // rtilde = r (shadow residual, kept constant)
  const rtilde = vecCopy(r);

  // Initialize scalars
  let rho = 1;
  let alpha = 1;
  let omega = 1;

  // v = zeros, p = zeros
  let v = new Float64Array(n);
  let p = new Float64Array(n);

  let iterations = 0;

  for (let k = 0; k < maxiter; k++) {
    iterations = k + 1;

    // rho_new = rtilde.dot(r)
    const rho_new = vecDot(rtilde, r);

    if (Math.abs(rho_new) < 1e-16) {
      // Breakdown
      return {
        x,
        info: -1,
        iterations,
        residualNorm: rnorm,
      };
    }

    // beta = (rho_new / rho) * (alpha / omega)
    const beta = (rho_new / rho) * (alpha / omega);
    rho = rho_new;

    // p = r + beta * (p - omega * v)
    p = vecAdd(r, vecScale(beta, vecSub(p, vecScale(omega, v))));

    // phat = M @ p (preconditioner)
    const phat = Mop.matvec(p);

    // v = A @ phat
    v = Aop.matvec(phat);

    // alpha = rho / rtilde.dot(v)
    const rtv = vecDot(rtilde, v);
    if (Math.abs(rtv) < 1e-16) {
      // Breakdown
      return {
        x,
        info: -2,
        iterations,
        residualNorm: rnorm,
      };
    }
    alpha = rho / rtv;

    // s = r - alpha * v
    const s = vecAxpy(-alpha, v, r);

    // Check for early convergence
    const snorm = vecNorm(s);
    if (snorm <= threshold) {
      // x = x + alpha * phat
      x = vecAxpy(alpha, phat, x);
      if (callback) callback(x);
      return {
        x,
        info: 0,
        iterations,
        residualNorm: snorm,
      };
    }

    // shat = M @ s (preconditioner)
    const shat = Mop.matvec(s);

    // t = A @ shat
    const t = Aop.matvec(shat);

    // omega = t.dot(s) / t.dot(t)
    const ts = vecDot(t, s);
    const tt = vecDot(t, t);
    if (Math.abs(tt) < 1e-16) {
      // Breakdown
      return {
        x,
        info: -3,
        iterations,
        residualNorm: rnorm,
      };
    }
    omega = ts / tt;

    // x = x + alpha * phat + omega * shat
    x = vecAxpy(alpha, phat, vecAxpy(omega, shat, x));

    // r = s - omega * t
    r = vecAxpy(-omega, t, s);

    // Check convergence
    rnorm = vecNorm(r);
    if (rnorm <= threshold) {
      if (callback) callback(x);
      return {
        x,
        info: 0,
        iterations,
        residualNorm: rnorm,
      };
    }

    if (Math.abs(omega) < 1e-16) {
      // Breakdown
      return {
        x,
        info: -4,
        iterations,
        residualNorm: rnorm,
      };
    }

    if (callback) callback(x);
  }

  // Did not converge
  return {
    x,
    info: maxiter,
    iterations,
    residualNorm: rnorm,
  };
}

/**
 * Compute Givens rotation parameters
 * Returns [c, s, r] such that:
 *   [c  s] [a]   [r]
 *   [-s c] [b] = [0]
 */
function givensRotation(a: number, b: number): [number, number, number] {
  if (b === 0) {
    return [1, 0, a];
  } else if (Math.abs(b) > Math.abs(a)) {
    const tau = -a / b;
    const s = 1 / Math.sqrt(1 + tau * tau);
    const c = s * tau;
    const r = b / s;
    return [c, s, r];
  } else {
    const tau = -b / a;
    const c = 1 / Math.sqrt(1 + tau * tau);
    const s = c * tau;
    const r = a / c;
    return [c, s, r];
  }
}

/**
 * Apply Givens rotation to two elements
 */
function applyGivens(c: number, s: number, x: number, y: number): [number, number] {
  return [c * x - s * y, s * x + c * y];
}

/**
 * GMRES (Generalized Minimal Residual) method for solving Ax = b
 *
 * Solves the linear system Ax = b for x using the restarted GMRES algorithm.
 * Works for general (possibly non-symmetric) matrices.
 *
 * @param A - Sparse matrix or LinearOperator
 * @param b - Right-hand side vector
 * @param options - Solver options including restart parameter
 * @returns Solver result with solution x and convergence info
 *
 * @example
 * ```typescript
 * const A = csr_matrix([
 *   [4, 1, 0],
 *   [1, 3, 1],
 *   [0, 1, 2]
 * ]);
 * const b = new Float64Array([1, 2, 3]);
 *
 * const result = gmres(A, b, { restart: 10 });
 * console.log(result.x);    // Solution
 * console.log(result.info); // 0 if converged
 * ```
 */
export function gmres(
  A: LinearOperatorLike,
  b: Float64Array,
  options: GMRESOptions = {}
): IterativeSolverResult {
  const Aop = aslinearoperator(A);
  const [m, n] = Aop.shape;

  if (m !== n) {
    throw new Error(`Matrix must be square, got ${m}x${n}`);
  }

  if (b.length !== n) {
    throw new Error(`Dimension mismatch: A is ${n}x${n}, b has length ${b.length}`);
  }

  // Solver parameters
  const tol = options.tol ?? 1e-5;
  const atol = options.atol ?? 0;
  const restart = options.restart ?? Math.min(20, n);
  const maxiter = options.maxiter ?? n;
  const callback = options.callback;

  // Preconditioner (default: identity)
  let Mop: LinearOperator;
  if (options.M) {
    Mop = aslinearoperator(options.M);
  } else {
    Mop = new IdentityOperator(n);
  }

  // Initial guess
  let x: Float64Array;
  if (options.x0) {
    x = vecCopy(options.x0);
  } else {
    x = new Float64Array(n);
  }

  // Convergence threshold
  const bnorm = vecNorm(b);
  const threshold = Math.max(tol * bnorm, atol);

  let totalIterations = 0;
  let rnorm = 0;

  // Outer loop (restarts)
  for (let outer = 0; outer < Math.ceil(maxiter / restart); outer++) {
    // Compute residual: r = b - A @ x
    const r = vecSub(b, Aop.matvec(x));

    // Apply preconditioner: r = M @ r
    const r_precond = Mop.matvec(r);

    rnorm = vecNorm(r_precond);

    // Check for convergence
    if (rnorm <= threshold) {
      return {
        x,
        info: 0,
        iterations: totalIterations,
        residualNorm: rnorm,
      };
    }

    // Normalize: v1 = r / ||r||
    const V: Float64Array[] = [vecScale(1 / rnorm, r_precond)];

    // Hessenberg matrix H (upper Hessenberg)
    const H: number[][] = [];
    for (let i = 0; i <= restart; i++) {
      H.push(new Array(restart).fill(0));
    }

    // RHS for least squares
    const g = new Float64Array(restart + 1);
    g[0] = rnorm;

    // Givens rotation parameters
    const cs: number[] = [];
    const sn: number[] = [];

    // Arnoldi iteration
    for (let j = 0; j < restart && totalIterations < maxiter; j++) {
      totalIterations++;

      // w = A @ M @ V[j]
      const Mvj = Mop.matvec(V[j]);
      let w = Aop.matvec(Mvj);

      // Modified Gram-Schmidt orthogonalization
      for (let i = 0; i <= j; i++) {
        H[i][j] = vecDot(w, V[i]);
        w = vecAxpy(-H[i][j], V[i], w);
      }

      H[j + 1][j] = vecNorm(w);

      if (Math.abs(H[j + 1][j]) < 1e-16) {
        // Lucky breakdown - exact solution found
        // Solve the least squares problem
        const y = solveUpperTriangular(H, g, j + 1);
        for (let i = 0; i <= j; i++) {
          const scaled = vecScale(y[i], V[i]);
          const precond = Mop.matvec(scaled);
          x = vecAdd(x, precond);
        }
        return {
          x,
          info: 0,
          iterations: totalIterations,
          residualNorm: Math.abs(g[j + 1]),
        };
      }

      // Normalize
      V.push(vecScale(1 / H[j + 1][j], w));

      // Apply previous Givens rotations to the j-th column of H
      for (let i = 0; i < j; i++) {
        [H[i][j], H[i + 1][j]] = applyGivens(cs[i], sn[i], H[i][j], H[i + 1][j]);
      }

      // Compute new Givens rotation to eliminate H[j+1][j]
      const [c, s, r] = givensRotation(H[j][j], H[j + 1][j]);
      cs.push(c);
      sn.push(s);

      // Apply to H
      H[j][j] = r;
      H[j + 1][j] = 0;

      // Apply to RHS
      [g[j], g[j + 1]] = applyGivens(c, s, g[j], g[j + 1]);

      rnorm = Math.abs(g[j + 1]);

      if (callback) {
        callback(rnorm);
      }

      // Check convergence
      if (rnorm <= threshold) {
        // Solve the least squares problem
        const y = solveUpperTriangular(H, g, j + 1);
        for (let i = 0; i <= j; i++) {
          const scaled = vecScale(y[i], V[i]);
          const precond = Mop.matvec(scaled);
          x = vecAdd(x, precond);
        }
        return {
          x,
          info: 0,
          iterations: totalIterations,
          residualNorm: rnorm,
        };
      }
    }

    // End of inner loop - solve least squares and update x
    const innerIters = Math.min(restart, totalIterations - outer * restart);
    const y = solveUpperTriangular(H, g, innerIters);
    for (let i = 0; i < y.length; i++) {
      const scaled = vecScale(y[i], V[i]);
      const precond = Mop.matvec(scaled);
      x = vecAdd(x, precond);
    }
  }

  // Did not converge
  return {
    x,
    info: maxiter,
    iterations: totalIterations,
    residualNorm: rnorm,
  };
}

/**
 * Solve upper triangular system Ry = g by back substitution
 */
function solveUpperTriangular(H: number[][], g: Float64Array, k: number): number[] {
  const y = new Array(k).fill(0);

  for (let i = k - 1; i >= 0; i--) {
    let sum = g[i];
    for (let j = i + 1; j < k; j++) {
      sum -= H[i][j] * y[j];
    }
    y[i] = sum / H[i][i];
  }

  return y;
}
