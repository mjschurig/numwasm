/**
 * Operator Combinators
 *
 * Utility functions for composing and transforming matrix-vector product functions.
 * These enable building complex operators from simpler ones without forming
 * explicit matrices.
 */

import type { MatVecFunction } from '../high-level-types.js';

/**
 * Create a linear combination operator: y = (α*A + β*B) * x
 *
 * @param A - First matrix operator
 * @param B - Second matrix operator
 * @param alpha - Scalar multiplier for A (default 1)
 * @param beta - Scalar multiplier for B (default 1)
 * @returns Function computing y = α*A*x + β*B*x
 *
 * @example
 * ```ts
 * // Create (2*A - B) operator
 * const combined = addMatvec(Amatvec, Bmatvec, 2, -1);
 * const y = combined(x);
 * ```
 */
export function addMatvec(
  A: MatVecFunction,
  B: MatVecFunction,
  alpha: number = 1,
  beta: number = 1
): MatVecFunction {
  return (x: Float64Array): Float64Array => {
    const Ax = A(x);
    const Bx = B(x);
    const n = Ax.length;
    const y = new Float64Array(n);

    for (let i = 0; i < n; i++) {
      y[i] = alpha * Ax[i] + beta * Bx[i];
    }

    return y;
  };
}

/**
 * Create a product operator: y = A * B * x
 *
 * Computes the composition of two operators. Note: this computes B first,
 * then A, matching mathematical notation A*B*x = A*(B*x).
 *
 * @param A - Outer operator (applied second)
 * @param B - Inner operator (applied first)
 * @returns Function computing y = A*(B*x)
 *
 * @example
 * ```ts
 * // Create A*B operator
 * const AB = mulMatvec(Amatvec, Bmatvec);
 * const y = AB(x); // Computes A*(B*x)
 * ```
 */
export function mulMatvec(
  A: MatVecFunction,
  B: MatVecFunction
): MatVecFunction {
  return (x: Float64Array): Float64Array => {
    const Bx = B(x);
    // Ensure Float64Array for A
    const BxArr = Bx instanceof Float64Array ? Bx : new Float64Array(Bx);
    const result = A(BxArr);
    return result instanceof Float64Array ? result : new Float64Array(result);
  };
}

/**
 * Create a shifted operator: y = (A - σI) * x
 *
 * Useful for shift-invert mode where you need to work with A - σI.
 *
 * @param matvec - Matrix operator A
 * @param sigma - Shift value σ
 * @returns Function computing y = A*x - σ*x
 *
 * @example
 * ```ts
 * // Create (A - 5*I) operator for shift-invert near σ=5
 * const shifted = shiftMatvec(Amatvec, 5);
 * const y = shifted(x);
 * ```
 */
export function shiftMatvec(
  matvec: MatVecFunction,
  sigma: number
): MatVecFunction {
  return (x: Float64Array): Float64Array => {
    const Ax = matvec(x);
    const n = x.length;
    const y = new Float64Array(n);

    for (let i = 0; i < n; i++) {
      y[i] = Ax[i] - sigma * x[i];
    }

    return y;
  };
}

/**
 * Create a scaled operator: y = α * A * x
 *
 * @param matvec - Matrix operator A
 * @param alpha - Scalar multiplier
 * @returns Function computing y = α*A*x
 *
 * @example
 * ```ts
 * // Create 0.5*A operator
 * const scaled = scaleMatvec(Amatvec, 0.5);
 * const y = scaled(x);
 * ```
 */
export function scaleMatvec(
  matvec: MatVecFunction,
  alpha: number
): MatVecFunction {
  if (alpha === 1) {
    return matvec; // No-op optimization
  }

  return (x: Float64Array): Float64Array => {
    const Ax = matvec(x);
    const n = Ax.length;
    const y = new Float64Array(n);

    for (let i = 0; i < n; i++) {
      y[i] = alpha * Ax[i];
    }

    return y;
  };
}

/**
 * Create a transpose operator via explicit computation.
 *
 * This builds A^T by probing A with standard basis vectors.
 * WARNING: This is O(n²) and should only be used for small matrices
 * or when no other option exists.
 *
 * For large matrices, prefer providing an explicit transpose operator.
 *
 * @param matvec - Matrix operator A (m x n)
 * @param m - Number of rows in A
 * @param n - Number of columns in A
 * @returns Function computing y = A^T * x (n x m)
 *
 * @example
 * ```ts
 * // Create A^T from A (expensive for large matrices!)
 * const AT = transposeMatvec(Amatvec, 100, 50);
 * const y = AT(x); // x is length 100, y is length 50
 * ```
 */
export function transposeMatvec(
  matvec: MatVecFunction,
  m: number,
  n: number
): MatVecFunction {
  // Build explicit transpose matrix by probing with basis vectors
  // A^T[i,j] = A[j,i] = (A * e_i)[j]
  const AT = new Float64Array(n * m);

  for (let i = 0; i < n; i++) {
    const ei = new Float64Array(n);
    ei[i] = 1;
    const Aei = matvec(ei);

    // Column i of A^T = row i of A = A*e_i
    for (let j = 0; j < m; j++) {
      AT[j * n + i] = Aei[j];
    }
  }

  // Return matvec for the computed transpose
  return (x: Float64Array): Float64Array => {
    if (x.length !== m) {
      throw new Error(`Input vector length must be ${m}, got ${x.length}`);
    }

    const y = new Float64Array(n);

    for (let i = 0; i < n; i++) {
      let sum = 0;
      for (let j = 0; j < m; j++) {
        sum += AT[j * n + i] * x[j];
      }
      y[i] = sum;
    }

    return y;
  };
}

/**
 * Create a symmetrized operator: y = (A + A^T)/2 * x
 *
 * This makes a non-symmetric operator symmetric. Requires both A and A^T.
 *
 * @param matvec - Matrix operator A
 * @param matvecT - Transpose operator A^T
 * @returns Function computing y = (A*x + A^T*x) / 2
 *
 * @example
 * ```ts
 * // Symmetrize a matrix
 * const Asym = symmetrizeMatvec(Amatvec, ATmatvec);
 * const y = Asym(x);
 * ```
 */
export function symmetrizeMatvec(
  matvec: MatVecFunction,
  matvecT: MatVecFunction
): MatVecFunction {
  return (x: Float64Array): Float64Array => {
    const Ax = matvec(x);
    const ATx = matvecT(x);
    const n = Ax.length;
    const y = new Float64Array(n);

    for (let i = 0; i < n; i++) {
      y[i] = 0.5 * (Ax[i] + ATx[i]);
    }

    return y;
  };
}

/**
 * Create an identity operator: y = x
 *
 * @param n - Dimension
 * @returns Function returning a copy of the input
 */
export function identityMatvec(n: number): MatVecFunction {
  return (x: Float64Array): Float64Array => {
    if (x.length !== n) {
      throw new Error(`Input vector length must be ${n}, got ${x.length}`);
    }
    return new Float64Array(x);
  };
}

/**
 * Create a negated operator: y = -A * x
 *
 * @param matvec - Matrix operator A
 * @returns Function computing y = -A*x
 */
export function negateMatvec(matvec: MatVecFunction): MatVecFunction {
  return scaleMatvec(matvec, -1);
}

/**
 * Create a power operator: y = A^k * x
 *
 * Computes the k-th power of an operator by repeated application.
 *
 * @param matvec - Matrix operator A
 * @param k - Power (must be non-negative integer)
 * @returns Function computing y = A^k * x
 *
 * @example
 * ```ts
 * // Compute A³ * x
 * const A3 = powerMatvec(Amatvec, 3);
 * const y = A3(x);
 * ```
 */
export function powerMatvec(
  matvec: MatVecFunction,
  k: number
): MatVecFunction {
  if (k < 0 || !Number.isInteger(k)) {
    throw new Error('Power k must be a non-negative integer');
  }

  if (k === 0) {
    // A^0 = I
    return (x: Float64Array): Float64Array => new Float64Array(x);
  }

  if (k === 1) {
    return matvec;
  }

  return (x: Float64Array): Float64Array => {
    let result = new Float64Array(x);
    for (let i = 0; i < k; i++) {
      const next = matvec(result);
      result = next instanceof Float64Array ? next : new Float64Array(next);
    }
    return result;
  };
}

/**
 * Create a block diagonal operator from multiple sub-operators.
 *
 * Given operators A₁, A₂, ..., Aₖ, creates:
 * ```
 * [A₁  0   0  ]
 * [ 0  A₂  0  ]
 * [ 0   0  A₃ ]
 * ```
 *
 * @param operators - Array of matrix operators
 * @param sizes - Array of dimensions for each operator
 * @returns Function computing y = blkdiag(A₁, A₂, ...) * x
 */
export function blockDiagMatvec(
  operators: MatVecFunction[],
  sizes: number[]
): MatVecFunction {
  if (operators.length !== sizes.length) {
    throw new Error('operators and sizes must have same length');
  }

  const totalSize = sizes.reduce((a, b) => a + b, 0);
  const offsets = new Array(sizes.length);
  let offset = 0;
  for (let i = 0; i < sizes.length; i++) {
    offsets[i] = offset;
    offset += sizes[i];
  }

  return (x: Float64Array): Float64Array => {
    if (x.length !== totalSize) {
      throw new Error(`Input vector length must be ${totalSize}, got ${x.length}`);
    }

    const y = new Float64Array(totalSize);

    for (let k = 0; k < operators.length; k++) {
      const n = sizes[k];
      const off = offsets[k];

      // Extract sub-vector
      const xk = x.subarray(off, off + n);

      // Apply operator
      const yk = operators[k](xk);

      // Copy result
      for (let i = 0; i < n; i++) {
        y[off + i] = yk[i];
      }
    }

    return y;
  };
}
