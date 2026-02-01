/**
 * Matrix Power
 *
 * Compute A^p for a square matrix A and integer or real power p.
 */

import {
  prepareMatrix,
  getMatrixDimensions,
} from '../helpers.js';
import { expm } from './expm.js';
import { logm } from './logm.js';
import type { Matrix, PowmResult } from './types.js';

/**
 * Compute matrix power A^p.
 *
 * For integer p >= 0: uses binary exponentiation.
 * For integer p < 0: computes inverse then uses binary exponentiation.
 * For non-integer p: uses A^p = expm(p * logm(A)).
 *
 * @param A - Input square matrix (n × n)
 * @param p - Power (integer or real)
 * @returns Matrix power A^p
 *
 * @example
 * ```typescript
 * const A = [[1, 1], [0, 1]];
 *
 * // A^2
 * const { P: A2 } = powm(A, 2);
 * // A2 = [[1, 2], [0, 1]]
 *
 * // A^(-1)
 * const { P: Ainv } = powm(A, -1);
 * // Ainv = [[1, -1], [0, 1]]
 *
 * // A^0.5 (square root)
 * const { P: Ahalf } = powm(A, 0.5);
 * ```
 */
export function powm(A: Matrix, p: number): PowmResult {
  // Get dimensions
  const [m, n] = getMatrixDimensions(A);

  if (m !== n) {
    throw new Error(`Matrix must be square, got ${m}×${n}`);
  }

  if (n === 0) {
    return {
      P: new Float64Array(0),
      n: 0,
      p,
      success: true,
      message: 'Success',
    };
  }

  // Prepare matrix (column-major)
  const aData = prepareMatrix(A);

  // Identity matrix
  const identity = (): Float64Array => {
    const I = new Float64Array(n * n);
    for (let i = 0; i < n; i++) {
      I[i + i * n] = 1;
    }
    return I;
  };

  // Matrix multiplication
  const matmul = (X: Float64Array, Y: Float64Array): Float64Array => {
    const result = new Float64Array(n * n);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        let sum = 0;
        for (let k = 0; k < n; k++) {
          sum += X[i + k * n] * Y[k + j * n];
        }
        result[i + j * n] = sum;
      }
    }
    return result;
  };

  // Matrix inverse using Gauss-Jordan
  const matrixInverse = (M: Float64Array): Float64Array | null => {
    const aug = new Float64Array(n * 2 * n);

    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        aug[i + j * n] = M[i + j * n];
        aug[i + (j + n) * n] = i === j ? 1 : 0;
      }
    }

    for (let col = 0; col < n; col++) {
      let maxVal = Math.abs(aug[col + col * n]);
      let maxRow = col;
      for (let row = col + 1; row < n; row++) {
        if (Math.abs(aug[row + col * n]) > maxVal) {
          maxVal = Math.abs(aug[row + col * n]);
          maxRow = row;
        }
      }

      if (maxVal < 1e-15) {
        return null;
      }

      if (maxRow !== col) {
        for (let j = 0; j < 2 * n; j++) {
          const temp = aug[col + j * n];
          aug[col + j * n] = aug[maxRow + j * n];
          aug[maxRow + j * n] = temp;
        }
      }

      const pivot = aug[col + col * n];
      for (let j = 0; j < 2 * n; j++) {
        aug[col + j * n] /= pivot;
      }

      for (let row = 0; row < n; row++) {
        if (row !== col) {
          const factor = aug[row + col * n];
          for (let j = 0; j < 2 * n; j++) {
            aug[row + j * n] -= factor * aug[col + j * n];
          }
        }
      }
    }

    const inv = new Float64Array(n * n);
    for (let i = 0; i < n; i++) {
      for (let j = 0; j < n; j++) {
        inv[i + j * n] = aug[i + (j + n) * n];
      }
    }

    return inv;
  };

  // Handle special cases
  if (p === 0) {
    return {
      P: identity(),
      n,
      p,
      success: true,
      message: 'Success',
    };
  }

  if (p === 1) {
    const result = new Float64Array(n * n);
    for (let i = 0; i < n * n; i++) {
      result[i] = aData[i];
    }
    return {
      P: result,
      n,
      p,
      success: true,
      message: 'Success',
    };
  }

  // Check if p is an integer
  const isInteger = Number.isInteger(p);

  if (isInteger) {
    // Binary exponentiation for integer powers
    let exp = Math.abs(p);
    let base = new Float64Array(n * n);
    for (let i = 0; i < n * n; i++) {
      base[i] = aData[i];
    }

    if (p < 0) {
      const inv = matrixInverse(base);
      if (!inv) {
        return {
          P: new Float64Array(0),
          n,
          p,
          success: false,
          message: 'Matrix is singular, cannot compute negative power',
        };
      }
      base = inv;
    }

    let result = identity();

    while (exp > 0) {
      if (exp & 1) {
        result = matmul(result, base);
      }
      base = matmul(base, base);
      exp >>= 1;
    }

    return {
      P: result,
      n,
      p,
      success: true,
      message: 'Success',
    };
  } else {
    // Non-integer power: A^p = expm(p * logm(A))
    // Convert column-major data to 2D array for logm
    const A2D: number[][] = [];
    for (let i = 0; i < n; i++) {
      A2D.push([]);
      for (let j = 0; j < n; j++) {
        A2D[i].push(aData[i + j * n]);
      }
    }

    const logResult = logm(A2D);
    if (!logResult.success) {
      return {
        P: new Float64Array(0),
        n,
        p,
        success: false,
        message: logResult.message,
      };
    }

    // Multiply by p
    const pLogA = new Float64Array(n * n);
    for (let i = 0; i < n * n; i++) {
      pLogA[i] = p * logResult.L[i];
    }

    // Convert to 2D array for expm
    const pLogA2D: number[][] = [];
    for (let i = 0; i < n; i++) {
      pLogA2D.push([]);
      for (let j = 0; j < n; j++) {
        pLogA2D[i].push(pLogA[i + j * n]);
      }
    }

    const expResult = expm(pLogA2D);

    return {
      P: expResult.E,
      n,
      p,
      success: expResult.success,
      message: expResult.message,
    };
  }
}
