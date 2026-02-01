/**
 * Hermitian Rank-2k Update (Pure TypeScript)
 *
 * Compute C = alpha*A*B^H + conj(alpha)*B*A^H + beta*C or
 *         C = alpha*A^H*B + conj(alpha)*B^H*A + beta*C
 *
 * Note: ZHER2K is for complex matrices. This implementation provides
 * a real-valued version (equivalent to symmetric rank-2k update).
 */

import {
  prepareMatrix,
  getMatrixDimensions,
} from '../helpers.js';
import type { Matrix, Her2kOptions, Her2kResult } from './types.js';

/**
 * Compute Hermitian rank-2k update.
 *
 * For real matrices, this is equivalent to symmetric rank-2k update:
 * If trans='N': C = alpha*A*B^T + alpha*B*A^T + beta*C
 * If trans='C': C = alpha*A^T*B + alpha*B^T*A + beta*C
 *
 * The result C is Hermitian (symmetric for real matrices).
 *
 * @param A - First input matrix
 * @param B - Second input matrix (same dimensions as A)
 * @param options - Update options
 * @returns Hermitian result matrix C
 *
 * @example
 * ```typescript
 * const A = [[1, 2, 3], [4, 5, 6]];
 * const B = [[1, 0, 1], [0, 1, 0]];
 *
 * // C = A*B^H + B*A^H (2×2 result)
 * const { C } = her2k(A, B);
 *
 * // C = A^H*B + B^H*A (3×3 result)
 * const { C: C2 } = her2k(A, B, { trans: 'C' });
 * ```
 */
export function her2k(A: Matrix, B: Matrix, options: Her2kOptions = {}): Her2kResult {
  const {
    uplo = 'lower',
    trans = 'N',
    alpha = 1.0,
    beta = 0.0,
    C,
  } = options;

  // Get dimensions
  const [rowsA, colsA] = getMatrixDimensions(A);
  const [rowsB, colsB] = getMatrixDimensions(B);

  // Check A and B have same dimensions
  if (rowsA !== rowsB || colsA !== colsB) {
    throw new Error(
      `A (${rowsA}×${colsA}) and B (${rowsB}×${colsB}) must have same dimensions`
    );
  }

  // Determine n (dimension of result) and k (common dimension)
  const n = trans === 'N' ? rowsA : colsA;
  const k = trans === 'N' ? colsA : rowsA;

  // Prepare matrices (column-major)
  const aData = prepareMatrix(A);
  const bData = prepareMatrix(B);

  // Prepare or create C
  let cData: Float64Array;
  if (C !== undefined) {
    const [rowsC, colsC] = getMatrixDimensions(C);
    if (rowsC !== n || colsC !== n) {
      throw new Error(
        `C dimensions ${rowsC}×${colsC} incompatible with result ${n}×${n}`
      );
    }
    cData = prepareMatrix(C);
  } else {
    cData = new Float64Array(n * n);
  }

  const result = new Float64Array(n * n);

  // Apply beta to existing C
  if (beta !== 0) {
    for (let i = 0; i < n * n; i++) {
      result[i] = beta * cData[i];
    }
  }

  // For real matrices: C = alpha * (A * B^T + B * A^T) + beta * C
  // or C = alpha * (A^T * B + B^T * A) + beta * C

  if (trans === 'N') {
    // C = alpha * A * B^T + alpha * B * A^T
    // A, B are rowsA × colsA, result is rowsA × rowsA
    for (let i = 0; i < n; i++) {
      for (let j = 0; j <= i; j++) {
        let sum = 0;
        for (let l = 0; l < k; l++) {
          // A(i,l) * B(j,l) + B(i,l) * A(j,l)
          const aIL = aData[i + l * rowsA];
          const aJL = aData[j + l * rowsA];
          const bIL = bData[i + l * rowsA];
          const bJL = bData[j + l * rowsA];
          sum += aIL * bJL + bIL * aJL;
        }
        result[i + j * n] += alpha * sum;
        if (i !== j) {
          result[j + i * n] += alpha * sum;
        }
      }
    }
  } else {
    // C = alpha * A^T * B + alpha * B^T * A (or A^H * B + B^H * A)
    // A, B are rowsA × colsA, result is colsA × colsA
    for (let i = 0; i < n; i++) {
      for (let j = 0; j <= i; j++) {
        let sum = 0;
        for (let l = 0; l < k; l++) {
          // A(l,i) * B(l,j) + B(l,i) * A(l,j)
          const aLI = aData[l + i * rowsA];
          const aLJ = aData[l + j * rowsA];
          const bLI = bData[l + i * rowsA];
          const bLJ = bData[l + j * rowsA];
          sum += aLI * bLJ + bLI * aLJ;
        }
        result[i + j * n] += alpha * sum;
        if (i !== j) {
          result[j + i * n] += alpha * sum;
        }
      }
    }
  }

  // Ensure symmetry based on uplo preference
  if (uplo === 'upper') {
    for (let j = 0; j < n; j++) {
      for (let i = j + 1; i < n; i++) {
        result[j + i * n] = result[i + j * n];
      }
    }
  } else {
    for (let j = 0; j < n; j++) {
      for (let i = j + 1; i < n; i++) {
        result[i + j * n] = result[j + i * n];
      }
    }
  }

  return {
    C: result,
    n,
    success: true,
    message: 'Success',
  };
}
