/**
 * Selected Eigenvalues
 *
 * Compute selected eigenvalues and eigenvectors of a symmetric matrix.
 * Uses DSYEVR which can efficiently compute a subset of eigenvalues.
 */

import { getLAPACKModule } from '../loader.js';
import {
  allocateDoubles,
  allocateInts,
  readDoubles,
  readInt,
  freeAll,
  prepareMatrix,
  getMatrixDimensions,
  getLapackErrorMessage,
  CHAR,
} from '../helpers.js';
import type { Matrix, EigSelectOptions, EigSelectResult } from './types.js';

/**
 * Compute selected eigenvalues and eigenvectors of a real symmetric n×n matrix A.
 *
 * Uses the Relatively Robust Representations (RRR) algorithm which can efficiently
 * compute a subset of eigenvalues.
 *
 * @param A - Input symmetric matrix (n × n). Only the specified triangle is accessed.
 * @param options - Selection and computation options
 * @returns Selected eigenvalues (in ascending order) and optionally eigenvectors
 *
 * @example
 * ```typescript
 * // Compute all eigenvalues
 * const A = [[2, 1, 0], [1, 2, 1], [0, 1, 2]];
 * const { values, vectors } = eigSelect(A);
 * ```
 *
 * @example
 * ```typescript
 * // Compute eigenvalues in interval (0.5, 2.5]
 * const { values, m } = eigSelect(A, { range: 'value', vl: 0.5, vu: 2.5 });
 * // m is the number of eigenvalues found
 * ```
 *
 * @example
 * ```typescript
 * // Compute the 2nd and 3rd smallest eigenvalues
 * const { values } = eigSelect(A, { range: 'index', il: 2, iu: 3 });
 * ```
 */
export function eigSelect(A: Matrix, options: EigSelectOptions = {}): EigSelectResult {
  const Module = getLAPACKModule();

  const {
    range = 'all',
    vl = 0,
    vu = 0,
    il = 1,
    iu = 1,
    computeVectors = true,
    uplo = 'lower',
    abstol = 0,
  } = options;

  // Get dimensions
  const [m, n] = getMatrixDimensions(A);
  if (m !== n) {
    throw new Error(`Matrix A must be square, got ${m}×${n}`);
  }

  // Prepare matrix
  const aData = prepareMatrix(A);

  // Determine range character
  let rangeChar: number;
  if (range === 'all') {
    rangeChar = CHAR.A;
  } else if (range === 'value') {
    rangeChar = CHAR.V;
  } else {
    rangeChar = 73; // 'I' for index
  }

  // Maximum number of eigenvalues
  let maxEig = n;
  if (range === 'index') {
    maxEig = iu - il + 1;
  }

  // Allocate WASM memory
  const aPtr = allocateDoubles(Module, aData, n * n);
  const wPtr = allocateDoubles(Module, null, n);
  const zPtr = allocateDoubles(Module, null, computeVectors ? n * maxEig : 1);
  const isuppzPtr = allocateInts(Module, null, 2 * n);
  const mPtr = allocateInts(Module, [0]);
  const infoPtr = allocateInts(Module, [0]);

  // Parameter pointers
  const jobzPtr = allocateInts(Module, [computeVectors ? CHAR.V : CHAR.N]);
  const rangePtr = allocateInts(Module, [rangeChar]);
  const uploPtr = allocateInts(Module, [uplo === 'upper' ? CHAR.U : CHAR.L]);
  const nPtr = allocateInts(Module, [n]);
  const ldaPtr = allocateInts(Module, [n]);
  const vlPtr = allocateDoubles(Module, [vl]);
  const vuPtr = allocateDoubles(Module, [vu]);
  const ilPtr = allocateInts(Module, [il]);
  const iuPtr = allocateInts(Module, [iu]);
  const abstolPtr = allocateDoubles(Module, [abstol]);
  const ldzPtr = allocateInts(Module, [computeVectors ? n : 1]);

  // Workspace query pointers
  const lworkPtr = allocateInts(Module, [-1]);
  const liworkPtr = allocateInts(Module, [-1]);
  const workQueryPtr = allocateDoubles(Module, null, 1);
  const iworkQueryPtr = allocateInts(Module, [0]);

  try {
    // Workspace query
    Module._dsyevr_(
      jobzPtr, rangePtr, uploPtr, nPtr, aPtr, ldaPtr,
      vlPtr, vuPtr, ilPtr, iuPtr, abstolPtr, mPtr,
      wPtr, zPtr, ldzPtr, isuppzPtr,
      workQueryPtr, lworkPtr, iworkQueryPtr, liworkPtr, infoPtr
    );

    let info = readInt(Module, infoPtr);
    if (info !== 0) {
      return createErrorResult(n, info, 'DSYEVR');
    }

    const lwork = Math.ceil(Module.HEAPF64[workQueryPtr >> 3]);
    const liwork = Module.HEAP32[iworkQueryPtr >> 2];

    const workPtr = allocateDoubles(Module, null, lwork);
    const iworkPtr = allocateInts(Module, null, liwork);

    Module.HEAP32[lworkPtr >> 2] = lwork;
    Module.HEAP32[liworkPtr >> 2] = liwork;

    // Compute eigenvalues and eigenvectors
    Module._dsyevr_(
      jobzPtr, rangePtr, uploPtr, nPtr, aPtr, ldaPtr,
      vlPtr, vuPtr, ilPtr, iuPtr, abstolPtr, mPtr,
      wPtr, zPtr, ldzPtr, isuppzPtr,
      workPtr, lworkPtr, iworkPtr, liworkPtr, infoPtr
    );

    info = readInt(Module, infoPtr);
    freeAll(Module, [workPtr, iworkPtr]);

    if (info !== 0) {
      return createErrorResult(n, info, 'DSYEVR');
    }

    // Read number of eigenvalues found
    const numFound = readInt(Module, mPtr);

    // Read results
    const result: EigSelectResult = {
      values: readDoubles(Module, wPtr, numFound),
      m: numFound,
      n,
      info: 0,
      success: true,
      message: 'Success',
    };

    if (computeVectors && numFound > 0) {
      result.vectors = readDoubles(Module, zPtr, n * numFound);
    }

    return result;
  } finally {
    freeAll(Module, [
      aPtr, wPtr, zPtr, isuppzPtr, mPtr, infoPtr,
      jobzPtr, rangePtr, uploPtr, nPtr, ldaPtr,
      vlPtr, vuPtr, ilPtr, iuPtr, abstolPtr, ldzPtr,
      lworkPtr, liworkPtr, workQueryPtr, iworkQueryPtr
    ]);
  }
}

function createErrorResult(n: number, info: number, routine: string): EigSelectResult {
  return {
    values: new Float64Array(0),
    m: 0,
    n,
    info,
    success: false,
    message: getLapackErrorMessage(routine, info),
  };
}
