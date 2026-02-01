/**
 * ARWasm - ARPACK WebAssembly Module
 *
 * This package provides a WebAssembly build of ARPACK (ARnoldi PACKage),
 * a collection of Fortran77 subroutines designed to solve large scale
 * eigenvalue problems.
 *
 * ARPACK provides routines for:
 * - Real symmetric eigenvalue problems (dsaupd/dseupd, ssaupd/sseupd)
 * - Real nonsymmetric eigenvalue problems (dnaupd/dneupd, snaupd/sneupd)
 * - Complex eigenvalue problems (cnaupd/cneupd, znaupd/zneupd)
 *
 * ## High-Level API
 *
 * For most use cases, use the high-level functions that handle all the
 * complexity of the reverse communication interface:
 *
 * @example
 * ```typescript
 * import { eigs } from 'arwasm';
 *
 * // Find the 6 smallest eigenvalues of a symmetric matrix
 * const n = 100;
 * const matvec = (x: Float64Array): Float64Array => {
 *   const y = new Float64Array(n);
 *   for (let i = 0; i < n; i++) {
 *     y[i] = 2 * x[i];
 *     if (i > 0) y[i] -= x[i - 1];
 *     if (i < n - 1) y[i] -= x[i + 1];
 *   }
 *   return y;
 * };
 *
 * const result = await eigs(matvec, n, 6, { which: 'SM' });
 * console.log('Smallest eigenvalues:', result.eigenvalues);
 * ```
 *
 * ## Low-Level API
 *
 * For advanced use cases, you can use the raw WASM module interface:
 *
 * @example
 * ```typescript
 * import { loadARPACKModule, getARPACKModule } from 'arwasm';
 *
 * await loadARPACKModule();
 * const arpack = getARPACKModule();
 * // Use reverse communication interface directly
 * ```
 *
 * @packageDocumentation
 */

// ============================================================
// HIGH-LEVEL API (recommended for most users)
// ============================================================

// Core eigenvalue solvers (real)
export { eigs } from './ts/core/eigs.js';
export { eign } from './ts/core/eign.js';
export { eigsh, isEigsResult, isEignResult } from './ts/core/eigsh.js';

// Complex eigenvalue solvers
export { zeigs } from './ts/complex/zeigs.js';
export { zeigsh } from './ts/complex/zeigsh.js';

// Singular value decomposition
export { svds } from './ts/svd/svds.js';

// Generalized eigenvalue problem
export { geigs } from './ts/generalized/geigs.js';

// Shift-invert convenience functions
export { eigsNear } from './ts/generalized/eigsNear.js';
export { eignNear } from './ts/generalized/eignNear.js';

// Graph/network analysis
export { laplacianEigs } from './ts/graph/laplacianEigs.js';
export { pagerank, pagerankEigs } from './ts/graph/pagerank.js';

// Dimensionality reduction
export { spectralEmbedding } from './ts/dimred/spectralEmbedding.js';
export { truncatedPCA, truncatedPCAfromData } from './ts/dimred/truncatedPCA.js';

// Linear algebra utilities
export { spectralRadius } from './ts/linalg/spectralRadius.js';
export { spectralNorm } from './ts/linalg/spectralNorm.js';
export { condest } from './ts/linalg/condest.js';
export { nuclearNormApprox } from './ts/linalg/nuclearNorm.js';

// Sparse matrix helpers
export { csrMatvec, csrMatvecT } from './ts/sparse/csrMatvec.js';
export { cscMatvec, cscMatvecT } from './ts/sparse/cscMatvec.js';
export { cooMatvec, cooMatvecT } from './ts/sparse/cooMatvec.js';
export { denseMatvec, denseMatvecT } from './ts/sparse/denseMatvec.js';
export {
  diagMatvec,
  diagMatvecInv,
  diagMatvecSqrt,
  diagMatvecInvSqrt,
} from './ts/sparse/diagMatvec.js';
export {
  tridiagMatvec,
  symTridiagMatvec,
  toeplitzTridiagMatvec,
} from './ts/sparse/tridiagMatvec.js';
export {
  bandedMatvec,
  symBandedMatvec,
  toeplitzBandedMatvec,
} from './ts/sparse/bandedMatvec.js';

// Operator combinators
export {
  addMatvec,
  mulMatvec,
  shiftMatvec,
  scaleMatvec,
  transposeMatvec,
  symmetrizeMatvec,
  identityMatvec,
  negateMatvec,
  powerMatvec,
  blockDiagMatvec,
} from './ts/operators/combinators.js';

// Validation and diagnostics
export {
  verifyEigs,
  verifyEign,
  verifySvds,
  checkSymmetry,
  checkPositiveDefinite,
  checkNormalization,
} from './ts/validation/verify.js';
export type {
  VerifyOptions,
  VerifyEigsResult,
  VerifySvdsResult,
  SymmetryCheckResult,
  PositiveDefiniteCheckResult,
} from './ts/validation/verify.js';

// Advanced eigenvalue modes
export { bucklingEigs, criticalBucklingLoad } from './ts/advanced/bucklingEigs.js';
export type { BucklingEigsOptions, BucklingResult } from './ts/advanced/bucklingEigs.js';
export { cayleyEigs, eigsInInterval } from './ts/advanced/cayleyEigs.js';
export type { CayleyEigsOptions } from './ts/advanced/cayleyEigs.js';

// Matrix functions
export { expmv, expmvMultiple } from './ts/matfun/expmv.js';
export type { ExpmvOptions, ExpmvResult } from './ts/matfun/expmv.js';
export { sqrtmv, invsqrtmv, matpowv } from './ts/matfun/sqrtmv.js';
export type { SqrtmvOptions, SqrtmvResult } from './ts/matfun/sqrtmv.js';

// Continuation and deflation
export { eigsDeflated, eigsContinue, eignDeflated } from './ts/continuation/eigsDeflated.js';
export type { EigsDeflatedOptions, EignDeflatedOptions } from './ts/continuation/eigsDeflated.js';

// High-level type definitions
export type {
  RealArray,
  ComplexArray,
  MatVecFunction,
  OperatorFunction,
  BMatVecFunction,
  ComplexMatVecFunction,
  ComplexOperatorFunction,
  ComplexBMatVecFunction,
  EigOptionsBase,
  ComplexEigOptionsBase,
  EigsOptions,
  EignOptions,
  EigsResult,
  EignResult,
  EigshOptions,
  ProblemType,
  ZeigsOptions,
  ZeigshOptions,
  ZeigsResult,
  ZeigshResult,
  SvdsOptions,
  SvdsResult,
  GeigsOptions,
  Complex,
  EigsNearOptions,
  EignNearOptions,
  // Graph analysis types
  LaplacianEigsOptions,
  LaplacianEigsResult,
  PagerankOptions,
  PagerankResult,
  // Dimensionality reduction types
  SpectralEmbeddingOptions,
  SpectralEmbeddingResult,
  TruncatedPCAOptions,
  TruncatedPCAResult,
  // Linear algebra types
  SpectralRadiusOptions,
  SpectralRadiusResult,
  SpectralNormOptions,
  SpectralNormResult,
  CondestOptions,
  CondestResult,
  NuclearNormOptions,
  NuclearNormResult,
} from './ts/high-level-types.js';

// Helper utilities
export {
  dsaupdWorklSize,
  dnaupdWorklSize,
  dneupdWorkevSize,
  znaupdWorklSize,
  zneupdWorkevSize,
  znaupdRworkSize,
  defaultNcv,
  getDsaupdMessage,
  getDnaupdMessage,
  getDseupdMessage,
  getDneupdMessage,
  getZnaupdMessage,
  getZneupdMessage,
} from './ts/helpers.js';

// ============================================================
// LOW-LEVEL API (for advanced usage)
// ============================================================

// Module loader functions
export {
  loadARPACKModule,
  getARPACKModule,
  isARPACKLoaded,
  resetARPACKModule,
  configureARPACK,
  type ARPACKLoadConfig,
} from './ts/loader.js';

// Low-level module types
export type {
  ARPACKModule,
  ARPACKModuleFactory,
  ARPACKModuleOptions,
  WhichSymmetric,
  WhichNonSymmetric,
  WhichComplex,
} from './ts/types.js';

// Error code mappings and constants
export {
  DSAUPD_ERRORS,
  DSEUPD_ERRORS,
  DNAUPD_ERRORS,
  DNEUPD_ERRORS,
  ZNAUPD_ERRORS,
  ZNEUPD_ERRORS,
  WHICH_SYMMETRIC,
  WHICH_NONSYMMETRIC,
  WHICH_COMPLEX,
} from './ts/types.js';

// ============================================================
// NAMESPACE EXPORTS (for organized module access)
// ============================================================

/** Core eigenvalue solvers for real matrices */
export * as core from './ts/core/index.js';

/** Complex eigenvalue solvers */
export * as complex from './ts/complex/index.js';

/** Singular value decomposition */
export * as svd from './ts/svd/index.js';

/** Generalized eigenvalue problems */
export * as generalized from './ts/generalized/index.js';

/** Graph and network analysis */
export * as graph from './ts/graph/index.js';

/** Dimensionality reduction */
export * as dimred from './ts/dimred/index.js';

/** Linear algebra utilities */
export * as linalg from './ts/linalg/index.js';

/** Sparse matrix helpers */
export * as sparse from './ts/sparse/index.js';

/** Operator combinators */
export * as operators from './ts/operators/index.js';

/** Validation and diagnostics */
export * as validation from './ts/validation/index.js';

/** Advanced eigenvalue modes */
export * as advanced from './ts/advanced/index.js';

/** Matrix functions */
export * as matfun from './ts/matfun/index.js';

/** Continuation and deflation */
export * as continuation from './ts/continuation/index.js';
