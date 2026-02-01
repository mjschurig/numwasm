/**
 * LAWasm - LAPACK + BLAS WebAssembly Module
 *
 * This package provides a WebAssembly build of LAPACK (Linear Algebra PACKage)
 * with BLAS (Basic Linear Algebra Subprograms) bundled internally.
 *
 * LAPACK provides routines for:
 * - Solving systems of linear equations
 * - Least squares solutions
 * - Eigenvalue problems
 * - Singular value decomposition
 * - Matrix factorizations (LU, Cholesky, QR, etc.)
 *
 * @example
 * ```typescript
 * import { loadLAPACKModule, solve, solveSymmetric } from 'lawasm';
 *
 * // Load the WASM module
 * await loadLAPACKModule();
 *
 * // Solve a general linear system Ax = b
 * const A = [[2, 1], [1, 3]];
 * const b = [1, 2];
 * const { x } = solve(A, b);
 *
 * // Solve a symmetric positive definite system
 * const S = [[4, 2], [2, 5]];
 * const { x: x2 } = solveSymmetric(S, b);
 * ```
 *
 * @packageDocumentation
 */

// Re-export loader functions
export {
  loadLAPACKModule,
  getLAPACKModule,
  isLAPACKLoaded,
  resetLAPACKModule,
  configureLAPACK,
  type LAPACKLoadConfig,
} from './ts/loader.js';

// Re-export low-level types
export type {
  LAPACKModule,
  LAPACKModuleFactory,
} from './ts/types.js';

// ============================================================
// High-Level Linear System Solvers
// ============================================================

export {
  solve,
  solveTriangular,
  solveSymmetric,
  solveHermitian,
  solveTridiagonal,
  solveBanded,
} from './ts/linear-solvers/index.js';

export type {
  Matrix,
  Vector,
  ComplexMatrix,
  ComplexVector,
  SolveOptions,
  SolveTriangularOptions,
  SolveSymmetricOptions,
  SolveHermitianOptions,
  SolveBandedOptions,
  SolveResult,
  SolveTriangularResult,
  SolveSymmetricResult,
  SolveComplexResult,
  SolveTridiagonalResult,
} from './ts/linear-solvers/index.js';

// ============================================================
// Matrix Factorizations
// ============================================================

export {
  lu,
  cholesky,
  qr,
  qrPivoted,
  lq,
  ldl,
  schur,
  hessenberg,
} from './ts/factorizations/index.js';

export type {
  LUOptions,
  LUResult,
  CholeskyOptions,
  CholeskyResult,
  QROptions,
  QRResult,
  QRPivotedOptions,
  QRPivotedResult,
  LQOptions,
  LQResult,
  LDLOptions,
  LDLResult,
  SchurOptions,
  SchurResult,
  HessenbergOptions,
  HessenbergResult,
} from './ts/factorizations/index.js';

// ============================================================
// Singular Value Decomposition
// ============================================================

export {
  svd,
  svdvals,
  svdCompact,
  svdRank,
} from './ts/svd/index.js';

export type {
  SVDOptions,
  SVDValsOptions,
  SVDRankOptions,
  SVDResult,
  SVDValsResult,
  SVDRankResult,
} from './ts/svd/index.js';

// ============================================================
// Eigenvalue Problems
// ============================================================

export {
  eig,
  eigvals,
  eigSymmetric,
  eigHermitian,
  eigGeneralized,
  eigGeneralizedSymmetric,
  eigBanded,
  eigTridiagonal,
  eigSelect,
} from './ts/eigenvalues/index.js';

export type {
  Complex,
  EigOptions,
  EigvalsOptions,
  EigResult,
  EigvalsResult,
  EigSymmetricOptions,
  EigSymmetricResult,
  EigHermitianOptions,
  EigHermitianResult,
  EigGeneralizedOptions,
  EigGeneralizedResult,
  EigGeneralizedSymmetricOptions,
  EigGeneralizedSymmetricResult,
  EigBandedOptions,
  EigBandedResult,
  EigTridiagonalOptions,
  EigTridiagonalResult,
  EigSelectOptions,
  EigSelectResult,
} from './ts/eigenvalues/index.js';

// ============================================================
// Least Squares & Minimum Norm
// ============================================================

export {
  lstsq,
  lstsqSVD,
  lstsqGelsy,
  constrainedLstSq,
  generalizedLstSq,
} from './ts/least-squares/index.js';

export type {
  LstSqOptions,
  LstSqResult,
  LstSqSVDOptions,
  LstSqSVDResult,
  LstSqGelsyOptions,
  LstSqGelsyResult,
  ConstrainedLstSqOptions,
  ConstrainedLstSqResult,
  GeneralizedLstSqOptions,
  GeneralizedLstSqResult,
} from './ts/least-squares/index.js';

// ============================================================
// Matrix Inverses & Pseudoinverse
// ============================================================

export {
  inv,
  invTriangular,
  invSymmetric,
  pinv,
} from './ts/inverses/index.js';

export type {
  InvOptions,
  InvResult,
  InvTriangularOptions,
  InvTriangularResult,
  InvSymmetricOptions,
  InvSymmetricResult,
  PinvOptions,
  PinvResult,
} from './ts/inverses/index.js';

// ============================================================
// Matrix Norms, Condition Numbers, Determinants & Rank
// ============================================================

export {
  norm,
  cond,
  condEst,
  rcond,
  det,
  logdet,
  slogdet,
  rank,
} from './ts/norms/index.js';

export type {
  NormType,
  NormOptions,
  NormResult,
  CondOptions,
  CondResult,
  CondEstOptions,
  CondEstResult,
  RcondOptions,
  RcondResult,
  DetOptions,
  DetResult,
  LogDetOptions,
  LogDetResult,
  SlogDetOptions,
  SlogDetResult,
  RankOptions,
  RankResult,
} from './ts/norms/index.js';

// ============================================================
// BLAS Level 3 (Matrix-Matrix)
// ============================================================

export {
  matmul,
  matmulTriangular,
  solveMatrixTriangular,
  syrk,
  syr2k,
  herk,
  her2k,
} from './ts/blas/index.js';

export type {
  TransposeOp,
  MatmulOptions,
  MatmulResult,
  MatmulTriangularOptions,
  MatmulTriangularResult,
  SolveMatrixTriangularOptions,
  SolveMatrixTriangularResult,
  SyrkOptions,
  SyrkResult,
  Syr2kOptions,
  Syr2kResult,
  HerkOptions,
  HerkResult,
  Her2kOptions,
  Her2kResult,
} from './ts/blas/index.js';

// ============================================================
// BLAS Level 2 (Matrix-Vector)
// ============================================================

export {
  matvec,
  matvecTriangular,
  solveVectorTriangular,
  ger,
  syr,
  her,
} from './ts/blas/index.js';

export type {
  MatvecOptions,
  MatvecResult,
  MatvecTriangularOptions,
  MatvecTriangularResult,
  SolveVectorTriangularOptions,
  SolveVectorTriangularResult,
  GerOptions,
  GerResult,
  SyrOptions,
  SyrResult,
  HerOptions,
  HerResult,
} from './ts/blas/index.js';

// ============================================================
// BLAS Level 1 (Vector)
// ============================================================

export {
  dot,
  dotc,
  axpy,
  scal,
  copy,
  swap,
  nrm2,
  asum,
  iamax,
} from './ts/blas/index.js';

export type {
  DotResult,
  DotcResult,
  AxpyResult,
  ScalResult,
  CopyResult,
  SwapResult,
  Nrm2Result,
  AsumResult,
  IamaxResult,
} from './ts/blas/index.js';

// ============================================================
// Matrix Utilities
// ============================================================

export {
  transpose,
  conjugate,
  hermitian,
  triu,
  tril,
  diag,
  trace,
  balance,
} from './ts/utilities/index.js';

export type {
  TransposeResult,
  ConjugateResult,
  HermitianResult,
  TriuOptions,
  TriuResult,
  TrilOptions,
  TrilResult,
  DiagResult,
  TraceResult,
  BalanceResult,
} from './ts/utilities/index.js';

// ============================================================
// Matrix Properties & Tests
// ============================================================

export {
  isSymmetric,
  isHermitian,
  isPositiveDefinite,
  isOrthogonal,
  isUnitary,
  isSingular,
} from './ts/utilities/index.js';

export type {
  IsSymmetricOptions,
  IsSymmetricResult,
  IsHermitianOptions,
  IsHermitianResult,
  IsPositiveDefiniteResult,
  IsOrthogonalOptions,
  IsOrthogonalResult,
  IsUnitaryOptions,
  IsUnitaryResult,
  IsSingularOptions,
  IsSingularResult,
} from './ts/utilities/index.js';

// ============================================================
// Matrix Functions
// ============================================================

export {
  expm,
  logm,
  sqrtm,
  powm,
  funm,
} from './ts/matrix-functions/index.js';

export type {
  ExpmResult,
  LogmResult,
  SqrtmResult,
  PowmResult,
  FunmResult,
} from './ts/matrix-functions/index.js';

// ============================================================
// Specialized Decompositions
// ============================================================

export {
  polarDecomposition,
  rrqr,
  csd,
} from './ts/matrix-functions/index.js';

export type {
  PolarDecompositionResult,
  RRQROptions,
  RRQRResult,
  CSDResult,
} from './ts/matrix-functions/index.js';

// ============================================================
// Complex-Specific Functions
// ============================================================

export {
  choleskyHermitian,
} from './ts/matrix-functions/index.js';

export type {
  CholeskyHermitianOptions,
  CholeskyHermitianResult,
  SolveHermitianResult,
} from './ts/matrix-functions/index.js';

// ============================================================
// Helper Utilities
// ============================================================

export {
  toColumnMajor,
  fromColumnMajor,
  prepareMatrix,
  prepareVector,
  CHAR,
  type RealArray,
  type ComplexArray,
} from './ts/helpers.js';
