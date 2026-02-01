/**
 * Eigenvalue Problems
 *
 * High-level functions for computing eigenvalues and eigenvectors.
 */

// Export all eigenvalue functions
export { eig } from './eig.js';
export { eigvals } from './eigvals.js';
export { eigSymmetric } from './eig-symmetric.js';
export { eigHermitian } from './eig-hermitian.js';
export { eigGeneralized } from './eig-generalized.js';
export { eigGeneralizedSymmetric } from './eig-generalized-symmetric.js';
export { eigBanded } from './eig-banded.js';
export { eigTridiagonal } from './eig-tridiagonal.js';
export { eigSelect } from './eig-select.js';

// Export all types
export type {
  Matrix,
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
} from './types.js';
