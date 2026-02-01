/**
 * Matrix Inverses & Pseudoinverse
 *
 * High-level functions for computing matrix inverses.
 */

// Export all inverse functions
export { inv } from './inv.js';
export { invTriangular } from './inv-triangular.js';
export { invSymmetric } from './inv-symmetric.js';
export { pinv } from './pinv.js';

// Export all types
export type {
  Matrix,
  InvOptions,
  InvResult,
  InvTriangularOptions,
  InvTriangularResult,
  InvSymmetricOptions,
  InvSymmetricResult,
  PinvOptions,
  PinvResult,
} from './types.js';
