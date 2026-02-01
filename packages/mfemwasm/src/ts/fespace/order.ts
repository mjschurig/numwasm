/**
 * Finite element space polynomial order.
 *
 * @module fespace
 */

import { FiniteElementSpace } from './FiniteElementSpace.js';

/**
 * Gets the maximum polynomial order of the finite element space.
 *
 * For spaces with variable order (hp-adaptivity), this returns the
 * maximum order across all elements.
 *
 * @param fespace - The finite element space
 * @returns The polynomial order
 * @throws Error if the space has been destroyed
 *
 * @category FiniteElementSpace
 *
 * @example Check space order
 * ```typescript
 * const fespace = await createH1(mesh, 3);
 * console.log(`Order: ${order(fespace)}`); // 3
 * ```
 */
export function order(fespace: FiniteElementSpace): number {
  const module = fespace.getModule();
  const ptr = fespace.getPointer();
  return module._mfem_fespace_get_order(ptr);
}
