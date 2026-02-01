/**
 * Compute divergence at a point.
 *
 * @module gridfunction
 */

import { GridFunction } from './GridFunction.js';
import type { IntegrationPoint } from './getValue.js';

/**
 * Computes the divergence of a vector grid function at a point within an element.
 *
 * The divergence is computed in physical coordinates:
 * div(F) = ∂Fx/∂x + ∂Fy/∂y + ∂Fz/∂z
 *
 * This is typically used with H(div) finite element spaces (Raviart-Thomas elements).
 *
 * @param gf - Vector GridFunction (e.g., from H(div) space)
 * @param elemIdx - Index of the element containing the point (0-based)
 * @param ip - Integration point in reference coordinates
 * @returns Divergence value (scalar)
 * @throws Error if the grid function has been destroyed
 * @throws Error if divergence computation fails
 *
 * @category GridFunction
 *
 * @example Check incompressibility (div(v) = 0)
 * ```typescript
 * const div = getDivergence(velocity, elemIdx, { x: 0.5, y: 0.5 });
 * if (Math.abs(div) < 1e-10) {
 *   console.log('Flow is approximately incompressible');
 * }
 * ```
 *
 * @example Compute charge density from electric field (ρ = ε₀ div(E))
 * ```typescript
 * const epsilon0 = 8.854e-12;
 * const divE = getDivergence(electricField, elemIdx, ip);
 * const chargeDensity = epsilon0 * divE;
 * ```
 */
export function getDivergence(
  gf: GridFunction,
  elemIdx: number,
  ip: IntegrationPoint
): number {
  const module = gf.getModule();
  const ptr = gf.getPointer();

  // Determine dimension from ip
  const dim = ip.z !== undefined ? 3 : ip.y !== undefined ? 2 : 1;

  // Allocate memory for coordinates
  const xiPtr = module._malloc(dim * 8);
  if (xiPtr === 0) {
    throw new Error('Failed to allocate memory for coordinates');
  }

  try {
    // Copy coordinates to WASM memory
    module.HEAPF64[xiPtr >> 3] = ip.x;
    if (dim > 1) {
      module.HEAPF64[(xiPtr >> 3) + 1] = ip.y ?? 0;
    }
    if (dim > 2) {
      module.HEAPF64[(xiPtr >> 3) + 2] = ip.z ?? 0;
    }

    const value = module._mfem_gridfunc_get_divergence(ptr, elemIdx, xiPtr, dim);

    if (Number.isNaN(value)) {
      throw new Error(`Failed to compute divergence at element ${elemIdx}`);
    }

    return value;
  } finally {
    module._free(xiPtr);
  }
}
