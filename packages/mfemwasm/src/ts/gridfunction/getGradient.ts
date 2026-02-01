/**
 * Compute gradient at a point.
 *
 * @module gridfunction
 */

import { GridFunction } from './GridFunction.js';
import type { IntegrationPoint } from './getValue.js';

/**
 * Computes the gradient of a scalar grid function at a point within an element.
 *
 * The gradient is computed in physical (not reference) coordinates. The result
 * is a vector with components [∂f/∂x, ∂f/∂y, ∂f/∂z] (depending on dimension).
 *
 * @param gf - Scalar GridFunction to differentiate
 * @param elemIdx - Index of the element containing the point (0-based)
 * @param ip - Integration point in reference coordinates
 * @returns Gradient vector in physical coordinates
 * @throws Error if the grid function has been destroyed
 * @throws Error if gradient computation fails
 *
 * @category GridFunction
 *
 * @example Compute gradient at element center
 * ```typescript
 * const grad = getGradient(temperature, 0, { x: 0.5, y: 0.5 });
 * console.log(`Temperature gradient: [${grad.join(', ')}]`);
 * ```
 *
 * @example Compute heat flux (q = -k * grad(T))
 * ```typescript
 * const k = 1.5; // thermal conductivity
 * const grad = getGradient(temperature, elemIdx, ip);
 * const flux = grad.map(g => -k * g);
 * ```
 */
export function getGradient(
  gf: GridFunction,
  elemIdx: number,
  ip: IntegrationPoint
): number[] {
  const module = gf.getModule();
  const ptr = gf.getPointer();

  // Determine dimension from ip
  const dim = ip.z !== undefined ? 3 : ip.y !== undefined ? 2 : 1;

  // Allocate memory for coordinates and gradient output
  const xiPtr = module._malloc(dim * 8);
  const gradPtr = module._malloc(3 * 8); // Max 3D gradient

  if (xiPtr === 0 || gradPtr === 0) {
    if (xiPtr) module._free(xiPtr);
    if (gradPtr) module._free(gradPtr);
    throw new Error('Failed to allocate memory');
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

    const result = module._mfem_gridfunc_get_gradient(ptr, elemIdx, xiPtr, dim, gradPtr);

    if (result < 0) {
      throw new Error(`Failed to compute gradient at element ${elemIdx}`);
    }

    // Read gradient from WASM memory
    const grad: number[] = [];
    for (let i = 0; i < dim; i++) {
      grad.push(module.HEAPF64[(gradPtr >> 3) + i]);
    }

    return grad;
  } finally {
    module._free(xiPtr);
    module._free(gradPtr);
  }
}
