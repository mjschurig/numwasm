/**
 * Compute curl at a point.
 *
 * @module gridfunction
 */

import { GridFunction } from './GridFunction.js';
import type { IntegrationPoint } from './getValue.js';

/**
 * Computes the curl of a vector grid function at a point within an element.
 *
 * The curl is computed in physical coordinates:
 * - In 3D: curl = [∂Fz/∂y - ∂Fy/∂z, ∂Fx/∂z - ∂Fz/∂x, ∂Fy/∂x - ∂Fx/∂y]
 * - In 2D: curl is a scalar (z-component of the 3D curl)
 *
 * This is typically used with H(curl) finite element spaces (Nedelec elements).
 *
 * @param gf - Vector GridFunction (e.g., from H(curl) space)
 * @param elemIdx - Index of the element containing the point (0-based)
 * @param ip - Integration point in reference coordinates
 * @returns Curl vector (or scalar in 2D represented as single-element array)
 * @throws Error if the grid function has been destroyed
 * @throws Error if curl computation fails
 *
 * @category GridFunction
 *
 * @example Compute magnetic field from vector potential (B = curl A)
 * ```typescript
 * const A = await GridFunction.create(hcurlSpace);
 * // ... initialize vector potential A ...
 * const B = getCurl(A, elemIdx, { x: 0.5, y: 0.5, z: 0.5 });
 * console.log(`Magnetic field: [${B.join(', ')}]`);
 * ```
 *
 * @example Check if field is curl-free
 * ```typescript
 * const curl = getCurl(electricField, 0, { x: 0.5, y: 0.5 });
 * const magnitude = Math.sqrt(curl.reduce((sum, c) => sum + c * c, 0));
 * if (magnitude < 1e-10) {
 *   console.log('Field is approximately curl-free');
 * }
 * ```
 */
export function getCurl(
  gf: GridFunction,
  elemIdx: number,
  ip: IntegrationPoint
): number[] {
  const module = gf.getModule();
  const ptr = gf.getPointer();

  // Determine dimension from ip
  const dim = ip.z !== undefined ? 3 : ip.y !== undefined ? 2 : 1;

  if (dim < 2) {
    throw new Error('Curl is only defined for 2D and 3D problems');
  }

  // Allocate memory for coordinates and curl output
  const xiPtr = module._malloc(dim * 8);
  const curlPtr = module._malloc(3 * 8); // Max 3D curl

  if (xiPtr === 0 || curlPtr === 0) {
    if (xiPtr) module._free(xiPtr);
    if (curlPtr) module._free(curlPtr);
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

    const result = module._mfem_gridfunc_get_curl(ptr, elemIdx, xiPtr, dim, curlPtr);

    if (result < 0) {
      throw new Error(`Failed to compute curl at element ${elemIdx}`);
    }

    // Read curl from WASM memory
    // In 2D, curl is a scalar (z-component); in 3D it's a 3-vector
    const curlDim = dim === 2 ? 1 : 3;
    const curl: number[] = [];
    for (let i = 0; i < curlDim; i++) {
      curl.push(module.HEAPF64[(curlPtr >> 3) + i]);
    }

    return curl;
  } finally {
    module._free(xiPtr);
    module._free(curlPtr);
  }
}
