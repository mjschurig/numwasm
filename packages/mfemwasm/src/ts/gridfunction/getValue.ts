/**
 * Evaluate grid function at a point.
 *
 * @module gridfunction
 */

import { GridFunction } from './GridFunction.js';

/**
 * Integration point coordinates in the reference element.
 *
 * @category GridFunction
 */
export interface IntegrationPoint {
  /** X coordinate in reference element [0, 1] */
  x: number;
  /** Y coordinate in reference element [0, 1] (for 2D/3D) */
  y?: number;
  /** Z coordinate in reference element [0, 1] (for 3D) */
  z?: number;
}

/**
 * Evaluates a grid function at a point within a specific element.
 *
 * The point is specified in reference (local) coordinates of the element.
 * For most elements, these coordinates range from 0 to 1.
 *
 * @param gf - GridFunction to evaluate
 * @param elemIdx - Index of the element containing the point (0-based)
 * @param ip - Integration point in reference coordinates
 * @returns The scalar value at the point
 * @throws Error if the grid function has been destroyed
 * @throws Error if evaluation fails
 *
 * @category GridFunction
 *
 * @example Evaluate at element center
 * ```typescript
 * // For a 2D element, center is at (0.5, 0.5) in reference coords
 * const value = getValue(gf, 0, { x: 0.5, y: 0.5 });
 * console.log(`Value at center of element 0: ${value}`);
 * ```
 *
 * @example Sample along an element edge
 * ```typescript
 * for (let t = 0; t <= 1; t += 0.1) {
 *   const v = getValue(gf, elemIdx, { x: t, y: 0 });
 *   console.log(`Value at (${t}, 0): ${v}`);
 * }
 * ```
 */
export function getValue(
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

    const value = module._mfem_gridfunc_get_value(ptr, elemIdx, xiPtr, dim);

    if (Number.isNaN(value)) {
      throw new Error(`Failed to evaluate grid function at element ${elemIdx}`);
    }

    return value;
  } finally {
    module._free(xiPtr);
  }
}
