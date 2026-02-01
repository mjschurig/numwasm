/**
 * Element volume computation.
 *
 * @module mesh
 */

import { getMFEMModule } from '../loader.js';
import { Mesh } from './Mesh.js';

/**
 * Computes the volume (or area/length) of each element in the mesh.
 *
 * Returns the measure of each element according to the mesh dimension:
 * - 1D: element lengths
 * - 2D: element areas
 * - 3D: element volumes
 *
 * This is useful for computing integrals, checking mesh quality,
 * or applying element-based operations.
 *
 * @param mesh - The mesh to compute element volumes for
 * @returns Float64Array where index i contains the volume of element i
 * @throws Error if the MFEM module is not loaded
 * @throws Error if memory allocation fails
 *
 * @category Mesh
 *
 * @example Get volumes of all elements
 * ```typescript
 * const mesh = await makeCartesian2D(5, 5);
 * const volumes = getElementVolumes(mesh);
 * console.log(volumes.length); // 25
 * console.log(volumes[0]); // 0.04 (each element is 0.2×0.2)
 * mesh.destroy();
 * ```
 *
 * @example Compute total domain volume
 * ```typescript
 * const mesh = await makeCartesian3D(4, 4, 4);
 * const volumes = getElementVolumes(mesh);
 * const totalVolume = volumes.reduce((a, b) => a + b, 0);
 * console.log(totalVolume); // ~1.0 (unit cube)
 * mesh.destroy();
 * ```
 *
 * @example Check mesh quality (volume variation)
 * ```typescript
 * const mesh = await makeDisk(100);
 * const volumes = getElementVolumes(mesh);
 * const minVol = Math.min(...volumes);
 * const maxVol = Math.max(...volumes);
 * console.log(`Volume ratio: ${maxVol / minVol}`);
 * // Lower ratio indicates more uniform mesh quality
 * mesh.destroy();
 * ```
 *
 * @example Identify small elements
 * ```typescript
 * const mesh = await makeCartesian2D(10, 10);
 * refineLocal(mesh, [0, 1, 2]); // Refine some elements
 * const volumes = getElementVolumes(mesh);
 * const threshold = 0.001;
 * const smallElements = volumes
 *   .map((v, i) => ({ index: i, volume: v }))
 *   .filter(e => e.volume < threshold);
 * console.log(`${smallElements.length} elements below threshold`);
 * mesh.destroy();
 * ```
 *
 * @see {@link getElementCentroids} - Get element center coordinates
 */
export function getElementVolumes(mesh: Mesh): Float64Array {
  const module = getMFEMModule();
  if (!module) {
    throw new Error('MFEM module not loaded');
  }

  const ne = mesh.numElements;
  const volumesPtr = module._malloc(ne * 8);

  if (volumesPtr === 0) {
    throw new Error('Failed to allocate memory for volumes');
  }

  try {
    const result = module._mfem_mesh_get_element_volumes(mesh.getPointer(), volumesPtr);
    if (result !== 0) {
      throw new Error('Failed to get element volumes');
    }

    const volumes = new Float64Array(ne);
    for (let i = 0; i < ne; i++) {
      volumes[i] = module.HEAPF64[(volumesPtr >> 3) + i];
    }
    return volumes;
  } finally {
    module._free(volumesPtr);
  }
}
