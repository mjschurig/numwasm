/**
 * MFEM format mesh export.
 *
 * @module mesh
 */

import { Mesh } from './Mesh.js';

/**
 * Exports a mesh to MFEM native format string.
 *
 * Generates a string representation of the mesh in MFEM's native format,
 * which can be saved to a file or used with other MFEM tools. This format
 * preserves all mesh information including high-order geometry nodes.
 *
 * @param mesh - The mesh to export
 * @returns MFEM format mesh string
 * @throws Error if memory allocation fails
 * @throws Error if export fails
 *
 * @category Mesh
 *
 * @example Export mesh to string
 * ```typescript
 * const mesh = await makeCartesian2D(4, 4);
 * const mfemStr = toMFEM(mesh);
 * console.log(mfemStr.substring(0, 100)); // Preview first 100 chars
 * mesh.destroy();
 * ```
 *
 * @example Save mesh to file (Node.js)
 * ```typescript
 * import { writeFileSync } from 'fs';
 *
 * const mesh = await makeCartesian3D(8, 8, 8);
 * refineUniform(mesh);
 * const mfemStr = toMFEM(mesh);
 * writeFileSync('output.mesh', mfemStr);
 * mesh.destroy();
 * ```
 *
 * @example Round-trip: export and reimport
 * ```typescript
 * const mesh1 = await makeCartesian2D(5, 5);
 * refineLocal(mesh1, [0, 1, 2]);
 *
 * // Export to string
 * const meshStr = toMFEM(mesh1);
 *
 * // Re-import (using Mesh.fromString)
 * const mesh2 = await Mesh.fromString(meshStr);
 *
 * console.log(mesh1.numElements === mesh2.numElements); // true
 * mesh1.destroy();
 * mesh2.destroy();
 * ```
 *
 * @example Export high-order mesh
 * ```typescript
 * const mesh = await makeDisk(50);
 * setCurvature(mesh, 3); // Cubic geometry
 * const mfemStr = toMFEM(mesh);
 * // The string includes the high-order nodes
 * mesh.destroy();
 * ```
 *
 * @see {@link toVTK} - Export to VTK format
 * @see {@link toGmsh} - Export to Gmsh format
 * @see {@link fromGmsh} - Import from Gmsh format
 */
export function toMFEM(mesh: Mesh): string {
  const module = mesh.getModule();
  const ptr = mesh.getPointer();

  // Allocate a large buffer for the output string
  const maxLen = 1024 * 1024; // 1MB should be sufficient for most meshes
  const outPtr = module._malloc(maxLen);

  if (outPtr === 0) {
    throw new Error('Failed to allocate memory for mesh export');
  }

  try {
    const result = module._mfem_mesh_to_mfem_string(ptr, outPtr, maxLen);

    if (result < 0) {
      throw new Error('Failed to export mesh to MFEM format');
    }

    // Read the null-terminated string from WASM memory
    let str = '';
    let i = 0;
    while (i < result) {
      const char = module.HEAPU8[outPtr + i];
      if (char === 0) break;
      str += String.fromCharCode(char);
      i++;
    }

    return str;
  } finally {
    module._free(outPtr);
  }
}
