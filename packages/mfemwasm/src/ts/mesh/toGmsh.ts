/**
 * Gmsh format mesh export.
 *
 * @module mesh
 */

import { Mesh } from './Mesh.js';

/**
 * Exports a mesh to Gmsh MSH 2.2 format string.
 *
 * Generates a string representation of the mesh in Gmsh MSH format,
 * which can be used with Gmsh or other tools that support this format.
 *
 * **Note:** This is a simplified implementation that exports node coordinates
 * and basic element information. Full element connectivity is not yet supported.
 * For complete mesh export, prefer {@link toVTK} or {@link toMFEM}.
 *
 * @param mesh - The mesh to export
 * @returns Gmsh MSH 2.2 format mesh string
 * @throws Error if memory allocation fails
 *
 * @category Mesh
 *
 * @example Export mesh to Gmsh format
 * ```typescript
 * const mesh = await makeCartesian2D(5, 5, 'tri');
 * const gmshStr = toGmsh(mesh);
 * console.log(gmshStr.substring(0, 300)); // Preview
 * mesh.destroy();
 * ```
 *
 * @example Save to file (Node.js)
 * ```typescript
 * import { writeFileSync } from 'fs';
 *
 * const mesh = await makeCartesian3D(4, 4, 4, 'tet');
 * const gmshStr = toGmsh(mesh);
 * writeFileSync('mesh.msh', gmshStr);
 * mesh.destroy();
 * ```
 *
 * @example Export with attributes
 * ```typescript
 * const mesh = await makeCartesian2D(8, 8);
 * // Elements have attributes that are preserved in export
 * const gmshStr = toGmsh(mesh);
 * // The $Elements section includes element attributes
 * mesh.destroy();
 * ```
 *
 * @see {@link toMFEM} - Export to MFEM format (recommended for MFEM workflows)
 * @see {@link toVTK} - Export to VTK format (recommended for visualization)
 * @see {@link fromGmsh} - Import from Gmsh format
 */
export function toGmsh(mesh: Mesh): string {
  const module = mesh.getModule();
  const ptr = mesh.getPointer();

  const dim = mesh.dimension;
  const nv = mesh.numVertices;
  const ne = mesh.numElements;
  const sdim = module._mfem_mesh_get_space_dimension(ptr);

  // Build Gmsh MSH 2.2 format output
  let output = '';

  // Header
  output += '$MeshFormat\n';
  output += '2.2 0 8\n';
  output += '$EndMeshFormat\n';

  // Nodes section
  output += '$Nodes\n';
  output += `${nv}\n`;

  const coordsPtr = module._malloc(sdim * 8);
  if (coordsPtr === 0) {
    throw new Error('Failed to allocate memory');
  }

  try {
    for (let i = 0; i < nv; i++) {
      module._mfem_mesh_get_vertex(ptr, i, coordsPtr);
      const x = module.HEAPF64[(coordsPtr >> 3)];
      const y = sdim > 1 ? module.HEAPF64[(coordsPtr >> 3) + 1] : 0;
      const z = sdim > 2 ? module.HEAPF64[(coordsPtr >> 3) + 2] : 0;
      output += `${i + 1} ${x} ${y} ${z}\n`;
    }
  } finally {
    module._free(coordsPtr);
  }

  output += '$EndNodes\n';

  // Elements section
  output += '$Elements\n';
  output += `${ne}\n`;

  // Element type mapping (MFEM geometry type to Gmsh element type)
  // Gmsh types: 1=line, 2=triangle, 3=quad, 4=tetrahedron, 5=hexahedron
  const elementVerticesPtr = module._malloc(8 * 4); // Max 8 vertices per element

  try {
    for (let i = 0; i < ne; i++) {
      // Get element vertices - we need to read from mesh
      // For now, output a placeholder based on dimension
      const attr = module._mfem_mesh_get_attribute(ptr, i);

      // Get element type based on dimension (simplified)
      let gmshType: number;
      if (dim === 1) {
        gmshType = 1; // Line
      } else if (dim === 2) {
        gmshType = 2; // Triangle (assuming triangular mesh)
      } else {
        gmshType = 4; // Tetrahedron (assuming tet mesh)
      }

      // For a proper implementation, we would need to get element vertices
      // This is a simplified version that outputs element info
      output += `${i + 1} ${gmshType} 2 ${attr} ${attr}`;

      // Note: A full implementation would read actual vertex indices
      // using mfem_mesh_get_element_vertices if available
      output += '\n';
    }
  } finally {
    module._free(elementVerticesPtr);
  }

  output += '$EndElements\n';

  return output;
}
