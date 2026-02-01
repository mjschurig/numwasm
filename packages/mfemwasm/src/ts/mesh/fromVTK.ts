/**
 * VTK mesh format loading.
 *
 * @module mesh
 */

import { Mesh } from './Mesh.js';

/**
 * Loads a mesh from a VTK format string.
 *
 * Parses mesh data in VTK legacy format (.vtk). This is useful for importing
 * meshes from visualization tools like ParaView or VTK-based applications.
 *
 * @param data - VTK legacy format string (.vtk file contents)
 * @returns A promise that resolves to the loaded Mesh instance
 * @throws Error if the mesh data is invalid or parsing fails
 *
 * @category Mesh
 *
 * @example Load mesh from VTK string
 * ```typescript
 * const vtkData = `# vtk DataFile Version 3.0
 * Simple mesh
 * ASCII
 * DATASET UNSTRUCTURED_GRID
 * POINTS 4 float
 * 0 0 0
 * 1 0 0
 * 1 1 0
 * 0 1 0
 * CELLS 1 5
 * 4 0 1 2 3
 * CELL_TYPES 1
 * 9`;
 *
 * const mesh = await fromVTK(vtkData);
 * console.log(mesh.numVertices); // 4
 * mesh.destroy();
 * ```
 *
 * @example Load from file (Node.js)
 * ```typescript
 * import { readFileSync } from 'fs';
 *
 * const data = readFileSync('mesh.vtk', 'utf8');
 * const mesh = await fromVTK(data);
 * ```
 *
 * @see {@link fromGmsh} for loading Gmsh format meshes
 * @see {@link toVTK} for exporting to VTK format
 */
export async function fromVTK(data: string): Promise<Mesh> {
  return Mesh.fromString(data);
}
