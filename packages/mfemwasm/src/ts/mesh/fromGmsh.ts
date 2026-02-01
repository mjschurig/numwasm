/**
 * Gmsh mesh format loading.
 *
 * @module mesh
 */

import { Mesh } from './Mesh.js';

/**
 * Loads a mesh from a Gmsh format string.
 *
 * Parses mesh data in Gmsh MSH format (versions 2.x and 4.x are supported).
 * This is useful for importing meshes created with the Gmsh mesh generator
 * or exported from other tools in Gmsh format.
 *
 * @param data - Gmsh mesh format string (MSH file contents)
 * @returns A promise that resolves to the loaded Mesh instance
 * @throws Error if the mesh data is invalid or parsing fails
 *
 * @category Mesh
 *
 * @example Load mesh from Gmsh string
 * ```typescript
 * const gmshData = `$MeshFormat
 * 2.2 0 8
 * $EndMeshFormat
 * $Nodes
 * 4
 * 1 0 0 0
 * 2 1 0 0
 * 3 1 1 0
 * 4 0 1 0
 * $EndNodes
 * $Elements
 * 2
 * 1 2 2 0 0 1 2 3
 * 2 2 2 0 0 1 3 4
 * $EndElements`;
 *
 * const mesh = await fromGmsh(gmshData);
 * console.log(mesh.numVertices); // 4
 * console.log(mesh.numElements); // 2
 * mesh.destroy();
 * ```
 *
 * @example Load from file (Node.js)
 * ```typescript
 * import { readFileSync } from 'fs';
 *
 * const data = readFileSync('mesh.msh', 'utf8');
 * const mesh = await fromGmsh(data);
 * ```
 *
 * @see {@link fromVTK} for loading VTK format meshes
 * @see {@link toGmsh} for exporting to Gmsh format
 */
export async function fromGmsh(data: string): Promise<Mesh> {
  return Mesh.fromString(data);
}
