/**
 * VTK format mesh export.
 *
 * @module mesh
 */

import { Mesh } from './Mesh.js';

/**
 * Exports a mesh to VTK legacy format string.
 *
 * Generates a string representation of the mesh in VTK legacy format (.vtk),
 * which can be visualized in ParaView, VisIt, or other VTK-compatible tools.
 *
 * @param mesh - The mesh to export
 * @returns VTK legacy format mesh string
 * @throws Error if memory allocation fails
 * @throws Error if export fails
 *
 * @category Mesh
 *
 * @example Export mesh for visualization
 * ```typescript
 * const mesh = await makeCartesian2D(10, 10);
 * refineUniform(mesh);
 * const vtkStr = toVTK(mesh);
 * console.log(vtkStr.substring(0, 200)); // Preview header
 * mesh.destroy();
 * ```
 *
 * @example Save to file for ParaView (Node.js)
 * ```typescript
 * import { writeFileSync } from 'fs';
 *
 * const mesh = await makeDisk(100);
 * setCurvature(mesh, 2);
 * const vtkStr = toVTK(mesh);
 * writeFileSync('disk.vtk', vtkStr);
 * // Open disk.vtk in ParaView
 * mesh.destroy();
 * ```
 *
 * @example Export 3D mesh
 * ```typescript
 * const mesh = await makeSphere(500);
 * setCurvature(mesh, 3);
 * const vtkStr = toVTK(mesh);
 * // vtkStr contains UNSTRUCTURED_GRID with tetrahedral cells
 * mesh.destroy();
 * ```
 *
 * @example Browser: download as file
 * ```typescript
 * const mesh = await makeCartesian3D(4, 4, 4);
 * const vtkStr = toVTK(mesh);
 *
 * const blob = new Blob([vtkStr], { type: 'text/plain' });
 * const url = URL.createObjectURL(blob);
 * const a = document.createElement('a');
 * a.href = url;
 * a.download = 'mesh.vtk';
 * a.click();
 * URL.revokeObjectURL(url);
 * mesh.destroy();
 * ```
 *
 * @see {@link toMFEM} - Export to MFEM format
 * @see {@link toGmsh} - Export to Gmsh format
 * @see {@link fromVTK} - Import from VTK format
 */
export function toVTK(mesh: Mesh): string {
  const module = mesh.getModule();
  const ptr = mesh.getPointer();

  // Allocate a large buffer for the output string
  const maxLen = 1024 * 1024; // 1MB should be sufficient for most meshes
  const outPtr = module._malloc(maxLen);

  if (outPtr === 0) {
    throw new Error('Failed to allocate memory for mesh export');
  }

  try {
    const result = module._mfem_mesh_to_vtk_string(ptr, outPtr, maxLen);

    if (result < 0) {
      throw new Error('Failed to export mesh to VTK format');
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
