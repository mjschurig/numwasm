/**
 * Export grid function to VTK format.
 *
 * @module gridfunction
 */

import { GridFunction } from './GridFunction.js';

/**
 * Options for VTK export.
 *
 * @category GridFunction
 */
export interface VTKExportOptions {
  /** Name of the field in the VTK file (default: "solution") */
  fieldName?: string;
  /** Refinement level for high-order visualization (default: 1) */
  refinement?: number;
}

/**
 * Exports a grid function to VTK format string for visualization.
 *
 * The returned string contains VTK data for the grid function values.
 * **Note:** This exports only the field data. For complete VTK output,
 * you also need to export the mesh using `mesh.toVTK()`.
 *
 * For visualization in ParaView or similar tools, combine the mesh
 * and field data into a complete VTK file.
 *
 * @param gf - GridFunction to export
 * @param options - Export options
 * @returns The grid function data in VTK format
 * @throws Error if the grid function has been destroyed
 * @throws Error if export fails
 *
 * @category GridFunction
 *
 * @example Export for ParaView
 * ```typescript
 * const meshVTK = mesh.toVTK();
 * const fieldVTK = toVTK(solution, { fieldName: 'temperature' });
 * // Combine for complete VTK file
 * ```
 *
 * @example Export with high-order refinement
 * ```typescript
 * const vtk = toVTK(solution, {
 *   fieldName: 'pressure',
 *   refinement: 3 // More points for curved elements
 * });
 * ```
 *
 * @see {@link save} for MFEM format output suitable for reloading
 */
export function toVTK(gf: GridFunction, options: VTKExportOptions = {}): string {
  const module = gf.getModule();
  const ptr = gf.getPointer();

  const fieldName = options.fieldName ?? 'solution';
  const ref = options.refinement ?? 1;

  // Allocate memory for field name
  const fieldNameBytes = new TextEncoder().encode(fieldName);
  const fieldNamePtr = module._malloc(fieldNameBytes.length + 1);

  if (fieldNamePtr === 0) {
    throw new Error('Failed to allocate memory for field name');
  }

  try {
    // Copy field name to WASM memory
    module.HEAPU8.set(fieldNameBytes, fieldNamePtr);
    module.HEAPU8[fieldNamePtr + fieldNameBytes.length] = 0;

    // Get required buffer size
    const length = module._mfem_gridfunc_get_vtk_string_length(ptr, fieldNamePtr, ref);

    if (length < 0) {
      throw new Error('Failed to get VTK string length');
    }

    // Allocate output buffer
    const bufferSize = length + 1;
    const outPtr = module._malloc(bufferSize);

    if (outPtr === 0) {
      throw new Error('Failed to allocate memory for VTK output');
    }

    try {
      const actualLen = module._mfem_gridfunc_to_vtk_string(
        ptr,
        fieldNamePtr,
        ref,
        outPtr,
        bufferSize
      );

      if (actualLen < 0) {
        throw new Error('Failed to export grid function to VTK format');
      }

      return module.UTF8ToString(outPtr);
    } finally {
      module._free(outPtr);
    }
  } finally {
    module._free(fieldNamePtr);
  }
}
