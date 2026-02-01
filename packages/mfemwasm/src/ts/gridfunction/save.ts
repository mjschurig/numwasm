/**
 * Save grid function to MFEM format string.
 *
 * @module gridfunction
 */

import { GridFunction } from './GridFunction.js';

/**
 * Exports a grid function to MFEM format string.
 *
 * The returned string can be saved to a file and later loaded back using
 * MFEM's GridFunction constructor that takes a mesh and input stream.
 *
 * The format includes the finite element collection type, ordering,
 * and all DOF values.
 *
 * @param gf - GridFunction to export
 * @returns The grid function data in MFEM format
 * @throws Error if the grid function has been destroyed
 * @throws Error if export fails
 *
 * @category GridFunction
 *
 * @example Save solution to string
 * ```typescript
 * const mfemData = save(solution);
 * // Can be saved to file in browser:
 * // const blob = new Blob([mfemData], { type: 'text/plain' });
 * // const url = URL.createObjectURL(blob);
 * ```
 *
 * @example Save for later analysis
 * ```typescript
 * const data = save(temperature);
 * localStorage.setItem('temperature_solution', data);
 * ```
 *
 * @see {@link toVTK} for VTK format output suitable for visualization
 */
export function save(gf: GridFunction): string {
  const module = gf.getModule();
  const ptr = gf.getPointer();

  // First, get the required buffer size
  const length = module._mfem_gridfunc_get_mfem_string_length(ptr);

  if (length < 0) {
    throw new Error('Failed to get grid function string length');
  }

  // Allocate buffer with extra space for null terminator
  const bufferSize = length + 1;
  const outPtr = module._malloc(bufferSize);

  if (outPtr === 0) {
    throw new Error('Failed to allocate memory for output');
  }

  try {
    const actualLen = module._mfem_gridfunc_to_mfem_string(ptr, outPtr, bufferSize);

    if (actualLen < 0) {
      throw new Error('Failed to export grid function to MFEM format');
    }

    // Read the string from WASM memory
    return module.UTF8ToString(outPtr);
  } finally {
    module._free(outPtr);
  }
}
