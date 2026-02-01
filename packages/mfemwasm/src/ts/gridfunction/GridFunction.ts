/**
 * GridFunction class for MFEM WASM.
 *
 * @module gridfunction
 */

import type { MFEMModule } from '../types.js';
import { loadMFEMModule } from '../loader.js';
import type { FiniteElementSpace } from '../fespace/FiniteElementSpace.js';

/**
 * High-level wrapper for MFEM GridFunction objects.
 *
 * A grid function represents a discrete function defined on a finite element space.
 * It stores the coefficients (degrees of freedom) that define the function and provides
 * methods for evaluation, projection, and data access.
 *
 * Grid functions are the primary objects for representing solutions to PDEs,
 * boundary conditions, coefficients, and other field quantities.
 *
 * **Important:** Resources must be explicitly freed by calling {@link GridFunction.destroy | destroy()}.
 *
 * @category GridFunction
 *
 * @example Creating and initializing a grid function
 * ```typescript
 * import { Mesh, FiniteElementSpace, GridFunction } from 'mfemwasm';
 *
 * const mesh = await Mesh.fromString(meshData);
 * const fespace = await FiniteElementSpace.createH1(mesh, 2);
 * const gf = await GridFunction.create(fespace);
 *
 * // Initialize to zero
 * gf.projectConstant(0.0);
 *
 * // Or set data directly
 * const data = new Float64Array(gf.size);
 * data.fill(1.0);
 * gf.setData(data);
 *
 * // Clean up
 * gf.destroy();
 * fespace.destroy();
 * mesh.destroy();
 * ```
 *
 * @example Accessing solution data
 * ```typescript
 * // After solving a PDE...
 * const solution = await GridFunction.create(fespace);
 * // ... solve and populate solution ...
 *
 * // Get all DOF values
 * const data = solution.getData();
 * console.log(`Min value: ${Math.min(...data)}`);
 * console.log(`Max value: ${Math.max(...data)}`);
 * ```
 *
 * @see {@link https://mfem.org/howto/gridfunc/ | MFEM GridFunction Documentation}
 */
export class GridFunction {
  private ptr: number;
  private module: MFEMModule;
  private destroyed = false;

  /**
   * Creates a GridFunction wrapper from an existing WASM pointer.
   *
   * **Note:** This constructor is primarily for internal use. Prefer using
   * the static factory method {@link GridFunction.create | create()}.
   *
   * @param module - The loaded MFEM WASM module
   * @param ptr - Pointer to an existing MFEM GridFunction object in WASM memory
   */
  constructor(module: MFEMModule, ptr: number) {
    this.module = module;
    this.ptr = ptr;
  }

  /**
   * Creates a new grid function on the specified finite element space.
   *
   * The grid function is initialized with unspecified values. Use
   * {@link GridFunction.projectConstant | projectConstant()} or
   * {@link GridFunction.setData | setData()} to initialize the values.
   *
   * @param fespace - The finite element space on which to define the function
   * @returns A promise that resolves to the new GridFunction
   * @throws Error if memory allocation fails
   *
   * @example
   * ```typescript
   * const gf = await GridFunction.create(fespace);
   * gf.projectConstant(0.0); // Initialize to zero
   * ```
   */
  static async create(fespace: FiniteElementSpace): Promise<GridFunction> {
    const module = await loadMFEMModule();
    const ptr = module._mfem_gridfunc_create(fespace.getPointer());
    if (ptr === 0) {
      throw new Error('Failed to create grid function');
    }
    return new GridFunction(module, ptr);
  }

  /**
   * The total number of degrees of freedom (DOFs) in this grid function.
   *
   * This equals the number of DOFs in the underlying finite element space.
   * For vector-valued spaces, this is the total number of scalar DOFs.
   *
   * @throws Error if the grid function has been destroyed
   */
  get size(): number {
    this.checkDestroyed();
    return this.module._mfem_gridfunc_get_size(this.ptr);
  }

  /**
   * Gets a copy of all DOF values as a Float64Array.
   *
   * The returned array is a copy; modifications will not affect the grid function.
   * Use {@link GridFunction.setData | setData()} to update the values.
   *
   * @returns A Float64Array containing all DOF values
   * @throws Error if the grid function has been destroyed
   * @throws Error if data access fails
   *
   * @example
   * ```typescript
   * const data = gf.getData();
   * console.log(`Number of DOFs: ${data.length}`);
   * console.log(`First value: ${data[0]}`);
   *
   * // Compute statistics
   * const sum = data.reduce((a, b) => a + b, 0);
   * const mean = sum / data.length;
   * ```
   */
  getData(): Float64Array {
    this.checkDestroyed();
    const dataPtr = this.module._mfem_gridfunc_get_data(this.ptr);
    const size = this.size;

    if (dataPtr === 0) {
      throw new Error('Failed to get grid function data');
    }

    const result = new Float64Array(size);
    for (let i = 0; i < size; i++) {
      result[i] = this.module.HEAPF64[(dataPtr >> 3) + i];
    }
    return result;
  }

  /**
   * Sets all DOF values from an array.
   *
   * The array length must exactly match the grid function size. Values are
   * copied into WASM memory; subsequent modifications to the input array
   * will not affect the grid function.
   *
   * @param data - Array of DOF values (Float64Array or number[])
   * @throws Error if the grid function has been destroyed
   * @throws Error if the array length doesn't match the grid function size
   * @throws Error if data access fails
   *
   * @example
   * ```typescript
   * // Set all values to 1.0
   * const data = new Float64Array(gf.size).fill(1.0);
   * gf.setData(data);
   *
   * // Set values from computation
   * const computed = new Float64Array(gf.size);
   * for (let i = 0; i < computed.length; i++) {
   *   computed[i] = Math.sin(i * 0.1);
   * }
   * gf.setData(computed);
   * ```
   */
  setData(data: Float64Array | number[]): void {
    this.checkDestroyed();
    const dataPtr = this.module._mfem_gridfunc_get_data(this.ptr);
    const size = this.size;

    if (dataPtr === 0) {
      throw new Error('Failed to get grid function data pointer');
    }

    if (data.length !== size) {
      throw new Error(
        `Data length ${data.length} does not match grid function size ${size}`
      );
    }

    for (let i = 0; i < size; i++) {
      this.module.HEAPF64[(dataPtr >> 3) + i] = data[i];
    }
  }

  /**
   * Projects a constant value onto the grid function.
   *
   * This sets all DOF values such that the grid function represents the
   * constant function f(x) = value everywhere in the domain.
   *
   * @param value - The constant value to project
   * @throws Error if the grid function has been destroyed
   * @throws Error if projection fails
   *
   * @example
   * ```typescript
   * // Initialize temperature field to 20 degrees
   * temperature.projectConstant(20.0);
   *
   * // Set initial concentration to zero
   * concentration.projectConstant(0.0);
   * ```
   */
  projectConstant(value: number): void {
    this.checkDestroyed();
    const result = this.module._mfem_gridfunc_project_coefficient(
      this.ptr,
      value
    );
    if (result !== 0) {
      throw new Error('Failed to project coefficient');
    }
  }

  /**
   * Gets the raw WASM pointer to the underlying MFEM GridFunction object.
   *
   * **Warning:** The pointer is only valid while this object exists and
   * has not been destroyed.
   *
   * @returns The WASM memory pointer to the MFEM GridFunction object
   * @throws Error if the grid function has been destroyed
   */
  getPointer(): number {
    this.checkDestroyed();
    return this.ptr;
  }

  /**
   * Gets the internal MFEM WASM module reference.
   *
   * @returns The MFEMModule instance
   * @internal
   */
  getModule(): MFEMModule {
    return this.module;
  }

  /**
   * Checks if this grid function has been destroyed.
   *
   * @returns `true` if {@link GridFunction.destroy | destroy()} has been called
   */
  isDestroyed(): boolean {
    return this.destroyed;
  }

  /**
   * Destroys this grid function and frees all associated WASM memory.
   *
   * After calling this method, the grid function object should not be used.
   * Any subsequent method calls will throw an error.
   *
   * **Important:** Always call this method when you are done with a grid function
   * to prevent memory leaks. The associated finite element space and mesh should
   * be destroyed separately.
   *
   * @example
   * ```typescript
   * const mesh = await Mesh.fromString(meshData);
   * const fespace = await FiniteElementSpace.createH1(mesh, 1);
   * const gf = await GridFunction.create(fespace);
   * try {
   *   // ... use gf
   * } finally {
   *   gf.destroy();
   *   fespace.destroy();
   *   mesh.destroy();
   * }
   * ```
   */
  destroy(): void {
    if (!this.destroyed && this.ptr !== 0) {
      this.module._mfem_gridfunc_destroy(this.ptr);
      this.ptr = 0;
      this.destroyed = true;
    }
  }

  private checkDestroyed(): void {
    if (this.destroyed) {
      throw new Error('GridFunction has been destroyed');
    }
  }
}
