/**
 * FiniteElementSpace class for MFEM WASM.
 *
 * @module fespace
 */

import type { MFEMModule } from '../types.js';
import { loadMFEMModule } from '../loader.js';
import type { Mesh } from '../mesh/Mesh.js';

/**
 * High-level wrapper for MFEM FiniteElementSpace objects.
 *
 * A finite element space defines the function space used for discretizing PDEs.
 * It combines a mesh with a choice of finite element basis functions (e.g., H1, L2,
 * H(curl), H(div)) and polynomial order.
 *
 * **Important:** Resources must be explicitly freed by calling {@link FiniteElementSpace.destroy | destroy()}.
 *
 * @category FiniteElementSpace
 *
 * @example Creating an H1 space for scalar problems
 * ```typescript
 * import { Mesh } from 'mfemwasm';
 * import { FiniteElementSpace } from 'mfemwasm/fespace';
 *
 * const mesh = await Mesh.fromString(meshData);
 *
 * // Create a linear (order 1) H1 space
 * const fespace = await FiniteElementSpace.createH1(mesh, 1);
 * console.log(`Degrees of freedom: ${fespace.ndofs}`);
 *
 * // Clean up
 * fespace.destroy();
 * mesh.destroy();
 * ```
 *
 * @example Creating a vector-valued H1 space
 * ```typescript
 * // For elasticity: 2D mesh with 2-component displacement field
 * const mesh = await makeCartesian2D(10, 10);
 * const fespace = await FiniteElementSpace.createH1(mesh, 2, 2); // order=2, vdim=2
 *
 * console.log(`Vector dimension: ${fespace.vdim}`);
 * console.log(`Total DOFs: ${fespace.ndofs}`);
 * ```
 *
 * @see {@link https://mfem.org/fem/ | MFEM Finite Element Spaces}
 */
export class FiniteElementSpace {
  private ptr: number;
  private module: MFEMModule;
  private destroyed = false;

  /**
   * Creates a FiniteElementSpace wrapper from an existing WASM pointer.
   *
   * **Note:** This constructor is primarily for internal use. Prefer using
   * the static factory methods like {@link FiniteElementSpace.createH1 | createH1()}.
   *
   * @param module - The loaded MFEM WASM module
   * @param ptr - Pointer to an existing MFEM FiniteElementSpace object in WASM memory
   */
  constructor(module: MFEMModule, ptr: number) {
    this.module = module;
    this.ptr = ptr;
  }

  /**
   * Creates an H1 (continuous Lagrange) finite element space.
   *
   * H1 spaces use continuous piecewise polynomial functions. They are the most
   * common choice for scalar problems like heat transfer, and for vector problems
   * like elasticity.
   *
   * **Basis functions:**
   * - Order 1: Linear (nodal values at vertices)
   * - Order 2: Quadratic (nodal values at vertices and edge midpoints)
   * - Order p: Polynomial of degree p
   *
   * @param mesh - The mesh on which to define the space
   * @param order - Polynomial order (1 = linear, 2 = quadratic, etc.)
   * @param vdim - Vector dimension (default 1 for scalar fields, 2 or 3 for vector fields)
   * @returns A promise that resolves to the new FiniteElementSpace
   * @throws Error if space creation fails
   *
   * @example Scalar Poisson problem
   * ```typescript
   * const mesh = await makeCartesian2D(10, 10);
   * const fespace = await FiniteElementSpace.createH1(mesh, 2); // Quadratic H1
   * ```
   *
   * @example 2D Elasticity (vector-valued)
   * ```typescript
   * const mesh = await makeCartesian2D(10, 10);
   * const fespace = await FiniteElementSpace.createH1(mesh, 1, 2); // Linear, 2 components
   * ```
   */
  static async createH1(
    mesh: Mesh,
    order: number,
    vdim = 1
  ): Promise<FiniteElementSpace> {
    const module = await loadMFEMModule();
    const ptr = module._mfem_fespace_create_h1(mesh.getPointer(), order, vdim);
    if (ptr === 0) {
      throw new Error('Failed to create finite element space');
    }
    return new FiniteElementSpace(module, ptr);
  }

  /**
   * The total number of degrees of freedom (DOFs) in the space.
   *
   * For vector-valued spaces, this is the total number of scalar DOFs
   * (i.e., number of nodes × vector dimension).
   *
   * @throws Error if the space has been destroyed
   */
  get ndofs(): number {
    this.checkDestroyed();
    return this.module._mfem_fespace_get_ndofs(this.ptr);
  }

  /**
   * The vector dimension of the space.
   *
   * - 1 for scalar fields (temperature, pressure, etc.)
   * - 2 for 2D vector fields (velocity, displacement in 2D)
   * - 3 for 3D vector fields (velocity, displacement in 3D)
   *
   * @throws Error if the space has been destroyed
   */
  get vdim(): number {
    this.checkDestroyed();
    return this.module._mfem_fespace_get_vdim(this.ptr);
  }

  /**
   * Gets the raw WASM pointer to the underlying MFEM FiniteElementSpace object.
   *
   * **Warning:** The pointer is only valid while this object exists and
   * has not been destroyed.
   *
   * @returns The WASM memory pointer to the MFEM FiniteElementSpace object
   * @throws Error if the space has been destroyed
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
   * Checks if this finite element space has been destroyed.
   *
   * @returns `true` if {@link FiniteElementSpace.destroy | destroy()} has been called
   */
  isDestroyed(): boolean {
    return this.destroyed;
  }

  /**
   * Destroys this finite element space and frees all associated WASM memory.
   *
   * After calling this method, the space object should not be used. Any
   * subsequent method calls will throw an error.
   *
   * **Important:** Always call this method when you are done with a space
   * to prevent memory leaks. The associated mesh should be destroyed separately.
   *
   * @example
   * ```typescript
   * const mesh = await Mesh.fromString(meshData);
   * const fespace = await FiniteElementSpace.createH1(mesh, 1);
   * try {
   *   // ... use fespace
   * } finally {
   *   fespace.destroy();
   *   mesh.destroy();
   * }
   * ```
   */
  destroy(): void {
    if (!this.destroyed && this.ptr !== 0) {
      this.module._mfem_fespace_destroy(this.ptr);
      this.ptr = 0;
      this.destroyed = true;
    }
  }

  private checkDestroyed(): void {
    if (this.destroyed) {
      throw new Error('FiniteElementSpace has been destroyed');
    }
  }
}
