/**
 * High-level Mesh wrapper for MFEM WASM
 *
 * Provides a JavaScript-friendly API for working with MFEM meshes.
 */

import type { MFEMModule } from './types.js';
import { loadMFEMModule, getMFEMModule } from './loader.js';

/**
 * Bounding box representation
 */
export interface BoundingBox {
  min: number[];
  max: number[];
}

/**
 * High-level wrapper for MFEM Mesh objects.
 *
 * Provides a JavaScript-friendly API for mesh operations.
 * Resources must be explicitly freed by calling destroy().
 *
 * @example
 * ```typescript
 * // Load a mesh from string
 * const mesh = await Mesh.fromString(meshData);
 * console.log(`Dimension: ${mesh.dimension}`);
 * console.log(`Elements: ${mesh.numElements}`);
 *
 * // Refine the mesh
 * mesh.refineUniform();
 *
 * // Clean up
 * mesh.destroy();
 * ```
 */
export class Mesh {
  private ptr: number;
  private module: MFEMModule;
  private destroyed = false;

  /**
   * Create a Mesh wrapper from an existing pointer.
   * Typically you should use the static factory methods instead.
   */
  constructor(module: MFEMModule, ptr: number) {
    this.module = module;
    this.ptr = ptr;
  }

  /**
   * Create an empty mesh.
   *
   * @returns Promise resolving to a new empty Mesh
   */
  static async create(): Promise<Mesh> {
    const module = await loadMFEMModule();
    const ptr = module._mfem_mesh_create();
    if (ptr === 0) {
      throw new Error('Failed to create mesh');
    }
    return new Mesh(module, ptr);
  }

  /**
   * Load a mesh from a string containing mesh data.
   *
   * Supports MFEM mesh format and other formats.
   *
   * @param meshData String containing mesh data
   * @returns Promise resolving to the loaded Mesh
   * @throws Error if loading fails
   *
   * @example
   * ```typescript
   * const meshData = `
   * MFEM mesh v1.0
   *
   * dimension
   * 2
   *
   * elements
   * 1
   * 1 3 0 1 2
   *
   * boundary
   * 3
   * 1 1 0 1
   * 1 1 1 2
   * 1 1 2 0
   *
   * vertices
   * 3
   * 2
   * 0 0
   * 1 0
   * 0 1
   * `;
   * const mesh = await Mesh.fromString(meshData);
   * ```
   */
  static async fromString(meshData: string): Promise<Mesh> {
    const module = await loadMFEMModule();
    const ptr = module._mfem_mesh_create();
    if (ptr === 0) {
      throw new Error('Failed to create mesh');
    }

    const mesh = new Mesh(module, ptr);

    try {
      mesh.loadFromString(meshData);
      return mesh;
    } catch (error) {
      mesh.destroy();
      throw error;
    }
  }

  /**
   * Load mesh data from a string into this mesh object.
   *
   * @param meshData String containing mesh data
   * @throws Error if loading fails
   */
  loadFromString(meshData: string): void {
    this.checkDestroyed();

    // Allocate memory for the string
    const encoder = new TextEncoder();
    const bytes = encoder.encode(meshData);
    const dataPtr = this.module._malloc(bytes.length + 1);

    if (dataPtr === 0) {
      throw new Error('Failed to allocate memory for mesh data');
    }

    try {
      // Copy string to WASM memory
      this.module.HEAPU8.set(bytes, dataPtr);
      this.module.HEAPU8[dataPtr + bytes.length] = 0; // null terminator

      // Load the mesh
      const result = this.module._mfem_mesh_load_from_string(
        this.ptr,
        dataPtr,
        bytes.length
      );

      if (result !== 0) {
        throw new Error('Failed to parse mesh data');
      }
    } finally {
      this.module._free(dataPtr);
    }
  }

  /**
   * Get the spatial dimension of the mesh (1, 2, or 3).
   */
  get dimension(): number {
    this.checkDestroyed();
    return this.module._mfem_mesh_get_dimension(this.ptr);
  }

  /**
   * Get the number of elements in the mesh.
   */
  get numElements(): number {
    this.checkDestroyed();
    return this.module._mfem_mesh_get_num_elements(this.ptr);
  }

  /**
   * Get the number of vertices in the mesh.
   */
  get numVertices(): number {
    this.checkDestroyed();
    return this.module._mfem_mesh_get_num_vertices(this.ptr);
  }

  /**
   * Get the number of boundary elements in the mesh.
   */
  get numBoundaryElements(): number {
    this.checkDestroyed();
    return this.module._mfem_mesh_get_num_boundary_elements(this.ptr);
  }

  /**
   * Uniformly refine the mesh once.
   * Each element is subdivided into smaller elements.
   */
  refineUniform(): void {
    this.checkDestroyed();
    const result = this.module._mfem_mesh_refine_uniform(this.ptr);
    if (result !== 0) {
      throw new Error('Failed to refine mesh');
    }
  }

  /**
   * Get the bounding box of the mesh.
   *
   * @returns Object with min and max coordinate arrays
   */
  getBoundingBox(): BoundingBox {
    this.checkDestroyed();

    const dim = this.dimension;
    const minPtr = this.module._malloc(dim * 8); // 8 bytes per double
    const maxPtr = this.module._malloc(dim * 8);

    if (minPtr === 0 || maxPtr === 0) {
      if (minPtr) this.module._free(minPtr);
      if (maxPtr) this.module._free(maxPtr);
      throw new Error('Failed to allocate memory for bounding box');
    }

    try {
      const result = this.module._mfem_mesh_get_bounding_box(
        this.ptr,
        minPtr,
        maxPtr
      );

      if (result !== 0) {
        throw new Error('Failed to get bounding box');
      }

      const min: number[] = [];
      const max: number[] = [];

      for (let i = 0; i < dim; i++) {
        min.push(this.module.HEAPF64[(minPtr >> 3) + i]);
        max.push(this.module.HEAPF64[(maxPtr >> 3) + i]);
      }

      return { min, max };
    } finally {
      this.module._free(minPtr);
      this.module._free(maxPtr);
    }
  }

  /**
   * Get the raw pointer to the underlying MFEM Mesh object.
   * Use with caution - the pointer is only valid while this object exists.
   */
  getPointer(): number {
    this.checkDestroyed();
    return this.ptr;
  }

  /**
   * Check if this mesh has been destroyed.
   */
  isDestroyed(): boolean {
    return this.destroyed;
  }

  /**
   * Destroy this mesh and free associated resources.
   * After calling destroy(), this object should not be used.
   */
  destroy(): void {
    if (!this.destroyed && this.ptr !== 0) {
      this.module._mfem_mesh_destroy(this.ptr);
      this.ptr = 0;
      this.destroyed = true;
    }
  }

  private checkDestroyed(): void {
    if (this.destroyed) {
      throw new Error('Mesh has been destroyed');
    }
  }
}

/**
 * High-level wrapper for MFEM FiniteElementSpace objects.
 */
export class FiniteElementSpace {
  private ptr: number;
  private module: MFEMModule;
  private destroyed = false;

  constructor(module: MFEMModule, ptr: number) {
    this.module = module;
    this.ptr = ptr;
  }

  /**
   * Create an H1 (continuous Lagrange) finite element space.
   *
   * @param mesh The mesh to build the space on
   * @param order Polynomial order (1 = linear, 2 = quadratic, etc.)
   * @param vdim Vector dimension (default 1 = scalar field)
   * @returns Promise resolving to the new space
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
   * Get the number of degrees of freedom.
   */
  get ndofs(): number {
    this.checkDestroyed();
    return this.module._mfem_fespace_get_ndofs(this.ptr);
  }

  /**
   * Get the vector dimension.
   */
  get vdim(): number {
    this.checkDestroyed();
    return this.module._mfem_fespace_get_vdim(this.ptr);
  }

  /**
   * Get the raw pointer to the underlying MFEM FiniteElementSpace object.
   */
  getPointer(): number {
    this.checkDestroyed();
    return this.ptr;
  }

  isDestroyed(): boolean {
    return this.destroyed;
  }

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

/**
 * High-level wrapper for MFEM GridFunction objects.
 */
export class GridFunction {
  private ptr: number;
  private module: MFEMModule;
  private destroyed = false;

  constructor(module: MFEMModule, ptr: number) {
    this.module = module;
    this.ptr = ptr;
  }

  /**
   * Create a grid function on a finite element space.
   *
   * @param fespace The finite element space
   * @returns Promise resolving to the new grid function
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
   * Get the size of the data array.
   */
  get size(): number {
    this.checkDestroyed();
    return this.module._mfem_gridfunc_get_size(this.ptr);
  }

  /**
   * Get a copy of the grid function data as a Float64Array.
   */
  getData(): Float64Array {
    this.checkDestroyed();
    const dataPtr = this.module._mfem_gridfunc_get_data(this.ptr);
    const size = this.size;

    if (dataPtr === 0) {
      throw new Error('Failed to get grid function data');
    }

    // Copy data from WASM memory
    const result = new Float64Array(size);
    for (let i = 0; i < size; i++) {
      result[i] = this.module.HEAPF64[(dataPtr >> 3) + i];
    }
    return result;
  }

  /**
   * Set the grid function data from a Float64Array.
   *
   * @param data Array of values (must match size)
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

    // Copy data to WASM memory
    for (let i = 0; i < size; i++) {
      this.module.HEAPF64[(dataPtr >> 3) + i] = data[i];
    }
  }

  /**
   * Project a constant value onto the grid function.
   *
   * @param value Constant value to project
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
   * Get the raw pointer to the underlying MFEM GridFunction object.
   */
  getPointer(): number {
    this.checkDestroyed();
    return this.ptr;
  }

  isDestroyed(): boolean {
    return this.destroyed;
  }

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
