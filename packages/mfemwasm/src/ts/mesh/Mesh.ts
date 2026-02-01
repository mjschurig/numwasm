/**
 * Mesh class and related types for MFEM WASM.
 *
 * @module mesh
 */

import type { MFEMModule } from '../types.js';
import { loadMFEMModule } from '../loader.js';

/**
 * Represents an axis-aligned bounding box in N-dimensional space.
 *
 * @category Mesh
 */
export interface BoundingBox {
  /** Minimum coordinates for each dimension */
  min: number[];
  /** Maximum coordinates for each dimension */
  max: number[];
}

/**
 * High-level wrapper for MFEM Mesh objects.
 *
 * The Mesh class provides a JavaScript-friendly API for working with finite element meshes.
 * It wraps the underlying MFEM C++ Mesh class and handles memory management automatically.
 *
 * **Important:** Resources must be explicitly freed by calling {@link Mesh.destroy | destroy()}.
 * Failure to do so will result in memory leaks.
 *
 * @category Mesh
 *
 * @example Creating a mesh from MFEM format string
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
 *
 * const mesh = await Mesh.fromString(meshData);
 * console.log(`Dimension: ${mesh.dimension}`);
 * console.log(`Elements: ${mesh.numElements}`);
 * console.log(`Vertices: ${mesh.numVertices}`);
 *
 * // Always clean up when done
 * mesh.destroy();
 * ```
 *
 * @example Using factory functions to create meshes
 * ```typescript
 * import { makeCartesian2D, makeSphere } from 'mfemwasm/mesh';
 *
 * // Create a 10x10 quadrilateral mesh
 * const quadMesh = await makeCartesian2D(10, 10, 'quad');
 *
 * // Create a 3D ball mesh
 * const sphereMesh = await makeSphere(100, 1.0);
 *
 * // Clean up
 * quadMesh.destroy();
 * sphereMesh.destroy();
 * ```
 *
 * @example Mesh refinement
 * ```typescript
 * const mesh = await makeCartesian2D(4, 4);
 * console.log(`Before: ${mesh.numElements} elements`);
 *
 * mesh.refineUniform();
 * console.log(`After: ${mesh.numElements} elements`);
 *
 * mesh.destroy();
 * ```
 */
export class Mesh {
  private ptr: number;
  private module: MFEMModule;
  private destroyed = false;

  /**
   * Creates a Mesh wrapper from an existing WASM pointer.
   *
   * **Note:** This constructor is primarily for internal use. Prefer using the static
   * factory methods {@link Mesh.create | create()} or {@link Mesh.fromString | fromString()}
   * instead.
   *
   * @param module - The loaded MFEM WASM module
   * @param ptr - Pointer to an existing MFEM Mesh object in WASM memory
   */
  constructor(module: MFEMModule, ptr: number) {
    this.module = module;
    this.ptr = ptr;
  }

  /**
   * Creates an empty mesh.
   *
   * The returned mesh has no elements, vertices, or boundaries. Use
   * {@link Mesh.loadFromString | loadFromString()} to populate it with mesh data.
   *
   * @returns A promise that resolves to a new empty Mesh instance
   * @throws Error if memory allocation fails
   *
   * @example
   * ```typescript
   * const mesh = await Mesh.create();
   * mesh.loadFromString(meshData);
   * // ... use mesh
   * mesh.destroy();
   * ```
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
   * Loads a mesh from a string containing mesh data in MFEM format.
   *
   * Supports the MFEM mesh format (versions 1.0 and 1.1). The mesh format
   * includes sections for dimension, elements, boundary elements, and vertices.
   *
   * @param meshData - String containing the mesh data in MFEM format
   * @returns A promise that resolves to the loaded Mesh instance
   * @throws Error if the mesh data is invalid or cannot be parsed
   *
   * @example
   * ```typescript
   * const meshData = await fetch('/meshes/square.mesh').then(r => r.text());
   * const mesh = await Mesh.fromString(meshData);
   * ```
   *
   * @see {@link https://mfem.org/mesh-formats/ | MFEM Mesh Formats}
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
   * Loads mesh data from a string into this mesh object.
   *
   * This method replaces any existing mesh data. The string must be in
   * MFEM mesh format.
   *
   * @param meshData - String containing the mesh data in MFEM format
   * @throws Error if the mesh has been destroyed
   * @throws Error if memory allocation fails
   * @throws Error if the mesh data cannot be parsed
   */
  loadFromString(meshData: string): void {
    this.checkDestroyed();

    const encoder = new TextEncoder();
    const bytes = encoder.encode(meshData);
    const dataPtr = this.module._malloc(bytes.length + 1);

    if (dataPtr === 0) {
      throw new Error('Failed to allocate memory for mesh data');
    }

    try {
      this.module.HEAPU8.set(bytes, dataPtr);
      this.module.HEAPU8[dataPtr + bytes.length] = 0;

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
   * The topological dimension of the mesh.
   *
   * - 1 for curve meshes (line segments)
   * - 2 for surface meshes (triangles, quadrilaterals)
   * - 3 for volume meshes (tetrahedra, hexahedra)
   *
   * @throws Error if the mesh has been destroyed
   */
  get dimension(): number {
    this.checkDestroyed();
    return this.module._mfem_mesh_get_dimension(this.ptr);
  }

  /**
   * The space (embedding) dimension of the mesh.
   *
   * This is the dimension of the coordinate space in which the mesh vertices
   * are defined. For example, a 2D surface mesh embedded in 3D space would have
   * dimension=2 but spaceDimension=3.
   *
   * @throws Error if the mesh has been destroyed
   */
  get spaceDimension(): number {
    this.checkDestroyed();
    return this.module._mfem_mesh_get_space_dimension(this.ptr);
  }

  /**
   * The total number of elements in the mesh.
   *
   * Elements are the primary geometric entities of the mesh (e.g., triangles
   * in 2D, tetrahedra in 3D).
   *
   * @throws Error if the mesh has been destroyed
   */
  get numElements(): number {
    this.checkDestroyed();
    return this.module._mfem_mesh_get_num_elements(this.ptr);
  }

  /**
   * The total number of vertices (nodes) in the mesh.
   *
   * @throws Error if the mesh has been destroyed
   */
  get numVertices(): number {
    this.checkDestroyed();
    return this.module._mfem_mesh_get_num_vertices(this.ptr);
  }

  /**
   * The total number of edges in the mesh.
   *
   * @throws Error if the mesh has been destroyed
   */
  get numEdges(): number {
    this.checkDestroyed();
    return this.module._mfem_mesh_get_num_edges(this.ptr);
  }

  /**
   * The total number of faces in the mesh.
   *
   * For 2D meshes, faces are the same as edges. For 3D meshes, faces are
   * the 2D surfaces bounding each element.
   *
   * @throws Error if the mesh has been destroyed
   */
  get numFaces(): number {
    this.checkDestroyed();
    return this.module._mfem_mesh_get_num_faces(this.ptr);
  }

  /**
   * The total number of boundary elements in the mesh.
   *
   * Boundary elements are the lower-dimensional entities that form the
   * boundary of the domain (e.g., edges in 2D, faces in 3D).
   *
   * @throws Error if the mesh has been destroyed
   */
  get numBoundaryElements(): number {
    this.checkDestroyed();
    return this.module._mfem_mesh_get_num_boundary_elements(this.ptr);
  }

  /**
   * Gets the material attribute of an element.
   *
   * Attributes are integer labels used to identify different regions or
   * materials in the mesh. They are used when setting up boundary conditions
   * and material properties.
   *
   * @param elemIndex - Zero-based index of the element
   * @returns The attribute value for the specified element
   * @throws Error if the mesh has been destroyed
   * @throws Error if the element index is out of bounds
   */
  getAttribute(elemIndex: number): number {
    this.checkDestroyed();
    return this.module._mfem_mesh_get_attribute(this.ptr, elemIndex);
  }

  /**
   * Gets the boundary attribute of a boundary element.
   *
   * Boundary attributes are integer labels used to identify different parts
   * of the domain boundary for applying boundary conditions.
   *
   * @param bdrIndex - Zero-based index of the boundary element
   * @returns The attribute value for the specified boundary element
   * @throws Error if the mesh has been destroyed
   * @throws Error if the boundary element index is out of bounds
   */
  getBdrAttribute(bdrIndex: number): number {
    this.checkDestroyed();
    return this.module._mfem_mesh_get_bdr_attribute(this.ptr, bdrIndex);
  }

  /**
   * Uniformly refines the mesh by subdividing each element.
   *
   * Each element is split into smaller elements of the same type:
   * - Triangles are split into 4 triangles
   * - Quadrilaterals are split into 4 quadrilaterals
   * - Tetrahedra are split into 8 tetrahedra
   * - Hexahedra are split into 8 hexahedra
   *
   * @throws Error if the mesh has been destroyed
   * @throws Error if refinement fails
   *
   * @example
   * ```typescript
   * const mesh = await makeCartesian2D(2, 2);
   * console.log(mesh.numElements); // 4
   *
   * mesh.refineUniform();
   * console.log(mesh.numElements); // 16
   *
   * mesh.refineUniform();
   * console.log(mesh.numElements); // 64
   * ```
   */
  refineUniform(): void {
    this.checkDestroyed();
    const result = this.module._mfem_mesh_refine_uniform(this.ptr);
    if (result !== 0) {
      throw new Error('Failed to refine mesh');
    }
  }

  /**
   * Computes the axis-aligned bounding box of the mesh.
   *
   * The bounding box is the smallest box (aligned with coordinate axes)
   * that contains all vertices of the mesh.
   *
   * @returns An object containing the minimum and maximum coordinates
   * @throws Error if the mesh has been destroyed
   * @throws Error if memory allocation fails
   *
   * @example
   * ```typescript
   * const mesh = await makeCartesian2D(10, 10, 'quad', 2.0, 3.0);
   * const bbox = mesh.getBoundingBox();
   * console.log(bbox.min); // [0, 0]
   * console.log(bbox.max); // [2, 3]
   * ```
   */
  getBoundingBox(): BoundingBox {
    this.checkDestroyed();

    const dim = this.dimension;
    const minPtr = this.module._malloc(dim * 8);
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
   * Whether this mesh is a NURBS mesh.
   *
   * NURBS (Non-Uniform Rational B-Splines) meshes are used for isogeometric
   * analysis and provide exact geometry representation from CAD data.
   *
   * @throws Error if the mesh has been destroyed
   *
   * @example Check NURBS support before creating NURBS space
   * ```typescript
   * if (mesh.isNURBS) {
   *   const fespace = createNURBS(mesh, mesh.nurbsOrder);
   * }
   * ```
   */
  get isNURBS(): boolean {
    this.checkDestroyed();
    return this.module._mfem_mesh_is_nurbs(this.ptr) === 1;
  }

  /**
   * The NURBS order of this mesh.
   *
   * Returns the polynomial order of the NURBS patches if this is a NURBS mesh,
   * or -1 if this is not a NURBS mesh.
   *
   * @throws Error if the mesh has been destroyed
   */
  get nurbsOrder(): number {
    this.checkDestroyed();
    return this.module._mfem_mesh_get_nurbs_order(this.ptr);
  }

  /**
   * Whether this mesh has NCMesh (non-conforming mesh) support.
   *
   * NCMesh support enables local refinement with hanging nodes and
   * derefinement (coarsening) operations. Use the `enableNCMesh()` function to
   * enable NCMesh support before refinement.
   *
   * @throws Error if the mesh has been destroyed
   */
  get isNCMesh(): boolean {
    this.checkDestroyed();
    return this.module._mfem_mesh_is_ncmesh(this.ptr) === 1;
  }

  /**
   * Gets the raw WASM pointer to the underlying MFEM Mesh object.
   *
   * **Warning:** The pointer is only valid while this Mesh object exists and
   * has not been destroyed. Use with caution for advanced interop scenarios.
   *
   * @returns The WASM memory pointer to the MFEM Mesh object
   * @throws Error if the mesh has been destroyed
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
   * Checks if this mesh has been destroyed.
   *
   * @returns `true` if {@link Mesh.destroy | destroy()} has been called, `false` otherwise
   */
  isDestroyed(): boolean {
    return this.destroyed;
  }

  /**
   * Destroys this mesh and frees all associated WASM memory.
   *
   * After calling this method, the mesh object should not be used. Any
   * subsequent method calls will throw an error.
   *
   * **Important:** Always call this method when you are done with a mesh
   * to prevent memory leaks.
   *
   * @example
   * ```typescript
   * const mesh = await Mesh.fromString(meshData);
   * try {
   *   // ... use mesh
   * } finally {
   *   mesh.destroy();
   * }
   * ```
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
