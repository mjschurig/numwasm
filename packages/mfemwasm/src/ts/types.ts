/**
 * MFEM WASM Module Type Definitions
 *
 * MFEM is a finite element library for solving partial differential equations.
 * This module provides WebAssembly bindings for core MFEM functionality including:
 * - Mesh loading and manipulation
 * - Finite element spaces
 * - Grid functions
 */

/**
 * MFEM WebAssembly Module Interface
 *
 * All C API functions are prefixed with underscore and operate on opaque pointers.
 * Memory must be managed explicitly using malloc/free.
 */
export interface MFEMModule {
  // ============================================================
  // MESH API
  // ============================================================

  /**
   * Create an empty mesh object.
   * @returns Pointer to new Mesh, or 0 (null) on failure
   */
  _mfem_mesh_create(): number;

  /**
   * Destroy a mesh object and free its memory.
   * @param meshPtr Pointer to Mesh object
   */
  _mfem_mesh_destroy(meshPtr: number): void;

  /**
   * Load a mesh from a string containing mesh data.
   * Supports MFEM mesh format, VTK, and other formats.
   *
   * @param meshPtr Pointer to Mesh object
   * @param dataPtr Pointer to mesh data string
   * @param len Length of data string
   * @returns 0 on success, -1 on failure
   */
  _mfem_mesh_load_from_string(
    meshPtr: number,
    dataPtr: number,
    len: number
  ): number;

  /**
   * Get the spatial dimension of the mesh.
   * @param meshPtr Pointer to Mesh object
   * @returns Dimension (1, 2, or 3), or -1 on error
   */
  _mfem_mesh_get_dimension(meshPtr: number): number;

  /**
   * Get the number of elements in the mesh.
   * @param meshPtr Pointer to Mesh object
   * @returns Number of elements, or -1 on error
   */
  _mfem_mesh_get_num_elements(meshPtr: number): number;

  /**
   * Get the number of vertices in the mesh.
   * @param meshPtr Pointer to Mesh object
   * @returns Number of vertices, or -1 on error
   */
  _mfem_mesh_get_num_vertices(meshPtr: number): number;

  /**
   * Get the number of boundary elements in the mesh.
   * @param meshPtr Pointer to Mesh object
   * @returns Number of boundary elements, or -1 on error
   */
  _mfem_mesh_get_num_boundary_elements(meshPtr: number): number;

  /**
   * Uniformly refine the mesh once.
   * @param meshPtr Pointer to Mesh object
   * @returns 0 on success, -1 on failure
   */
  _mfem_mesh_refine_uniform(meshPtr: number): number;

  /**
   * Get the bounding box of the mesh.
   * @param meshPtr Pointer to Mesh object
   * @param minCoordsPtr Pointer to output array for minimum coordinates
   * @param maxCoordsPtr Pointer to output array for maximum coordinates
   * @returns 0 on success, -1 on failure
   */
  _mfem_mesh_get_bounding_box(
    meshPtr: number,
    minCoordsPtr: number,
    maxCoordsPtr: number
  ): number;

  // ============================================================
  // FINITE ELEMENT SPACE API
  // ============================================================

  /**
   * Create an H1 (continuous Lagrange) finite element space.
   * @param meshPtr Pointer to Mesh object
   * @param order Polynomial order (1 = linear, 2 = quadratic, etc.)
   * @param vdim Vector dimension (1 = scalar, >1 = vector field)
   * @returns Pointer to new FiniteElementSpace, or 0 on failure
   */
  _mfem_fespace_create_h1(
    meshPtr: number,
    order: number,
    vdim: number
  ): number;

  /**
   * Destroy a finite element space and free its memory.
   * @param fesPtr Pointer to FiniteElementSpace object
   */
  _mfem_fespace_destroy(fesPtr: number): void;

  /**
   * Get the number of degrees of freedom in the space.
   * @param fesPtr Pointer to FiniteElementSpace object
   * @returns Number of DOFs, or -1 on error
   */
  _mfem_fespace_get_ndofs(fesPtr: number): number;

  /**
   * Get the vector dimension of the space.
   * @param fesPtr Pointer to FiniteElementSpace object
   * @returns Vector dimension, or -1 on error
   */
  _mfem_fespace_get_vdim(fesPtr: number): number;

  // ============================================================
  // GRID FUNCTION API
  // ============================================================

  /**
   * Create a grid function on a finite element space.
   * @param fesPtr Pointer to FiniteElementSpace object
   * @returns Pointer to new GridFunction, or 0 on failure
   */
  _mfem_gridfunc_create(fesPtr: number): number;

  /**
   * Destroy a grid function and free its memory.
   * @param gfPtr Pointer to GridFunction object
   */
  _mfem_gridfunc_destroy(gfPtr: number): void;

  /**
   * Get a pointer to the grid function's data array.
   * @param gfPtr Pointer to GridFunction object
   * @returns Pointer to double array, or 0 on error
   */
  _mfem_gridfunc_get_data(gfPtr: number): number;

  /**
   * Get the size of the grid function's data array.
   * @param gfPtr Pointer to GridFunction object
   * @returns Size, or -1 on error
   */
  _mfem_gridfunc_get_size(gfPtr: number): number;

  /**
   * Project a constant value onto the grid function.
   * @param gfPtr Pointer to GridFunction object
   * @param value Constant value to project
   * @returns 0 on success, -1 on failure
   */
  _mfem_gridfunc_project_coefficient(gfPtr: number, value: number): number;

  // ============================================================
  // MEMORY MANAGEMENT
  // ============================================================

  /**
   * Allocate memory in the WASM heap.
   * @param size Number of bytes to allocate
   * @returns Pointer to allocated memory, or 0 on failure
   */
  _malloc(size: number): number;

  /**
   * Free previously allocated memory.
   * @param ptr Pointer returned by _malloc
   */
  _free(ptr: number): void;

  // ============================================================
  // HEAP VIEWS - Direct Access to WASM Memory
  // ============================================================

  /** Float64 view of WASM heap. Index = byteOffset / 8 for doubles. */
  HEAPF64: Float64Array;

  /** Float32 view of WASM heap. Index = byteOffset / 4 for floats. */
  HEAPF32: Float32Array;

  /** Int32 view of WASM heap. Index = byteOffset / 4 for 32-bit integers. */
  HEAP32: Int32Array;

  /** Int8 view of WASM heap. Index = byteOffset directly. */
  HEAP8: Int8Array;

  /** Uint8 view of WASM heap. Index = byteOffset directly. */
  HEAPU8: Uint8Array;

  // ============================================================
  // EMSCRIPTEN RUNTIME UTILITIES
  // ============================================================

  /**
   * Read a value from WASM memory at the specified pointer.
   * @param ptr Memory address (byte offset)
   * @param type Value type: 'i8', 'i16', 'i32', 'i64', 'float', 'double'
   * @returns The value at that address
   */
  getValue(ptr: number, type: string): number;

  /**
   * Write a value to WASM memory at the specified pointer.
   * @param ptr Memory address (byte offset)
   * @param value Value to write
   * @param type Value type: 'i8', 'i16', 'i32', 'i64', 'float', 'double'
   */
  setValue(ptr: number, value: number, type: string): void;

  /**
   * Convert a C string pointer to a JavaScript string.
   * @param ptr Pointer to null-terminated C string
   * @returns JavaScript string
   */
  UTF8ToString(ptr: number): string;

  /**
   * Write a JavaScript string to WASM memory as UTF-8.
   * @param str JavaScript string to write
   * @param outPtr Pointer to output buffer
   * @param maxBytesToWrite Maximum bytes to write (including null terminator)
   */
  stringToUTF8(str: string, outPtr: number, maxBytesToWrite: number): void;

  /**
   * Get the number of bytes needed to encode a string as UTF-8.
   * @param str JavaScript string
   * @returns Number of bytes (not including null terminator)
   */
  lengthBytesUTF8(str: string): number;
}

/**
 * Options for configuring the MFEM WASM module loader.
 */
export interface MFEMModuleOptions {
  /**
   * Custom function to locate WASM files.
   * @param path Filename being requested (e.g., 'mfem.wasm')
   * @param scriptDirectory Directory of the JS loader script
   * @returns Full URL or path to the file
   */
  locateFile?: (path: string, scriptDirectory: string) => string;
}

/**
 * Factory function to create an MFEM WASM module instance.
 *
 * @param options Optional configuration
 * @returns Promise resolving to initialized MFEM module
 */
export type MFEMModuleFactory = (
  options?: MFEMModuleOptions
) => Promise<MFEMModule>;
