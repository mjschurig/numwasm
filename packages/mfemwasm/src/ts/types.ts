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

  /**
   * Get the space dimension (embedding dimension) of the mesh.
   * @param meshPtr Pointer to Mesh object
   * @returns Space dimension, or -1 on error
   */
  _mfem_mesh_get_space_dimension(meshPtr: number): number;

  /**
   * Get the number of edges in the mesh.
   * @param meshPtr Pointer to Mesh object
   * @returns Number of edges, or -1 on error
   */
  _mfem_mesh_get_num_edges(meshPtr: number): number;

  /**
   * Get the number of faces in the mesh.
   * @param meshPtr Pointer to Mesh object
   * @returns Number of faces, or -1 on error
   */
  _mfem_mesh_get_num_faces(meshPtr: number): number;

  /**
   * Get element attribute.
   * @param meshPtr Pointer to Mesh object
   * @param elemIdx Element index
   * @returns Element attribute, or -1 on error
   */
  _mfem_mesh_get_attribute(meshPtr: number, elemIdx: number): number;

  /**
   * Get boundary element attribute.
   * @param meshPtr Pointer to Mesh object
   * @param bdrIdx Boundary element index
   * @returns Boundary attribute, or -1 on error
   */
  _mfem_mesh_get_bdr_attribute(meshPtr: number, bdrIdx: number): number;

  /**
   * Get element volumes (or areas in 2D, lengths in 1D).
   * @param meshPtr Pointer to Mesh object
   * @param volumesPtr Output array for volumes
   * @returns 0 on success, -1 on failure
   */
  _mfem_mesh_get_element_volumes(meshPtr: number, volumesPtr: number): number;

  /**
   * Get element centroids.
   * @param meshPtr Pointer to Mesh object
   * @param centroidsPtr Output array for centroids (size >= num_elements * space_dim)
   * @returns 0 on success, -1 on failure
   */
  _mfem_mesh_get_element_centroids(meshPtr: number, centroidsPtr: number): number;

  /**
   * Create a 1D Cartesian mesh.
   * @param n Number of elements
   * @param sx Domain size
   * @returns Pointer to new Mesh, or 0 on failure
   */
  _mfem_mesh_make_cartesian_1d(n: number, sx: number): number;

  /**
   * Create a 2D Cartesian mesh.
   * @param nx Number of elements in x
   * @param ny Number of elements in y
   * @param type Element type: 0 = quad, 1 = triangle
   * @param sx Domain size in x
   * @param sy Domain size in y
   * @param sfcOrdering Use space-filling curve ordering
   * @returns Pointer to new Mesh, or 0 on failure
   */
  _mfem_mesh_make_cartesian_2d(
    nx: number,
    ny: number,
    type: number,
    sx: number,
    sy: number,
    sfcOrdering: number
  ): number;

  /**
   * Create a 3D Cartesian mesh.
   * @param nx Number of elements in x
   * @param ny Number of elements in y
   * @param nz Number of elements in z
   * @param type Element type: 0 = hex, 1 = tet, 2 = wedge
   * @param sx Domain size in x
   * @param sy Domain size in y
   * @param sz Domain size in z
   * @param sfcOrdering Use space-filling curve ordering
   * @returns Pointer to new Mesh, or 0 on failure
   */
  _mfem_mesh_make_cartesian_3d(
    nx: number,
    ny: number,
    nz: number,
    type: number,
    sx: number,
    sy: number,
    sz: number,
    sfcOrdering: number
  ): number;

  /**
   * Local refinement of specified elements.
   * @param meshPtr Pointer to Mesh object
   * @param elementsPtr Array of element indices
   * @param numElements Number of elements
   * @returns 0 on success, -1 on failure
   */
  _mfem_mesh_refine_local(meshPtr: number, elementsPtr: number, numElements: number): number;

  /**
   * Scale mesh coordinates.
   * @param meshPtr Pointer to Mesh object
   * @param sx Scale in x
   * @param sy Scale in y
   * @param sz Scale in z
   * @returns 0 on success, -1 on failure
   */
  _mfem_mesh_scale(meshPtr: number, sx: number, sy: number, sz: number): number;

  /**
   * Translate mesh coordinates.
   * @param meshPtr Pointer to Mesh object
   * @param dx Translation in x
   * @param dy Translation in y
   * @param dz Translation in z
   * @returns 0 on success, -1 on failure
   */
  _mfem_mesh_translate(meshPtr: number, dx: number, dy: number, dz: number): number;

  /**
   * Rotate mesh coordinates (2D rotation).
   * @param meshPtr Pointer to Mesh object
   * @param angle Rotation angle in radians
   * @returns 0 on success, -1 on failure
   */
  _mfem_mesh_rotate_2d(meshPtr: number, angle: number): number;

  /**
   * Set mesh curvature order.
   * @param meshPtr Pointer to Mesh object
   * @param order Polynomial order
   * @returns 0 on success, -1 on failure
   */
  _mfem_mesh_set_curvature(meshPtr: number, order: number): number;

  /**
   * Export mesh to MFEM format string.
   * @param meshPtr Pointer to Mesh object
   * @param outPtr Output buffer
   * @param maxLen Maximum length
   * @returns Actual length, or -1 on failure
   */
  _mfem_mesh_to_mfem_string(meshPtr: number, outPtr: number, maxLen: number): number;

  /**
   * Get length of mesh in MFEM format.
   * @param meshPtr Pointer to Mesh object
   * @returns String length, or -1 on failure
   */
  _mfem_mesh_get_mfem_string_length(meshPtr: number): number;

  /**
   * Export mesh to VTK format string.
   * @param meshPtr Pointer to Mesh object
   * @param outPtr Output buffer
   * @param maxLen Maximum length
   * @returns Actual length, or -1 on failure
   */
  _mfem_mesh_to_vtk_string(meshPtr: number, outPtr: number, maxLen: number): number;

  /**
   * Get length of mesh in VTK format.
   * @param meshPtr Pointer to Mesh object
   * @returns String length, or -1 on failure
   */
  _mfem_mesh_get_vtk_string_length(meshPtr: number): number;

  /**
   * Get vertex coordinates.
   * @param meshPtr Pointer to Mesh object
   * @param vertexIdx Vertex index
   * @param coordsPtr Output array for coordinates
   * @returns 0 on success, -1 on failure
   */
  _mfem_mesh_get_vertex(meshPtr: number, vertexIdx: number, coordsPtr: number): number;

  /**
   * Set vertex coordinates.
   * @param meshPtr Pointer to Mesh object
   * @param vertexIdx Vertex index
   * @param coordsPtr New coordinates
   * @returns 0 on success, -1 on failure
   */
  _mfem_mesh_set_vertex(meshPtr: number, vertexIdx: number, coordsPtr: number): number;

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

  /**
   * Get the polynomial order of the space.
   * @param fesPtr Pointer to FiniteElementSpace object
   * @returns Polynomial order, or -1 on error
   */
  _mfem_fespace_get_order(fesPtr: number): number;

  /**
   * Get the true (global) vector size.
   * @param fesPtr Pointer to FiniteElementSpace object
   * @returns True vector size, or -1 on error
   */
  _mfem_fespace_get_true_vsize(fesPtr: number): number;

  /**
   * Create an H1 space with positive (Bernstein) basis.
   * @param meshPtr Pointer to Mesh object
   * @param order Polynomial order
   * @param vdim Vector dimension
   * @returns Pointer to new FiniteElementSpace, or 0 on failure
   */
  _mfem_fespace_create_h1_positive(meshPtr: number, order: number, vdim: number): number;

  /**
   * Create an L2 (discontinuous) finite element space.
   * @param meshPtr Pointer to Mesh object
   * @param order Polynomial order
   * @param vdim Vector dimension
   * @returns Pointer to new FiniteElementSpace, or 0 on failure
   */
  _mfem_fespace_create_l2(meshPtr: number, order: number, vdim: number): number;

  /**
   * Create an H(curl) (Nedelec) finite element space.
   * @param meshPtr Pointer to Mesh object
   * @param order Polynomial order
   * @returns Pointer to new FiniteElementSpace, or 0 on failure
   */
  _mfem_fespace_create_nd(meshPtr: number, order: number): number;

  /**
   * Create an H(div) (Raviart-Thomas) finite element space.
   * @param meshPtr Pointer to Mesh object
   * @param order Polynomial order
   * @returns Pointer to new FiniteElementSpace, or 0 on failure
   */
  _mfem_fespace_create_rt(meshPtr: number, order: number): number;

  /**
   * Get the DOF indices for an element.
   * @param fesPtr Pointer to FiniteElementSpace object
   * @param elemIdx Element index
   * @param dofsPtr Output array for DOF indices
   * @returns Number of DOFs, or -1 on error
   */
  _mfem_fespace_get_element_dofs(fesPtr: number, elemIdx: number, dofsPtr: number): number;

  /**
   * Get boundary DOFs for a given boundary attribute.
   * @param fesPtr Pointer to FiniteElementSpace object
   * @param bdrAttr Boundary attribute (1-based)
   * @param dofsPtr Output array for DOF indices
   * @param maxDofs Maximum DOFs to return
   * @returns Number of DOFs found, or -1 on error
   */
  _mfem_fespace_get_boundary_dofs(fesPtr: number, bdrAttr: number, dofsPtr: number, maxDofs: number): number;

  /**
   * Get essential (Dirichlet) DOFs for given boundary attributes.
   * @param fesPtr Pointer to FiniteElementSpace object
   * @param bdrAttrsPtr Array of boundary attributes
   * @param numAttrs Number of attributes
   * @param dofsPtr Output array for DOF indices
   * @param maxDofs Maximum DOFs to return
   * @returns Number of DOFs found, or -1 on error
   */
  _mfem_fespace_get_essential_dofs(fesPtr: number, bdrAttrsPtr: number, numAttrs: number, dofsPtr: number, maxDofs: number): number;

  /**
   * Create a NURBS finite element space.
   * @param meshPtr Pointer to Mesh object (must have NURBS extension)
   * @param order Polynomial order (-1 for variable order from mesh)
   * @param vdim Vector dimension
   * @returns Pointer to new FiniteElementSpace, or 0 on failure
   */
  _mfem_fespace_create_nurbs(meshPtr: number, order: number, vdim: number): number;

  // ============================================================
  // NURBS AND NCMESH API
  // ============================================================

  /**
   * Check if mesh has NURBS extension.
   * @param meshPtr Pointer to Mesh object
   * @returns 1 if NURBS mesh, 0 otherwise
   */
  _mfem_mesh_is_nurbs(meshPtr: number): number;

  /**
   * Get the NURBS order from a NURBS mesh.
   * @param meshPtr Pointer to Mesh object
   * @returns NURBS order, or -1 if not a NURBS mesh
   */
  _mfem_mesh_get_nurbs_order(meshPtr: number): number;

  /**
   * Check if mesh has NCMesh (non-conforming mesh) extension.
   * @param meshPtr Pointer to Mesh object
   * @returns 1 if NCMesh, 0 otherwise
   */
  _mfem_mesh_is_ncmesh(meshPtr: number): number;

  /**
   * Enable non-conforming mesh mode.
   * Must be called before local refinement if derefinement is desired later.
   * @param meshPtr Pointer to Mesh object
   * @returns 0 on success, -1 on failure
   */
  _mfem_mesh_enable_ncmesh(meshPtr: number): number;

  /**
   * Derefine mesh elements (coarsening).
   * The mesh must have been refined with NCMesh enabled.
   * @param meshPtr Pointer to Mesh object
   * @param derefsPtr Array of derefinement table rows
   * @param numDerefs Number of derefinements
   * @returns 0 on success, -1 on failure
   */
  _mfem_mesh_derefine(meshPtr: number, derefsPtr: number, numDerefs: number): number;

  /**
   * Get derefinement table size for NCMesh.
   * @param meshPtr Pointer to Mesh object
   * @returns Number of possible derefinements, or -1 on error
   */
  _mfem_mesh_get_deref_table_size(meshPtr: number): number;

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

  /**
   * Get the value of a grid function at a point within an element.
   * @param gfPtr Pointer to GridFunction object
   * @param elemIdx Element index
   * @param xiPtr Pointer to local coordinates array
   * @param dim Dimension of coordinates
   * @returns Value at point, or NaN on error
   */
  _mfem_gridfunc_get_value(gfPtr: number, elemIdx: number, xiPtr: number, dim: number): number;

  /**
   * Get the gradient of a grid function at a point.
   * @param gfPtr Pointer to GridFunction object
   * @param elemIdx Element index
   * @param xiPtr Pointer to local coordinates
   * @param dim Dimension of coordinates
   * @param gradPtr Output pointer for gradient
   * @returns 0 on success, -1 on failure
   */
  _mfem_gridfunc_get_gradient(gfPtr: number, elemIdx: number, xiPtr: number, dim: number, gradPtr: number): number;

  /**
   * Get the curl of a vector grid function at a point.
   * @param gfPtr Pointer to GridFunction object
   * @param elemIdx Element index
   * @param xiPtr Pointer to local coordinates
   * @param dim Dimension of coordinates
   * @param curlPtr Output pointer for curl
   * @returns 0 on success, -1 on failure
   */
  _mfem_gridfunc_get_curl(gfPtr: number, elemIdx: number, xiPtr: number, dim: number, curlPtr: number): number;

  /**
   * Get the divergence of a vector grid function at a point.
   * @param gfPtr Pointer to GridFunction object
   * @param elemIdx Element index
   * @param xiPtr Pointer to local coordinates
   * @param dim Dimension of coordinates
   * @returns Divergence value, or NaN on error
   */
  _mfem_gridfunc_get_divergence(gfPtr: number, elemIdx: number, xiPtr: number, dim: number): number;

  /**
   * Compute the L2 norm of a grid function.
   * @param gfPtr Pointer to GridFunction object
   * @returns L2 norm, or -1 on error
   */
  _mfem_gridfunc_norm_l2(gfPtr: number): number;

  /**
   * Compute the L-infinity norm of a grid function.
   * @param gfPtr Pointer to GridFunction object
   * @returns L-infinity norm, or -1 on error
   */
  _mfem_gridfunc_norm_linf(gfPtr: number): number;

  /**
   * Compute the H1 seminorm of a grid function.
   * @param gfPtr Pointer to GridFunction object
   * @returns H1 seminorm, or -1 on error
   */
  _mfem_gridfunc_norm_h1_semi(gfPtr: number): number;

  /**
   * Compute L2 error against a constant exact solution.
   * @param gfPtr Pointer to GridFunction object
   * @param exactValue Exact constant value
   * @returns L2 error, or -1 on error
   */
  _mfem_gridfunc_compute_l2_error_const(gfPtr: number, exactValue: number): number;

  /**
   * Project a function coefficient onto boundary DOFs.
   * @param gfPtr Pointer to GridFunction object
   * @param bdrAttr Boundary attribute
   * @param value Constant value to project
   * @returns 0 on success, -1 on failure
   */
  _mfem_gridfunc_project_bdr_coefficient(gfPtr: number, bdrAttr: number, value: number): number;

  /**
   * Save grid function to MFEM format string.
   * @param gfPtr Pointer to GridFunction object
   * @param outPtr Output buffer
   * @param maxLen Maximum length
   * @returns Actual length, or -1 on failure
   */
  _mfem_gridfunc_to_mfem_string(gfPtr: number, outPtr: number, maxLen: number): number;

  /**
   * Get length of grid function in MFEM format.
   * @param gfPtr Pointer to GridFunction object
   * @returns String length, or -1 on failure
   */
  _mfem_gridfunc_get_mfem_string_length(gfPtr: number): number;

  /**
   * Save grid function to VTK format string.
   * @param gfPtr Pointer to GridFunction object
   * @param fieldNamePtr Pointer to field name string
   * @param ref Refinement level
   * @param outPtr Output buffer
   * @param maxLen Maximum length
   * @returns Actual length, or -1 on failure
   */
  _mfem_gridfunc_to_vtk_string(gfPtr: number, fieldNamePtr: number, ref: number, outPtr: number, maxLen: number): number;

  /**
   * Get length of grid function in VTK format.
   * @param gfPtr Pointer to GridFunction object
   * @param fieldNamePtr Pointer to field name string
   * @param ref Refinement level
   * @returns String length, or -1 on failure
   */
  _mfem_gridfunc_get_vtk_string_length(gfPtr: number, fieldNamePtr: number, ref: number): number;

  /**
   * Get the finite element space from a grid function.
   * @param gfPtr Pointer to GridFunction object
   * @returns Pointer to FiniteElementSpace, or 0 on error
   */
  _mfem_gridfunc_get_fespace(gfPtr: number): number;

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
