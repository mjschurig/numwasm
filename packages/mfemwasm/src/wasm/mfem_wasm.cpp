/**
 * MFEM WebAssembly C API Wrapper
 *
 * This file provides a C-compatible API for MFEM functionality,
 * designed to be exported via Emscripten for use in JavaScript/TypeScript.
 *
 * All functions use extern "C" linkage and operate on opaque pointers
 * to MFEM objects.
 */

#include "mfem.hpp"
#include <cstring>
#include <sstream>
#include <cstdlib>

using namespace mfem;

extern "C" {

// ============================================================================
// Mesh API
// ============================================================================

/**
 * Create an empty mesh object.
 * @returns Pointer to new Mesh, or nullptr on failure
 */
void* mfem_mesh_create() {
    try {
        return new Mesh();
    } catch (...) {
        return nullptr;
    }
}

/**
 * Destroy a mesh object and free its memory.
 * @param mesh_ptr Pointer to Mesh object
 */
void mfem_mesh_destroy(void* mesh_ptr) {
    if (mesh_ptr) {
        delete static_cast<Mesh*>(mesh_ptr);
    }
}

/**
 * Load a mesh from a string containing mesh data.
 * Supports MFEM mesh format, VTK, and other formats.
 *
 * @param mesh_ptr Pointer to Mesh object
 * @param data Mesh data as a null-terminated string
 * @param len Length of data string
 * @returns 0 on success, -1 on failure
 */
int mfem_mesh_load_from_string(void* mesh_ptr, const char* data, int len) {
    if (!mesh_ptr || !data) return -1;

    try {
        Mesh* mesh = static_cast<Mesh*>(mesh_ptr);
        std::string mesh_data(data, len);
        std::istringstream mesh_stream(mesh_data);
        *mesh = Mesh(mesh_stream);
        return 0;
    } catch (...) {
        return -1;
    }
}

/**
 * Get the spatial dimension of the mesh.
 * @param mesh_ptr Pointer to Mesh object
 * @returns Dimension (1, 2, or 3), or -1 on error
 */
int mfem_mesh_get_dimension(void* mesh_ptr) {
    if (!mesh_ptr) return -1;
    return static_cast<Mesh*>(mesh_ptr)->Dimension();
}

/**
 * Get the number of elements in the mesh.
 * @param mesh_ptr Pointer to Mesh object
 * @returns Number of elements, or -1 on error
 */
int mfem_mesh_get_num_elements(void* mesh_ptr) {
    if (!mesh_ptr) return -1;
    return static_cast<Mesh*>(mesh_ptr)->GetNE();
}

/**
 * Get the number of vertices in the mesh.
 * @param mesh_ptr Pointer to Mesh object
 * @returns Number of vertices, or -1 on error
 */
int mfem_mesh_get_num_vertices(void* mesh_ptr) {
    if (!mesh_ptr) return -1;
    return static_cast<Mesh*>(mesh_ptr)->GetNV();
}

/**
 * Get the number of boundary elements in the mesh.
 * @param mesh_ptr Pointer to Mesh object
 * @returns Number of boundary elements, or -1 on error
 */
int mfem_mesh_get_num_boundary_elements(void* mesh_ptr) {
    if (!mesh_ptr) return -1;
    return static_cast<Mesh*>(mesh_ptr)->GetNBE();
}

/**
 * Uniformly refine the mesh once.
 * @param mesh_ptr Pointer to Mesh object
 * @returns 0 on success, -1 on failure
 */
int mfem_mesh_refine_uniform(void* mesh_ptr) {
    if (!mesh_ptr) return -1;
    try {
        static_cast<Mesh*>(mesh_ptr)->UniformRefinement();
        return 0;
    } catch (...) {
        return -1;
    }
}

/**
 * Get the bounding box of the mesh.
 * @param mesh_ptr Pointer to Mesh object
 * @param min_coords Output array for minimum coordinates (size >= dimension)
 * @param max_coords Output array for maximum coordinates (size >= dimension)
 * @returns 0 on success, -1 on failure
 */
int mfem_mesh_get_bounding_box(void* mesh_ptr, double* min_coords, double* max_coords) {
    if (!mesh_ptr || !min_coords || !max_coords) return -1;

    try {
        Mesh* mesh = static_cast<Mesh*>(mesh_ptr);
        int dim = mesh->Dimension();
        Vector min_v(dim), max_v(dim);
        mesh->GetBoundingBox(min_v, max_v);

        for (int i = 0; i < dim; i++) {
            min_coords[i] = min_v(i);
            max_coords[i] = max_v(i);
        }
        return 0;
    } catch (...) {
        return -1;
    }
}

/**
 * Get the space dimension (embedding dimension) of the mesh.
 * @param mesh_ptr Pointer to Mesh object
 * @returns Space dimension, or -1 on error
 */
int mfem_mesh_get_space_dimension(void* mesh_ptr) {
    if (!mesh_ptr) return -1;
    return static_cast<Mesh*>(mesh_ptr)->SpaceDimension();
}

/**
 * Get the number of edges in the mesh.
 * @param mesh_ptr Pointer to Mesh object
 * @returns Number of edges, or -1 on error
 */
int mfem_mesh_get_num_edges(void* mesh_ptr) {
    if (!mesh_ptr) return -1;
    return static_cast<Mesh*>(mesh_ptr)->GetNEdges();
}

/**
 * Get the number of faces in the mesh.
 * @param mesh_ptr Pointer to Mesh object
 * @returns Number of faces, or -1 on error
 */
int mfem_mesh_get_num_faces(void* mesh_ptr) {
    if (!mesh_ptr) return -1;
    return static_cast<Mesh*>(mesh_ptr)->GetNFaces();
}

/**
 * Get element attribute.
 * @param mesh_ptr Pointer to Mesh object
 * @param elem_idx Element index
 * @returns Element attribute, or -1 on error
 */
int mfem_mesh_get_attribute(void* mesh_ptr, int elem_idx) {
    if (!mesh_ptr) return -1;
    Mesh* mesh = static_cast<Mesh*>(mesh_ptr);
    if (elem_idx < 0 || elem_idx >= mesh->GetNE()) return -1;
    return mesh->GetAttribute(elem_idx);
}

/**
 * Get boundary element attribute.
 * @param mesh_ptr Pointer to Mesh object
 * @param bdr_idx Boundary element index
 * @returns Boundary attribute, or -1 on error
 */
int mfem_mesh_get_bdr_attribute(void* mesh_ptr, int bdr_idx) {
    if (!mesh_ptr) return -1;
    Mesh* mesh = static_cast<Mesh*>(mesh_ptr);
    if (bdr_idx < 0 || bdr_idx >= mesh->GetNBE()) return -1;
    return mesh->GetBdrAttribute(bdr_idx);
}

/**
 * Get element volumes (or areas in 2D, lengths in 1D).
 * @param mesh_ptr Pointer to Mesh object
 * @param volumes Output array for volumes (size >= num_elements)
 * @returns 0 on success, -1 on failure
 */
int mfem_mesh_get_element_volumes(void* mesh_ptr, double* volumes) {
    if (!mesh_ptr || !volumes) return -1;
    try {
        Mesh* mesh = static_cast<Mesh*>(mesh_ptr);
        int ne = mesh->GetNE();
        for (int i = 0; i < ne; i++) {
            volumes[i] = mesh->GetElementVolume(i);
        }
        return 0;
    } catch (...) {
        return -1;
    }
}

/**
 * Get element centroids.
 * @param mesh_ptr Pointer to Mesh object
 * @param centroids Output array for centroids (size >= num_elements * space_dim)
 * @returns 0 on success, -1 on failure
 */
int mfem_mesh_get_element_centroids(void* mesh_ptr, double* centroids) {
    if (!mesh_ptr || !centroids) return -1;
    try {
        Mesh* mesh = static_cast<Mesh*>(mesh_ptr);
        int ne = mesh->GetNE();
        int sdim = mesh->SpaceDimension();
        Vector center(sdim);
        for (int i = 0; i < ne; i++) {
            mesh->GetElementCenter(i, center);
            for (int d = 0; d < sdim; d++) {
                centroids[i * sdim + d] = center(d);
            }
        }
        return 0;
    } catch (...) {
        return -1;
    }
}

/**
 * Create a 1D Cartesian mesh.
 * @param n Number of elements
 * @param sx Domain size (length)
 * @returns Pointer to new Mesh, or nullptr on failure
 */
void* mfem_mesh_make_cartesian_1d(int n, double sx) {
    if (n < 1) return nullptr;
    try {
        return new Mesh(Mesh::MakeCartesian1D(n, sx));
    } catch (...) {
        return nullptr;
    }
}

/**
 * Create a 2D Cartesian mesh.
 * @param nx Number of elements in x
 * @param ny Number of elements in y
 * @param type Element type: 0 = quadrilateral, 1 = triangle
 * @param sx Domain size in x
 * @param sy Domain size in y
 * @param sfc_ordering Use space-filling curve ordering
 * @returns Pointer to new Mesh, or nullptr on failure
 */
void* mfem_mesh_make_cartesian_2d(int nx, int ny, int type, double sx, double sy, int sfc_ordering) {
    if (nx < 1 || ny < 1) return nullptr;
    try {
        Element::Type elem_type = (type == 0) ? Element::QUADRILATERAL : Element::TRIANGLE;
        return new Mesh(Mesh::MakeCartesian2D(nx, ny, elem_type, false, sx, sy, sfc_ordering != 0));
    } catch (...) {
        return nullptr;
    }
}

/**
 * Create a 3D Cartesian mesh.
 * @param nx Number of elements in x
 * @param ny Number of elements in y
 * @param nz Number of elements in z
 * @param type Element type: 0 = hexahedron, 1 = tetrahedron, 2 = wedge
 * @param sx Domain size in x
 * @param sy Domain size in y
 * @param sz Domain size in z
 * @param sfc_ordering Use space-filling curve ordering
 * @returns Pointer to new Mesh, or nullptr on failure
 */
void* mfem_mesh_make_cartesian_3d(int nx, int ny, int nz, int type, double sx, double sy, double sz, int sfc_ordering) {
    if (nx < 1 || ny < 1 || nz < 1) return nullptr;
    try {
        Element::Type elem_type;
        switch (type) {
            case 1: elem_type = Element::TETRAHEDRON; break;
            case 2: elem_type = Element::WEDGE; break;
            default: elem_type = Element::HEXAHEDRON; break;
        }
        return new Mesh(Mesh::MakeCartesian3D(nx, ny, nz, elem_type, sx, sy, sz, sfc_ordering != 0));
    } catch (...) {
        return nullptr;
    }
}

/**
 * Local refinement of specified elements.
 * @param mesh_ptr Pointer to Mesh object
 * @param elements Array of element indices to refine
 * @param num_elements Number of elements to refine
 * @returns 0 on success, -1 on failure
 */
int mfem_mesh_refine_local(void* mesh_ptr, const int* elements, int num_elements) {
    if (!mesh_ptr || !elements || num_elements < 1) return -1;
    try {
        Mesh* mesh = static_cast<Mesh*>(mesh_ptr);
        Array<int> el_to_refine(num_elements);
        for (int i = 0; i < num_elements; i++) {
            el_to_refine[i] = elements[i];
        }
        mesh->GeneralRefinement(el_to_refine);
        return 0;
    } catch (...) {
        return -1;
    }
}

/**
 * Transform mesh coordinates using a vector function.
 * Note: This requires a callback mechanism; for now we provide specific transforms.
 */

/**
 * Scale mesh coordinates.
 * @param mesh_ptr Pointer to Mesh object
 * @param sx Scale factor in x
 * @param sy Scale factor in y (ignored for 1D)
 * @param sz Scale factor in z (ignored for 1D/2D)
 * @returns 0 on success, -1 on failure
 */
int mfem_mesh_scale(void* mesh_ptr, double sx, double sy, double sz) {
    if (!mesh_ptr) return -1;
    try {
        Mesh* mesh = static_cast<Mesh*>(mesh_ptr);
        int sdim = mesh->SpaceDimension();
        double* coords = mesh->GetVertex(0);
        int nv = mesh->GetNV();

        for (int i = 0; i < nv; i++) {
            double* v = mesh->GetVertex(i);
            v[0] *= sx;
            if (sdim > 1) v[1] *= sy;
            if (sdim > 2) v[2] *= sz;
        }
        return 0;
    } catch (...) {
        return -1;
    }
}

/**
 * Translate mesh coordinates.
 * @param mesh_ptr Pointer to Mesh object
 * @param dx Translation in x
 * @param dy Translation in y (ignored for 1D)
 * @param dz Translation in z (ignored for 1D/2D)
 * @returns 0 on success, -1 on failure
 */
int mfem_mesh_translate(void* mesh_ptr, double dx, double dy, double dz) {
    if (!mesh_ptr) return -1;
    try {
        Mesh* mesh = static_cast<Mesh*>(mesh_ptr);
        int sdim = mesh->SpaceDimension();
        int nv = mesh->GetNV();

        for (int i = 0; i < nv; i++) {
            double* v = mesh->GetVertex(i);
            v[0] += dx;
            if (sdim > 1) v[1] += dy;
            if (sdim > 2) v[2] += dz;
        }
        return 0;
    } catch (...) {
        return -1;
    }
}

/**
 * Rotate mesh coordinates (2D rotation around origin, or 3D around z-axis).
 * @param mesh_ptr Pointer to Mesh object
 * @param angle Rotation angle in radians
 * @returns 0 on success, -1 on failure
 */
int mfem_mesh_rotate_2d(void* mesh_ptr, double angle) {
    if (!mesh_ptr) return -1;
    try {
        Mesh* mesh = static_cast<Mesh*>(mesh_ptr);
        int sdim = mesh->SpaceDimension();
        if (sdim < 2) return -1;

        int nv = mesh->GetNV();
        double c = cos(angle);
        double s = sin(angle);

        for (int i = 0; i < nv; i++) {
            double* v = mesh->GetVertex(i);
            double x = v[0], y = v[1];
            v[0] = c * x - s * y;
            v[1] = s * x + c * y;
        }
        return 0;
    } catch (...) {
        return -1;
    }
}

/**
 * Set mesh curvature order (convert to high-order mesh).
 * @param mesh_ptr Pointer to Mesh object
 * @param order Polynomial order for curved elements
 * @returns 0 on success, -1 on failure
 */
int mfem_mesh_set_curvature(void* mesh_ptr, int order) {
    if (!mesh_ptr || order < 1) return -1;
    try {
        Mesh* mesh = static_cast<Mesh*>(mesh_ptr);
        mesh->SetCurvature(order);
        return 0;
    } catch (...) {
        return -1;
    }
}

/**
 * Export mesh to MFEM format string.
 * @param mesh_ptr Pointer to Mesh object
 * @param out_ptr Output buffer pointer
 * @param max_len Maximum output length
 * @returns Actual length written, or -1 on failure
 */
int mfem_mesh_to_mfem_string(void* mesh_ptr, char* out_ptr, int max_len) {
    if (!mesh_ptr || !out_ptr || max_len < 1) return -1;
    try {
        Mesh* mesh = static_cast<Mesh*>(mesh_ptr);
        std::ostringstream oss;
        mesh->Print(oss);
        std::string result = oss.str();
        int len = std::min((int)result.size(), max_len - 1);
        memcpy(out_ptr, result.c_str(), len);
        out_ptr[len] = '\0';
        return len;
    } catch (...) {
        return -1;
    }
}

/**
 * Get the length of mesh when exported to MFEM format.
 * @param mesh_ptr Pointer to Mesh object
 * @returns String length, or -1 on failure
 */
int mfem_mesh_get_mfem_string_length(void* mesh_ptr) {
    if (!mesh_ptr) return -1;
    try {
        Mesh* mesh = static_cast<Mesh*>(mesh_ptr);
        std::ostringstream oss;
        mesh->Print(oss);
        return oss.str().size();
    } catch (...) {
        return -1;
    }
}

/**
 * Export mesh to VTK format string.
 * @param mesh_ptr Pointer to Mesh object
 * @param out_ptr Output buffer pointer
 * @param max_len Maximum output length
 * @returns Actual length written, or -1 on failure
 */
int mfem_mesh_to_vtk_string(void* mesh_ptr, char* out_ptr, int max_len) {
    if (!mesh_ptr || !out_ptr || max_len < 1) return -1;
    try {
        Mesh* mesh = static_cast<Mesh*>(mesh_ptr);
        std::ostringstream oss;
        mesh->PrintVTK(oss);
        std::string result = oss.str();
        int len = std::min((int)result.size(), max_len - 1);
        memcpy(out_ptr, result.c_str(), len);
        out_ptr[len] = '\0';
        return len;
    } catch (...) {
        return -1;
    }
}

/**
 * Get the length of mesh when exported to VTK format.
 * @param mesh_ptr Pointer to Mesh object
 * @returns String length, or -1 on failure
 */
int mfem_mesh_get_vtk_string_length(void* mesh_ptr) {
    if (!mesh_ptr) return -1;
    try {
        Mesh* mesh = static_cast<Mesh*>(mesh_ptr);
        std::ostringstream oss;
        mesh->PrintVTK(oss);
        return oss.str().size();
    } catch (...) {
        return -1;
    }
}

/**
 * Get vertex coordinates.
 * @param mesh_ptr Pointer to Mesh object
 * @param vertex_idx Vertex index
 * @param coords Output array for coordinates (size >= space_dim)
 * @returns 0 on success, -1 on failure
 */
int mfem_mesh_get_vertex(void* mesh_ptr, int vertex_idx, double* coords) {
    if (!mesh_ptr || !coords) return -1;
    try {
        Mesh* mesh = static_cast<Mesh*>(mesh_ptr);
        if (vertex_idx < 0 || vertex_idx >= mesh->GetNV()) return -1;
        double* v = mesh->GetVertex(vertex_idx);
        int sdim = mesh->SpaceDimension();
        for (int d = 0; d < sdim; d++) {
            coords[d] = v[d];
        }
        return 0;
    } catch (...) {
        return -1;
    }
}

/**
 * Set vertex coordinates.
 * @param mesh_ptr Pointer to Mesh object
 * @param vertex_idx Vertex index
 * @param coords New coordinates (size >= space_dim)
 * @returns 0 on success, -1 on failure
 */
int mfem_mesh_set_vertex(void* mesh_ptr, int vertex_idx, const double* coords) {
    if (!mesh_ptr || !coords) return -1;
    try {
        Mesh* mesh = static_cast<Mesh*>(mesh_ptr);
        if (vertex_idx < 0 || vertex_idx >= mesh->GetNV()) return -1;
        double* v = mesh->GetVertex(vertex_idx);
        int sdim = mesh->SpaceDimension();
        for (int d = 0; d < sdim; d++) {
            v[d] = coords[d];
        }
        return 0;
    } catch (...) {
        return -1;
    }
}

// ============================================================================
// FiniteElementSpace API
// ============================================================================

/**
 * Create an H1 (continuous Lagrange) finite element space.
 * @param mesh_ptr Pointer to Mesh object
 * @param order Polynomial order (1 = linear, 2 = quadratic, etc.)
 * @param vdim Vector dimension (1 = scalar, >1 = vector field)
 * @returns Pointer to new FiniteElementSpace, or nullptr on failure
 */
void* mfem_fespace_create_h1(void* mesh_ptr, int order, int vdim) {
    if (!mesh_ptr || order < 1) return nullptr;

    try {
        Mesh* mesh = static_cast<Mesh*>(mesh_ptr);
        FiniteElementCollection* fec = new H1_FECollection(order, mesh->Dimension());
        FiniteElementSpace* fes = new FiniteElementSpace(mesh, fec, vdim);
        // Note: fes takes ownership of fec, but we store pointer to allow cleanup
        return fes;
    } catch (...) {
        return nullptr;
    }
}

/**
 * Destroy a finite element space and free its memory.
 * @param fes_ptr Pointer to FiniteElementSpace object
 */
void mfem_fespace_destroy(void* fes_ptr) {
    if (fes_ptr) {
        FiniteElementSpace* fes = static_cast<FiniteElementSpace*>(fes_ptr);
        const FiniteElementCollection* fec = fes->FEColl();
        delete fes;
        delete fec;
    }
}

/**
 * Get the number of degrees of freedom in the space.
 * @param fes_ptr Pointer to FiniteElementSpace object
 * @returns Number of DOFs, or -1 on error
 */
int mfem_fespace_get_ndofs(void* fes_ptr) {
    if (!fes_ptr) return -1;
    return static_cast<FiniteElementSpace*>(fes_ptr)->GetNDofs();
}

/**
 * Get the vector dimension of the space.
 * @param fes_ptr Pointer to FiniteElementSpace object
 * @returns Vector dimension, or -1 on error
 */
int mfem_fespace_get_vdim(void* fes_ptr) {
    if (!fes_ptr) return -1;
    return static_cast<FiniteElementSpace*>(fes_ptr)->GetVDim();
}

/**
 * Get the polynomial order of the space.
 * @param fes_ptr Pointer to FiniteElementSpace object
 * @returns Polynomial order, or -1 on error
 */
int mfem_fespace_get_order(void* fes_ptr) {
    if (!fes_ptr) return -1;
    return static_cast<FiniteElementSpace*>(fes_ptr)->GetMaxElementOrder();
}

/**
 * Get the true (global) number of DOFs.
 * @param fes_ptr Pointer to FiniteElementSpace object
 * @returns True DOF count, or -1 on error
 */
int mfem_fespace_get_true_vsize(void* fes_ptr) {
    if (!fes_ptr) return -1;
    return static_cast<FiniteElementSpace*>(fes_ptr)->GetTrueVSize();
}

/**
 * Create an H1 space with positive (Bernstein) basis.
 * @param mesh_ptr Pointer to Mesh object
 * @param order Polynomial order
 * @param vdim Vector dimension
 * @returns Pointer to new FiniteElementSpace, or nullptr on failure
 */
void* mfem_fespace_create_h1_positive(void* mesh_ptr, int order, int vdim) {
    if (!mesh_ptr || order < 1) return nullptr;

    try {
        Mesh* mesh = static_cast<Mesh*>(mesh_ptr);
        FiniteElementCollection* fec = new H1Pos_FECollection(order, mesh->Dimension());
        FiniteElementSpace* fes = new FiniteElementSpace(mesh, fec, vdim);
        return fes;
    } catch (...) {
        return nullptr;
    }
}

/**
 * Create an L2 (discontinuous) finite element space.
 * @param mesh_ptr Pointer to Mesh object
 * @param order Polynomial order
 * @param vdim Vector dimension
 * @returns Pointer to new FiniteElementSpace, or nullptr on failure
 */
void* mfem_fespace_create_l2(void* mesh_ptr, int order, int vdim) {
    if (!mesh_ptr || order < 0) return nullptr;

    try {
        Mesh* mesh = static_cast<Mesh*>(mesh_ptr);
        FiniteElementCollection* fec = new L2_FECollection(order, mesh->Dimension());
        FiniteElementSpace* fes = new FiniteElementSpace(mesh, fec, vdim);
        return fes;
    } catch (...) {
        return nullptr;
    }
}

/**
 * Create an H(curl) (Nedelec) finite element space.
 * @param mesh_ptr Pointer to Mesh object
 * @param order Polynomial order
 * @returns Pointer to new FiniteElementSpace, or nullptr on failure
 */
void* mfem_fespace_create_nd(void* mesh_ptr, int order) {
    if (!mesh_ptr || order < 1) return nullptr;

    try {
        Mesh* mesh = static_cast<Mesh*>(mesh_ptr);
        FiniteElementCollection* fec = new ND_FECollection(order, mesh->Dimension());
        FiniteElementSpace* fes = new FiniteElementSpace(mesh, fec);
        return fes;
    } catch (...) {
        return nullptr;
    }
}

/**
 * Create an H(div) (Raviart-Thomas) finite element space.
 * @param mesh_ptr Pointer to Mesh object
 * @param order Polynomial order
 * @returns Pointer to new FiniteElementSpace, or nullptr on failure
 */
void* mfem_fespace_create_rt(void* mesh_ptr, int order) {
    if (!mesh_ptr || order < 0) return nullptr;

    try {
        Mesh* mesh = static_cast<Mesh*>(mesh_ptr);
        FiniteElementCollection* fec = new RT_FECollection(order, mesh->Dimension());
        FiniteElementSpace* fes = new FiniteElementSpace(mesh, fec);
        return fes;
    } catch (...) {
        return nullptr;
    }
}

/**
 * Get the DOF map for an element.
 * @param fes_ptr Pointer to FiniteElementSpace object
 * @param elem_idx Element index
 * @param dofs Output array for DOF indices (must be large enough)
 * @returns Number of DOFs for this element, or -1 on error
 */
int mfem_fespace_get_element_dofs(void* fes_ptr, int elem_idx, int* dofs) {
    if (!fes_ptr || !dofs) return -1;
    try {
        FiniteElementSpace* fes = static_cast<FiniteElementSpace*>(fes_ptr);
        Array<int> dof_array;
        fes->GetElementDofs(elem_idx, dof_array);
        for (int i = 0; i < dof_array.Size(); i++) {
            dofs[i] = dof_array[i];
        }
        return dof_array.Size();
    } catch (...) {
        return -1;
    }
}

/**
 * Get boundary DOFs for a given boundary attribute.
 * @param fes_ptr Pointer to FiniteElementSpace object
 * @param bdr_attr Boundary attribute (1-based)
 * @param dofs Output array for DOF indices
 * @param max_dofs Maximum number of DOFs to return
 * @returns Number of boundary DOFs found, or -1 on error
 */
int mfem_fespace_get_boundary_dofs(void* fes_ptr, int bdr_attr, int* dofs, int max_dofs) {
    if (!fes_ptr || !dofs) return -1;
    try {
        FiniteElementSpace* fes = static_cast<FiniteElementSpace*>(fes_ptr);
        Mesh* mesh = fes->GetMesh();

        Array<int> bdr_marker(mesh->bdr_attributes.Max());
        bdr_marker = 0;
        if (bdr_attr > 0 && bdr_attr <= bdr_marker.Size()) {
            bdr_marker[bdr_attr - 1] = 1;
        }

        Array<int> ess_tdof_list;
        fes->GetEssentialTrueDofs(bdr_marker, ess_tdof_list);

        int count = std::min(ess_tdof_list.Size(), max_dofs);
        for (int i = 0; i < count; i++) {
            dofs[i] = ess_tdof_list[i];
        }
        return ess_tdof_list.Size();
    } catch (...) {
        return -1;
    }
}

/**
 * Get essential (Dirichlet) DOFs for given boundary attributes.
 * @param fes_ptr Pointer to FiniteElementSpace object
 * @param bdr_attrs Array of boundary attributes (1-based)
 * @param num_attrs Number of attributes
 * @param dofs Output array for DOF indices
 * @param max_dofs Maximum number of DOFs to return
 * @returns Number of essential DOFs found, or -1 on error
 */
int mfem_fespace_get_essential_dofs(void* fes_ptr, const int* bdr_attrs, int num_attrs, int* dofs, int max_dofs) {
    if (!fes_ptr || !dofs) return -1;
    try {
        FiniteElementSpace* fes = static_cast<FiniteElementSpace*>(fes_ptr);
        Mesh* mesh = fes->GetMesh();

        Array<int> bdr_marker(mesh->bdr_attributes.Max());
        bdr_marker = 0;
        for (int i = 0; i < num_attrs; i++) {
            int attr = bdr_attrs[i];
            if (attr > 0 && attr <= bdr_marker.Size()) {
                bdr_marker[attr - 1] = 1;
            }
        }

        Array<int> ess_tdof_list;
        fes->GetEssentialTrueDofs(bdr_marker, ess_tdof_list);

        int count = std::min(ess_tdof_list.Size(), max_dofs);
        for (int i = 0; i < count; i++) {
            dofs[i] = ess_tdof_list[i];
        }
        return ess_tdof_list.Size();
    } catch (...) {
        return -1;
    }
}

/**
 * Check if mesh has NURBS extension.
 * @param mesh_ptr Pointer to Mesh object
 * @returns 1 if NURBS mesh, 0 otherwise
 */
int mfem_mesh_is_nurbs(void* mesh_ptr) {
    if (!mesh_ptr) return 0;
    Mesh* mesh = static_cast<Mesh*>(mesh_ptr);
    return mesh->NURBSext != nullptr ? 1 : 0;
}

/**
 * Get the NURBS order from a NURBS mesh.
 * @param mesh_ptr Pointer to Mesh object
 * @returns NURBS order, or -1 if not a NURBS mesh
 */
int mfem_mesh_get_nurbs_order(void* mesh_ptr) {
    if (!mesh_ptr) return -1;
    Mesh* mesh = static_cast<Mesh*>(mesh_ptr);
    if (!mesh->NURBSext) return -1;
    return mesh->NURBSext->GetOrder();
}

/**
 * Create a NURBS finite element space.
 * @param mesh_ptr Pointer to Mesh object (must be a NURBS mesh)
 * @param order Polynomial order (or -1 for variable order from mesh)
 * @param vdim Vector dimension
 * @returns Pointer to new FiniteElementSpace, or nullptr on failure
 */
void* mfem_fespace_create_nurbs(void* mesh_ptr, int order, int vdim) {
    if (!mesh_ptr) return nullptr;

    try {
        Mesh* mesh = static_cast<Mesh*>(mesh_ptr);
        if (!mesh->NURBSext) {
            // Not a NURBS mesh
            return nullptr;
        }

        // Use NURBSFECollection::VariableOrder (-1) to use the mesh's native order
        int nurbs_order = (order > 0) ? order : NURBSFECollection::VariableOrder;
        FiniteElementCollection* fec = new NURBSFECollection(nurbs_order);
        FiniteElementSpace* fes = new FiniteElementSpace(mesh, fec, vdim);
        return fes;
    } catch (...) {
        return nullptr;
    }
}

/**
 * Check if mesh has NCMesh (non-conforming mesh) extension.
 * @param mesh_ptr Pointer to Mesh object
 * @returns 1 if NCMesh, 0 otherwise
 */
int mfem_mesh_is_ncmesh(void* mesh_ptr) {
    if (!mesh_ptr) return 0;
    Mesh* mesh = static_cast<Mesh*>(mesh_ptr);
    return mesh->ncmesh != nullptr ? 1 : 0;
}

/**
 * Enable non-conforming mesh mode.
 * This must be called before local refinement if derefinement is desired later.
 * @param mesh_ptr Pointer to Mesh object
 * @returns 0 on success, -1 on failure
 */
int mfem_mesh_enable_ncmesh(void* mesh_ptr) {
    if (!mesh_ptr) return -1;
    try {
        Mesh* mesh = static_cast<Mesh*>(mesh_ptr);
        mesh->EnsureNCMesh();
        return 0;
    } catch (...) {
        return -1;
    }
}

/**
 * Derefine mesh elements (coarsening).
 * The mesh must have been refined with NCMesh enabled.
 * @param mesh_ptr Pointer to Mesh object
 * @param derefinements Array of derefinement table rows (from NCMesh)
 * @param num_derefs Number of derefinements
 * @returns 0 on success, -1 on failure
 */
int mfem_mesh_derefine(void* mesh_ptr, const int* derefinements, int num_derefs) {
    if (!mesh_ptr) return -1;
    try {
        Mesh* mesh = static_cast<Mesh*>(mesh_ptr);
        if (!mesh->ncmesh) {
            // Not a non-conforming mesh
            return -1;
        }

        Array<int> derefs(num_derefs);
        for (int i = 0; i < num_derefs; i++) {
            derefs[i] = derefinements[i];
        }

        mesh->Derefine(derefs);
        return 0;
    } catch (...) {
        return -1;
    }
}

/**
 * Get derefinement table size for NCMesh.
 * @param mesh_ptr Pointer to Mesh object
 * @returns Number of possible derefinements, or -1 on error
 */
int mfem_mesh_get_deref_table_size(void* mesh_ptr) {
    if (!mesh_ptr) return -1;
    try {
        Mesh* mesh = static_cast<Mesh*>(mesh_ptr);
        if (!mesh->ncmesh) return -1;
        return mesh->ncmesh->GetDerefinementTable().Size();
    } catch (...) {
        return -1;
    }
}

// ============================================================================
// GridFunction API
// ============================================================================

/**
 * Create a grid function on a finite element space.
 * @param fes_ptr Pointer to FiniteElementSpace object
 * @returns Pointer to new GridFunction, or nullptr on failure
 */
void* mfem_gridfunc_create(void* fes_ptr) {
    if (!fes_ptr) return nullptr;

    try {
        FiniteElementSpace* fes = static_cast<FiniteElementSpace*>(fes_ptr);
        return new GridFunction(fes);
    } catch (...) {
        return nullptr;
    }
}

/**
 * Destroy a grid function and free its memory.
 * @param gf_ptr Pointer to GridFunction object
 */
void mfem_gridfunc_destroy(void* gf_ptr) {
    if (gf_ptr) {
        delete static_cast<GridFunction*>(gf_ptr);
    }
}

/**
 * Get a pointer to the grid function's data array.
 * @param gf_ptr Pointer to GridFunction object
 * @returns Pointer to double array, or nullptr on error
 */
double* mfem_gridfunc_get_data(void* gf_ptr) {
    if (!gf_ptr) return nullptr;
    return static_cast<GridFunction*>(gf_ptr)->GetData();
}

/**
 * Get the size (number of elements) of the grid function's data array.
 * @param gf_ptr Pointer to GridFunction object
 * @returns Size, or -1 on error
 */
int mfem_gridfunc_get_size(void* gf_ptr) {
    if (!gf_ptr) return -1;
    return static_cast<GridFunction*>(gf_ptr)->Size();
}

/**
 * Project a constant value onto the grid function.
 * @param gf_ptr Pointer to GridFunction object
 * @param value Constant value to project
 * @returns 0 on success, -1 on failure
 */
int mfem_gridfunc_project_coefficient(void* gf_ptr, double value) {
    if (!gf_ptr) return -1;

    try {
        GridFunction* gf = static_cast<GridFunction*>(gf_ptr);
        ConstantCoefficient coef(value);
        gf->ProjectCoefficient(coef);
        return 0;
    } catch (...) {
        return -1;
    }
}

/**
 * Get the value of a grid function at a point within an element.
 * @param gf_ptr Pointer to GridFunction object
 * @param elem_idx Element index
 * @param xi Local coordinates (reference element) array
 * @param dim Dimension of xi array
 * @returns Value at the point, or NaN on error
 */
double mfem_gridfunc_get_value(void* gf_ptr, int elem_idx, const double* xi, int dim) {
    if (!gf_ptr || !xi) return std::nan("");

    try {
        GridFunction* gf = static_cast<GridFunction*>(gf_ptr);
        IntegrationPoint ip;
        ip.x = (dim > 0) ? xi[0] : 0.0;
        ip.y = (dim > 1) ? xi[1] : 0.0;
        ip.z = (dim > 2) ? xi[2] : 0.0;
        return gf->GetValue(elem_idx, ip);
    } catch (...) {
        return std::nan("");
    }
}

/**
 * Get the gradient of a grid function at a point within an element.
 * @param gf_ptr Pointer to GridFunction object
 * @param elem_idx Element index
 * @param xi Local coordinates (reference element) array
 * @param dim Dimension of xi array
 * @param grad_out Output array for gradient (size = space dimension)
 * @returns 0 on success, -1 on failure
 */
int mfem_gridfunc_get_gradient(void* gf_ptr, int elem_idx, const double* xi, int dim, double* grad_out) {
    if (!gf_ptr || !xi || !grad_out) return -1;

    try {
        GridFunction* gf = static_cast<GridFunction*>(gf_ptr);
        Mesh* mesh = gf->FESpace()->GetMesh();

        IntegrationPoint ip;
        ip.x = (dim > 0) ? xi[0] : 0.0;
        ip.y = (dim > 1) ? xi[1] : 0.0;
        ip.z = (dim > 2) ? xi[2] : 0.0;

        ElementTransformation* tr = mesh->GetElementTransformation(elem_idx);
        tr->SetIntPoint(&ip);

        Vector grad;
        gf->GetGradient(*tr, grad);

        for (int i = 0; i < grad.Size(); i++) {
            grad_out[i] = grad(i);
        }
        return 0;
    } catch (...) {
        return -1;
    }
}

/**
 * Get the curl of a vector grid function at a point within an element.
 * @param gf_ptr Pointer to GridFunction object
 * @param elem_idx Element index
 * @param xi Local coordinates array
 * @param dim Dimension of xi array
 * @param curl_out Output array for curl
 * @returns 0 on success, -1 on failure
 */
int mfem_gridfunc_get_curl(void* gf_ptr, int elem_idx, const double* xi, int dim, double* curl_out) {
    if (!gf_ptr || !xi || !curl_out) return -1;

    try {
        GridFunction* gf = static_cast<GridFunction*>(gf_ptr);
        Mesh* mesh = gf->FESpace()->GetMesh();

        IntegrationPoint ip;
        ip.x = (dim > 0) ? xi[0] : 0.0;
        ip.y = (dim > 1) ? xi[1] : 0.0;
        ip.z = (dim > 2) ? xi[2] : 0.0;

        ElementTransformation* tr = mesh->GetElementTransformation(elem_idx);
        tr->SetIntPoint(&ip);

        Vector curl;
        gf->GetCurl(*tr, curl);

        for (int i = 0; i < curl.Size(); i++) {
            curl_out[i] = curl(i);
        }
        return 0;
    } catch (...) {
        return -1;
    }
}

/**
 * Get the divergence of a vector grid function at a point within an element.
 * @param gf_ptr Pointer to GridFunction object
 * @param elem_idx Element index
 * @param xi Local coordinates array
 * @param dim Dimension of xi array
 * @returns Divergence value, or NaN on error
 */
double mfem_gridfunc_get_divergence(void* gf_ptr, int elem_idx, const double* xi, int dim) {
    if (!gf_ptr || !xi) return std::nan("");

    try {
        GridFunction* gf = static_cast<GridFunction*>(gf_ptr);
        Mesh* mesh = gf->FESpace()->GetMesh();

        IntegrationPoint ip;
        ip.x = (dim > 0) ? xi[0] : 0.0;
        ip.y = (dim > 1) ? xi[1] : 0.0;
        ip.z = (dim > 2) ? xi[2] : 0.0;

        ElementTransformation* tr = mesh->GetElementTransformation(elem_idx);
        tr->SetIntPoint(&ip);

        return gf->GetDivergence(*tr);
    } catch (...) {
        return std::nan("");
    }
}

/**
 * Compute the L2 norm of a grid function.
 * @param gf_ptr Pointer to GridFunction object
 * @returns L2 norm, or -1.0 on error
 */
double mfem_gridfunc_norm_l2(void* gf_ptr) {
    if (!gf_ptr) return -1.0;

    try {
        GridFunction* gf = static_cast<GridFunction*>(gf_ptr);
        ConstantCoefficient zero(0.0);
        return gf->ComputeL2Error(zero);
    } catch (...) {
        return -1.0;
    }
}

/**
 * Compute the L-infinity (max) norm of a grid function.
 * @param gf_ptr Pointer to GridFunction object
 * @returns L-infinity norm, or -1.0 on error
 */
double mfem_gridfunc_norm_linf(void* gf_ptr) {
    if (!gf_ptr) return -1.0;

    try {
        GridFunction* gf = static_cast<GridFunction*>(gf_ptr);
        return gf->Normlinf();
    } catch (...) {
        return -1.0;
    }
}

/**
 * Compute the H1 seminorm of a grid function (gradient norm).
 * @param gf_ptr Pointer to GridFunction object
 * @returns H1 seminorm, or -1.0 on error
 */
double mfem_gridfunc_norm_h1_semi(void* gf_ptr) {
    if (!gf_ptr) return -1.0;

    try {
        GridFunction* gf = static_cast<GridFunction*>(gf_ptr);
        ConstantCoefficient zero(0.0);
        VectorFunctionCoefficient zero_vec(gf->FESpace()->GetMesh()->SpaceDimension(),
            [](const Vector &x, Vector &v) { v = 0.0; });
        // H1 error with zero exact solution gives us the H1 norm
        return gf->ComputeH1Error(&zero, &zero_vec);
    } catch (...) {
        return -1.0;
    }
}

/**
 * Compute L2 error against a constant exact solution.
 * @param gf_ptr Pointer to GridFunction object
 * @param exact_value Exact constant value
 * @returns L2 error, or -1.0 on error
 */
double mfem_gridfunc_compute_l2_error_const(void* gf_ptr, double exact_value) {
    if (!gf_ptr) return -1.0;

    try {
        GridFunction* gf = static_cast<GridFunction*>(gf_ptr);
        ConstantCoefficient exact(exact_value);
        return gf->ComputeL2Error(exact);
    } catch (...) {
        return -1.0;
    }
}

/**
 * Project a function coefficient onto the grid function.
 * The function is evaluated at DOF locations for nodal interpolation.
 * @param gf_ptr Pointer to GridFunction object
 * @param func_ptr Function pointer that takes (x, y, z) and returns value
 * @returns 0 on success, -1 on failure
 */
typedef double (*ScalarFuncPtr)(double x, double y, double z);

int mfem_gridfunc_project_function(void* gf_ptr, ScalarFuncPtr func) {
    if (!gf_ptr || !func) return -1;

    try {
        GridFunction* gf = static_cast<GridFunction*>(gf_ptr);
        FunctionCoefficient coef([func](const Vector &x) {
            double px = x.Size() > 0 ? x(0) : 0.0;
            double py = x.Size() > 1 ? x(1) : 0.0;
            double pz = x.Size() > 2 ? x(2) : 0.0;
            return func(px, py, pz);
        });
        gf->ProjectCoefficient(coef);
        return 0;
    } catch (...) {
        return -1;
    }
}

/**
 * Project a coefficient onto boundary DOFs only.
 * @param gf_ptr Pointer to GridFunction object
 * @param bdr_attr Boundary attribute to project on
 * @param value Constant value to project
 * @returns 0 on success, -1 on failure
 */
int mfem_gridfunc_project_bdr_coefficient(void* gf_ptr, int bdr_attr, double value) {
    if (!gf_ptr) return -1;

    try {
        GridFunction* gf = static_cast<GridFunction*>(gf_ptr);
        FiniteElementSpace* fes = gf->FESpace();
        Mesh* mesh = fes->GetMesh();

        // Get boundary DOFs for this attribute
        Array<int> bdr_dofs;
        Array<int> bdr_attr_marker(mesh->bdr_attributes.Max());
        bdr_attr_marker = 0;
        if (bdr_attr > 0 && bdr_attr <= mesh->bdr_attributes.Max()) {
            bdr_attr_marker[bdr_attr - 1] = 1;
        }
        fes->GetEssentialTrueDofs(bdr_attr_marker, bdr_dofs);

        // Set the value at boundary DOFs
        for (int i = 0; i < bdr_dofs.Size(); i++) {
            (*gf)(bdr_dofs[i]) = value;
        }
        return 0;
    } catch (...) {
        return -1;
    }
}

/**
 * Save grid function to MFEM format string.
 * @param gf_ptr Pointer to GridFunction object
 * @param out_ptr Output buffer
 * @param max_len Maximum buffer length
 * @returns Actual length, or -1 on failure
 */
int mfem_gridfunc_to_mfem_string(void* gf_ptr, char* out_ptr, int max_len) {
    if (!gf_ptr || !out_ptr || max_len <= 0) return -1;

    try {
        GridFunction* gf = static_cast<GridFunction*>(gf_ptr);
        std::ostringstream oss;
        gf->Save(oss);
        std::string result = oss.str();

        int len = static_cast<int>(result.length());
        if (len >= max_len) {
            return len + 1; // Return required size
        }

        std::memcpy(out_ptr, result.c_str(), len + 1);
        return len;
    } catch (...) {
        return -1;
    }
}

/**
 * Get the length of grid function in MFEM format.
 * @param gf_ptr Pointer to GridFunction object
 * @returns String length, or -1 on failure
 */
int mfem_gridfunc_get_mfem_string_length(void* gf_ptr) {
    if (!gf_ptr) return -1;

    try {
        GridFunction* gf = static_cast<GridFunction*>(gf_ptr);
        std::ostringstream oss;
        gf->Save(oss);
        return static_cast<int>(oss.str().length());
    } catch (...) {
        return -1;
    }
}

/**
 * Save grid function to VTK format string (data only, mesh must be saved separately).
 * @param gf_ptr Pointer to GridFunction object
 * @param field_name Name of the field in VTK output
 * @param ref Refinement level for output
 * @param out_ptr Output buffer
 * @param max_len Maximum buffer length
 * @returns Actual length, or -1 on failure
 */
int mfem_gridfunc_to_vtk_string(void* gf_ptr, const char* field_name, int ref,
                                 char* out_ptr, int max_len) {
    if (!gf_ptr || !field_name || !out_ptr || max_len <= 0) return -1;

    try {
        GridFunction* gf = static_cast<GridFunction*>(gf_ptr);
        std::ostringstream oss;
        gf->SaveVTK(oss, field_name, ref);
        std::string result = oss.str();

        int len = static_cast<int>(result.length());
        if (len >= max_len) {
            return len + 1;
        }

        std::memcpy(out_ptr, result.c_str(), len + 1);
        return len;
    } catch (...) {
        return -1;
    }
}

/**
 * Get the length of grid function in VTK format.
 * @param gf_ptr Pointer to GridFunction object
 * @param field_name Name of the field
 * @param ref Refinement level
 * @returns String length, or -1 on failure
 */
int mfem_gridfunc_get_vtk_string_length(void* gf_ptr, const char* field_name, int ref) {
    if (!gf_ptr || !field_name) return -1;

    try {
        GridFunction* gf = static_cast<GridFunction*>(gf_ptr);
        std::ostringstream oss;
        gf->SaveVTK(oss, field_name, ref);
        return static_cast<int>(oss.str().length());
    } catch (...) {
        return -1;
    }
}

/**
 * Get the finite element space pointer from a grid function.
 * @param gf_ptr Pointer to GridFunction object
 * @returns Pointer to FiniteElementSpace, or nullptr on error
 */
void* mfem_gridfunc_get_fespace(void* gf_ptr) {
    if (!gf_ptr) return nullptr;
    return static_cast<GridFunction*>(gf_ptr)->FESpace();
}

} // extern "C"
