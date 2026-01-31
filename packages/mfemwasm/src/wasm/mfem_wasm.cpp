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

} // extern "C"
