/**
 * Eigen WebAssembly C API Wrapper
 *
 * This file provides a C-compatible API for Eigen functionality,
 * designed to be exported via Emscripten for use in JavaScript/TypeScript.
 *
 * All functions use extern "C" linkage and operate on opaque pointers
 * to Eigen objects.
 */

#include "Eigen/Dense"
#include "Eigen/Sparse"
#include "Eigen/Eigenvalues"
#include <cstring>
#include <cstdlib>
#include <vector>

using namespace Eigen;

// Type aliases for clarity
using MatrixXdPtr = MatrixXd*;
using VectorXdPtr = VectorXd*;
using SparseMatrixPtr = SparseMatrix<double>*;

extern "C" {

// ============================================================================
// Dense Matrix API
// ============================================================================

/**
 * Create a new matrix with given dimensions.
 * @param rows Number of rows
 * @param cols Number of columns
 * @returns Pointer to new MatrixXd, or nullptr on failure
 */
void* eigen_matrix_create(int rows, int cols) {
    if (rows <= 0 || cols <= 0) return nullptr;
    try {
        return new MatrixXd(rows, cols);
    } catch (...) {
        return nullptr;
    }
}

/**
 * Destroy a matrix object and free its memory.
 * @param mat_ptr Pointer to MatrixXd object
 */
void eigen_matrix_destroy(void* mat_ptr) {
    if (mat_ptr) {
        delete static_cast<MatrixXdPtr>(mat_ptr);
    }
}

/**
 * Create a matrix from existing data (column-major order).
 * @param rows Number of rows
 * @param cols Number of columns
 * @param data Pointer to data array (column-major, size = rows * cols)
 * @returns Pointer to new MatrixXd, or nullptr on failure
 */
void* eigen_matrix_create_from_data(int rows, int cols, const double* data) {
    if (rows <= 0 || cols <= 0 || !data) return nullptr;
    try {
        MatrixXd* mat = new MatrixXd(rows, cols);
        memcpy(mat->data(), data, rows * cols * sizeof(double));
        return mat;
    } catch (...) {
        return nullptr;
    }
}

/**
 * Get the number of rows in a matrix.
 */
int eigen_matrix_get_rows(void* mat_ptr) {
    if (!mat_ptr) return 0;
    return static_cast<MatrixXdPtr>(mat_ptr)->rows();
}

/**
 * Get the number of columns in a matrix.
 */
int eigen_matrix_get_cols(void* mat_ptr) {
    if (!mat_ptr) return 0;
    return static_cast<MatrixXdPtr>(mat_ptr)->cols();
}

/**
 * Get pointer to the matrix data (column-major order).
 */
double* eigen_matrix_get_data(void* mat_ptr) {
    if (!mat_ptr) return nullptr;
    return static_cast<MatrixXdPtr>(mat_ptr)->data();
}

/**
 * Set a single element in the matrix.
 */
void eigen_matrix_set_element(void* mat_ptr, int row, int col, double value) {
    if (!mat_ptr) return;
    MatrixXd* mat = static_cast<MatrixXdPtr>(mat_ptr);
    if (row >= 0 && row < mat->rows() && col >= 0 && col < mat->cols()) {
        (*mat)(row, col) = value;
    }
}

/**
 * Get a single element from the matrix.
 */
double eigen_matrix_get_element(void* mat_ptr, int row, int col) {
    if (!mat_ptr) return 0.0;
    MatrixXd* mat = static_cast<MatrixXdPtr>(mat_ptr);
    if (row >= 0 && row < mat->rows() && col >= 0 && col < mat->cols()) {
        return (*mat)(row, col);
    }
    return 0.0;
}

/**
 * Set matrix to identity.
 */
void eigen_matrix_set_identity(void* mat_ptr) {
    if (!mat_ptr) return;
    static_cast<MatrixXdPtr>(mat_ptr)->setIdentity();
}

/**
 * Set all elements to zero.
 */
void eigen_matrix_set_zero(void* mat_ptr) {
    if (!mat_ptr) return;
    static_cast<MatrixXdPtr>(mat_ptr)->setZero();
}

/**
 * Set all elements to one.
 */
void eigen_matrix_set_ones(void* mat_ptr) {
    if (!mat_ptr) return;
    static_cast<MatrixXdPtr>(mat_ptr)->setOnes();
}

/**
 * Set elements to random values in [-1, 1].
 */
void eigen_matrix_set_random(void* mat_ptr) {
    if (!mat_ptr) return;
    static_cast<MatrixXdPtr>(mat_ptr)->setRandom();
}

/**
 * Copy a matrix.
 * @returns Pointer to new MatrixXd copy, or nullptr on failure
 */
void* eigen_matrix_copy(void* mat_ptr) {
    if (!mat_ptr) return nullptr;
    try {
        return new MatrixXd(*static_cast<MatrixXdPtr>(mat_ptr));
    } catch (...) {
        return nullptr;
    }
}

/**
 * Transpose a matrix.
 * @returns Pointer to new transposed MatrixXd, or nullptr on failure
 */
void* eigen_matrix_transpose(void* mat_ptr) {
    if (!mat_ptr) return nullptr;
    try {
        MatrixXd* mat = static_cast<MatrixXdPtr>(mat_ptr);
        return new MatrixXd(mat->transpose());
    } catch (...) {
        return nullptr;
    }
}

/**
 * Complex conjugate of a matrix (no-op for real matrices).
 */
void* eigen_matrix_conjugate(void* mat_ptr) {
    if (!mat_ptr) return nullptr;
    try {
        return new MatrixXd(*static_cast<MatrixXdPtr>(mat_ptr));
    } catch (...) {
        return nullptr;
    }
}

/**
 * Adjoint (conjugate transpose) of a matrix.
 */
void* eigen_matrix_adjoint(void* mat_ptr) {
    if (!mat_ptr) return nullptr;
    try {
        MatrixXd* mat = static_cast<MatrixXdPtr>(mat_ptr);
        return new MatrixXd(mat->adjoint());
    } catch (...) {
        return nullptr;
    }
}

// ============================================================================
// Matrix Arithmetic
// ============================================================================

/**
 * Add two matrices: result = a + b
 */
void* eigen_matrix_add(void* a_ptr, void* b_ptr) {
    if (!a_ptr || !b_ptr) return nullptr;
    try {
        MatrixXd* a = static_cast<MatrixXdPtr>(a_ptr);
        MatrixXd* b = static_cast<MatrixXdPtr>(b_ptr);
        if (a->rows() != b->rows() || a->cols() != b->cols()) return nullptr;
        return new MatrixXd(*a + *b);
    } catch (...) {
        return nullptr;
    }
}

/**
 * Subtract two matrices: result = a - b
 */
void* eigen_matrix_subtract(void* a_ptr, void* b_ptr) {
    if (!a_ptr || !b_ptr) return nullptr;
    try {
        MatrixXd* a = static_cast<MatrixXdPtr>(a_ptr);
        MatrixXd* b = static_cast<MatrixXdPtr>(b_ptr);
        if (a->rows() != b->rows() || a->cols() != b->cols()) return nullptr;
        return new MatrixXd(*a - *b);
    } catch (...) {
        return nullptr;
    }
}

/**
 * Multiply two matrices: result = a * b
 */
void* eigen_matrix_multiply(void* a_ptr, void* b_ptr) {
    if (!a_ptr || !b_ptr) return nullptr;
    try {
        MatrixXd* a = static_cast<MatrixXdPtr>(a_ptr);
        MatrixXd* b = static_cast<MatrixXdPtr>(b_ptr);
        if (a->cols() != b->rows()) return nullptr;
        return new MatrixXd(*a * *b);
    } catch (...) {
        return nullptr;
    }
}

/**
 * Scalar multiplication: result = mat * scalar
 */
void* eigen_matrix_scalar_multiply(void* mat_ptr, double scalar) {
    if (!mat_ptr) return nullptr;
    try {
        MatrixXd* mat = static_cast<MatrixXdPtr>(mat_ptr);
        return new MatrixXd(*mat * scalar);
    } catch (...) {
        return nullptr;
    }
}

/**
 * Scalar addition: result = mat + scalar (adds to each element)
 */
void* eigen_matrix_scalar_add(void* mat_ptr, double scalar) {
    if (!mat_ptr) return nullptr;
    try {
        MatrixXd* mat = static_cast<MatrixXdPtr>(mat_ptr);
        return new MatrixXd(mat->array() + scalar);
    } catch (...) {
        return nullptr;
    }
}

/**
 * Element-wise multiplication: result = a .* b
 */
void* eigen_matrix_elementwise_multiply(void* a_ptr, void* b_ptr) {
    if (!a_ptr || !b_ptr) return nullptr;
    try {
        MatrixXd* a = static_cast<MatrixXdPtr>(a_ptr);
        MatrixXd* b = static_cast<MatrixXdPtr>(b_ptr);
        if (a->rows() != b->rows() || a->cols() != b->cols()) return nullptr;
        return new MatrixXd(a->array() * b->array());
    } catch (...) {
        return nullptr;
    }
}

/**
 * Element-wise division: result = a ./ b
 */
void* eigen_matrix_elementwise_divide(void* a_ptr, void* b_ptr) {
    if (!a_ptr || !b_ptr) return nullptr;
    try {
        MatrixXd* a = static_cast<MatrixXdPtr>(a_ptr);
        MatrixXd* b = static_cast<MatrixXdPtr>(b_ptr);
        if (a->rows() != b->rows() || a->cols() != b->cols()) return nullptr;
        return new MatrixXd(a->array() / b->array());
    } catch (...) {
        return nullptr;
    }
}

// ============================================================================
// Matrix Reductions
// ============================================================================

/**
 * Compute the Frobenius norm of a matrix.
 */
double eigen_matrix_norm(void* mat_ptr) {
    if (!mat_ptr) return 0.0;
    return static_cast<MatrixXdPtr>(mat_ptr)->norm();
}

/**
 * Compute the squared Frobenius norm of a matrix.
 */
double eigen_matrix_squared_norm(void* mat_ptr) {
    if (!mat_ptr) return 0.0;
    return static_cast<MatrixXdPtr>(mat_ptr)->squaredNorm();
}

/**
 * Normalize a matrix in-place (divide by norm).
 */
void eigen_matrix_normalize(void* mat_ptr) {
    if (!mat_ptr) return;
    static_cast<MatrixXdPtr>(mat_ptr)->normalize();
}

/**
 * Sum of all elements.
 */
double eigen_matrix_sum(void* mat_ptr) {
    if (!mat_ptr) return 0.0;
    return static_cast<MatrixXdPtr>(mat_ptr)->sum();
}

/**
 * Product of all elements.
 */
double eigen_matrix_prod(void* mat_ptr) {
    if (!mat_ptr) return 0.0;
    return static_cast<MatrixXdPtr>(mat_ptr)->prod();
}

/**
 * Mean of all elements.
 */
double eigen_matrix_mean(void* mat_ptr) {
    if (!mat_ptr) return 0.0;
    return static_cast<MatrixXdPtr>(mat_ptr)->mean();
}

/**
 * Minimum coefficient.
 */
double eigen_matrix_min_coeff(void* mat_ptr) {
    if (!mat_ptr) return 0.0;
    return static_cast<MatrixXdPtr>(mat_ptr)->minCoeff();
}

/**
 * Maximum coefficient.
 */
double eigen_matrix_max_coeff(void* mat_ptr) {
    if (!mat_ptr) return 0.0;
    return static_cast<MatrixXdPtr>(mat_ptr)->maxCoeff();
}

/**
 * Trace (sum of diagonal elements).
 */
double eigen_matrix_trace(void* mat_ptr) {
    if (!mat_ptr) return 0.0;
    return static_cast<MatrixXdPtr>(mat_ptr)->trace();
}

/**
 * Determinant of a square matrix.
 */
double eigen_matrix_determinant(void* mat_ptr) {
    if (!mat_ptr) return 0.0;
    MatrixXd* mat = static_cast<MatrixXdPtr>(mat_ptr);
    if (mat->rows() != mat->cols()) return 0.0;
    return mat->determinant();
}

/**
 * Inverse of a square matrix.
 * @returns Pointer to new inverse MatrixXd, or nullptr on failure
 */
void* eigen_matrix_inverse(void* mat_ptr) {
    if (!mat_ptr) return nullptr;
    try {
        MatrixXd* mat = static_cast<MatrixXdPtr>(mat_ptr);
        if (mat->rows() != mat->cols()) return nullptr;
        return new MatrixXd(mat->inverse());
    } catch (...) {
        return nullptr;
    }
}

// ============================================================================
// Matrix Decompositions
// ============================================================================

/**
 * LU decomposition with partial pivoting.
 * @param mat_ptr Input matrix
 * @param l_out Output: lower triangular matrix L
 * @param u_out Output: upper triangular matrix U
 * @param p_out Output: permutation matrix P (as indices)
 * @returns 0 on success, -1 on failure
 */
int eigen_matrix_lu_decompose(void* mat_ptr, void** l_out, void** u_out, int* p_out) {
    if (!mat_ptr || !l_out || !u_out) return -1;
    try {
        MatrixXd* mat = static_cast<MatrixXdPtr>(mat_ptr);
        PartialPivLU<MatrixXd> lu(*mat);

        *l_out = new MatrixXd(lu.matrixLU().triangularView<StrictlyLower>());
        static_cast<MatrixXdPtr>(*l_out)->diagonal().setOnes();
        *u_out = new MatrixXd(lu.matrixLU().triangularView<Upper>());

        if (p_out) {
            auto perm = lu.permutationP();
            for (int i = 0; i < perm.size(); ++i) {
                p_out[i] = perm.indices()(i);
            }
        }
        return 0;
    } catch (...) {
        return -1;
    }
}

/**
 * Solve linear system Ax = b using LU decomposition.
 * @param a_ptr Coefficient matrix A
 * @param b_ptr Right-hand side b
 * @returns Solution vector x, or nullptr on failure
 */
void* eigen_matrix_lu_solve(void* a_ptr, void* b_ptr) {
    if (!a_ptr || !b_ptr) return nullptr;
    try {
        MatrixXd* A = static_cast<MatrixXdPtr>(a_ptr);
        MatrixXd* b = static_cast<MatrixXdPtr>(b_ptr);
        return new MatrixXd(A->partialPivLu().solve(*b));
    } catch (...) {
        return nullptr;
    }
}

/**
 * QR decomposition.
 * @param mat_ptr Input matrix
 * @param q_out Output: orthogonal matrix Q
 * @param r_out Output: upper triangular matrix R
 * @returns 0 on success, -1 on failure
 */
int eigen_matrix_qr_decompose(void* mat_ptr, void** q_out, void** r_out) {
    if (!mat_ptr || !q_out || !r_out) return -1;
    try {
        MatrixXd* mat = static_cast<MatrixXdPtr>(mat_ptr);
        HouseholderQR<MatrixXd> qr(*mat);

        int m = mat->rows();
        int n = mat->cols();
        int k = std::min(m, n);

        *q_out = new MatrixXd(qr.householderQ() * MatrixXd::Identity(m, k));
        *r_out = new MatrixXd(qr.matrixQR().topRows(k).triangularView<Upper>());
        return 0;
    } catch (...) {
        return -1;
    }
}

/**
 * Solve linear system Ax = b using QR decomposition.
 */
void* eigen_matrix_qr_solve(void* a_ptr, void* b_ptr) {
    if (!a_ptr || !b_ptr) return nullptr;
    try {
        MatrixXd* A = static_cast<MatrixXdPtr>(a_ptr);
        MatrixXd* b = static_cast<MatrixXdPtr>(b_ptr);
        return new MatrixXd(A->householderQr().solve(*b));
    } catch (...) {
        return nullptr;
    }
}

/**
 * Cholesky decomposition (for symmetric positive definite matrices).
 * @param mat_ptr Input matrix (must be symmetric positive definite)
 * @param l_out Output: lower triangular matrix L such that A = L * L^T
 * @returns 0 on success, -1 on failure
 */
int eigen_matrix_cholesky_decompose(void* mat_ptr, void** l_out) {
    if (!mat_ptr || !l_out) return -1;
    try {
        MatrixXd* mat = static_cast<MatrixXdPtr>(mat_ptr);
        LLT<MatrixXd> llt(*mat);
        if (llt.info() != Success) return -1;
        *l_out = new MatrixXd(llt.matrixL());
        return 0;
    } catch (...) {
        return -1;
    }
}

/**
 * Solve linear system Ax = b using Cholesky decomposition.
 */
void* eigen_matrix_cholesky_solve(void* a_ptr, void* b_ptr) {
    if (!a_ptr || !b_ptr) return nullptr;
    try {
        MatrixXd* A = static_cast<MatrixXdPtr>(a_ptr);
        MatrixXd* b = static_cast<MatrixXdPtr>(b_ptr);
        LLT<MatrixXd> llt(*A);
        if (llt.info() != Success) return nullptr;
        return new MatrixXd(llt.solve(*b));
    } catch (...) {
        return nullptr;
    }
}

/**
 * Singular Value Decomposition (SVD).
 * @param mat_ptr Input matrix
 * @param u_out Output: left singular vectors U
 * @param s_out Output: singular values (as diagonal)
 * @param v_out Output: right singular vectors V
 * @returns 0 on success, -1 on failure
 */
int eigen_matrix_svd(void* mat_ptr, void** u_out, double* s_out, void** v_out) {
    if (!mat_ptr || !u_out || !s_out || !v_out) return -1;
    try {
        MatrixXd* mat = static_cast<MatrixXdPtr>(mat_ptr);
        JacobiSVD<MatrixXd> svd(*mat, ComputeThinU | ComputeThinV);

        *u_out = new MatrixXd(svd.matrixU());
        *v_out = new MatrixXd(svd.matrixV());

        VectorXd sv = svd.singularValues();
        for (int i = 0; i < sv.size(); ++i) {
            s_out[i] = sv(i);
        }
        return 0;
    } catch (...) {
        return -1;
    }
}

/**
 * Eigenvalue decomposition for symmetric matrices.
 * @param mat_ptr Input symmetric matrix
 * @param eigenvalues_out Output: eigenvalues (real)
 * @param eigenvectors_out Output: eigenvectors matrix
 * @returns 0 on success, -1 on failure
 */
int eigen_matrix_eigenvalues_symmetric(void* mat_ptr, double* eigenvalues_out, void** eigenvectors_out) {
    if (!mat_ptr || !eigenvalues_out || !eigenvectors_out) return -1;
    try {
        MatrixXd* mat = static_cast<MatrixXdPtr>(mat_ptr);
        SelfAdjointEigenSolver<MatrixXd> solver(*mat);
        if (solver.info() != Success) return -1;

        VectorXd ev = solver.eigenvalues();
        for (int i = 0; i < ev.size(); ++i) {
            eigenvalues_out[i] = ev(i);
        }
        *eigenvectors_out = new MatrixXd(solver.eigenvectors());
        return 0;
    } catch (...) {
        return -1;
    }
}

/**
 * Eigenvalue decomposition for general matrices.
 * @param mat_ptr Input matrix
 * @param eigenvalues_real_out Output: real parts of eigenvalues
 * @param eigenvalues_imag_out Output: imaginary parts of eigenvalues
 * @param eigenvectors_out Output: eigenvectors matrix (may be complex)
 * @returns 0 on success, -1 on failure
 */
int eigen_matrix_eigenvalues_general(void* mat_ptr, double* eigenvalues_real_out,
                                      double* eigenvalues_imag_out, void** eigenvectors_out) {
    if (!mat_ptr || !eigenvalues_real_out || !eigenvalues_imag_out) return -1;
    try {
        MatrixXd* mat = static_cast<MatrixXdPtr>(mat_ptr);
        EigenSolver<MatrixXd> solver(*mat);
        if (solver.info() != Success) return -1;

        auto ev = solver.eigenvalues();
        for (int i = 0; i < ev.size(); ++i) {
            eigenvalues_real_out[i] = ev(i).real();
            eigenvalues_imag_out[i] = ev(i).imag();
        }

        if (eigenvectors_out) {
            // Return real part of eigenvectors for real eigenvalues
            *eigenvectors_out = new MatrixXd(solver.eigenvectors().real());
        }
        return 0;
    } catch (...) {
        return -1;
    }
}

// ============================================================================
// Vector API
// ============================================================================

/**
 * Create a new vector with given size.
 */
void* eigen_vector_create(int size) {
    if (size <= 0) return nullptr;
    try {
        return new VectorXd(size);
    } catch (...) {
        return nullptr;
    }
}

/**
 * Destroy a vector object.
 */
void eigen_vector_destroy(void* vec_ptr) {
    if (vec_ptr) {
        delete static_cast<VectorXdPtr>(vec_ptr);
    }
}

/**
 * Create a vector from existing data.
 */
void* eigen_vector_create_from_data(int size, const double* data) {
    if (size <= 0 || !data) return nullptr;
    try {
        VectorXd* vec = new VectorXd(size);
        memcpy(vec->data(), data, size * sizeof(double));
        return vec;
    } catch (...) {
        return nullptr;
    }
}

/**
 * Get vector size.
 */
int eigen_vector_get_size(void* vec_ptr) {
    if (!vec_ptr) return 0;
    return static_cast<VectorXdPtr>(vec_ptr)->size();
}

/**
 * Get pointer to vector data.
 */
double* eigen_vector_get_data(void* vec_ptr) {
    if (!vec_ptr) return nullptr;
    return static_cast<VectorXdPtr>(vec_ptr)->data();
}

/**
 * Set a single element.
 */
void eigen_vector_set_element(void* vec_ptr, int index, double value) {
    if (!vec_ptr) return;
    VectorXd* vec = static_cast<VectorXdPtr>(vec_ptr);
    if (index >= 0 && index < vec->size()) {
        (*vec)(index) = value;
    }
}

/**
 * Get a single element.
 */
double eigen_vector_get_element(void* vec_ptr, int index) {
    if (!vec_ptr) return 0.0;
    VectorXd* vec = static_cast<VectorXdPtr>(vec_ptr);
    if (index >= 0 && index < vec->size()) {
        return (*vec)(index);
    }
    return 0.0;
}

/**
 * Dot product of two vectors.
 */
double eigen_vector_dot(void* a_ptr, void* b_ptr) {
    if (!a_ptr || !b_ptr) return 0.0;
    VectorXd* a = static_cast<VectorXdPtr>(a_ptr);
    VectorXd* b = static_cast<VectorXdPtr>(b_ptr);
    if (a->size() != b->size()) return 0.0;
    return a->dot(*b);
}

/**
 * Cross product of two 3D vectors.
 */
void* eigen_vector_cross(void* a_ptr, void* b_ptr) {
    if (!a_ptr || !b_ptr) return nullptr;
    try {
        VectorXd* a = static_cast<VectorXdPtr>(a_ptr);
        VectorXd* b = static_cast<VectorXdPtr>(b_ptr);
        if (a->size() != 3 || b->size() != 3) return nullptr;

        Vector3d av((*a)(0), (*a)(1), (*a)(2));
        Vector3d bv((*b)(0), (*b)(1), (*b)(2));
        Vector3d result = av.cross(bv);

        VectorXd* out = new VectorXd(3);
        (*out)(0) = result(0);
        (*out)(1) = result(1);
        (*out)(2) = result(2);
        return out;
    } catch (...) {
        return nullptr;
    }
}

/**
 * Euclidean norm of a vector.
 */
double eigen_vector_norm(void* vec_ptr) {
    if (!vec_ptr) return 0.0;
    return static_cast<VectorXdPtr>(vec_ptr)->norm();
}

/**
 * Normalize a vector in-place.
 */
void eigen_vector_normalize(void* vec_ptr) {
    if (!vec_ptr) return;
    static_cast<VectorXdPtr>(vec_ptr)->normalize();
}

// ============================================================================
// Sparse Matrix API
// ============================================================================

/**
 * Create a new sparse matrix with given dimensions.
 */
void* eigen_sparse_matrix_create(int rows, int cols) {
    if (rows <= 0 || cols <= 0) return nullptr;
    try {
        return new SparseMatrix<double>(rows, cols);
    } catch (...) {
        return nullptr;
    }
}

/**
 * Destroy a sparse matrix.
 */
void eigen_sparse_matrix_destroy(void* sp_ptr) {
    if (sp_ptr) {
        delete static_cast<SparseMatrixPtr>(sp_ptr);
    }
}

/**
 * Set sparse matrix from triplets (row, col, value).
 * @param sp_ptr Sparse matrix pointer
 * @param rows Array of row indices
 * @param cols Array of column indices
 * @param vals Array of values
 * @param nnz Number of non-zeros
 * @returns 0 on success, -1 on failure
 */
int eigen_sparse_matrix_set_from_triplets(void* sp_ptr, const int* rows,
                                           const int* cols, const double* vals, int nnz) {
    if (!sp_ptr || !rows || !cols || !vals || nnz < 0) return -1;
    try {
        SparseMatrix<double>* sp = static_cast<SparseMatrixPtr>(sp_ptr);
        std::vector<Triplet<double>> triplets;
        triplets.reserve(nnz);
        for (int i = 0; i < nnz; ++i) {
            triplets.emplace_back(rows[i], cols[i], vals[i]);
        }
        sp->setFromTriplets(triplets.begin(), triplets.end());
        return 0;
    } catch (...) {
        return -1;
    }
}

/**
 * Get number of non-zeros in sparse matrix.
 */
int eigen_sparse_matrix_get_nnz(void* sp_ptr) {
    if (!sp_ptr) return 0;
    return static_cast<SparseMatrixPtr>(sp_ptr)->nonZeros();
}

/**
 * Sparse matrix-vector multiplication: result = A * x
 */
void* eigen_sparse_matrix_multiply_vector(void* sp_ptr, void* vec_ptr) {
    if (!sp_ptr || !vec_ptr) return nullptr;
    try {
        SparseMatrix<double>* A = static_cast<SparseMatrixPtr>(sp_ptr);
        VectorXd* x = static_cast<VectorXdPtr>(vec_ptr);
        if (A->cols() != x->size()) return nullptr;
        return new VectorXd(*A * *x);
    } catch (...) {
        return nullptr;
    }
}

/**
 * Solve sparse linear system Ax = b using SparseLU.
 */
void* eigen_sparse_solve_lu(void* sp_ptr, void* b_ptr) {
    if (!sp_ptr || !b_ptr) return nullptr;
    try {
        SparseMatrix<double>* A = static_cast<SparseMatrixPtr>(sp_ptr);
        VectorXd* b = static_cast<VectorXdPtr>(b_ptr);

        SparseLU<SparseMatrix<double>> solver;
        solver.compute(*A);
        if (solver.info() != Success) return nullptr;

        return new VectorXd(solver.solve(*b));
    } catch (...) {
        return nullptr;
    }
}

/**
 * Solve sparse linear system Ax = b using Conjugate Gradient (for SPD matrices).
 */
void* eigen_sparse_solve_cg(void* sp_ptr, void* b_ptr) {
    if (!sp_ptr || !b_ptr) return nullptr;
    try {
        SparseMatrix<double>* A = static_cast<SparseMatrixPtr>(sp_ptr);
        VectorXd* b = static_cast<VectorXdPtr>(b_ptr);

        ConjugateGradient<SparseMatrix<double>> solver;
        solver.compute(*A);
        if (solver.info() != Success) return nullptr;

        return new VectorXd(solver.solve(*b));
    } catch (...) {
        return nullptr;
    }
}

/**
 * Solve sparse linear system Ax = b using BiCGSTAB (for general matrices).
 */
void* eigen_sparse_solve_bicgstab(void* sp_ptr, void* b_ptr) {
    if (!sp_ptr || !b_ptr) return nullptr;
    try {
        SparseMatrix<double>* A = static_cast<SparseMatrixPtr>(sp_ptr);
        VectorXd* b = static_cast<VectorXdPtr>(b_ptr);

        BiCGSTAB<SparseMatrix<double>> solver;
        solver.compute(*A);
        if (solver.info() != Success) return nullptr;

        return new VectorXd(solver.solve(*b));
    } catch (...) {
        return nullptr;
    }
}

} // extern "C"
