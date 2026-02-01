#!/bin/bash
#
# Build Eigen WebAssembly module
# Compiles C++ wrapper around Eigen library using Emscripten
#
# Eigen is a header-only C++ template library, so we compile
# our wrapper code that exposes key functionality via a C API.
#
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
SRC_DIR="$PROJECT_DIR/src/wasm"
OUT_DIR="$PROJECT_DIR/dist/wasm"
OBJ_DIR="$PROJECT_DIR/.build-eigen"

echo "Eigen WASM Build"
echo "================"
echo "Source: $SRC_DIR"
echo "Output: $OUT_DIR"
echo ""

# Create output/temp directories
mkdir -p "$OUT_DIR" "$OBJ_DIR"

# Check for Emscripten
if ! command -v emcc &> /dev/null; then
    echo "Error: emcc (Emscripten) not found in PATH"
    echo "Please install Emscripten and activate it with 'source emsdk_env.sh'"
    exit 1
fi

echo "Using Emscripten: $(emcc --version | head -1)"
echo ""

# Check for source files
if [ ! -f "$SRC_DIR/eigen_wasm.cpp" ]; then
    echo "Error: eigen_wasm.cpp not found at $SRC_DIR"
    exit 1
fi

if [ ! -d "$SRC_DIR/Eigen" ]; then
    echo "Error: Eigen headers not found at $SRC_DIR/Eigen"
    exit 1
fi

# Compiler flags
CXX_FLAGS=(
    -O2
    -std=c++17
    -I"$SRC_DIR"
    -DNDEBUG
    # Eigen-specific flags
    -DEIGEN_DONT_VECTORIZE          # Disable SIMD (basic build)
    -DEIGEN_DISABLE_UNALIGNED_ARRAY_ASSERT
    # Warning suppressions
    -Wno-deprecated-declarations
    -Wno-unused-parameter
    -Wno-sign-compare
)

# Optional: Enable WASM SIMD for better performance (requires modern browsers)
# Uncomment these lines to enable SIMD support:
# CXX_FLAGS+=(
#     -msimd128
#     -DEIGEN_WASM_SIMD128
# )
# Remove -DEIGEN_DONT_VECTORIZE if enabling SIMD

EIGEN_OBJS=()

echo "Compiling Eigen WASM wrapper..."
emcc -c "$SRC_DIR/eigen_wasm.cpp" -o "$OBJ_DIR/eigen_wasm.o" "${CXX_FLAGS[@]}"
EIGEN_OBJS+=("$OBJ_DIR/eigen_wasm.o")

echo "  Total object files: ${#EIGEN_OBJS[@]}"
echo ""

# Exported functions - all C functions from eigen_wasm.cpp
EXPORTED_FUNCTIONS='[
    "_eigen_matrix_create",
    "_eigen_matrix_destroy",
    "_eigen_matrix_create_from_data",
    "_eigen_matrix_get_rows",
    "_eigen_matrix_get_cols",
    "_eigen_matrix_get_data",
    "_eigen_matrix_set_element",
    "_eigen_matrix_get_element",
    "_eigen_matrix_set_identity",
    "_eigen_matrix_set_zero",
    "_eigen_matrix_set_ones",
    "_eigen_matrix_set_random",
    "_eigen_matrix_copy",
    "_eigen_matrix_transpose",
    "_eigen_matrix_conjugate",
    "_eigen_matrix_adjoint",
    "_eigen_matrix_add",
    "_eigen_matrix_subtract",
    "_eigen_matrix_multiply",
    "_eigen_matrix_scalar_multiply",
    "_eigen_matrix_scalar_add",
    "_eigen_matrix_elementwise_multiply",
    "_eigen_matrix_elementwise_divide",
    "_eigen_matrix_norm",
    "_eigen_matrix_squared_norm",
    "_eigen_matrix_normalize",
    "_eigen_matrix_sum",
    "_eigen_matrix_prod",
    "_eigen_matrix_mean",
    "_eigen_matrix_min_coeff",
    "_eigen_matrix_max_coeff",
    "_eigen_matrix_trace",
    "_eigen_matrix_determinant",
    "_eigen_matrix_inverse",
    "_eigen_matrix_lu_decompose",
    "_eigen_matrix_lu_solve",
    "_eigen_matrix_qr_decompose",
    "_eigen_matrix_qr_solve",
    "_eigen_matrix_cholesky_decompose",
    "_eigen_matrix_cholesky_solve",
    "_eigen_matrix_svd",
    "_eigen_matrix_eigenvalues_symmetric",
    "_eigen_matrix_eigenvalues_general",
    "_eigen_vector_create",
    "_eigen_vector_destroy",
    "_eigen_vector_create_from_data",
    "_eigen_vector_get_size",
    "_eigen_vector_get_data",
    "_eigen_vector_set_element",
    "_eigen_vector_get_element",
    "_eigen_vector_dot",
    "_eigen_vector_cross",
    "_eigen_vector_norm",
    "_eigen_vector_normalize",
    "_eigen_sparse_matrix_create",
    "_eigen_sparse_matrix_destroy",
    "_eigen_sparse_matrix_set_from_triplets",
    "_eigen_sparse_matrix_get_nnz",
    "_eigen_sparse_matrix_multiply_vector",
    "_eigen_sparse_solve_lu",
    "_eigen_sparse_solve_cg",
    "_eigen_sparse_solve_bicgstab",
    "_malloc",
    "_free"
]'

EXPORTED_RUNTIME='["ccall", "cwrap", "getValue", "setValue", "HEAPF64", "HEAPF32", "HEAP32", "HEAPU8"]'

LINK_FLAGS=(
    -s WASM=1
    -s MODULARIZE=1
    -s EXPORTED_FUNCTIONS="$EXPORTED_FUNCTIONS"
    -s EXPORTED_RUNTIME_METHODS="$EXPORTED_RUNTIME"
    -s ALLOW_MEMORY_GROWTH=1
    -s INITIAL_MEMORY=16777216
    -s STACK_SIZE=1048576
    -O2
)

# Build CJS version
echo "Linking CJS module..."
emcc \
    "${EIGEN_OBJS[@]}" \
    "${LINK_FLAGS[@]}" \
    -s EXPORT_NAME="createEigenModule" \
    -o "$OUT_DIR/eigen.cjs"

# Build ESM version
echo "Linking ESM module..."
emcc \
    "${EIGEN_OBJS[@]}" \
    "${LINK_FLAGS[@]}" \
    -s EXPORT_NAME="createEigenModule" \
    -s EXPORT_ES6=1 \
    -o "$OUT_DIR/eigen.mjs"

# Clean up object files
rm -rf "$OBJ_DIR"

echo ""
echo "Build complete!"
echo "  CJS: $OUT_DIR/eigen.cjs"
echo "  ESM: $OUT_DIR/eigen.mjs"
echo "  WASM: $OUT_DIR/eigen.wasm"
echo ""

if [ -f "$OUT_DIR/eigen.wasm" ]; then
    echo "File sizes:"
    ls -lh "$OUT_DIR/eigen.cjs" "$OUT_DIR/eigen.mjs" "$OUT_DIR/eigen.wasm"
fi
