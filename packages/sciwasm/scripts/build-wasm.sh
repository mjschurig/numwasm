#!/bin/bash
#
# Build SciWASM WebAssembly module using Emscripten
#
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
SRC_DIR="$PROJECT_DIR/src/wasm"
OUT_DIR="$PROJECT_DIR/dist/wasm"
OBJ_DIR="$PROJECT_DIR/.build-wasm"

echo "SciWASM Build"
echo "============="
echo "Source: $SRC_DIR"
echo "Output: $OUT_DIR"
echo ""

# Create output/temp directories
mkdir -p "$OUT_DIR" "$OBJ_DIR"

# Check for Emscripten
if ! command -v emcc &> /dev/null; then
    echo "Error: emcc (Emscripten) not found in PATH"
    echo "Please install Emscripten or rebuild the devcontainer"
    exit 1
fi

echo "Using Emscripten: $(emcc --version | head -1)"
echo ""

# Step 1: Compile C source files to .o
C_OBJS=()
for cfile in "$SRC_DIR/optimize/nelder_mead.c" "$SRC_DIR/optimize/bfgs.c" "$SRC_DIR/optimize/lbfgsb.c" "$SRC_DIR/optimize/blas_lite.c" "$SRC_DIR/quadpack.c" "$SRC_DIR/quadpack_wasm.c" "$SRC_DIR/special/comb.c"; do
    if [ -f "$cfile" ]; then
        base=$(basename "$cfile" .c)
        obj="$OBJ_DIR/${base}.o"
        echo "Compiling C: $cfile"
        emcc -c "$cfile" -o "$obj" -O2
        C_OBJS+=("$obj")
    fi
done

# Step 2: Compile C++ source files to .o
CPP_OBJS=()
for cppfile in "$SRC_DIR/sparsetools.cpp" "$SRC_DIR/spatial/build.cxx" "$SRC_DIR/spatial/query.cxx" "$SRC_DIR/spatial/query_ball_point.cxx" "$SRC_DIR/spatial/kdtree_wrapper.cpp" "$SRC_DIR/special/gamma.cpp"; do
    if [ -f "$cppfile" ]; then
        base=$(basename "$cppfile" | sed 's/\.[^.]*$//')
        obj="$OBJ_DIR/${base}.o"
        echo "Compiling C++: $cppfile"
        emcc -c "$cppfile" -o "$obj" -std=c++17 -I"$SRC_DIR" -I"$SRC_DIR/spatial" -I"$SRC_DIR/xsf" -O2
        CPP_OBJS+=("$obj")
    fi
done

ALL_OBJS=("${C_OBJS[@]}" "${CPP_OBJS[@]}")

echo ""

# All exported functions
EXPORTED_FUNCTIONS='[
    "_nelder_mead_minimize",
    "_bfgs_minimize",
    "_setulb",
    "_wasm_dqagse",
    "_wasm_dqagie",
    "_wasm_gamma",
    "_wasm_gammaln",
    "_wasm_rgamma",
    "_wasm_binom",
    "_wasm_binom_exact",
    "_wasm_poch",
    "_wasm_perm_exact",
    "_sp_csr_matvec_f64",
    "_sp_csr_matvecs_f64",
    "_sp_csr_tocsc_f64",
    "_sp_csr_todense_f64",
    "_sp_csr_diagonal_f64",
    "_sp_csr_sort_indices_f64",
    "_sp_csr_has_sorted_indices",
    "_sp_csr_has_canonical_format",
    "_sp_csr_sum_duplicates_f64",
    "_sp_csr_eliminate_zeros_f64",
    "_sp_csr_matmat_maxnnz",
    "_sp_csr_matmat_f64",
    "_sp_csr_plus_csr_f64",
    "_sp_csr_minus_csr_f64",
    "_sp_csr_elmul_csr_f64",
    "_sp_csr_eldiv_csr_f64",
    "_sp_csr_scale_rows_f64",
    "_sp_csr_scale_columns_f64",
    "_sp_expandptr",
    "_sp_csr_row_index_f64",
    "_sp_csr_row_slice_f64",
    "_sp_csr_column_index1",
    "_sp_csr_column_index2_f64",
    "_sp_csr_sample_offsets",
    "_sp_csr_sample_values_f64",
    "_sp_get_csr_submatrix_nnz",
    "_sp_get_csr_submatrix_f64",
    "_sp_csc_matvec_f64",
    "_sp_csc_matvecs_f64",
    "_sp_coo_tocsr_f64",
    "_sp_coo_todense_f64",
    "_sp_coo_matvec_f64",
    "_sp_malloc",
    "_sp_free",
    "_kdtree_build",
    "_kdtree_query_knn",
    "_kdtree_query_ball_point",
    "_kdtree_free",
    "_malloc",
    "_free"
]'

EXPORTED_RUNTIME='["ccall", "cwrap", "getValue", "setValue", "HEAPF64", "HEAPF32", "HEAP32", "HEAP8", "HEAPU8", "HEAPU32", "addFunction", "removeFunction"]'

# Link flags (shared between CJS and ESM)
LINK_FLAGS=(
    -s WASM=1
    -s MODULARIZE=1
    -s EXPORTED_FUNCTIONS="$EXPORTED_FUNCTIONS"
    -s EXPORTED_RUNTIME_METHODS="$EXPORTED_RUNTIME"
    -s ALLOW_MEMORY_GROWTH=1
    -s ALLOW_TABLE_GROWTH=1
    -s INITIAL_MEMORY=16777216
    -s STACK_SIZE=1048576
    -O2
)

# Build CJS version (Node.js)
echo "Linking CJS module for Node.js..."
emcc \
    "${ALL_OBJS[@]}" \
    "${LINK_FLAGS[@]}" \
    -s EXPORT_NAME="createSciWASMModule" \
    -o "$OUT_DIR/sciwasm.cjs"

echo ""
echo "CJS build complete!"
echo ""

# Build ESM version (Browser)
echo "Linking ESM module for browser support..."
emcc \
    "${ALL_OBJS[@]}" \
    "${LINK_FLAGS[@]}" \
    -s EXPORT_NAME="createSciWASMModule" \
    -s EXPORT_ES6=1 \
    -o "$OUT_DIR/sciwasm.mjs"

# Clean up object files
rm -rf "$OBJ_DIR"

echo ""
echo "Build complete!"
echo "  CJS Module:  $OUT_DIR/sciwasm.cjs (Node.js)"
echo "  ESM Module:  $OUT_DIR/sciwasm.mjs (Browser)"
echo "  WASM:        $OUT_DIR/sciwasm.wasm"

# Show file sizes
echo ""
echo "File sizes:"
ls -lh "$OUT_DIR/sciwasm.cjs" "$OUT_DIR/sciwasm.mjs" "$OUT_DIR/sciwasm.wasm" 2>/dev/null || true
