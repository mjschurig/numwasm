#!/bin/bash
#
# Build MFEM WebAssembly module
# Serial-only build with minimal dependencies
#
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
SRC_DIR="$PROJECT_DIR/src/wasm"
OUT_DIR="$PROJECT_DIR/dist/wasm"
OBJ_DIR="$PROJECT_DIR/.build-mfem"

echo "MFEM WASM Build"
echo "==============="
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

# Common compiler flags
COMMON_FLAGS=(
    -O2
    -DNDEBUG
    -I"$SRC_DIR"
    -I"$SRC_DIR/config"
    -I"$SRC_DIR/general"
    -I"$SRC_DIR/linalg"
    -I"$SRC_DIR/mesh"
    -I"$SRC_DIR/fem"
    -I"$SRC_DIR/include"
)

# C++ specific flags
# Note: MFEM uses dynamic_cast so we need RTTI enabled
CXX_FLAGS=(
    "${COMMON_FLAGS[@]}"
    -std=c++17
    -frtti
    -fexceptions
    -Wno-deprecated-declarations
    -Wno-unused-parameter
    -Wno-sign-compare
)

# C specific flags
C_FLAGS=(
    "${COMMON_FLAGS[@]}"
    -Wno-implicit-function-declaration
)

MFEM_OBJS=()

compile_cpp() {
    local src="$1"
    local module="$2"
    local base=$(basename "$src" .cpp)
    local obj="$OBJ_DIR/${module}_${base}.o"

    emcc -c "$src" -o "$obj" "${CXX_FLAGS[@]}"
    MFEM_OBJS+=("$obj")
}

compile_c() {
    local src="$1"
    local name="$2"
    local obj="$OBJ_DIR/${name}.o"

    emcc -c "$src" -o "$obj" "${C_FLAGS[@]}"
    MFEM_OBJS+=("$obj")
}

# ==============================================================================
# Compile BLAS (from arwasm)
# ==============================================================================
echo "Compiling BLAS..."
for src in "$SRC_DIR/BLAS"/*.c; do
    [ -f "$src" ] || continue
    base=$(basename "$src" .c)
    obj="$OBJ_DIR/blas_${base}.o"
    emcc -c "$src" -o "$obj" -O2 -I"$SRC_DIR/include" \
        -Wno-implicit-function-declaration \
        -Wno-incompatible-pointer-types \
        -Wno-parentheses
    MFEM_OBJS+=("$obj")
done
echo "  $(ls "$SRC_DIR/BLAS"/*.c | wc -l) BLAS files compiled"

# ==============================================================================
# Compile LAPACK (from arwasm)
# ==============================================================================
echo "Compiling LAPACK..."
for src in "$SRC_DIR/LAPACK"/*.c; do
    [ -f "$src" ] || continue
    base=$(basename "$src" .c)
    obj="$OBJ_DIR/lapack_${base}.o"
    emcc -c "$src" -o "$obj" -O2 -I"$SRC_DIR/include" \
        -Wno-implicit-function-declaration \
        -Wno-incompatible-pointer-types \
        -Wno-parentheses \
        -Wno-logical-op-parentheses \
        -Wno-sometimes-uninitialized
    MFEM_OBJS+=("$obj")
done
echo "  $(ls "$SRC_DIR/LAPACK"/*.c 2>/dev/null | wc -l) LAPACK files compiled"

# Compile f2c runtime
if [ -f "$SRC_DIR/f2c_runtime.c" ]; then
    echo "Compiling f2c runtime..."
    emcc -c "$SRC_DIR/f2c_runtime.c" -o "$OBJ_DIR/f2c_runtime.o" -O2 -I"$SRC_DIR/include"
    MFEM_OBJS+=("$OBJ_DIR/f2c_runtime.o")
fi

# ==============================================================================
# Compile general module
# ==============================================================================
echo "Compiling general module..."
for src in "$SRC_DIR/general"/*.cpp; do
    [ -f "$src" ] || continue
    base=$(basename "$src" .cpp)

    # Skip socket-related files (won't work in WASM)
    case "$base" in
        isockstream|osockstream|socketstream) continue ;;
        communication|adios2stream) continue ;;
    esac

    echo "  $base.cpp"
    compile_cpp "$src" "general"
done

# ==============================================================================
# Compile linalg module
# ==============================================================================
echo "Compiling linalg module..."
for src in "$SRC_DIR/linalg"/*.cpp; do
    [ -f "$src" ] || continue
    base=$(basename "$src" .cpp)
    echo "  $base.cpp"
    compile_cpp "$src" "linalg"
done

# Compile batched subdirectory
if [ -d "$SRC_DIR/linalg/batched" ]; then
    echo "  Compiling batched..."
    for src in "$SRC_DIR/linalg/batched"/*.cpp; do
        [ -f "$src" ] || continue
        base=$(basename "$src" .cpp)
        echo "    $base.cpp"
        compile_cpp "$src" "linalg_batched"
    done
fi

# ==============================================================================
# Compile mesh module
# ==============================================================================
echo "Compiling mesh module..."
for src in "$SRC_DIR/mesh"/*.cpp; do
    [ -f "$src" ] || continue
    base=$(basename "$src" .cpp)
    echo "  $base.cpp"
    compile_cpp "$src" "mesh"
done

# Compile submesh subdirectory
if [ -d "$SRC_DIR/mesh/submesh" ]; then
    echo "  Compiling submesh..."
    for src in "$SRC_DIR/mesh/submesh"/*.cpp; do
        [ -f "$src" ] || continue
        base=$(basename "$src" .cpp)
        echo "    $base.cpp"
        compile_cpp "$src" "mesh_submesh"
    done
fi

# ==============================================================================
# Compile fem module
# ==============================================================================
echo "Compiling fem module..."
for src in "$SRC_DIR/fem"/*.cpp; do
    [ -f "$src" ] || continue
    base=$(basename "$src" .cpp)
    echo "  $base.cpp"
    compile_cpp "$src" "fem"
done

# Compile fem subdirectories
for subdir in fe integ lor qinterp eltrans dfem; do
    if [ -d "$SRC_DIR/fem/$subdir" ]; then
        echo "  Compiling $subdir..."
        for src in "$SRC_DIR/fem/$subdir"/*.cpp; do
            [ -f "$src" ] || continue
            base=$(basename "$src" .cpp)
            echo "    $base.cpp"
            compile_cpp "$src" "fem_${subdir}"
        done
    fi
done

# ==============================================================================
# Compile WASM wrapper
# ==============================================================================
if [ -f "$SRC_DIR/mfem_wasm.cpp" ]; then
    echo "Compiling WASM wrapper..."
    compile_cpp "$SRC_DIR/mfem_wasm.cpp" "wasm"
fi

echo ""
echo "Total object files: ${#MFEM_OBJS[@]}"
echo ""

# ==============================================================================
# Link
# ==============================================================================

# Exported functions
EXPORTED_FUNCTIONS='[
    "_mfem_mesh_create",
    "_mfem_mesh_destroy",
    "_mfem_mesh_load_from_string",
    "_mfem_mesh_get_dimension",
    "_mfem_mesh_get_num_elements",
    "_mfem_mesh_get_num_vertices",
    "_mfem_mesh_get_num_boundary_elements",
    "_mfem_mesh_refine_uniform",
    "_mfem_mesh_get_bounding_box",
    "_mfem_fespace_create_h1",
    "_mfem_fespace_destroy",
    "_mfem_fespace_get_ndofs",
    "_mfem_fespace_get_vdim",
    "_mfem_gridfunc_create",
    "_mfem_gridfunc_destroy",
    "_mfem_gridfunc_get_data",
    "_mfem_gridfunc_get_size",
    "_mfem_gridfunc_project_coefficient",
    "_malloc",
    "_free"
]'

EXPORTED_RUNTIME='["ccall", "cwrap", "getValue", "setValue", "HEAPF64", "HEAPF32", "HEAP32", "HEAPU8", "UTF8ToString", "stringToUTF8", "lengthBytesUTF8"]'

LINK_FLAGS=(
    -s WASM=1
    -s MODULARIZE=1
    -s EXPORTED_FUNCTIONS="$EXPORTED_FUNCTIONS"
    -s EXPORTED_RUNTIME_METHODS="$EXPORTED_RUNTIME"
    -s ALLOW_MEMORY_GROWTH=1
    -s INITIAL_MEMORY=33554432
    -s STACK_SIZE=2097152
    -s DISABLE_EXCEPTION_CATCHING=0
    -O2
    -fexceptions
)

# Build CJS version
echo "Linking CJS module..."
emcc \
    "${MFEM_OBJS[@]}" \
    "${LINK_FLAGS[@]}" \
    -s EXPORT_NAME="createMFEMModule" \
    -o "$OUT_DIR/mfem.cjs"

# Build ESM version
echo "Linking ESM module..."
emcc \
    "${MFEM_OBJS[@]}" \
    "${LINK_FLAGS[@]}" \
    -s EXPORT_NAME="createMFEMModule" \
    -s EXPORT_ES6=1 \
    -o "$OUT_DIR/mfem.mjs"

# Clean up object files
rm -rf "$OBJ_DIR"

echo ""
echo "Build complete!"
echo "  CJS: $OUT_DIR/mfem.cjs"
echo "  ESM: $OUT_DIR/mfem.mjs"
echo "  WASM: $OUT_DIR/mfem.wasm"
echo ""

if [ -f "$OUT_DIR/mfem.wasm" ]; then
    echo "File sizes:"
    ls -lh "$OUT_DIR/mfem.cjs" "$OUT_DIR/mfem.mjs" "$OUT_DIR/mfem.wasm"
fi
