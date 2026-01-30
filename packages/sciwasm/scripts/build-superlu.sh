#!/bin/bash
#
# Build SuperLU WebAssembly module (separate from sciwasm)
# This creates a standalone superlu.wasm for sparse direct solvers
#
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
SRC_DIR="$PROJECT_DIR/src/wasm"
OUT_DIR="$PROJECT_DIR/dist/wasm"
OBJ_DIR="$PROJECT_DIR/.build-superlu"

echo "SuperLU WASM Build"
echo "=================="
echo "Source: $SRC_DIR/superlu"
echo "Output: $OUT_DIR"
echo ""

# Create output/temp directories
mkdir -p "$OUT_DIR" "$OBJ_DIR"

# Check for Emscripten
if ! command -v emcc &> /dev/null; then
    echo "Error: emcc (Emscripten) not found in PATH"
    exit 1
fi

echo "Using Emscripten: $(emcc --version | head -1)"
echo ""

SUPERLU_DIR="$SRC_DIR/superlu"
SUPERLU_SRC="$SUPERLU_DIR/SRC"
SUPERLU_CBLAS="$SUPERLU_DIR/CBLAS"

if [ ! -d "$SUPERLU_SRC" ]; then
    echo "Error: SuperLU source not found at $SUPERLU_SRC"
    exit 1
fi

SUPERLU_OBJS=()

# Compile all SuperLU SRC files
echo "Compiling SuperLU SRC..."
for cfile in "$SUPERLU_SRC"/*.c; do
    if [ -f "$cfile" ]; then
        base=$(basename "$cfile" .c)
        obj="$OBJ_DIR/superlu_${base}.o"
        emcc -c "$cfile" -o "$obj" -O2 \
            -I"$SUPERLU_SRC" \
            -I"$SUPERLU_DIR" \
            -Wno-implicit-function-declaration \
            -Wno-incompatible-pointer-types \
            -Wno-parentheses
        SUPERLU_OBJS+=("$obj")
    fi
done

# Compile all CBLAS files (SuperLU has its own complete BLAS)
echo ""
echo "Compiling CBLAS..."
for cfile in "$SUPERLU_CBLAS"/*.c; do
    if [ -f "$cfile" ]; then
        base=$(basename "$cfile" .c)
        obj="$OBJ_DIR/cblas_${base}.o"
        emcc -c "$cfile" -o "$obj" -O2 \
            -I"$SUPERLU_SRC" \
            -I"$SUPERLU_DIR" \
            -Wno-implicit-function-declaration \
            -Wno-incompatible-pointer-types \
            -Wno-parentheses
        SUPERLU_OBJS+=("$obj")
    fi
done

echo ""
echo "  Total object files: ${#SUPERLU_OBJS[@]}"
echo ""

# SuperLU exported functions - all precision variants
# s = single, d = double, c = complex, z = double complex
EXPORTED_FUNCTIONS='[
    "_sgssv", "_sgstrf", "_sgstrs", "_sgssvx", "_sgsisx", "_sgsitrf",
    "_sgscon", "_sgsequ", "_sgsrfs", "_slaqgs",
    "_sCreate_CompCol_Matrix", "_sCreate_CompRow_Matrix", "_sCreate_Dense_Matrix",
    "_sCopy_CompCol_Matrix", "_sCopy_Dense_Matrix", "_sCompRow_to_CompCol",
    "_sallocateA", "_sPrint_CompCol_Matrix", "_sPrint_Dense_Matrix",
    "_sGenXtrue", "_sFillRHS", "_sinf_norm_error",
    "_sQuerySpace", "_smemory_usage",
    "_dgssv", "_dgstrf", "_dgstrs", "_dgssvx", "_dgsisx", "_dgsitrf",
    "_dgscon", "_dgsequ", "_dgsrfs", "_dlaqgs",
    "_dCreate_CompCol_Matrix", "_dCreate_CompRow_Matrix", "_dCreate_Dense_Matrix",
    "_dCopy_CompCol_Matrix", "_dCopy_Dense_Matrix", "_dCompRow_to_CompCol",
    "_dallocateA", "_dPrint_CompCol_Matrix", "_dPrint_Dense_Matrix",
    "_dGenXtrue", "_dFillRHS", "_dinf_norm_error",
    "_dQuerySpace", "_dmemory_usage",
    "_cgssv", "_cgstrf", "_cgstrs", "_cgssvx", "_cgsisx", "_cgsitrf",
    "_cgscon", "_cgsequ", "_cgsrfs", "_claqgs",
    "_cCreate_CompCol_Matrix", "_cCreate_CompRow_Matrix", "_cCreate_Dense_Matrix",
    "_cCopy_CompCol_Matrix", "_cCopy_Dense_Matrix", "_cCompRow_to_CompCol",
    "_callocateA", "_cPrint_CompCol_Matrix", "_cPrint_Dense_Matrix",
    "_cGenXtrue", "_cFillRHS", "_cinf_norm_error",
    "_cQuerySpace", "_cmemory_usage",
    "_zgssv", "_zgstrf", "_zgstrs", "_zgssvx", "_zgsisx", "_zgsitrf",
    "_zgscon", "_zgsequ", "_zgsrfs", "_zlaqgs",
    "_zCreate_CompCol_Matrix", "_zCreate_CompRow_Matrix", "_zCreate_Dense_Matrix",
    "_zCopy_CompCol_Matrix", "_zCopy_Dense_Matrix", "_zCompRow_to_CompCol",
    "_zallocateA", "_zPrint_CompCol_Matrix", "_zPrint_Dense_Matrix",
    "_zGenXtrue", "_zFillRHS", "_zinf_norm_error",
    "_zQuerySpace", "_zmemory_usage",
    "_set_default_options", "_StatInit", "_StatFree",
    "_Destroy_CompCol_Matrix", "_Destroy_CompRow_Matrix", "_Destroy_Dense_Matrix",
    "_Destroy_SuperNode_Matrix", "_Destroy_SuperMatrix_Store",
    "_get_perm_c", "_sp_preorder", "_sp_ienv",
    "_input_error", "_superlu_abort_and_exit",
    "_intMalloc", "_intCalloc",
    "_floatMalloc", "_floatCalloc",
    "_doubleMalloc", "_doubleCalloc",
    "_singlecomplexMalloc", "_singlecomplexCalloc",
    "_doublecomplexMalloc", "_doublecomplexCalloc",
    "_superlu_malloc", "_superlu_free",
    "_malloc", "_free"
]'

EXPORTED_RUNTIME='["ccall", "cwrap", "getValue", "setValue", "HEAPF64", "HEAPF32", "HEAP32", "HEAPU8", "addFunction", "removeFunction"]'

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

# Build CJS version
echo "Linking CJS module..."
emcc \
    "${SUPERLU_OBJS[@]}" \
    "${LINK_FLAGS[@]}" \
    -s EXPORT_NAME="createSuperLUModule" \
    -o "$OUT_DIR/superlu.cjs"

# Build ESM version
echo "Linking ESM module..."
emcc \
    "${SUPERLU_OBJS[@]}" \
    "${LINK_FLAGS[@]}" \
    -s EXPORT_NAME="createSuperLUModule" \
    -s EXPORT_ES6=1 \
    -o "$OUT_DIR/superlu.mjs"

# Clean up
rm -rf "$OBJ_DIR"

echo ""
echo "Build complete!"
echo "  CJS: $OUT_DIR/superlu.cjs"
echo "  ESM: $OUT_DIR/superlu.mjs"
echo "  WASM: $OUT_DIR/superlu.wasm"
echo ""
echo "File sizes:"
ls -lh "$OUT_DIR/superlu.cjs" "$OUT_DIR/superlu.mjs" "$OUT_DIR/superlu.wasm" 2>/dev/null || true
