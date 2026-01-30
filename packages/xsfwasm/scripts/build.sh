#!/bin/bash
#
# Build XSF (Special Functions) WebAssembly module (separate from sciwasm)
# This creates a standalone xsf.wasm for special mathematical functions
#
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
SRC_DIR="$PROJECT_DIR/src/wasm"
OUT_DIR="$PROJECT_DIR/dist/wasm"
OBJ_DIR="$PROJECT_DIR/.build-xsf"

XSF_DIR="$SRC_DIR"
SPECIAL_DIR="$SRC_DIR/special"

echo "XSF WASM Build"
echo "=============="
echo "Source: $XSF_DIR"
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

XSF_OBJS=()

# Compile C source files
echo "Compiling C sources..."
for cfile in "$SPECIAL_DIR/comb.c"; do
    if [ -f "$cfile" ]; then
        base=$(basename "$cfile" .c)
        obj="$OBJ_DIR/${base}.o"
        echo "  Compiling: $(basename "$cfile")"
        emcc -c "$cfile" -o "$obj" -O2
        XSF_OBJS+=("$obj")
    fi
done

# Compile C++ source files
# NOTE: All special functions are consolidated in gamma.cpp to avoid
# duplicate symbol errors from header-only xsf library
echo ""
echo "Compiling C++ sources..."
for cppfile in "$SPECIAL_DIR/gamma.cpp"; do
    if [ -f "$cppfile" ]; then
        base=$(basename "$cppfile" .cpp)
        obj="$OBJ_DIR/${base}.o"
        echo "  Compiling: $(basename "$cppfile")"
        emcc -c "$cppfile" -o "$obj" -std=c++17 -I"$XSF_DIR" -O2
        XSF_OBJS+=("$obj")
    fi
done

echo ""
echo "  Total object files: ${#XSF_OBJS[@]}"
echo ""

# XSF exported functions
EXPORTED_FUNCTIONS='[
    "_wasm_gamma",
    "_wasm_gammaln",
    "_wasm_rgamma",
    "_wasm_binom",
    "_wasm_binom_exact",
    "_wasm_poch",
    "_wasm_perm_exact",
    "_wasm_beta",
    "_wasm_betaln",
    "_wasm_erf",
    "_wasm_erfc",
    "_wasm_erfcx",
    "_wasm_erfi",
    "_wasm_j0",
    "_wasm_j1",
    "_wasm_jv",
    "_wasm_y0",
    "_wasm_y1",
    "_wasm_yv",
    "_wasm_i0",
    "_wasm_i1",
    "_wasm_iv",
    "_wasm_k0",
    "_wasm_k1",
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
    "${XSF_OBJS[@]}" \
    "${LINK_FLAGS[@]}" \
    -s EXPORT_NAME="createXSFModule" \
    -o "$OUT_DIR/xsf.cjs"

# Build ESM version
echo "Linking ESM module..."
emcc \
    "${XSF_OBJS[@]}" \
    "${LINK_FLAGS[@]}" \
    -s EXPORT_NAME="createXSFModule" \
    -s EXPORT_ES6=1 \
    -o "$OUT_DIR/xsf.mjs"

# Clean up
rm -rf "$OBJ_DIR"

echo ""
echo "Build complete!"
echo "  CJS: $OUT_DIR/xsf.cjs"
echo "  ESM: $OUT_DIR/xsf.mjs"
echo "  WASM: $OUT_DIR/xsf.wasm"
echo ""
echo "File sizes:"
ls -lh "$OUT_DIR/xsf.cjs" "$OUT_DIR/xsf.mjs" "$OUT_DIR/xsf.wasm" 2>/dev/null || true
