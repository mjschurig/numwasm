#!/bin/bash
#
# Build QUADPACK WebAssembly module (separate from sciwasm)
# This creates a standalone quadpack.wasm for numerical integration
#
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
SRC_DIR="$PROJECT_DIR/src/wasm"
OUT_DIR="$PROJECT_DIR/dist/wasm"
OBJ_DIR="$PROJECT_DIR/.build-quadpack"

echo "QUADPACK WASM Build"
echo "==================="
echo "Source: $SRC_DIR/quadpack"
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

QUADPACK_DIR="$SRC_DIR/quadpack"

if [ ! -d "$QUADPACK_DIR" ]; then
    echo "Error: QUADPACK source not found at $QUADPACK_DIR"
    echo "Run scripts/convert-quadpack.sh first"
    exit 1
fi

QUADPACK_OBJS=()

echo "Compiling QUADPACK..."

# Compile all QUADPACK C files (excluding f2c_runtime.c which we handle separately)
for cfile in "$QUADPACK_DIR"/*.c; do
    if [ -f "$cfile" ]; then
        base=$(basename "$cfile" .c)
        # Skip f2c_runtime.c - compile it separately
        if [ "$base" = "f2c_runtime" ]; then
            continue
        fi
        obj="$OBJ_DIR/quadpack_${base}.o"
        emcc -c "$cfile" -o "$obj" -O2 \
            -I"$QUADPACK_DIR/include" \
            -Wno-implicit-function-declaration \
            -Wno-incompatible-pointer-types \
            -Wno-parentheses \
            -Wno-logical-op-parentheses
        QUADPACK_OBJS+=("$obj")
    fi
done

# Compile f2c runtime
if [ -f "$QUADPACK_DIR/f2c_runtime.c" ]; then
    echo "  Compiling f2c runtime..."
    emcc -c "$QUADPACK_DIR/f2c_runtime.c" -o "$OBJ_DIR/f2c_runtime.o" -O2 \
        -I"$QUADPACK_DIR/include"
    QUADPACK_OBJS+=("$OBJ_DIR/f2c_runtime.o")
fi

echo "  Total object files: ${#QUADPACK_OBJS[@]}"
echo ""

# QUADPACK exported functions (double precision primary, single precision secondary)
# All function names have trailing underscore from f2c
EXPORTED_FUNCTIONS='[
    "_dqagse_", "_dqagie_",
    "_dqags_", "_dqagi_",
    "_dqag_", "_dqage_",
    "_dqng_",
    "_dqawoe_", "_dqawfe_", "_dqawo_", "_dqawf_",
    "_dqawse_", "_dqawce_", "_dqaws_", "_dqawc_",
    "_dqagpe_", "_dqagp_",
    "_dqk15_", "_dqk21_", "_dqk31_", "_dqk41_", "_dqk51_", "_dqk61_",
    "_dqk15i_", "_dqk15w_",
    "_qagse_", "_qagie_",
    "_qags_", "_qagi_",
    "_qag_", "_qage_",
    "_qng_",
    "_qawoe_", "_qawfe_", "_qawo_", "_qawf_",
    "_qawse_", "_qawce_", "_qaws_", "_qawc_",
    "_qagpe_", "_qagp_",
    "_qk15_", "_qk21_", "_qk31_", "_qk41_", "_qk51_", "_qk61_",
    "_qk15i_", "_qk15w_",
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
    "${QUADPACK_OBJS[@]}" \
    "${LINK_FLAGS[@]}" \
    -s EXPORT_NAME="createQUADPACKModule" \
    -o "$OUT_DIR/quadpack.cjs"

# Build ESM version
echo "Linking ESM module..."
emcc \
    "${QUADPACK_OBJS[@]}" \
    "${LINK_FLAGS[@]}" \
    -s EXPORT_NAME="createQUADPACKModule" \
    -s EXPORT_ES6=1 \
    -o "$OUT_DIR/quadpack.mjs"

# Clean up
rm -rf "$OBJ_DIR"

echo ""
echo "Build complete!"
echo "  CJS: $OUT_DIR/quadpack.cjs"
echo "  ESM: $OUT_DIR/quadpack.mjs"
echo "  WASM: $OUT_DIR/quadpack.wasm"
echo ""
echo "File sizes:"
ls -lh "$OUT_DIR/quadpack.cjs" "$OUT_DIR/quadpack.mjs" "$OUT_DIR/quadpack.wasm" 2>/dev/null || true
