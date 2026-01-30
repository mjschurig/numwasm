#!/bin/bash
#
# Build ARPACK WebAssembly module (separate from sciwasm)
# This creates a standalone arpack.wasm with its own bundled BLAS/LAPACK
#
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
SRC_DIR="$PROJECT_DIR/src/wasm"
OUT_DIR="$PROJECT_DIR/dist/wasm"
OBJ_DIR="$PROJECT_DIR/.build-arpack"

echo "ARPACK WASM Build"
echo "================="
echo "Source: $SRC_DIR/arpack"
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

ARPACK_DIR="$SRC_DIR"

if [ ! -d "$ARPACK_DIR/SRC" ]; then
    echo "Error: ARPACK source not found at $ARPACK_DIR"
    echo "Expected directories: SRC/, BLAS/, LAPACK/, UTIL/"
    exit 1
fi

ARPACK_OBJS=()

# UTIL functions that are provided by arpack_stubs.c (timing and output)
ARPACK_UTIL_SKIP="second second_NONE ivout dvout svout cvout zvout dmout smout cmout zmout"

# LAPACK functions to skip (duplicates or provided by f2c_runtime)
ARPACK_LAPACK_SKIP="zdscal xerbla"

echo "Compiling ARPACK..."

# Compile in dependency order: BLAS -> LAPACK -> UTIL -> SRC
for subdir in BLAS LAPACK UTIL SRC; do
    if [ -d "$ARPACK_DIR/$subdir" ]; then
        echo "  Compiling $subdir..."
        for cfile in "$ARPACK_DIR/$subdir"/*.c; do
            if [ -f "$cfile" ]; then
                base=$(basename "$cfile" .c)

                # Skip files that are provided by stubs
                skip=0
                case "$subdir" in
                    LAPACK)
                        for s in $ARPACK_LAPACK_SKIP; do
                            if [ "$base" = "$s" ]; then
                                skip=1
                                break
                            fi
                        done
                        ;;
                    UTIL)
                        for s in $ARPACK_UTIL_SKIP; do
                            if [ "$base" = "$s" ]; then
                                skip=1
                                break
                            fi
                        done
                        ;;
                esac
                if [ $skip -eq 1 ]; then
                    continue
                fi

                obj="$OBJ_DIR/arpack_${subdir}_${base}.o"
                emcc -c "$cfile" -o "$obj" -O2 \
                    -I"$ARPACK_DIR/include" \
                    -Wno-implicit-function-declaration \
                    -Wno-incompatible-pointer-types \
                    -Wno-parentheses \
                    -Wno-logical-op-parentheses
                ARPACK_OBJS+=("$obj")
            fi
        done
    fi
done

# Compile stubs
if [ -f "$ARPACK_DIR/arpack_stubs.c" ]; then
    echo "  Compiling stubs..."
    emcc -c "$ARPACK_DIR/arpack_stubs.c" -o "$OBJ_DIR/arpack_stubs.o" -O2 \
        -I"$ARPACK_DIR/include"
    ARPACK_OBJS+=("$OBJ_DIR/arpack_stubs.o")
fi

# Compile f2c runtime
if [ -f "$ARPACK_DIR/f2c_runtime.c" ]; then
    echo "  Compiling f2c runtime..."
    emcc -c "$ARPACK_DIR/f2c_runtime.c" -o "$OBJ_DIR/f2c_runtime.o" -O2 \
        -I"$ARPACK_DIR/include"
    ARPACK_OBJS+=("$OBJ_DIR/f2c_runtime.o")
fi

echo "  Total object files: ${#ARPACK_OBJS[@]}"
echo ""

# ARPACK exported functions (all precision variants)
EXPORTED_FUNCTIONS='[
    "_dsaupd_", "_dseupd_",
    "_dnaupd_", "_dneupd_",
    "_ssaupd_", "_sseupd_",
    "_snaupd_", "_sneupd_",
    "_cnaupd_", "_cneupd_",
    "_znaupd_", "_zneupd_",
    "_malloc", "_free"
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
    "${ARPACK_OBJS[@]}" \
    "${LINK_FLAGS[@]}" \
    -s EXPORT_NAME="createARPACKModule" \
    -o "$OUT_DIR/arpack.cjs"

# Build ESM version
echo "Linking ESM module..."
emcc \
    "${ARPACK_OBJS[@]}" \
    "${LINK_FLAGS[@]}" \
    -s EXPORT_NAME="createARPACKModule" \
    -s EXPORT_ES6=1 \
    -o "$OUT_DIR/arpack.mjs"

# Clean up
rm -rf "$OBJ_DIR"

echo ""
echo "Build complete!"
echo "  CJS: $OUT_DIR/arpack.cjs"
echo "  ESM: $OUT_DIR/arpack.mjs"
echo "  WASM: $OUT_DIR/arpack.wasm"
echo ""
echo "File sizes:"
ls -lh "$OUT_DIR/arpack.cjs" "$OUT_DIR/arpack.mjs" "$OUT_DIR/arpack.wasm" 2>/dev/null || true
