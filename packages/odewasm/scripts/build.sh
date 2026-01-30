#!/bin/bash
#
# Build ODE WASM module (DOPRI5, DOP853, RADAU5)
# Compiles f2c-converted Fortran ODE solvers to WebAssembly
#
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
SRC_DIR="$PROJECT_DIR/src/wasm"
OUT_DIR="$PROJECT_DIR/dist/wasm"
OBJ_DIR="$PROJECT_DIR/.build-ode"

echo "ODE WASM Build"
echo "=============="
echo "Source: $SRC_DIR"
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

# Check for source files
if [ ! -f "$SRC_DIR/dopri5.c" ]; then
    echo "Error: ODE source not found at $SRC_DIR"
    echo "Expected files: dopri5.c, dop853.c, radau5.c, decsol.c"
    exit 1
fi

ODE_OBJS=()

echo "Compiling ODE solvers..."

# Common compilation flags
CFLAGS=(
    -O2
    -I"$SRC_DIR/include"
    -Wno-implicit-function-declaration
    -Wno-incompatible-pointer-types
    -Wno-parentheses
    -Wno-logical-op-parentheses
    -Wno-sometimes-uninitialized
    -Wno-unused-variable
    -Wno-unused-but-set-variable
)

# Compile f2c runtime
echo "  Compiling f2c_runtime.c..."
emcc -c "$SRC_DIR/f2c_runtime.c" -o "$OBJ_DIR/f2c_runtime.o" "${CFLAGS[@]}"
ODE_OBJS+=("$OBJ_DIR/f2c_runtime.o")

# Compile DOPRI5
echo "  Compiling dopri5.c..."
emcc -c "$SRC_DIR/dopri5.c" -o "$OBJ_DIR/dopri5.o" "${CFLAGS[@]}"
ODE_OBJS+=("$OBJ_DIR/dopri5.o")

# Compile DOP853
echo "  Compiling dop853.c..."
emcc -c "$SRC_DIR/dop853.c" -o "$OBJ_DIR/dop853.o" "${CFLAGS[@]}"
ODE_OBJS+=("$OBJ_DIR/dop853.o")

# Compile RADAU5
echo "  Compiling radau5.c..."
emcc -c "$SRC_DIR/radau5.c" -o "$OBJ_DIR/radau5.o" "${CFLAGS[@]}"
ODE_OBJS+=("$OBJ_DIR/radau5.o")

# Compile DC_DECSOL (decomposition/solve routines for RADAU5)
echo "  Compiling dc_decsol.c..."
emcc -c "$SRC_DIR/dc_decsol.c" -o "$OBJ_DIR/dc_decsol.o" "${CFLAGS[@]}"
ODE_OBJS+=("$OBJ_DIR/dc_decsol.o")

# Compile DECSOL (dense linear solver for RADAU5)
echo "  Compiling decsol.c..."
emcc -c "$SRC_DIR/decsol.c" -o "$OBJ_DIR/decsol.o" "${CFLAGS[@]}"
ODE_OBJS+=("$OBJ_DIR/decsol.o")

# Compile WASM wrapper
echo "  Compiling ode_wasm.c..."
emcc -c "$SRC_DIR/ode_wasm.c" -o "$OBJ_DIR/ode_wasm.o" "${CFLAGS[@]}"
ODE_OBJS+=("$OBJ_DIR/ode_wasm.o")

echo "  Total object files: ${#ODE_OBJS[@]}"
echo ""

# Exported functions
EXPORTED_FUNCTIONS='[
    "_wasm_dopri5", "_wasm_dop853", "_wasm_radau5",
    "_wasm_contd5", "_wasm_contd8", "_wasm_contr5",
    "_wasm_set_fcn_callback", "_wasm_set_solout_callback", "_wasm_set_jac_callback",
    "_wasm_dopri5_work_size", "_wasm_dopri5_iwork_size",
    "_wasm_dop853_work_size", "_wasm_dop853_iwork_size",
    "_wasm_radau5_work_size", "_wasm_radau5_iwork_size",
    "_wasm_malloc", "_wasm_free",
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
    "${ODE_OBJS[@]}" \
    "${LINK_FLAGS[@]}" \
    -s EXPORT_NAME="createODEModule" \
    -o "$OUT_DIR/ode.cjs"

# Build ESM version
echo "Linking ESM module..."
emcc \
    "${ODE_OBJS[@]}" \
    "${LINK_FLAGS[@]}" \
    -s EXPORT_NAME="createODEModule" \
    -s EXPORT_ES6=1 \
    -o "$OUT_DIR/ode.mjs"

# Clean up
rm -rf "$OBJ_DIR"

echo ""
echo "Build complete!"
echo "  CJS: $OUT_DIR/ode.cjs"
echo "  ESM: $OUT_DIR/ode.mjs"
echo "  WASM: $OUT_DIR/ode.wasm"
echo ""
echo "File sizes:"
ls -lh "$OUT_DIR/ode.cjs" "$OUT_DIR/ode.mjs" "$OUT_DIR/ode.wasm" 2>/dev/null || true
