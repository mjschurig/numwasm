#!/bin/bash
#
# Build LINPACK WebAssembly module (separate from sciwasm)
# This creates a standalone linpack.wasm for linear algebra operations
#
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
SRC_DIR="$PROJECT_DIR/src/wasm"
OUT_DIR="$PROJECT_DIR/dist/wasm"
OBJ_DIR="$PROJECT_DIR/.build-linpack"

echo "LINPACK WASM Build"
echo "=================="
echo "Source: $SRC_DIR/linpack"
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

LINPACK_DIR="$SRC_DIR/linpack"

if [ ! -d "$LINPACK_DIR" ]; then
    echo "Error: LINPACK source not found at $LINPACK_DIR"
    echo "Run scripts/convert-linpack.sh first"
    exit 1
fi

LINPACK_OBJS=()

echo "Compiling LINPACK..."

# Compile all LINPACK/BLAS C files (excluding f2c_runtime.c which we handle separately)
for cfile in "$LINPACK_DIR"/*.c; do
    if [ -f "$cfile" ]; then
        base=$(basename "$cfile" .c)
        # Skip f2c_runtime.c - compile it separately
        if [ "$base" = "f2c_runtime" ]; then
            continue
        fi
        obj="$OBJ_DIR/linpack_${base}.o"
        emcc -c "$cfile" -o "$obj" -O2 \
            -I"$LINPACK_DIR/include" \
            -Wno-implicit-function-declaration \
            -Wno-incompatible-pointer-types \
            -Wno-parentheses \
            -Wno-logical-op-parentheses
        LINPACK_OBJS+=("$obj")
    fi
done

# Compile f2c runtime
if [ -f "$LINPACK_DIR/f2c_runtime.c" ]; then
    echo "  Compiling f2c runtime..."
    emcc -c "$LINPACK_DIR/f2c_runtime.c" -o "$OBJ_DIR/f2c_runtime.o" -O2 \
        -I"$LINPACK_DIR/include"
    LINPACK_OBJS+=("$OBJ_DIR/f2c_runtime.o")
fi

echo "  Total object files: ${#LINPACK_OBJS[@]}"
echo ""

# LINPACK exported functions (all have trailing underscore from f2c)
# All 176 LINPACK routines + BLAS routines
EXPORTED_FUNCTIONS='[
    "_dgefa_", "_dgesl_", "_dgeco_", "_dgedi_",
    "_dgbfa_", "_dgbsl_", "_dgbco_", "_dgbdi_",
    "_dpofa_", "_dposl_", "_dpoco_", "_dpodi_",
    "_dppfa_", "_dppsl_", "_dppco_", "_dppdi_",
    "_dpbfa_", "_dpbsl_", "_dpbco_", "_dpbdi_",
    "_dsifa_", "_dsisl_", "_dsico_", "_dsidi_",
    "_dspfa_", "_dspsl_", "_dspco_", "_dspdi_",
    "_dtrco_", "_dtrdi_", "_dtrsl_",
    "_dgtsl_", "_dptsl_",
    "_dqrdc_", "_dqrsl_",
    "_dsvdc_",
    "_dchdc_", "_dchdd_", "_dchex_", "_dchud_",
    "_sgefa_", "_sgesl_", "_sgeco_", "_sgedi_",
    "_sgbfa_", "_sgbsl_", "_sgbco_", "_sgbdi_",
    "_spofa_", "_sposl_", "_spoco_", "_spodi_",
    "_sppfa_", "_sppsl_", "_sppco_", "_sppdi_",
    "_spbfa_", "_spbsl_", "_spbco_", "_spbdi_",
    "_ssifa_", "_ssisl_", "_ssico_", "_ssidi_",
    "_sspfa_", "_sspsl_", "_sspco_", "_sspdi_",
    "_strco_", "_strdi_", "_strsl_",
    "_sgtsl_", "_sptsl_",
    "_sqrdc_", "_sqrsl_",
    "_ssvdc_",
    "_schdc_", "_schdd_", "_schex_", "_schud_",
    "_cgefa_", "_cgesl_", "_cgeco_", "_cgedi_",
    "_cgbfa_", "_cgbsl_", "_cgbco_", "_cgbdi_",
    "_cpofa_", "_cposl_", "_cpoco_", "_cpodi_",
    "_cppfa_", "_cppsl_", "_cppco_", "_cppdi_",
    "_cpbfa_", "_cpbsl_", "_cpbco_", "_cpbdi_",
    "_chifa_", "_chisl_", "_chico_", "_chidi_",
    "_chpfa_", "_chpsl_", "_chpco_", "_chpdi_",
    "_csifa_", "_csisl_", "_csico_", "_csidi_",
    "_cspfa_", "_cspsl_", "_cspco_", "_cspdi_",
    "_ctrco_", "_ctrdi_", "_ctrsl_",
    "_cgtsl_", "_cptsl_",
    "_cqrdc_", "_cqrsl_",
    "_csvdc_",
    "_cchdc_", "_cchdd_", "_cchex_", "_cchud_",
    "_zgefa_", "_zgesl_", "_zgeco_", "_zgedi_",
    "_zgbfa_", "_zgbsl_", "_zgbco_", "_zgbdi_",
    "_zpofa_", "_zposl_", "_zpoco_", "_zpodi_",
    "_zppfa_", "_zppsl_", "_zppco_", "_zppdi_",
    "_zpbfa_", "_zpbsl_", "_zpbco_", "_zpbdi_",
    "_zhifa_", "_zhisl_", "_zhico_", "_zhidi_",
    "_zhpfa_", "_zhpsl_", "_zhpco_", "_zhpdi_",
    "_zsifa_", "_zsisl_", "_zsico_", "_zsidi_",
    "_zspfa_", "_zspsl_", "_zspco_", "_zspdi_",
    "_ztrco_", "_ztrdi_", "_ztrsl_",
    "_zgtsl_", "_zptsl_",
    "_zqrdc_", "_zqrsl_",
    "_zsvdc_",
    "_zchdc_", "_zchdd_", "_zchex_", "_zchud_",
    "_daxpy_", "_dcopy_", "_dscal_", "_dswap_", "_ddot_", "_dasum_", "_drot_", "_dnrm2_", "_drotg_",
    "_saxpy_", "_scopy_", "_sscal_", "_sswap_", "_sdot_", "_sasum_", "_srot_", "_snrm2_", "_srotg_",
    "_caxpy_", "_ccopy_", "_cscal_", "_cswap_", "_cdotc_", "_cdotu_", "_csrot_", "_csscal_", "_scnrm2_", "_crotg_",
    "_zaxpy_", "_zcopy_", "_zscal_", "_zswap_", "_zdotc_", "_zdotu_", "_zdrot_", "_zdscal_", "_dznrm2_", "_zrotg_",
    "_idamax_", "_isamax_", "_icamax_", "_izamax_",
    "_dzasum_", "_scasum_", "_scabs1_", "_dcabs1_",
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
    "${LINPACK_OBJS[@]}" \
    "${LINK_FLAGS[@]}" \
    -s EXPORT_NAME="createLINPACKModule" \
    -o "$OUT_DIR/linpack.cjs"

# Build ESM version
echo "Linking ESM module..."
emcc \
    "${LINPACK_OBJS[@]}" \
    "${LINK_FLAGS[@]}" \
    -s EXPORT_NAME="createLINPACKModule" \
    -s EXPORT_ES6=1 \
    -o "$OUT_DIR/linpack.mjs"

# Clean up
rm -rf "$OBJ_DIR"

echo ""
echo "Build complete!"
echo "  CJS: $OUT_DIR/linpack.cjs"
echo "  ESM: $OUT_DIR/linpack.mjs"
echo "  WASM: $OUT_DIR/linpack.wasm"
echo ""
echo "File sizes:"
ls -lh "$OUT_DIR/linpack.cjs" "$OUT_DIR/linpack.mjs" "$OUT_DIR/linpack.wasm" 2>/dev/null || true
