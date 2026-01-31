#!/bin/bash
#
# Build ODE WASM module
# Compiles f2c-converted Fortran ODE solvers to WebAssembly
#
# Included solvers:
#   Hairer: DOPRI5, DOP853, RADAU5
#   Netlib: RKF45, DVERK, ODE, VODE, ZVODE, EPSODE, VODPK, RKSUITE, RKC
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
    echo "Expected files: dopri5.c, dop853.c, radau5.c, decsol.c, rkf45.c, etc."
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

# ============================================
# BLAS Level 1 routines
# ============================================

# Real BLAS
for blas_file in daxpy dcopy ddot dnrm2 dscal idamax; do
    if [ -f "$SRC_DIR/${blas_file}.c" ]; then
        echo "  Compiling ${blas_file}.c..."
        emcc -c "$SRC_DIR/${blas_file}.c" -o "$OBJ_DIR/${blas_file}.o" "${CFLAGS[@]}"
        ODE_OBJS+=("$OBJ_DIR/${blas_file}.o")
    fi
done

# Complex BLAS
for blas_file in zaxpy zcopy zdotc zdscal izamax zscal dcabs1; do
    if [ -f "$SRC_DIR/${blas_file}.c" ]; then
        echo "  Compiling ${blas_file}.c..."
        emcc -c "$SRC_DIR/${blas_file}.c" -o "$OBJ_DIR/${blas_file}.o" "${CFLAGS[@]}"
        ODE_OBJS+=("$OBJ_DIR/${blas_file}.o")
    fi
done

# ============================================
# LINPACK routines (used by VODE family)
# ============================================

# Compile DGEFA (general matrix factorization)
if [ -f "$SRC_DIR/dgefa.c" ]; then
    echo "  Compiling dgefa.c..."
    emcc -c "$SRC_DIR/dgefa.c" -o "$OBJ_DIR/dgefa.o" "${CFLAGS[@]}"
    ODE_OBJS+=("$OBJ_DIR/dgefa.o")
fi

# Compile DGESL (general matrix solve)
if [ -f "$SRC_DIR/dgesl.c" ]; then
    echo "  Compiling dgesl.c..."
    emcc -c "$SRC_DIR/dgesl.c" -o "$OBJ_DIR/dgesl.o" "${CFLAGS[@]}"
    ODE_OBJS+=("$OBJ_DIR/dgesl.o")
fi

# Compile DGBFA (banded matrix factorization)
if [ -f "$SRC_DIR/dgbfa.c" ]; then
    echo "  Compiling dgbfa.c..."
    emcc -c "$SRC_DIR/dgbfa.c" -o "$OBJ_DIR/dgbfa.o" "${CFLAGS[@]}"
    ODE_OBJS+=("$OBJ_DIR/dgbfa.o")
fi

# Compile DGBSL (banded matrix solve)
if [ -f "$SRC_DIR/dgbsl.c" ]; then
    echo "  Compiling dgbsl.c..."
    emcc -c "$SRC_DIR/dgbsl.c" -o "$OBJ_DIR/dgbsl.o" "${CFLAGS[@]}"
    ODE_OBJS+=("$OBJ_DIR/dgbsl.o")
fi

# ============================================
# Complex LINPACK routines (used by ZVODE)
# ============================================

# Compile ZGEFA (complex general matrix factorization)
if [ -f "$SRC_DIR/zgefa.c" ]; then
    echo "  Compiling zgefa.c..."
    emcc -c "$SRC_DIR/zgefa.c" -o "$OBJ_DIR/zgefa.o" "${CFLAGS[@]}"
    ODE_OBJS+=("$OBJ_DIR/zgefa.o")
fi

# Compile ZGESL (complex general matrix solve)
if [ -f "$SRC_DIR/zgesl.c" ]; then
    echo "  Compiling zgesl.c..."
    emcc -c "$SRC_DIR/zgesl.c" -o "$OBJ_DIR/zgesl.o" "${CFLAGS[@]}"
    ODE_OBJS+=("$OBJ_DIR/zgesl.o")
fi

# Compile ZGBFA (complex banded matrix factorization)
if [ -f "$SRC_DIR/zgbfa.c" ]; then
    echo "  Compiling zgbfa.c..."
    emcc -c "$SRC_DIR/zgbfa.c" -o "$OBJ_DIR/zgbfa.o" "${CFLAGS[@]}"
    ODE_OBJS+=("$OBJ_DIR/zgbfa.o")
fi

# Compile ZGBSL (complex banded matrix solve)
if [ -f "$SRC_DIR/zgbsl.c" ]; then
    echo "  Compiling zgbsl.c..."
    emcc -c "$SRC_DIR/zgbsl.c" -o "$OBJ_DIR/zgbsl.o" "${CFLAGS[@]}"
    ODE_OBJS+=("$OBJ_DIR/zgbsl.o")
fi

# ============================================
# Hairer Solvers
# ============================================

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

# ============================================
# Netlib Solvers
# ============================================

# Compile RKF45
if [ -f "$SRC_DIR/rkf45.c" ]; then
    echo "  Compiling rkf45.c..."
    emcc -c "$SRC_DIR/rkf45.c" -o "$OBJ_DIR/rkf45.o" "${CFLAGS[@]}"
    ODE_OBJS+=("$OBJ_DIR/rkf45.o")
fi

# Compile DVERK
if [ -f "$SRC_DIR/dverk.c" ]; then
    echo "  Compiling dverk.c..."
    emcc -c "$SRC_DIR/dverk.c" -o "$OBJ_DIR/dverk.o" "${CFLAGS[@]}"
    ODE_OBJS+=("$OBJ_DIR/dverk.o")
fi

# Compile ODE (Adams-Bashforth-Moulton)
if [ -f "$SRC_DIR/ode.c" ]; then
    echo "  Compiling ode.c..."
    emcc -c "$SRC_DIR/ode.c" -o "$OBJ_DIR/ode.o" "${CFLAGS[@]}"
    ODE_OBJS+=("$OBJ_DIR/ode.o")
fi

# Compile VODE
if [ -f "$SRC_DIR/vode.c" ]; then
    echo "  Compiling vode.c..."
    emcc -c "$SRC_DIR/vode.c" -o "$OBJ_DIR/vode.o" "${CFLAGS[@]}"
    ODE_OBJS+=("$OBJ_DIR/vode.o")
fi

# Compile ZVODE
if [ -f "$SRC_DIR/zvode.c" ]; then
    echo "  Compiling zvode.c..."
    emcc -c "$SRC_DIR/zvode.c" -o "$OBJ_DIR/zvode.o" "${CFLAGS[@]}"
    ODE_OBJS+=("$OBJ_DIR/zvode.o")
fi

# Compile EPSODE
if [ -f "$SRC_DIR/epsode.c" ]; then
    echo "  Compiling epsode.c..."
    emcc -c "$SRC_DIR/epsode.c" -o "$OBJ_DIR/epsode.o" "${CFLAGS[@]}"
    ODE_OBJS+=("$OBJ_DIR/epsode.o")
fi

# Compile VODPK
if [ -f "$SRC_DIR/vodpk.c" ]; then
    echo "  Compiling vodpk.c..."
    emcc -c "$SRC_DIR/vodpk.c" -o "$OBJ_DIR/vodpk.o" "${CFLAGS[@]}"
    ODE_OBJS+=("$OBJ_DIR/vodpk.o")
fi

# Compile RKSUITE
if [ -f "$SRC_DIR/rksuite.c" ]; then
    echo "  Compiling rksuite.c..."
    emcc -c "$SRC_DIR/rksuite.c" -o "$OBJ_DIR/rksuite.o" "${CFLAGS[@]}"
    ODE_OBJS+=("$OBJ_DIR/rksuite.o")
fi

# Compile RKC
if [ -f "$SRC_DIR/rkc.c" ]; then
    echo "  Compiling rkc.c..."
    emcc -c "$SRC_DIR/rkc.c" -o "$OBJ_DIR/rkc.o" "${CFLAGS[@]}"
    ODE_OBJS+=("$OBJ_DIR/rkc.o")
fi

echo "  Total object files: ${#ODE_OBJS[@]}"
echo ""

# Exported functions
EXPORTED_FUNCTIONS='[
    "_wasm_dopri5", "_wasm_dop853", "_wasm_radau5",
    "_wasm_rkf45", "_wasm_dverk", "_wasm_ode",
    "_wasm_vode", "_wasm_zvode", "_wasm_vodpk",
    "_wasm_rksuite_setup", "_wasm_rksuite_ut", "_wasm_rksuite_ct", "_wasm_rksuite_stat",
    "_wasm_rkc", "_wasm_rkc_int",
    "_wasm_contd5", "_wasm_contd8", "_wasm_contr5",
    "_wasm_set_fcn_callback", "_wasm_set_solout_callback", "_wasm_set_jac_callback",
    "_wasm_set_jac_vode_callback", "_wasm_set_psol_callback", "_wasm_set_spcrad_callback",
    "_wasm_dopri5_work_size", "_wasm_dopri5_iwork_size",
    "_wasm_dop853_work_size", "_wasm_dop853_iwork_size",
    "_wasm_radau5_work_size", "_wasm_radau5_iwork_size",
    "_wasm_rkf45_work_size", "_wasm_rkf45_iwork_size",
    "_wasm_dverk_c_size", "_wasm_dverk_work_size",
    "_wasm_ode_work_size", "_wasm_ode_iwork_size",
    "_wasm_vode_rwork_size", "_wasm_vode_iwork_size",
    "_wasm_zvode_zwork_size", "_wasm_zvode_rwork_size",
    "_wasm_vodpk_rwork_size", "_wasm_vodpk_iwork_size",
    "_wasm_rksuite_work_size",
    "_wasm_rkc_work_size",
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
    -Wl,--allow-multiple-definition
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
