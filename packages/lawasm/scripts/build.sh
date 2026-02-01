#!/bin/bash
#
# Build LAPACK + BLAS WebAssembly module using native Fortran compilation
# Uses the r-wasm/flang-wasm Docker container with patched LLVM Flang
#
# This compiles Fortran directly to WebAssembly - no f2c conversion needed!
#
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PROJECT_DIR="$(dirname "$SCRIPT_DIR")"
LAPACK_SRC="$PROJECT_DIR/src/wasm"
OUT_DIR="$PROJECT_DIR/dist/wasm"
BUILD_DIR="$PROJECT_DIR/.build-lapack-flang"

DOCKER_IMAGE="ghcr.io/r-wasm/flang-wasm:v20.1.4"

echo "=============================================="
echo "LAPACK + BLAS WASM Build (Native Fortran)"
echo "=============================================="
echo "Using Docker image: $DOCKER_IMAGE"
echo "Source: $LAPACK_SRC"
echo ""

# Check for Docker
if ! command -v docker &> /dev/null; then
    echo "Error: docker not found"
    exit 1
fi

# Check for LAPACK source files
if [ ! -d "$LAPACK_SRC/SRC" ] || [ ! -d "$LAPACK_SRC/BLAS" ]; then
    echo "Error: LAPACK source not found at $LAPACK_SRC"
    echo "Expected directories: SRC/, BLAS/, INSTALL/"
    exit 1
fi

# Create directories
mkdir -p "$OUT_DIR" "$BUILD_DIR"

# Create a build script to run inside the container
cat > "$BUILD_DIR/build-inside.sh" << 'INNERSCRIPT'
#!/bin/bash
set -e

FLANG="/opt/flang/host/bin/flang"
EMCC="/opt/emsdk/upstream/emscripten/emcc"
EMAR="/opt/emsdk/upstream/emscripten/emar"
EMRANLIB="/opt/emsdk/upstream/emscripten/emranlib"

LAPACK_SRC="/src/lapack"
BUILD="/build"
OUT="/out"

FFLAGS="-O2 -fPIC"

cd "$BUILD"

# ===========================================
# Build BLAS
# ===========================================
echo ""
echo "Building BLAS..."
echo "================"

mkdir -p blas_obj
cd blas_obj

count=0
# Compile .f files (F77)
for f in "$LAPACK_SRC/BLAS"/*.f; do
    if [ -f "$f" ]; then
        base=$(basename "$f")
        base_noext="${base%.*}"
        echo "  Compiling: $base"
        $FLANG $FFLAGS -c "$f" -o "${base_noext}.o" 2>&1 || {
            echo "  Warning: Failed to compile $base"
            continue
        }
        ((count++)) || true
    fi
done

# Compile .f90 files (Fortran 90) - includes dnrm2, snrm2, dznrm2, scnrm2, etc.
for f in "$LAPACK_SRC/BLAS"/*.f90; do
    if [ -f "$f" ]; then
        base=$(basename "$f")
        base_noext="${base%.*}"
        echo "  Compiling: $base"
        $FLANG $FFLAGS -c "$f" -o "${base_noext}.o" 2>&1 || {
            echo "  Warning: Failed to compile $base"
            continue
        }
        ((count++)) || true
    fi
done

echo "  Creating libblas.a ($count objects)"
$EMAR rcs libblas.a *.o
$EMRANLIB libblas.a

cp libblas.a "$BUILD/"
cd "$BUILD"

# ===========================================
# Build LAPACK utility/install files first
# ===========================================
echo ""
echo "Building LAPACK INSTALL utilities..."
echo "====================================="

mkdir -p install_obj
cd install_obj

# Compile machine-dependent files from INSTALL directory
for f in dlamch slamch droundup_lwork sroundup_lwork ilaver; do
    src="$LAPACK_SRC/INSTALL/${f}.f"
    if [ -f "$src" ]; then
        echo "  Compiling: ${f}.f"
        $FLANG $FFLAGS -c "$src" -o "${f}.o" 2>&1 || echo "  Warning: Failed to compile ${f}.f"
    fi
done

# Also compile second/dsecnd (timer functions)
for f in second_INT_CPU_TIME dsecnd_INT_CPU_TIME; do
    src="$LAPACK_SRC/INSTALL/${f}.f"
    if [ -f "$src" ]; then
        echo "  Compiling: ${f}.f"
        $FLANG $FFLAGS -c "$src" -o "${f}.o" 2>&1 || echo "  Warning: Failed"
    fi
done

cd "$BUILD"

# ===========================================
# Build Fortran 90 Modules first (dependency order)
# ===========================================
echo ""
echo "Building Fortran 90 Modules..."
echo "=============================="

mkdir -p f90_modules
cd f90_modules

# 1. First compile LA_CONSTANTS module (no dependencies)
echo "  Compiling: la_constants.f90"
$FLANG $FFLAGS -c "$LAPACK_SRC/SRC/la_constants.f90" -o la_constants.o 2>&1 || {
    echo "  ERROR: Failed to compile la_constants.f90"
    exit 1
}

# 2. Compile LA_XISNAN module (depends on LA_CONSTANTS)
#    The .F90 file needs preprocessing. Use -cpp flag or just compile directly.
#    We'll use the fallback mode (no USE_IEEE_INTRINSIC or USE_ISNAN defined)
echo "  Compiling: la_xisnan.F90"
$FLANG $FFLAGS -c "$LAPACK_SRC/SRC/la_xisnan.F90" -o la_xisnan.o 2>&1 || {
    echo "  ERROR: Failed to compile la_xisnan.F90"
    exit 1
}

# 3. Compile all remaining .f90 files from SRC (includes dlartg, slartg, etc.)
for f in "$LAPACK_SRC/SRC"/*.f90; do
    if [ -f "$f" ]; then
        base=$(basename "$f")
        base_noext="${base%.*}"
        # Skip already compiled modules
        if [ "$base_noext" = "la_constants" ]; then
            continue
        fi
        echo "  Compiling: $base"
        $FLANG $FFLAGS -c "$f" -o "${base_noext}.o" 2>&1 || {
            echo "  Warning: Failed to compile $base"
            continue
        }
    fi
done

# List compiled module objects
echo "  F90 module objects: $(ls -1 *.o 2>/dev/null | wc -l)"

cd "$BUILD"

# ===========================================
# Build LAPACK (F77 files)
# ===========================================
echo ""
echo "Building LAPACK..."
echo "=================="

mkdir -p lapack_obj
cd lapack_obj

count=0
# Compile .f files (F77)
for f in "$LAPACK_SRC/SRC"/*.f; do
    if [ -f "$f" ]; then
        base=$(basename "$f")
        base_noext="${base%.*}"
        # Skip deprecated files
        if [[ "$f" == *"/DEPRECATED/"* ]]; then
            continue
        fi
        echo "  Compiling: $base"
        $FLANG $FFLAGS -c "$f" -o "${base_noext}.o" 2>&1 || {
            echo "  Warning: Failed to compile $base"
            continue
        }
        ((count++)) || true
    fi
done

# Copy install objects
cp "$BUILD/install_obj"/*.o . 2>/dev/null || true

# Copy F90 module objects
cp "$BUILD/f90_modules"/*.o . 2>/dev/null || true

echo "  Creating liblapack.a ($count+ F77 objects + F90 modules)"
$EMAR rcs liblapack.a *.o
$EMRANLIB liblapack.a

cp liblapack.a "$BUILD/"
cd "$BUILD"

# ===========================================
# Link into WebAssembly module
# ===========================================
echo ""
echo "Linking WASM module..."
echo "======================"

# Create a minimal main to satisfy emcc
cat > empty_main.c << 'EOF'
// Empty main - LAPACK is a library
int main() { return 0; }
EOF

$EMCC -c empty_main.c -o empty_main.o

# Export the most important LAPACK/BLAS functions
EXPORTED_FUNCTIONS='[
    "_dgetrf_","_dgetrs_","_dgetri_","_dgesv_","_dgecon_",
    "_dpotrf_","_dpotrs_","_dpotri_","_dposv_","_dpocon_",
    "_dgeev_","_dsyev_","_dsyevd_","_dsyevr_",
    "_dgesvd_","_dgesdd_",
    "_dgeqrf_","_dgeqp3_","_dorgqr_","_dormqr_",
    "_dgels_","_dgelsd_","_dgelss_","_dgelsy_",
    "_dtrtrs_","_dtrtri_",
    "_dlamch_","_dlange_","_dlansy_","_dlaswp_","_dlassq_",
    "_sgetrf_","_sgetrs_","_sgetri_","_sgesv_",
    "_spotrf_","_spotrs_","_spotri_","_sposv_",
    "_sgeev_","_ssyev_","_ssyevd_","_ssyevr_",
    "_sgesvd_","_sgesdd_",
    "_sgeqrf_","_sgeqp3_","_sorgqr_","_sormqr_",
    "_sgels_","_sgelsd_","_sgelss_","_sgelsy_",
    "_strtrs_","_strtri_",
    "_slamch_","_slange_","_slansy_","_slaswp_","_slassq_",
    "_zgetrf_","_zgetrs_","_zgetri_","_zgesv_",
    "_zpotrf_","_zpotrs_","_zpotri_","_zposv_",
    "_zgeev_","_zheev_","_zheevd_","_zheevr_",
    "_zgesvd_","_zgesdd_",
    "_zgeqrf_","_zgeqp3_","_zungqr_","_zunmqr_",
    "_zgels_","_zgelsd_","_zgelss_","_zgelsy_",
    "_ztrtrs_","_ztrtri_","_zlassq_",
    "_cgetrf_","_cgetrs_","_cgetri_","_cgesv_",
    "_cpotrf_","_cpotrs_","_cpotri_","_cposv_",
    "_cgeev_","_cheev_","_cheevd_","_cheevr_",
    "_cgesvd_","_cgesdd_",
    "_cgeqrf_","_cgeqp3_","_cungqr_","_cunmqr_",
    "_cgels_","_cgelsd_","_cgelss_","_cgelsy_",
    "_ctrtrs_","_ctrtri_","_classq_",
    "_daxpy_","_dcopy_","_dscal_","_dswap_","_ddot_","_dnrm2_","_dasum_","_idamax_",
    "_dgemv_","_dtrsv_","_dtrmv_","_dger_","_dsyr_","_dsyr2_","_dsymv_",
    "_dgemm_","_dtrsm_","_dtrmm_","_dsyrk_","_dsyr2k_",
    "_saxpy_","_scopy_","_sscal_","_sswap_","_sdot_","_snrm2_","_sasum_","_isamax_",
    "_sgemv_","_strsv_","_strmv_","_sger_","_ssyr_","_ssyr2_","_ssymv_",
    "_sgemm_","_strsm_","_strmm_","_ssyrk_","_ssyr2k_",
    "_zaxpy_","_zcopy_","_zscal_","_zswap_","_zdotc_","_zdotu_","_dznrm2_","_dzasum_","_izamax_",
    "_zgemv_","_ztrsv_","_ztrmv_","_zgerc_","_zgeru_","_zher_","_zher2_","_zhemv_",
    "_zgemm_","_ztrsm_","_ztrmm_","_zherk_","_zher2k_",
    "_caxpy_","_ccopy_","_cscal_","_cswap_","_cdotc_","_cdotu_","_scnrm2_","_scasum_","_icamax_",
    "_cgemv_","_ctrsv_","_ctrmv_","_cgerc_","_cgeru_","_cher_","_cher2_","_chemv_",
    "_cgemm_","_ctrsm_","_ctrmm_","_cherk_","_cher2k_",
    "_lsame_",
    "_malloc","_free"
]'

EXPORTED_RUNTIME='["ccall","cwrap","getValue","setValue","HEAPF64","HEAPF32","HEAP32","HEAPU8"]'

# Build CJS module
echo "  Creating CJS module..."
$EMCC empty_main.o \
    -L. -llapack -lblas \
    -L/opt/flang/wasm/lib -lFortranRuntime \
    -s WASM=1 \
    -s MODULARIZE=1 \
    -s EXPORT_NAME="createLAPACKModule" \
    -s EXPORTED_FUNCTIONS="$EXPORTED_FUNCTIONS" \
    -s EXPORTED_RUNTIME_METHODS="$EXPORTED_RUNTIME" \
    -s ALLOW_MEMORY_GROWTH=1 \
    -s INITIAL_MEMORY=16777216 \
    -s STACK_SIZE=1048576 \
    -O2 \
    -o "$OUT/lapack.cjs"

# Build ESM module
echo "  Creating ESM module..."
$EMCC empty_main.o \
    -L. -llapack -lblas \
    -L/opt/flang/wasm/lib -lFortranRuntime \
    -s WASM=1 \
    -s MODULARIZE=1 \
    -s EXPORT_NAME="createLAPACKModule" \
    -s EXPORT_ES6=1 \
    -s EXPORTED_FUNCTIONS="$EXPORTED_FUNCTIONS" \
    -s EXPORTED_RUNTIME_METHODS="$EXPORTED_RUNTIME" \
    -s ALLOW_MEMORY_GROWTH=1 \
    -s INITIAL_MEMORY=16777216 \
    -s STACK_SIZE=1048576 \
    -O2 \
    -o "$OUT/lapack.mjs"

echo ""
echo "Build complete!"
ls -lh "$OUT"/lapack.*

INNERSCRIPT

chmod +x "$BUILD_DIR/build-inside.sh"

# Run the build inside Docker
echo "Starting Docker build..."
docker run --rm \
    -v "$LAPACK_SRC:/src/lapack:ro" \
    -v "$BUILD_DIR:/build" \
    -v "$OUT_DIR:/out" \
    "$DOCKER_IMAGE" \
    /bin/bash /build/build-inside.sh

echo ""
echo "=============================================="
echo "Build Complete!"
echo "=============================================="
echo ""
echo "Output files:"
ls -lh "$OUT_DIR/lapack."* 2>/dev/null || echo "  (no output files found)"
