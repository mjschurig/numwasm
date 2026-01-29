#!/bin/bash
set -e

echo "Building SymWASM module..."

# Directories
BUILD_DIR=".build-wasm"
OUT_DIR="dist/wasm"
SRC_DIR="src/wasm/symengine"

# Create build directories
mkdir -p "$BUILD_DIR"
mkdir -p "$OUT_DIR"

# Source files to compile - ALL SymEngine core files
# Excludes: MPFR, MPC, ARB, FLINT, Piranha, Boost, LLVM (optional deps)
# Excludes: tests, utilities/matchpycpp, utilities/catch, utilities/teuchos
SOURCES=(
  # Core
  "add.cpp"
  # "as_real_imag.cpp"  # duplicates symbols in solve.cpp
  "assumptions.cpp"
  "basic.cpp"
  "complex.cpp"
  "complex_double.cpp"
  "constants.cpp"
  "cse.cpp"
  "cwrapper.cpp"
  "dense_matrix.cpp"
  "derivative.cpp"
  "dict.cpp"
  "diophantine.cpp"
  "eval.cpp"
  "eval_double.cpp"
  "expand.cpp"
  "expression.cpp"
  "fields.cpp"
  "finitediff.cpp"
  "functions.cpp"
  "infinity.cpp"
  "integer.cpp"
  "logic.cpp"
  "matrix.cpp"
  "monomials.cpp"
  "mp_wrapper.cpp"
  "mul.cpp"
  "nan.cpp"
  "ntheory.cpp"
  "ntheory_funcs.cpp"
  "number.cpp"
  "numer_denom.cpp"
  "pow.cpp"
  "prime_sieve.cpp"
  "rational.cpp"
  "real_double.cpp"
  "refine.cpp"
  "rewrite.cpp"
  "rings.cpp"
  "series.cpp"
  "series_generic.cpp"
  "set_funcs.cpp"
  "sets.cpp"
  "simplify.cpp"
  "solve.cpp"
  "sparse_matrix.cpp"
  "symbol.cpp"
  "symengine_rcp.cpp"
  "tuple.cpp"
  "visitor.cpp"
  # Printers
  "printers/codegen.cpp"
  "printers/latex.cpp"
  "printers/mathml.cpp"
  "printers/sbml.cpp"
  "printers/stringbox.cpp"
  "printers/strprinter.cpp"
  "printers/unicode.cpp"
  # Polynomials
  "polys/basic_conversions.cpp"
  "polys/msymenginepoly.cpp"
  "polys/uexprpoly.cpp"
  "polys/uintpoly.cpp"
  "polys/uratpoly.cpp"
  # Matrices
  "matrices/conjugate_matrix.cpp"
  "matrices/diagonal_matrix.cpp"
  "matrices/hadamard_product.cpp"
  "matrices/identity_matrix.cpp"
  "matrices/immutable_dense_matrix.cpp"
  "matrices/is_diagonal.cpp"
  "matrices/is_lower.cpp"
  "matrices/is_real.cpp"
  "matrices/is_square.cpp"
  "matrices/is_symmetric.cpp"
  "matrices/is_toeplitz.cpp"
  "matrices/is_upper.cpp"
  "matrices/is_zero.cpp"
  "matrices/matrix_add.cpp"
  "matrices/matrix_mul.cpp"
  "matrices/matrix_symbol.cpp"
  "matrices/size.cpp"
  "matrices/trace.cpp"
  "matrices/transpose.cpp"
  "matrices/zero_matrix.cpp"
  # Parser - excluded (requires fast_float dependency)
  # "parser/parser.cpp"
  # "parser/tokenizer.cpp"
)

# GMP paths
GMP_PREFIX=".gmp-build"

# Compiler flags
CFLAGS=(
  "-O2"
  "-std=c++11"
  "-I" "src/wasm"
  "-I" "src/wasm/symengine"
  "-I" "src/wasm/symengine/utilities"
  "-I" "$GMP_PREFIX/include"
  "-D" "WITH_SYMENGINE_RCP"
  "-D" "WITH_SYMENGINE_THREAD_SAFE"
  "-D" "HAVE_SYMENGINE_GMP"
  "-D" "symengine_EXPORTS"
)

# Compile each source file
echo "Compiling C++ sources to object files..."
OBJECTS=()
for src in "${SOURCES[@]}"; do
  src_file="$SRC_DIR/$src"
  obj_file="$BUILD_DIR/$(basename ${src%.cpp}.o)"

  if [ -f "$src_file" ]; then
    echo "  Compiling $(basename $src)..."
    emcc -c "$src_file" -o "$obj_file" "${CFLAGS[@]}"
    OBJECTS+=("$obj_file")
  else
    echo "  ⚠ Warning: $src_file not found, skipping"
  fi
done

# Exported functions (C API from cwrapper.h) - JSON array format
EXPORTED_FUNCTIONS='[
    "_malloc",
    "_free",
    "_basic_new_stack",
    "_basic_free_stack",
    "_basic_new_heap",
    "_basic_free_heap",
    "_basic_assign",
    "_symbol_set",
    "_integer_set_si",
    "_integer_set_ui",
    "_rational_set",
    "_rational_set_si",
    "_rational_set_ui",
    "_real_double_set_d",
    "_basic_add",
    "_basic_sub",
    "_basic_mul",
    "_basic_div",
    "_basic_pow",
    "_basic_str",
    "_basic_get_type",
    "_basic_eq",
    "_basic_hash",
    "_basic_const_zero",
    "_basic_const_one",
    "_basic_const_minus_one",
    "_basic_const_I",
    "_basic_const_pi",
    "_basic_const_E",
    "_basic_free_symbols",
    "_setbasic_new",
    "_setbasic_free",
    "_setbasic_get",
    "_setbasic_size",
    "_complex_set",
    "_basic_neg",
    "_basic_get_args",
    "_vecbasic_new",
    "_vecbasic_free",
    "_vecbasic_size",
    "_vecbasic_get",
    "_vecbasic_push_back",
    "_basic_const_infinity",
    "_basic_const_neginfinity",
    "_basic_const_complex_infinity",
    "_basic_const_EulerGamma",
    "_basic_const_Catalan",
    "_basic_const_GoldenRatio",
    "_basic_const_nan",
    "_basic_subs2",
    "_basic_subs",
    "_mapbasicbasic_new",
    "_mapbasicbasic_free",
    "_mapbasicbasic_insert",
    "_mapbasicbasic_size",
    "_basic_evalf",
    "_real_double_get_d"
]'

EXPORTED_RUNTIME='["ccall", "cwrap", "getValue", "setValue", "UTF8ToString", "stringToUTF8", "lengthBytesUTF8", "HEAPF64", "HEAP32", "HEAPU8"]'

# Link flags (following sciwasm pattern - bare flags, no quotes around -s)
LINK_FLAGS=(
    -s WASM=1
    -s MODULARIZE=1
    -s EXPORTED_FUNCTIONS="$EXPORTED_FUNCTIONS"
    -s EXPORTED_RUNTIME_METHODS="$EXPORTED_RUNTIME"
    -s ALLOW_MEMORY_GROWTH=1
    -s INITIAL_MEMORY=16777216
    -s MAXIMUM_MEMORY=2147483648
    -O2
    -L"$GMP_PREFIX/lib"
    -lgmp
)

# Link to CommonJS module (for Node.js)
echo "Linking WASM module (CJS)..."
emcc \
    "${OBJECTS[@]}" \
    "${LINK_FLAGS[@]}" \
    -s EXPORT_NAME="createSymwasmModule" \
    -o "$OUT_DIR/symwasm.cjs"
echo "  ✓ $OUT_DIR/symwasm.cjs"
echo "  ✓ $OUT_DIR/symwasm.wasm"

# Link to ES module (for browsers)
echo "Linking WASM module (ESM)..."
emcc \
    "${OBJECTS[@]}" \
    "${LINK_FLAGS[@]}" \
    -s EXPORT_NAME="createSymwasmModule" \
    -s EXPORT_ES6=1 \
    -o "$OUT_DIR/symwasm.mjs"
echo "  ✓ $OUT_DIR/symwasm.mjs"

echo "✓ Build complete!"
echo ""
echo "Output files:"
ls -lh "$OUT_DIR"
