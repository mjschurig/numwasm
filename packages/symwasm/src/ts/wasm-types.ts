/**
 * Type definitions for the SymWASM module
 * These interfaces define the WASM module's API and memory access
 */

export interface EmscriptenModule {
  // Memory views
  HEAPF64: Float64Array;
  HEAP32: Int32Array;
  HEAPU8: Uint8Array;
  HEAPF32: Float32Array;
  HEAP16: Int16Array;
  HEAPU16: Uint16Array;
  HEAP8: Int8Array;
  HEAPU32: Uint32Array;

  // Memory management
  _malloc(size: number): number;
  _free(ptr: number): void;

  // Utility functions
  ccall(
    name: string,
    returnType: string | null,
    argTypes: string[],
    args: any[]
  ): any;
  cwrap(
    name: string,
    returnType: string | null,
    argTypes: string[]
  ): (...args: any[]) => any;
  UTF8ToString(ptr: number, maxBytesToRead?: number): string;
  stringToUTF8(
    str: string,
    outPtr: number,
    maxBytesToWrite: number
  ): void;
  lengthBytesUTF8(str: string): number;
  getValue(ptr: number, type: string): number;
  setValue(ptr: number, value: number, type: string): void;
}

/**
 * SymWASM module interface
 * Extends Emscripten module with SymEngine C API functions
 */
export interface SymwasmModule extends EmscriptenModule {
  // === Memory Management ===
  _basic_new_stack(): number; // For stack-allocated basic (requires pre-allocated memory)
  _basic_free_stack(ptr: number): void;
  _basic_new_heap(): number; // For heap-allocated basic (allocates and returns pointer)
  _basic_free_heap(ptr: number): void;

  // === Symbol Creation ===
  _symbol_set(ptr: number, namePtr: number): number; // namePtr is a C string pointer

  // === Number Creation ===
  _integer_set_si(ptr: number, value: number): number;
  _integer_set_ui(ptr: number, value: number): number;
  _rational_set(ptr: number, num: number, den: number): number;
  _rational_set_si(ptr: number, num: number, den: number): number;
  _rational_set_ui(ptr: number, num: number, den: number): number;
  _real_double_set_d(ptr: number, value: number): number;
  _complex_set(s: number, re: number, im: number): number; // re, im are Basic pointers

  // === Basic Arithmetic ===
  _basic_add(result: number, a: number, b: number): number;
  _basic_sub(result: number, a: number, b: number): number;
  _basic_mul(result: number, a: number, b: number): number;
  _basic_div(result: number, a: number, b: number): number;
  _basic_pow(result: number, base: number, exp: number): number;
  _basic_neg(result: number, a: number): number;

  // === Expression Arguments ===
  _basic_get_args(self: number, args: number): number; // args is CVecBasic*

  // === Vector Container Operations (CVecBasic) ===
  _vecbasic_new(): number; // Returns CVecBasic*
  _vecbasic_free(self: number): void;
  _vecbasic_size(self: number): number; // Returns size_t
  _vecbasic_get(self: number, n: number, result: number): number; // Gets nth element
  _vecbasic_push_back(self: number, value: number): number; // Appends element

  // === Constants ===
  _basic_const_zero(ptr: number): number;
  _basic_const_one(ptr: number): number;
  _basic_const_minus_one(ptr: number): number;
  _basic_const_I(ptr: number): number;
  _basic_const_pi(ptr: number): number;
  _basic_const_E(ptr: number): number;
  _basic_const_infinity(ptr: number): number;
  _basic_const_neginfinity(ptr: number): number;
  _basic_const_complex_infinity(ptr: number): number;
  _basic_const_EulerGamma(ptr: number): number;
  _basic_const_Catalan(ptr: number): number;
  _basic_const_GoldenRatio(ptr: number): number;
  _basic_const_nan(ptr: number): number;

  // === Substitution ===
  _basic_subs2(result: number, expr: number, old: number, new_: number): number;
  _basic_subs(result: number, expr: number, mapbb: number): number;

  // === Map Container Operations (CMapBasicBasic) ===
  _mapbasicbasic_new(): number;
  _mapbasicbasic_free(self: number): void;
  _mapbasicbasic_insert(self: number, key: number, mapped: number): void;
  _mapbasicbasic_size(self: number): number;

  // === Numerical Evaluation ===
  _basic_evalf(s: number, b: number, bits: number, real: number): number;
  _real_double_get_d(s: number): number;

  // === Elementary Functions (1-argument) ===
  _basic_sin(result: number, arg: number): number;
  _basic_cos(result: number, arg: number): number;
  _basic_tan(result: number, arg: number): number;
  _basic_cot(result: number, arg: number): number;
  _basic_sec(result: number, arg: number): number;
  _basic_csc(result: number, arg: number): number;
  _basic_asin(result: number, arg: number): number;
  _basic_acos(result: number, arg: number): number;
  _basic_atan(result: number, arg: number): number;
  _basic_acot(result: number, arg: number): number;
  _basic_asec(result: number, arg: number): number;
  _basic_acsc(result: number, arg: number): number;
  _basic_sinh(result: number, arg: number): number;
  _basic_cosh(result: number, arg: number): number;
  _basic_tanh(result: number, arg: number): number;
  _basic_coth(result: number, arg: number): number;
  _basic_sech(result: number, arg: number): number;
  _basic_csch(result: number, arg: number): number;
  _basic_asinh(result: number, arg: number): number;
  _basic_acosh(result: number, arg: number): number;
  _basic_atanh(result: number, arg: number): number;
  _basic_acoth(result: number, arg: number): number;
  _basic_asech(result: number, arg: number): number;
  _basic_acsch(result: number, arg: number): number;
  _basic_exp(result: number, arg: number): number;
  _basic_log(result: number, arg: number): number;
  _basic_sqrt(result: number, arg: number): number;
  _basic_cbrt(result: number, arg: number): number;
  _basic_lambertw(result: number, arg: number): number;
  _basic_abs(result: number, arg: number): number;
  _basic_sign(result: number, arg: number): number;
  _basic_floor(result: number, arg: number): number;
  _basic_ceiling(result: number, arg: number): number;
  _basic_gamma(result: number, arg: number): number;
  _basic_loggamma(result: number, arg: number): number;
  _basic_erf(result: number, arg: number): number;
  _basic_erfc(result: number, arg: number): number;
  _basic_zeta(result: number, arg: number): number;
  _basic_dirichlet_eta(result: number, arg: number): number;

  // === Elementary Functions (2-argument) ===
  _basic_atan2(result: number, y: number, x: number): number;
  _basic_beta(result: number, a: number, b: number): number;
  _basic_lowergamma(result: number, s: number, x: number): number;
  _basic_uppergamma(result: number, s: number, x: number): number;
  _basic_polygamma(result: number, n: number, x: number): number;
  _basic_kronecker_delta(result: number, i: number, j: number): number;

  // === Additional Functions (Phase 2.1b) ===
  _basic_max(result: number, args: number): number; // args is CVecBasic*
  _basic_min(result: number, args: number): number; // args is CVecBasic*
  _complex_base_real_part(result: number, com: number): number;
  _complex_base_imaginary_part(result: number, com: number): number;
  _basic_digamma(result: number, arg: number): number;
  _basic_conjugate(result: number, arg: number): number;

  // === Calculus (Phase 2.2) ===
  _basic_diff(result: number, expr: number, symbol: number): number;

  // === Series Expansion (Phase 2.3) ===
  _basic_series(result: number, expr: number, symbol: number, prec: number): number;

  // === Simplification (Phase 2.4) ===
  _basic_expand(result: number, a: number): number;
  _basic_simplify(result: number, a: number): number;
  _basic_as_numer_denom(numer: number, denom: number, x: number): number;
  _basic_rewrite_as_exp(result: number, a: number): number;
  _basic_rewrite_as_sin(result: number, a: number): number;
  _basic_rewrite_as_cos(result: number, a: number): number;
  _basic_as_real_imag(real: number, imag: number, x: number): number;

  // === Type Information ===
  _basic_get_type(ptr: number): number;

  // === Comparison ===
  _basic_eq(a: number, b: number): number;
  _basic_hash(ptr: number): number;

  // === String Conversion ===
  _basic_str(ptr: number): number; // Returns char* (must be freed)

  // === Free Symbols Extraction ===
  _basic_free_symbols(self: number, symbols: number): number; // Returns exception code

  // === Set Container Operations (CSetBasic) ===
  _setbasic_new(): number; // Returns CSetBasic*
  _setbasic_free(self: number): void;
  _setbasic_get(self: number, n: number, result: number): void; // Gets nth element into result
  _setbasic_size(self: number): number; // Returns size_t

  // === Dense Matrix Operations (Phase 3.1) ===
  _dense_matrix_new(): number; // Returns CDenseMatrix*
  _dense_matrix_new_rows_cols(r: number, c: number): number; // Returns CDenseMatrix*
  _dense_matrix_new_vec(rows: number, cols: number, l: number): number; // Returns CDenseMatrix*, l is CVecBasic*
  _dense_matrix_free(self: number): void;
  _dense_matrix_rows(s: number): number;
  _dense_matrix_cols(s: number): number;
  _dense_matrix_str(s: number): number; // Returns char* (must be freed)
  _dense_matrix_eq(lhs: number, rhs: number): number; // Returns 1 if equal
  _dense_matrix_get_basic(s: number, mat: number, r: number, c: number): number;
  _dense_matrix_set_basic(mat: number, r: number, c: number, s: number): number;
  _dense_matrix_eye(s: number, N: number, M: number, k: number): number;
  _dense_matrix_zeros(s: number, r: number, c: number): number;
  _dense_matrix_ones(s: number, r: number, c: number): number;
  _dense_matrix_diag(s: number, d: number, k: number): number; // d is CVecBasic*

  // === Dense Matrix — Basic Operations ===
  _dense_matrix_det(s: number, mat: number): number; // s is Basic*, mat is CDenseMatrix*
  _dense_matrix_inv(s: number, mat: number): number; // s is CDenseMatrix*, mat is CDenseMatrix*
  _dense_matrix_transpose(s: number, mat: number): number; // s is CDenseMatrix*, mat is CDenseMatrix*
  _dense_matrix_add_matrix(s: number, matA: number, matB: number): number; // s = matA + matB
  _dense_matrix_mul_matrix(s: number, matA: number, matB: number): number; // s = matA * matB
  _dense_matrix_add_scalar(s: number, matA: number, b: number): number; // s = matA + b (b is Basic*)
  _dense_matrix_mul_scalar(s: number, matA: number, b: number): number; // s = matA * b (b is Basic*)

  // === Dense Matrix — Submatrix Operations ===
  _dense_matrix_submatrix(
    s: number,
    mat: number,
    r1: number,
    c1: number,
    r2: number,
    c2: number,
    r: number,
    c: number
  ): number;
  _dense_matrix_row_join(A: number, B: number): number;
  _dense_matrix_col_join(A: number, B: number): number;
  _dense_matrix_row_del(C: number, k: number): number;
  _dense_matrix_col_del(C: number, k: number): number;

  // === Dense Matrix — Factorizations ===
  _dense_matrix_LU(l: number, u: number, mat: number): number; // LU factorization: mat = L*U
  _dense_matrix_LDL(l: number, d: number, mat: number): number; // LDL factorization: mat = L*D*L^T
  _dense_matrix_FFLU(lu: number, mat: number): number; // Fraction-free LU factorization
  _dense_matrix_FFLDU(l: number, d: number, u: number, mat: number): number; // Fraction-free LDU factorization
  _dense_matrix_LU_solve(x: number, A: number, b: number): number; // Solve A*x = b using LU

  // === Common Subexpression Elimination ===
  _basic_cse(
    replacement_syms: number,
    replacement_exprs: number,
    reduced_exprs: number,
    exprs: number
  ): number; // replacement_syms, replacement_exprs, reduced_exprs, exprs are all CVecBasic*

  // === Dense Matrix — Calculus ===
  _dense_matrix_diff(result: number, A: number, x: number): number; // Elementwise derivative of A w.r.t. x
  _dense_matrix_jacobian(result: number, A: number, x: number): number; // Jacobian of A w.r.t. x (A is column vector, x is column vector of symbols)

  // === Sparse Matrix Operations ===
  _sparse_matrix_new(): number; // Returns CSparseMatrix*
  _sparse_matrix_free(self: number): void;
  _sparse_matrix_init(s: number): void;
  _sparse_matrix_rows_cols(s: number, r: number, c: number): void;
  _sparse_matrix_str(s: number): number; // Returns char* (must be freed)
  _sparse_matrix_get_basic(s: number, mat: number, r: number, c: number): number;
  _sparse_matrix_set_basic(mat: number, r: number, c: number, s: number): number;
  _sparse_matrix_eq(lhs: number, rhs: number): number; // Returns 1 if equal
}

/**
 * SymWASM module factory function type
 */
export interface SymwasmModuleFactory {
  (options?: Partial<EmscriptenModule>): Promise<SymwasmModule>;
}

/**
 * SymEngine type IDs (from type_codes.inc)
 * These correspond to the TypeID enum in SymEngine
 * Note: Values depend on compile-time configuration (MPFR, MPC, PIRANHA, FLINT)
 * These values are empirically determined from the WASM build
 */
export enum SymEngineTypeID {
  SYMENGINE_INTEGER = 0,
  SYMENGINE_RATIONAL = 1,
  SYMENGINE_COMPLEX = 2,
  SYMENGINE_COMPLEX_DOUBLE = 3,
  SYMENGINE_REAL_MPFR = 4,
  SYMENGINE_COMPLEX_MPC = 5,
  SYMENGINE_REAL_DOUBLE = 6,
  SYMENGINE_INFTY = 7,
  SYMENGINE_NOT_A_NUMBER = 8,
  SYMENGINE_URATPSERIESPIRANHA = 9,
  SYMENGINE_UPSERIESPIRANHA = 10,
  SYMENGINE_URATPSERIESFLINT = 11,
  SYMENGINE_NUMBER_WRAPPER = 12,
  // NUMBER_WRAPPER marks the end of Number subclasses
  SYMENGINE_SYMBOL = 13,
  SYMENGINE_DUMMY = 14,
  SYMENGINE_MUL = 15,
  SYMENGINE_ADD = 16,
  SYMENGINE_POW = 17,
  // ... more types follow
  SYMENGINE_CONSTANT = 31,
  // Add more as needed
}

/**
 * SymEngine exception codes (from symengine_exception.h)
 */
export enum SymEngineException {
  NO_EXCEPTION = 0,
  RUNTIME_ERROR = 1,
  DIV_BY_ZERO = 2,
  NOT_IMPLEMENTED = 3,
  DOMAIN_ERROR = 4,
  PARSE_ERROR = 5,
}

/**
 * Domain for numerical evaluation (from eval.h)
 */
export enum EvalfDomain {
  Complex = 0,
  Real = 1,
  Symbolic = 2,
}
