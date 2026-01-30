/**
 * WASM entry points for special functions from xsf library.
 *
 * All special function wrappers are in this single file to avoid
 * duplicate symbol issues from header-only templates.
 */

#include "../gamma.h"
#include "../loggamma.h"
#include "../beta.h"
#include "../erf.h"
#include "../cephes/j0.h"
#include "../cephes/j1.h"
#include "../cephes/jv.h"
#include "../cephes/yv.h"
#include "../cephes/i0.h"
#include "../cephes/i1.h"
#include "../cephes/k0.h"
#include "../cephes/k1.h"
#include "../cephes/scipy_iv.h"

extern "C" {

// ============================================================
// Gamma functions
// ============================================================

double wasm_gamma(double x) {
    return xsf::gamma(x);
}

double wasm_gammaln(double x) {
    return xsf::gammaln(x);
}

double wasm_rgamma(double x) {
    return xsf::cephes::rgamma(x);
}

// ============================================================
// Beta functions
// ============================================================

double wasm_beta(double a, double b) {
    return xsf::beta(a, b);
}

double wasm_betaln(double a, double b) {
    return xsf::betaln(a, b);
}

// ============================================================
// Error functions
// ============================================================

double wasm_erf(double x) {
    return xsf::erf(x);
}

double wasm_erfc(double x) {
    return xsf::erfc(x);
}

double wasm_erfcx(double x) {
    return xsf::erfcx(x);
}

double wasm_erfi(double x) {
    return xsf::erfi(x);
}

// ============================================================
// Bessel functions of the first kind (J)
// ============================================================

double wasm_j0(double x) {
    return xsf::cephes::j0(x);
}

double wasm_j1(double x) {
    return xsf::cephes::j1(x);
}

double wasm_jv(double v, double x) {
    return xsf::cephes::jv(v, x);
}

// ============================================================
// Bessel functions of the second kind (Y)
// ============================================================

double wasm_y0(double x) {
    return xsf::cephes::y0(x);
}

double wasm_y1(double x) {
    return xsf::cephes::y1(x);
}

double wasm_yv(double v, double x) {
    return xsf::cephes::yv(v, x);
}

// ============================================================
// Modified Bessel functions of the first kind (I)
// ============================================================

double wasm_i0(double x) {
    return xsf::cephes::i0(x);
}

double wasm_i1(double x) {
    return xsf::cephes::i1(x);
}

double wasm_iv(double v, double x) {
    return xsf::cephes::iv(v, x);
}

// ============================================================
// Modified Bessel functions of the second kind (K)
// ============================================================

double wasm_k0(double x) {
    return xsf::cephes::k0(x);
}

double wasm_k1(double x) {
    return xsf::cephes::k1(x);
}

} // extern "C"
