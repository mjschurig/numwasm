/**
 * WASM entry points for gamma functions from xsf library.
 *
 * These thin wrappers expose the xsf::gamma and related functions
 * for calling from JavaScript/TypeScript via Emscripten.
 */

#include "../xsf/gamma.h"
#include "../xsf/loggamma.h"

extern "C" {

/**
 * Gamma function: Γ(x)
 * Returns the gamma function of x.
 */
double wasm_gamma(double x) {
    return xsf::gamma(x);
}

/**
 * Log gamma function: ln(|Γ(x)|)
 * Returns the natural logarithm of the absolute value of the gamma function.
 */
double wasm_gammaln(double x) {
    return xsf::gammaln(x);
}

/**
 * Reciprocal gamma function: 1/Γ(x)
 * Returns the reciprocal of the gamma function.
 * More numerically stable than computing 1/gamma(x) for large |x|.
 */
double wasm_rgamma(double x) {
    return xsf::cephes::rgamma(x);
}

} // extern "C"
