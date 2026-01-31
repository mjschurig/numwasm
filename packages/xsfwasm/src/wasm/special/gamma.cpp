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

// Additional headers for new functions
#include "../airy.h"
#include "../digamma.h"
#include "../ellip.h"
#include "../expint.h"
#include "../fresnel.h"
#include "../hyp2f1.h"
#include "../kelvin.h"
#include "../lambertw.h"
#include "../cephes/sici.h"
#include "../cephes/shichi.h"
#include "../bessel.h"  // Must be before sph_bessel.h
#include "../sph_bessel.h"
#include "../struve.h"
#include "../zeta.h"
#include "../legendre.h"

// Statistics
#include "../cephes/bdtr.h"
#include "../cephes/chdtr.h"
#include "../cephes/fdtr.h"
#include "../cephes/gdtr.h"
#include "../cephes/kolmogorov.h"
#include "../cephes/nbdtr.h"
#include "../cephes/ndtr.h"
#include "../cephes/ndtri.h"
#include "../cephes/owens_t.h"
#include "../cephes/pdtr.h"

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

// ============================================================
// Airy functions
// ============================================================

void wasm_airy(double x, double* ai, double* aip, double* bi, double* bip) {
    xsf::cephes::airy(x, ai, aip, bi, bip);
}

// ============================================================
// Digamma (psi) function
// ============================================================

double wasm_digamma(double x) {
    return xsf::digamma(x);
}

// ============================================================
// Elliptic integrals
// ============================================================

double wasm_ellipk(double m) {
    return xsf::ellipk(m);
}

double wasm_ellipe(double m) {
    return xsf::ellipe(m);
}

double wasm_ellipkinc(double phi, double m) {
    return xsf::ellipkinc(phi, m);
}

double wasm_ellipeinc(double phi, double m) {
    return xsf::ellipeinc(phi, m);
}

void wasm_ellipj(double u, double m, double* sn, double* cn, double* dn, double* ph) {
    xsf::ellipj(u, m, *sn, *cn, *dn, *ph);
}

// ============================================================
// Exponential integrals
// ============================================================

double wasm_exp1(double x) {
    return xsf::exp1(x);
}

double wasm_expi(double x) {
    return xsf::expi(x);
}

// ============================================================
// Fresnel integrals
// ============================================================

void wasm_fresnel(double x, double* s, double* c) {
    xsf::cephes::fresnl(x, s, c);
}

// ============================================================
// Hypergeometric function
// ============================================================

double wasm_hyp2f1(double a, double b, double c, double x) {
    return xsf::hyp2f1(a, b, c, x);
}

// ============================================================
// Kelvin functions
// ============================================================

double wasm_ber(double x) {
    return xsf::ber(x);
}

double wasm_bei(double x) {
    return xsf::bei(x);
}

double wasm_ker(double x) {
    return xsf::ker(x);
}

double wasm_kei(double x) {
    return xsf::kei(x);
}

double wasm_berp(double x) {
    return xsf::berp(x);
}

double wasm_beip(double x) {
    return xsf::beip(x);
}

double wasm_kerp(double x) {
    return xsf::kerp(x);
}

double wasm_keip(double x) {
    return xsf::keip(x);
}

// ============================================================
// Lambert W function
// ============================================================

double wasm_lambertw(double x, int k) {
    std::complex<double> z(x, 0.0);
    std::complex<double> result = xsf::lambertw(z, static_cast<long>(k), 1e-8);
    return result.real();
}

// ============================================================
// Sine and cosine integrals
// ============================================================

void wasm_sici(double x, double* si, double* ci) {
    xsf::cephes::sici(x, *si, *ci);
}

void wasm_shichi(double x, double* shi, double* chi) {
    xsf::cephes::shichi(x, *shi, *chi);
}

// ============================================================
// Spherical Bessel functions
// ============================================================

double wasm_spherical_jn(int n, double x) {
    return xsf::sph_bessel_j(static_cast<long>(n), x);
}

double wasm_spherical_yn(int n, double x) {
    return xsf::sph_bessel_y(static_cast<long>(n), x);
}

double wasm_spherical_in(int n, double x) {
    return xsf::sph_bessel_i(static_cast<long>(n), x);
}

double wasm_spherical_kn(int n, double x) {
    return xsf::sph_bessel_k(static_cast<long>(n), x);
}

// ============================================================
// Struve functions
// ============================================================

double wasm_struve_h(double v, double x) {
    return xsf::cephes::struve_h(v, x);
}

double wasm_struve_l(double v, double x) {
    return xsf::cephes::struve_l(v, x);
}

// ============================================================
// Riemann zeta function
// ============================================================

double wasm_zeta(double x) {
    return xsf::cephes::zeta(x, 1.0);
}

double wasm_zetac(double x) {
    return xsf::cephes::zetac(x);
}

// ============================================================
// Legendre polynomials
// ============================================================

double wasm_legendre_p(int n, double x) {
    return xsf::legendre_p(n, x);
}

// ============================================================
// Statistical distributions
// ============================================================

// Normal distribution
double wasm_ndtr(double x) {
    return xsf::cephes::ndtr(x);
}

double wasm_ndtri(double p) {
    return xsf::cephes::ndtri(p);
}

double wasm_log_ndtr(double x) {
    double t = x * 0.7071067811865476;  // M_SQRT1_2
    if (x < -1.0) {
        return log(xsf::erfcx(-t) / 2) - t * t;
    } else {
        return log1p(-xsf::erfc(t) / 2);
    }
}

// Chi-square distribution
double wasm_chdtr(double df, double x) {
    return xsf::cephes::chdtr(df, x);
}

double wasm_chdtrc(double df, double x) {
    return xsf::cephes::chdtrc(df, x);
}

double wasm_chdtri(double df, double p) {
    return xsf::cephes::chdtri(df, p);
}

// F distribution
double wasm_fdtr(double dfn, double dfd, double x) {
    return xsf::cephes::fdtr(dfn, dfd, x);
}

double wasm_fdtrc(double dfn, double dfd, double x) {
    return xsf::cephes::fdtrc(dfn, dfd, x);
}

double wasm_fdtri(double dfn, double dfd, double p) {
    return xsf::cephes::fdtri(dfn, dfd, p);
}

// Gamma distribution
double wasm_gdtr(double a, double b, double x) {
    return xsf::cephes::gdtr(a, b, x);
}

double wasm_gdtrc(double a, double b, double x) {
    return xsf::cephes::gdtrc(a, b, x);
}

// Poisson distribution
double wasm_pdtr(double k, double m) {
    return xsf::cephes::pdtr(k, m);
}

double wasm_pdtrc(double k, double m) {
    return xsf::cephes::pdtrc(k, m);
}

double wasm_pdtri(int k, double p) {
    return xsf::cephes::pdtri(k, p);
}

// Binomial distribution
double wasm_bdtr(double k, int n, double p) {
    return xsf::cephes::bdtr(k, n, p);
}

double wasm_bdtrc(double k, int n, double p) {
    return xsf::cephes::bdtrc(k, n, p);
}

double wasm_bdtri(double k, int n, double y) {
    return xsf::cephes::bdtri(k, n, y);
}

// Negative binomial distribution
double wasm_nbdtr(int k, int n, double p) {
    return xsf::cephes::nbdtr(k, n, p);
}

double wasm_nbdtrc(int k, int n, double p) {
    return xsf::cephes::nbdtrc(k, n, p);
}

double wasm_nbdtri(int k, int n, double p) {
    return xsf::cephes::nbdtri(k, n, p);
}

// Kolmogorov distribution
double wasm_kolmogorov(double x) {
    return xsf::cephes::kolmogorov(x);
}

double wasm_kolmogi(double p) {
    return xsf::cephes::kolmogi(p);
}

// Smirnov distribution
double wasm_smirnov(int n, double x) {
    return xsf::cephes::smirnov(n, x);
}

double wasm_smirnovi(int n, double p) {
    return xsf::cephes::smirnovi(n, p);
}

// Owen's T function
double wasm_owens_t(double h, double a) {
    return xsf::cephes::owens_t(h, a);
}

} // extern "C"
