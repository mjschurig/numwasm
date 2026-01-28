/**
 * WASM entry points for QUADPACK adaptive quadrature.
 *
 * These thin wrappers expose the QUADPACK dqagse and dqagie routines
 * for calling from JavaScript/TypeScript via Emscripten.
 */

#include "quadpack.h"

void wasm_dqagse(double(*fcn)(double*), double a, double b,
                 double epsabs, double epsrel, int limit,
                 double* result, double* abserr, int* neval, int* ier,
                 double* alist, double* blist, double* rlist, double* elist,
                 int* iord, int* last) {
    dqagse(fcn, a, b, epsabs, epsrel, limit, result, abserr, neval, ier,
           alist, blist, rlist, elist, iord, last);
}

void wasm_dqagie(double(*fcn)(double*), double bound, int inf,
                 double epsabs, double epsrel, int limit,
                 double* result, double* abserr, int* neval, int* ier,
                 double* alist, double* blist, double* rlist, double* elist,
                 int* iord, int* last) {
    dqagie(fcn, bound, inf, epsabs, epsrel, limit, result, abserr, neval, ier,
           alist, blist, rlist, elist, iord, last);
}
