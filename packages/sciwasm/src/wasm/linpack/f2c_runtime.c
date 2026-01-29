/*
 * f2c Runtime Library - Functions for f2c-generated LINPACK code
 * Includes complex math functions for c* and z* routines.
 */

#include "f2c.h"
#include <math.h>
#include <string.h>

/* String copy */
int s_copy(char *a, char *b, ftnlen la, ftnlen lb) {
    char *aend = a + la;
    if (la <= lb) {
        while (a < aend) *a++ = *b++;
    } else {
        char *bend = b + lb;
        while (b < bend) *a++ = *b++;
        while (a < aend) *a++ = ' ';
    }
    return 0;
}

/* String compare */
integer s_cmp(char *a, char *b, ftnlen la, ftnlen lb) {
    char *aend = a + la, *bend = b + lb;
    while (a < aend && b < bend) {
        if (*a != *b) return (*a - *b);
        ++a; ++b;
    }
    while (a < aend) { if (*a != ' ') return (*a - ' '); ++a; }
    while (b < bend) { if (*b != ' ') return (' ' - *b); ++b; }
    return 0;
}

/* Absolute value functions */
double d_abs(doublereal *x) { return *x >= 0 ? *x : -*x; }
double r_abs(real *x) { return *x >= 0 ? *x : -*x; }

/* Sign functions */
double d_sign(doublereal *a, doublereal *b) {
    double x = (*a >= 0 ? *a : -*a);
    return (*b >= 0 ? x : -x);
}

double r_sign(real *a, real *b) {
    double x = (*a >= 0 ? *a : -*a);
    return (*b >= 0 ? x : -x);
}

/* Log base 10 */
double d_lg10(doublereal *x) {
    return log10(*x);
}

double r_lg10(real *x) {
    return log10((double)*x);
}

/* Power functions */
double pow_dd(doublereal *ap, doublereal *bp) { return pow(*ap, *bp); }

double pow_di(doublereal *ap, integer *bp) {
    double r = 1, x = *ap;
    integer n = *bp;
    if (n == 0) return 1;
    if (n < 0) { n = -n; x = 1/x; }
    for (;;) {
        if (n & 1) r *= x;
        if ((n >>= 1) == 0) return r;
        x *= x;
    }
}

double pow_ri(real *ap, integer *bp) {
    double r = 1, x = *ap;
    integer n = *bp;
    if (n == 0) return 1;
    if (n < 0) { n = -n; x = 1/x; }
    for (;;) {
        if (n & 1) r *= x;
        if ((n >>= 1) == 0) return r;
        x *= x;
    }
}

integer pow_ii(integer *ap, integer *bp) {
    integer r = 1, x = *ap, n = *bp;
    if (n <= 0) {
        if (n == 0 || x == 1) return 1;
        if (x != -1) return x == 0 ? 1/x : 0;
        n = -n;
    }
    for (;;) {
        if (n & 1) r *= x;
        if ((n >>= 1) == 0) return r;
        x *= x;
    }
}

/* Nearest integer */
integer i_nint(real *x) {
    return (integer)(*x >= 0 ? floor(*x + 0.5) : -floor(0.5 - *x));
}

integer i_dnnt(doublereal *x) {
    return (integer)(*x >= 0 ? floor(*x + 0.5) : -floor(0.5 - *x));
}

/* Machine constants for IEEE 754 double precision */
doublereal d1mach_(integer *i) {
    switch (*i) {
        case 1: return 2.2250738585072014e-308;  /* smallest positive number */
        case 2: return 1.7976931348623157e+308;  /* largest number */
        case 3: return 1.1102230246251565e-16;   /* smallest relative spacing (eps/2) */
        case 4: return 2.2204460492503131e-16;   /* largest relative spacing (eps) */
        case 5: return 0.30102999566398120;      /* log10(2) */
        default: return 0.0;
    }
}

/* Single precision machine constants */
doublereal r1mach_(integer *i) {
    switch (*i) {
        case 1: return 1.1754944e-38;   /* smallest positive number */
        case 2: return 3.4028235e+38;   /* largest number */
        case 3: return 5.9604645e-08;   /* smallest relative spacing (eps/2) */
        case 4: return 1.1920929e-07;   /* largest relative spacing (eps) */
        case 5: return 0.30103;         /* log10(2) */
        default: return 0.0;
    }
}

/* Error handler stub - WASM doesn't have stderr */
int xerror_(char *mesg, integer *nmesg, integer *nerr, integer *level,
            ftnlen mesg_len) {
    (void)mesg; (void)nmesg; (void)nerr; (void)level; (void)mesg_len;
    return 0;
}

/* ========================================
 * Complex math functions for LINPACK
 * ======================================== */

/* Double complex absolute value: |z| = sqrt(r^2 + i^2) */
double z_abs(doublecomplex *z) {
    double ar = z->r >= 0 ? z->r : -z->r;
    double ai = z->i >= 0 ? z->i : -z->i;
    double s, t;
    if (ar == 0 && ai == 0) return 0;
    if (ar > ai) {
        t = ai / ar;
        return ar * sqrt(1.0 + t*t);
    }
    t = ar / ai;
    return ai * sqrt(1.0 + t*t);
}

/* Single complex absolute value */
double c_abs(complex *c) {
    float ar = c->r >= 0 ? c->r : -c->r;
    float ai = c->i >= 0 ? c->i : -c->i;
    float s, t;
    if (ar == 0 && ai == 0) return 0;
    if (ar > ai) {
        t = ai / ar;
        return ar * sqrt(1.0 + t*t);
    }
    t = ar / ai;
    return ai * sqrt(1.0 + t*t);
}

/* Double complex division: r = a / b */
void z_div(doublecomplex *r, doublecomplex *a, doublecomplex *b) {
    double ratio, den;
    double abr, abi;

    abr = b->r >= 0 ? b->r : -b->r;
    abi = b->i >= 0 ? b->i : -b->i;

    if (abr <= abi) {
        ratio = b->r / b->i;
        den = b->i * (1 + ratio*ratio);
        r->r = (a->r * ratio + a->i) / den;
        r->i = (a->i * ratio - a->r) / den;
    } else {
        ratio = b->i / b->r;
        den = b->r * (1 + ratio*ratio);
        r->r = (a->r + a->i * ratio) / den;
        r->i = (a->i - a->r * ratio) / den;
    }
}

/* Single complex division */
void c_div(complex *r, complex *a, complex *b) {
    float ratio, den;
    float abr, abi;

    abr = b->r >= 0 ? b->r : -b->r;
    abi = b->i >= 0 ? b->i : -b->i;

    if (abr <= abi) {
        ratio = b->r / b->i;
        den = b->i * (1 + ratio*ratio);
        r->r = (a->r * ratio + a->i) / den;
        r->i = (a->i * ratio - a->r) / den;
    } else {
        ratio = b->i / b->r;
        den = b->r * (1 + ratio*ratio);
        r->r = (a->r + a->i * ratio) / den;
        r->i = (a->i - a->r * ratio) / den;
    }
}

/* Double complex square root */
void z_sqrt(doublecomplex *r, doublecomplex *z) {
    double mag = z_abs(z);
    double t;

    if (mag == 0) {
        r->r = r->i = 0;
    } else if (z->r > 0) {
        t = sqrt(0.5 * (mag + z->r));
        r->r = t;
        r->i = z->i / (2 * t);
    } else {
        t = sqrt(0.5 * (mag - z->r));
        if (z->i < 0) t = -t;
        r->r = z->i / (2 * t);
        r->i = t;
    }
}

/* Single complex square root */
void c_sqrt(complex *r, complex *c) {
    doublecomplex zc, zr;
    zc.r = c->r;
    zc.i = c->i;
    z_sqrt(&zr, &zc);
    r->r = (float)zr.r;
    r->i = (float)zr.i;
}

/* Double complex conjugate: r = conj(z) */
void d_cnjg(doublecomplex *r, doublecomplex *z) {
    r->r = z->r;
    r->i = -z->i;
}

/* Single complex conjugate */
void r_cnjg(complex *r, complex *c) {
    r->r = c->r;
    r->i = -c->i;
}

/* Imaginary part of double complex */
double d_imag(doublecomplex *z) {
    return z->i;
}

/* Imaginary part of single complex */
double r_imag(complex *c) {
    return (double)c->i;
}

/* Double complex exponential: r = exp(z) = exp(r) * (cos(i) + i*sin(i)) */
void z_exp(doublecomplex *r, doublecomplex *z) {
    double expx = exp(z->r);
    r->r = expx * cos(z->i);
    r->i = expx * sin(z->i);
}

/* Single complex exponential */
void c_exp(complex *r, complex *c) {
    double expx = exp((double)c->r);
    r->r = (float)(expx * cos((double)c->i));
    r->i = (float)(expx * sin((double)c->i));
}

/* Double complex logarithm: r = log(z) */
void z_log(doublecomplex *r, doublecomplex *z) {
    r->r = log(z_abs(z));
    r->i = atan2(z->i, z->r);
}

/* Single complex logarithm */
void c_log(complex *r, complex *c) {
    doublecomplex zc, zr;
    zc.r = c->r;
    zc.i = c->i;
    z_log(&zr, &zc);
    r->r = (float)zr.r;
    r->i = (float)zr.i;
}

/* Double complex cosine */
void z_cos(doublecomplex *r, doublecomplex *z) {
    r->r = cos(z->r) * cosh(z->i);
    r->i = -sin(z->r) * sinh(z->i);
}

/* Double complex sine */
void z_sin(doublecomplex *r, doublecomplex *z) {
    r->r = sin(z->r) * cosh(z->i);
    r->i = cos(z->r) * sinh(z->i);
}

/* ========================================
 * BLAS Level 1 routines (Fortran 90 in modern BLAS)
 * crotg and zrotg are not in ARPACK, so we implement them in C.
 * (dnrm2, snrm2, scnrm2, dznrm2, drotg, srotg come from ARPACK/BLAS)
 * ======================================== */

/* CROTG: Construct Givens plane rotation (complex single precision) */
int crotg_(complex *a, complex *b, real *c, complex *s) {
    real d, f1, f2, g2, h2, p, u, uu, v, vv, w;
    complex fs, gs, r;

    f1 = fabsf(a->r) + fabsf(a->i);  /* |a|_1 */
    g2 = fabsf(b->r) + fabsf(b->i);  /* |b|_1 */

    if (g2 == 0.0f) {
        *c = 1.0f;
        s->r = 0.0f;
        s->i = 0.0f;
        r = *a;
    } else if (f1 == 0.0f) {
        *c = 0.0f;
        /* s = conj(b) / |b| */
        d = sqrtf(b->r * b->r + b->i * b->i);
        s->r = b->r / d;
        s->i = -b->i / d;
        r.r = d;
        r.i = 0.0f;
    } else {
        f2 = a->r * a->r + a->i * a->i;  /* |a|^2 */
        g2 = b->r * b->r + b->i * b->i;  /* |b|^2 */
        h2 = f2 + g2;
        d = sqrtf(h2);
        *c = sqrtf(f2) / d;
        /* s = conj(b) * (a / |a|) / d */
        p = sqrtf(f2);
        s->r = (b->r * a->r + b->i * a->i) / (p * d);
        s->i = (b->r * a->i - b->i * a->r) / (p * d);
        /* r = a * d / |a| */
        r.r = a->r * d / p;
        r.i = a->i * d / p;
    }
    *a = r;
    return 0;
}

/* ZROTG: Construct Givens plane rotation (complex double precision) */
int zrotg_(doublecomplex *a, doublecomplex *b, doublereal *c, doublecomplex *s) {
    doublereal d, f1, f2, g2, h2, p;
    doublecomplex r;

    f1 = fabs(a->r) + fabs(a->i);
    g2 = fabs(b->r) + fabs(b->i);

    if (g2 == 0.0) {
        *c = 1.0;
        s->r = 0.0;
        s->i = 0.0;
        r = *a;
    } else if (f1 == 0.0) {
        *c = 0.0;
        d = sqrt(b->r * b->r + b->i * b->i);
        s->r = b->r / d;
        s->i = -b->i / d;
        r.r = d;
        r.i = 0.0;
    } else {
        f2 = a->r * a->r + a->i * a->i;
        g2 = b->r * b->r + b->i * b->i;
        h2 = f2 + g2;
        d = sqrt(h2);
        *c = sqrt(f2) / d;
        p = sqrt(f2);
        s->r = (b->r * a->r + b->i * a->i) / (p * d);
        s->i = (b->r * a->i - b->i * a->r) / (p * d);
        r.r = a->r * d / p;
        r.i = a->i * d / p;
    }
    *a = r;
    return 0;
}
