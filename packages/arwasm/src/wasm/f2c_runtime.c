/*
 * f2c Runtime Library - Essential functions for f2c-generated code
 * Standalone implementations for WebAssembly compatibility.
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

/* Complex absolute value */
double z_abs(doublecomplex *z) { return sqrt(z->r * z->r + z->i * z->i); }
double c_abs(complex *z) { return sqrt(z->r * z->r + z->i * z->i); }

/* Complex imaginary part extraction */
double r_imag(complex *z) { return (double)z->i; }
double d_imag(doublecomplex *z) { return z->i; }

/* Complex division */
void z_div(doublecomplex *c, doublecomplex *a, doublecomplex *b) {
    double abr = b->r >= 0 ? b->r : -b->r;
    double abi = b->i >= 0 ? b->i : -b->i;
    double ratio, den;
    if (abr <= abi) {
        if (abi == 0) { c->r = c->i = abr / abi; return; }
        ratio = b->r / b->i;
        den = b->i * (1 + ratio * ratio);
        c->r = (a->r * ratio + a->i) / den;
        c->i = (a->i * ratio - a->r) / den;
    } else {
        ratio = b->i / b->r;
        den = b->r * (1 + ratio * ratio);
        c->r = (a->r + a->i * ratio) / den;
        c->i = (a->i - a->r * ratio) / den;
    }
}

void c_div(complex *c, complex *a, complex *b) {
    float abr = b->r >= 0 ? b->r : -b->r;
    float abi = b->i >= 0 ? b->i : -b->i;
    float ratio, den;
    if (abr <= abi) {
        if (abi == 0) { c->r = c->i = abr / abi; return; }
        ratio = b->r / b->i;
        den = b->i * (1 + ratio * ratio);
        c->r = (a->r * ratio + a->i) / den;
        c->i = (a->i * ratio - a->r) / den;
    } else {
        ratio = b->i / b->r;
        den = b->r * (1 + ratio * ratio);
        c->r = (a->r + a->i * ratio) / den;
        c->i = (a->i - a->r * ratio) / den;
    }
}

/* Complex square root */
void z_sqrt(doublecomplex *r, doublecomplex *z) {
    double mag = z_abs(z), t;
    if (mag == 0) { r->r = r->i = 0; return; }
    if (z->r > 0) {
        r->r = sqrt(0.5 * (mag + z->r));
        r->i = 0.5 * z->i / r->r;
    } else {
        t = sqrt(0.5 * (mag - z->r));
        if (z->i < 0) t = -t;
        r->i = t;
        r->r = 0.5 * z->i / t;
    }
}

void c_sqrt(complex *r, complex *z) {
    float mag = c_abs(z), t;
    if (mag == 0) { r->r = r->i = 0; return; }
    if (z->r > 0) {
        r->r = sqrtf(0.5f * (mag + z->r));
        r->i = 0.5f * z->i / r->r;
    } else {
        t = sqrtf(0.5f * (mag - z->r));
        if (z->i < 0) t = -t;
        r->i = t;
        r->r = 0.5f * z->i / t;
    }
}

/* Complex conjugate */
void d_cnjg(doublecomplex *r, doublecomplex *z) { r->r = z->r; r->i = -z->i; }
void r_cnjg(complex *r, complex *z) { r->r = z->r; r->i = -z->i; }

/* Complex exponential: e^(a+bi) = e^a * (cos(b) + i*sin(b)) */
void z_exp(doublecomplex *r, doublecomplex *z) {
    double ea = exp(z->r);
    r->r = ea * cos(z->i);
    r->i = ea * sin(z->i);
}

void c_exp(complex *r, complex *z) {
    float ea = expf(z->r);
    r->r = ea * cosf(z->i);
    r->i = ea * sinf(z->i);
}

/* Fortran I/O stubs - not supported in WASM */
int s_wsfe(cilist *a) { (void)a; return 0; }
int e_wsfe(void) { return 0; }
int do_fio(ftnint *number, char *ptr, ftnlen len) {
    (void)number; (void)ptr; (void)len; return 0;
}
int s_wsle(cilist *a) { (void)a; return 0; }
int e_wsle(void) { return 0; }
int s_stop(char *s, ftnlen n) { (void)s; (void)n; return 0; }

/* LAPACK error handler stub */
int xerbla_(char *srname, integer *info, ftnlen srname_len) {
    (void)srname; (void)info; (void)srname_len;
    return 0;
}
