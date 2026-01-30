/*
 * f2c Runtime Library - Essential functions for f2c-generated ODE solver code
 * Minimal implementation for WebAssembly compatibility.
 */

/* Include system headers BEFORE f2c.h to avoid macro conflicts */
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include "include/f2c.h"

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

/* String length */
integer i_len(char *s, ftnlen n) {
    (void)s;
    return n;
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

/* Modulo */
doublereal d_mod(doublereal *x, doublereal *y) {
    return fmod(*x, *y);
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

/* Formatted I/O stubs - ODE solvers use these for error messages */
int s_wsfe(cilist *a) { (void)a; return 0; }
int e_wsfe(void) { return 0; }
int do_fio(integer *n, char *ptr, ftnlen len) {
    (void)n; (void)ptr; (void)len;
    return 0;
}

/* List-directed I/O stubs */
int s_wsle(cilist *a) { (void)a; return 0; }
int e_wsle(void) { return 0; }
int do_lio(integer *type, integer *num, char *ptr, ftnlen len) {
    (void)type; (void)num; (void)ptr; (void)len;
    return 0;
}

/* Stop statement - just return in WASM */
int s_stop(char *s, ftnlen n) {
    (void)s; (void)n;
    return 0;
}
