#!/bin/bash
# LINPACK + BLAS (pure Fortran 77) to C Conversion Script
# Converts LINPACK and required BLAS files using f2c for WebAssembly compilation

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCIWASM_DIR="$(dirname "$SCRIPT_DIR")"
REFERENCE_DIR="$SCIWASM_DIR/reference"
WASM_DIR="$SCIWASM_DIR/src/wasm"
LINPACK_OUT="$WASM_DIR/linpack"

echo "========================================"
echo "LINPACK Fortran to C Conversion"
echo "========================================"
echo ""

# Check for f2c
if ! command -v f2c &> /dev/null; then
    echo "ERROR: f2c is not installed"
    echo "Install with: sudo apt-get install f2c libf2c2-dev"
    exit 1
fi

# Check for LINPACK reference
LINPACK_REF="$REFERENCE_DIR/LINPACK"
if [ ! -d "$LINPACK_REF" ]; then
    echo "ERROR: LINPACK not found at $LINPACK_REF"
    exit 1
fi

# Check for BLAS reference
BLAS_REF="$REFERENCE_DIR/BLAS/SRC"
if [ ! -d "$BLAS_REF" ]; then
    echo "ERROR: BLAS not found at $BLAS_REF"
    exit 1
fi

# Clean previous output
if [ -d "$LINPACK_OUT" ]; then
    echo "Cleaning previous output..."
    rm -rf "$LINPACK_OUT"
fi

# Create output directory structure
echo "Creating directory structure..."
mkdir -p "$LINPACK_OUT/include"

# Copy LINPACK files
echo ""
echo "Copying LINPACK files..."
cp "$LINPACK_REF"/*.f "$LINPACK_OUT/"
linpack_count=$(ls "$LINPACK_OUT"/*.f 2>/dev/null | wc -l)
echo "  LINPACK files: $linpack_count"

# Copy only needed BLAS files (dependencies for LINPACK)
echo ""
echo "Copying required BLAS files..."

# Fortran 77 BLAS routines from reference/BLAS/SRC
BLAS_NEEDED_F77="caxpy ccopy cdotc cdotu cscal csrot csscal cswap \
                 dasum daxpy dcopy ddot drot dscal dswap dzasum \
                 icamax idamax isamax izamax \
                 sasum saxpy scasum scopy sdot srot sscal sswap \
                 zaxpy zcopy zdotc zdotu zdrot zdscal zscal zswap \
                 scabs1 dcabs1"

# Fortran 90 BLAS routines in reference/BLAS/SRC - use ARPACK's F77 versions instead
# ARPACK has F77 versions: dnrm2, snrm2, scnrm2, dznrm2, drotg, srotg
BLAS_NEEDED_FROM_ARPACK="dnrm2 snrm2 scnrm2 dznrm2 drotg srotg"
# These need C implementation: crotg, zrotg (not in ARPACK)

ARPACK_BLAS_REF="$REFERENCE_DIR/ARPACK/BLAS"

blas_count=0
for routine in $BLAS_NEEDED_F77; do
    if [ -f "$BLAS_REF/${routine}.f" ]; then
        cp "$BLAS_REF/${routine}.f" "$LINPACK_OUT/"
        ((blas_count++)) || true
    else
        echo "  Warning: BLAS routine $routine not found in reference/BLAS/SRC"
    fi
done
echo "  BLAS F77 files from BLAS/SRC: $blas_count"

# Copy F77 versions from ARPACK for F90 routines in modern BLAS
arpack_count=0
for routine in $BLAS_NEEDED_FROM_ARPACK; do
    if [ -f "$ARPACK_BLAS_REF/${routine}.f" ]; then
        cp "$ARPACK_BLAS_REF/${routine}.f" "$LINPACK_OUT/"
        ((arpack_count++)) || true
    else
        echo "  Warning: BLAS routine $routine not found in ARPACK/BLAS"
    fi
done
echo "  BLAS F77 files from ARPACK/BLAS: $arpack_count"

# Create f2c.h header (same as QUADPACK)
echo ""
echo "Creating f2c.h..."
cat > "$LINPACK_OUT/include/f2c.h" << 'HEADEREOF'
/* f2c.h  --  Standard Fortran to C header file */
/* Modified for WebAssembly/Emscripten compatibility */

#ifndef F2C_INCLUDE
#define F2C_INCLUDE

#include <stdint.h>

typedef int32_t integer;
typedef uint32_t uinteger;
typedef char *address;
typedef int16_t shortint;
typedef float real;
typedef double doublereal;
typedef struct { real r, i; } complex;
typedef struct { doublereal r, i; } doublecomplex;
typedef int32_t logical;
typedef int16_t shortlogical;
typedef char logical1;
typedef char integer1;

#define TRUE_ (1)
#define FALSE_ (0)

#ifndef Extern
#define Extern extern
#endif

#ifdef f2c_i2
typedef int16_t flag;
typedef int16_t ftnlen;
typedef int16_t ftnint;
#else
typedef int32_t flag;
typedef int32_t ftnlen;
typedef int32_t ftnint;
#endif

typedef struct { flag cierr; ftnint ciunit; flag ciend; char *cifmt; ftnint cirec; } cilist;
typedef struct { flag icierr; char *iciunit; flag iciend; char *icifmt; ftnint icirlen; ftnint icirnum; } icilist;
typedef struct { flag oerr; ftnint ounit; char *ofnm; ftnlen ofnmlen; char *osta; char *oacc; char *ofm; ftnint orl; char *oblnk; } olist;
typedef struct { flag cerr; ftnint cunit; char *csta; } cllist;
typedef struct { flag aerr; ftnint aunit; } alist;
typedef struct {
    flag inerr; ftnint inunit; char *infile; ftnlen infilen;
    ftnint *inex, *inopen, *innum, *innamed; char *inname; ftnlen innamlen;
    char *inacc; ftnlen inacclen; char *inseq; ftnlen inseqlen;
    char *indir; ftnlen indirlen; char *infmt; ftnlen infmtlen;
    char *inform; ftnint informlen; char *inunf; ftnlen inunflen;
    ftnint *inrecl, *innrec; char *inblank; ftnlen inblanklen;
} inlist;

#define VOID void

union Multitype { integer1 g; shortint h; integer i; real r; doublereal d; complex c; doublecomplex z; };
typedef union Multitype Multitype;
struct Vardesc { char *name; char *addr; ftnlen *dims; int type; };
typedef struct Vardesc Vardesc;
struct Namelist { char *name; Vardesc **vars; int nvars; };
typedef struct Namelist Namelist;

#define abs(x) ((x) >= 0 ? (x) : -(x))
#define dabs(x) (doublereal)abs(x)
#define min(a,b) ((a) <= (b) ? (a) : (b))
#define max(a,b) ((a) >= (b) ? (a) : (b))
#define dmin(a,b) (doublereal)min(a,b)
#define dmax(a,b) (doublereal)max(a,b)
#define bit_test(a,b) ((a) >> (b) & 1)
#define bit_clear(a,b) ((a) & ~((uinteger)1 << (b)))
#define bit_set(a,b) ((a) | ((uinteger)1 << (b)))

#define F2C_proc_par_types 1
#ifdef __cplusplus
typedef int (*U_fp)(...);
typedef shortint (*J_fp)(...);
typedef integer (*I_fp)(...);
typedef real (*R_fp)(...);
typedef doublereal (*D_fp)(...), (*E_fp)(...);
typedef VOID (*C_fp)(...);
typedef VOID (*Z_fp)(...);
typedef logical (*L_fp)(...);
typedef shortlogical (*K_fp)(...);
typedef VOID (*H_fp)(...);
typedef int (*S_fp)(...);
#else
typedef int (*U_fp)();
typedef shortint (*J_fp)();
typedef integer (*I_fp)();
typedef real (*R_fp)();
typedef doublereal (*D_fp)(), (*E_fp)();
typedef VOID (*C_fp)();
typedef VOID (*Z_fp)();
typedef logical (*L_fp)();
typedef shortlogical (*K_fp)();
typedef VOID (*H_fp)();
typedef int (*S_fp)();
#endif

typedef VOID C_f;
typedef VOID H_f;
typedef VOID Z_f;
typedef doublereal E_f;

/* f2c runtime function declarations */
extern int s_copy(char *a, char *b, ftnlen la, ftnlen lb);
extern integer s_cmp(char *a, char *b, ftnlen la, ftnlen lb);
extern double d_abs(doublereal *x);
extern double r_abs(real *x);
extern double d_sign(doublereal *a, doublereal *b);
extern double r_sign(real *a, real *b);
extern double pow_dd(doublereal *ap, doublereal *bp);
extern double pow_di(doublereal *ap, integer *bp);
extern double pow_ri(real *ap, integer *bp);
extern integer pow_ii(integer *ap, integer *bp);
extern integer i_nint(real *x);
extern integer i_dnnt(doublereal *x);
extern double d_lg10(doublereal *x);
extern double r_lg10(real *x);
extern doublereal d1mach_(integer *i);
extern doublereal r1mach_(integer *i);
extern int xerror_(char *mesg, integer *nmesg, integer *nerr, integer *level, ftnlen mesg_len);

/* Complex math functions for LINPACK */
extern double z_abs(doublecomplex *z);
extern double c_abs(complex *c);
extern void z_div(doublecomplex *r, doublecomplex *a, doublecomplex *b);
extern void c_div(complex *r, complex *a, complex *b);
extern void z_sqrt(doublecomplex *r, doublecomplex *z);
extern void c_sqrt(complex *r, complex *c);
extern void d_cnjg(doublecomplex *r, doublecomplex *z);
extern void r_cnjg(complex *r, complex *c);
extern double d_imag(doublecomplex *z);
extern double r_imag(complex *c);

#ifndef Skip_f2c_Undefs
#undef cray
#undef gcos
#undef mc68010
#undef mc68020
#undef mips
#undef pdp11
#undef sgi
#undef sparc
#undef sun
#undef sun2
#undef sun3
#undef sun4
#undef u370
#undef u3b
#undef u3b2
#undef u3b5
#undef unix
#undef vax
#endif

#endif /* F2C_INCLUDE */
HEADEREOF

# Convert Fortran files to C
echo ""
echo "Converting Fortran files to C..."

cd "$LINPACK_OUT"
count=0
errors=0

for f in *.f; do
    if [ -f "$f" ]; then
        # f2c flags:
        # -a: Make local variables automatic (not static) - important for reentrancy
        # -A: Produce ANSI C (function prototypes)
        if f2c -a -A "$f" 2>/dev/null; then
            ((count++)) || true
        else
            echo "  Warning: Failed to convert $f"
            ((errors++)) || true
        fi
    fi
done

# Clean up Fortran files (keep only C files)
rm -f *.f

echo "  Converted: $count, Errors: $errors"

# Create f2c runtime library with complex math support
echo ""
echo "Creating f2c runtime library..."
cat > "$LINPACK_OUT/f2c_runtime.c" << 'EOF'
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
EOF

# Summary
echo ""
echo "========================================"
echo "Conversion Complete!"
echo "========================================"
echo ""
echo "Output directory: $LINPACK_OUT"
echo ""
echo "Files created:"
echo "  C files:  $(ls "$LINPACK_OUT"/*.c 2>/dev/null | wc -l)"
echo "  include/: f2c.h"
echo "  f2c_runtime.c"
echo ""
echo "Next steps:"
echo "  1. Run: bash scripts/build-linpack.sh"
