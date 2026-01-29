#!/bin/bash
# QUADPACK (pure Fortran 77) to C Conversion Script
# Converts all QUADPACK files using f2c for WebAssembly compilation with Emscripten

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCIWASM_DIR="$(dirname "$SCRIPT_DIR")"
REFERENCE_DIR="$SCIWASM_DIR/reference"
WASM_DIR="$SCIWASM_DIR/src/wasm"
QUADPACK_OUT="$WASM_DIR/quadpack"

echo "========================================"
echo "QUADPACK Fortran to C Conversion"
echo "========================================"
echo ""

# Check for f2c
if ! command -v f2c &> /dev/null; then
    echo "ERROR: f2c is not installed"
    echo "Install with: sudo apt-get install f2c libf2c2-dev"
    exit 1
fi

# Check for QUADPACK reference
QUADPACK_REF="$REFERENCE_DIR/QUADPACK"
if [ ! -d "$QUADPACK_REF" ]; then
    echo "ERROR: QUADPACK not found at $QUADPACK_REF"
    exit 1
fi

# Clean previous output
if [ -d "$QUADPACK_OUT" ]; then
    echo "Cleaning previous output..."
    rm -rf "$QUADPACK_OUT"
fi

# Create output directory structure
echo "Creating directory structure..."
mkdir -p "$QUADPACK_OUT/include"

# Create standalone f2c.h for WebAssembly (uses stdint.h instead of inttypes.h)
echo "Creating f2c.h..."
cat > "$QUADPACK_OUT/include/f2c.h" << 'HEADEREOF'
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

# Convert Fortran files
echo ""
echo "Converting QUADPACK Fortran files..."

# Copy all QUADPACK Fortran files to output directory
cp "$QUADPACK_REF"/*.f "$QUADPACK_OUT/"

# Also copy required LINPACK files (tridiagonal solvers needed by qc25f)
LINPACK_REF="$REFERENCE_DIR/LINPACK"
if [ -d "$LINPACK_REF" ]; then
    echo "Adding LINPACK dependencies (dgtsl, sgtsl)..."
    cp "$LINPACK_REF/dgtsl.f" "$QUADPACK_OUT/" 2>/dev/null || true
    cp "$LINPACK_REF/sgtsl.f" "$QUADPACK_OUT/" 2>/dev/null || true
fi

# Convert each file with f2c
cd "$QUADPACK_OUT"
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

# Create f2c runtime library (minimal for QUADPACK - no complex math needed)
echo ""
echo "Creating f2c runtime library..."
cat > "$QUADPACK_OUT/f2c_runtime.c" << 'EOF'
/*
 * f2c Runtime Library - Essential functions for f2c-generated QUADPACK code
 * Minimal implementation for WebAssembly compatibility.
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
EOF

# Summary
echo ""
echo "========================================"
echo "Conversion Complete!"
echo "========================================"
echo ""
echo "Output directory: $QUADPACK_OUT"
echo ""
echo "Files created:"
echo "  C files:  $(ls "$QUADPACK_OUT"/*.c 2>/dev/null | wc -l)"
echo "  include/: f2c.h"
echo "  f2c_runtime.c"
echo ""
echo "Next steps:"
echo "  1. Run: bash scripts/build-quadpack.sh"
