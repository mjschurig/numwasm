#!/bin/bash
# Original ARPACK (pure Fortran 77) to C Conversion Script
# Converts all ARPACK files including bundled LAPACK and BLAS using f2c
# for WebAssembly compilation with Emscripten

set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SCIWASM_DIR="$(dirname "$SCRIPT_DIR")"
REFERENCE_DIR="$SCIWASM_DIR/reference"
WASM_DIR="$SCIWASM_DIR/src/wasm"
ARPACK_OUT="$WASM_DIR/arpack"

echo "========================================"
echo "Original ARPACK Fortran to C Conversion"
echo "========================================"
echo ""

# Check for f2c
if ! command -v f2c &> /dev/null; then
    echo "ERROR: f2c is not installed"
    echo "Install with: sudo apt-get install f2c libf2c2-dev"
    exit 1
fi

# Check for ARPACK reference (uppercase - original ARPACK)
ARPACK_REF="$REFERENCE_DIR/ARPACK"
if [ ! -d "$ARPACK_REF" ]; then
    echo "ERROR: Original ARPACK not found at $ARPACK_REF"
    exit 1
fi

# Clean previous output
if [ -d "$ARPACK_OUT" ]; then
    echo "Cleaning previous output..."
    rm -rf "$ARPACK_OUT"
fi

# Create output directory structure
echo "Creating directory structure..."
mkdir -p "$ARPACK_OUT/SRC"
mkdir -p "$ARPACK_OUT/UTIL"
mkdir -p "$ARPACK_OUT/LAPACK"
mkdir -p "$ARPACK_OUT/BLAS"
mkdir -p "$ARPACK_OUT/include"

# Create standalone f2c.h for WebAssembly (uses stdint.h instead of inttypes.h)
echo "Creating f2c.h..."
cat > "$ARPACK_OUT/include/f2c.h" << 'HEADEREOF'
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
extern double d_sign(doublereal *a, doublereal *b);
extern double r_sign(real *a, real *b);
extern double pow_dd(doublereal *ap, doublereal *bp);
extern double pow_di(doublereal *ap, integer *bp);
extern double pow_ri(real *ap, integer *bp);
extern integer pow_ii(integer *ap, integer *bp);
extern integer i_nint(real *x);
extern integer i_dnnt(doublereal *x);
extern void z_div(doublecomplex *c, doublecomplex *a, doublecomplex *b);
extern void c_div(complex *c, complex *a, complex *b);
extern double z_abs(doublecomplex *z);
extern double c_abs(complex *z);
extern double r_imag(complex *z);
extern double d_imag(doublecomplex *z);
extern void z_sqrt(doublecomplex *r, doublecomplex *z);
extern void c_sqrt(complex *r, complex *z);
extern void d_cnjg(doublecomplex *r, doublecomplex *z);
extern void r_cnjg(complex *r, complex *z);

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

# Function to convert Fortran files in a directory
convert_directory() {
    local src_dir="$1"
    local out_dir="$2"
    local name="$3"
    local include_dir="$4"  # Optional directory with include files

    echo ""
    echo "Converting $name..."

    # Count source files
    local src_count=$(ls "$src_dir"/*.f 2>/dev/null | wc -l)
    echo "  Source files: $src_count"

    if [ "$src_count" -eq 0 ]; then
        echo "  No Fortran files found in $src_dir"
        return
    fi

    # Copy Fortran files to output directory
    cp "$src_dir"/*.f "$out_dir/"

    # Copy include files if specified (needed for ARPACK SRC files)
    if [ -n "$include_dir" ] && [ -d "$include_dir" ]; then
        cp "$include_dir"/*.h "$out_dir/" 2>/dev/null || true
    fi

    # Convert each file with f2c
    cd "$out_dir"
    local count=0
    local errors=0

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

    # Clean up Fortran files and include files (keep only C files)
    rm -f *.f *.h

    echo "  Converted: $count, Errors: $errors"
}

# Convert all directories (order doesn't matter for conversion, only for compilation)
# BLAS -> LAPACK -> UTIL -> SRC (dependency order for compilation)

convert_directory "$ARPACK_REF/BLAS" "$ARPACK_OUT/BLAS" "BLAS (76 files)"
convert_directory "$ARPACK_REF/LAPACK" "$ARPACK_OUT/LAPACK" "LAPACK (162 files)"
convert_directory "$ARPACK_REF/UTIL" "$ARPACK_OUT/UTIL" "UTIL (14 files)"
# SRC needs include files (debug.h, stat.h) for COMMON block declarations
convert_directory "$ARPACK_REF/SRC" "$ARPACK_OUT/SRC" "SRC (70 files)" "$ARPACK_REF/SRC"

# Post-process: Fix COMMON block definitions to be extern declarations
# f2c generates full struct definitions in every file that uses a COMMON block
# We need to make them extern since we define them once in arpack_stubs.c
echo ""
echo "Post-processing COMMON block declarations..."
for cfile in "$ARPACK_OUT/SRC"/*.c "$ARPACK_OUT/UTIL"/*.c; do
    if [ -f "$cfile" ]; then
        # Replace struct definitions with extern declarations for debug_ and timing_
        # Pattern: "struct { ... } debug_;" -> "extern struct { ... } debug_;"
        sed -i 's/^struct {$/extern struct {/' "$cfile"
    fi
done
echo "  Done."

# Create arpack_config.h
echo ""
echo "Creating configuration header..."
cat > "$ARPACK_OUT/include/arpack_config.h" << 'EOF'
#ifndef ARPACK_CONFIG_H
#define ARPACK_CONFIG_H

/*
 * ARPACK configuration for WebAssembly/Emscripten build
 * Converted from original ARPACK (pure Fortran 77) using f2c
 */

/* f2c type definitions */
#include "f2c.h"

/* Debug COMMON block - controls debug output levels
 * logfil: logical unit for output (6 = stdout)
 * ndigit: number of digits for output (-3 = suppress)
 * mgetv0, msaupd, etc.: debug levels for each routine (0 = none)
 */
extern struct {
    integer logfil, ndigit, mgetv0;
    integer msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd;
    integer mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd;
    integer mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd;
} debug_;

/* Timing COMMON block - performance statistics
 * nopx, nbx: operation counts
 * t*: timing for each phase
 */
extern struct {
    integer nopx, nbx, nrorth, nitref, nrstrt;
    real tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv;
    real tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv;
    real tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv;
    real tmvopx, tmvbx, tgetv0, titref, trvec;
} timing_;

#endif /* ARPACK_CONFIG_H */
EOF

# Create stubs for timing and I/O functions
echo "Creating stub functions..."
cat > "$ARPACK_OUT/arpack_stubs.c" << 'EOF'
/*
 * ARPACK stub functions for WebAssembly
 * These replace system-dependent functions that don't work in WASM
 */

#include "f2c.h"

/* ============================================================
 * Timing stubs - WASM doesn't have etime()
 * These are called by ARPACK to measure performance
 * ============================================================ */

int second_(real *t) {
    *t = 0.0f;
    return 0;
}

/* Double precision version used by some LAPACK routines */
int dsecnd_(doublereal *t) {
    *t = 0.0;
    return 0;
}

/* Alternative name used in some ARPACK versions */
int arscnd_(real *t) {
    *t = 0.0f;
    return 0;
}

/* ============================================================
 * COMMON block definitions
 * These are shared between all ARPACK routines
 * ============================================================ */

/* Debug COMMON block - initialized to suppress debug output */
struct {
    integer logfil, ndigit, mgetv0;
    integer msaupd, msaup2, msaitr, mseigt, msapps, msgets, mseupd;
    integer mnaupd, mnaup2, mnaitr, mneigh, mnapps, mngets, mneupd;
    integer mcaupd, mcaup2, mcaitr, mceigh, mcapps, mcgets, mceupd;
} debug_ = {
    6,    /* logfil: stdout */
    -3,   /* ndigit: suppress output */
    0,    /* mgetv0 */
    0, 0, 0, 0, 0, 0, 0,  /* symmetric routines */
    0, 0, 0, 0, 0, 0, 0,  /* non-symmetric routines */
    0, 0, 0, 0, 0, 0, 0   /* complex routines */
};

/* Timing COMMON block - all zeros (timing disabled) */
struct {
    integer nopx, nbx, nrorth, nitref, nrstrt;
    real tsaupd, tsaup2, tsaitr, tseigt, tsgets, tsapps, tsconv;
    real tnaupd, tnaup2, tnaitr, tneigh, tngets, tnapps, tnconv;
    real tcaupd, tcaup2, tcaitr, tceigh, tcgets, tcapps, tcconv;
    real tmvopx, tmvbx, tgetv0, titref, trvec;
} timing_ = {0};

/* ============================================================
 * Output stubs - these would normally write to files
 * In WASM we just ignore the output
 * ============================================================ */

/* Integer vector output */
int ivout_(integer *lout, integer *n, integer *ix, integer *idigit,
           char *ifmt, ftnlen ifmt_len) {
    (void)lout; (void)n; (void)ix; (void)idigit; (void)ifmt; (void)ifmt_len;
    return 0;
}

/* Double precision vector output */
int dvout_(integer *lout, integer *n, doublereal *sx, integer *idigit,
           char *ifmt, ftnlen ifmt_len) {
    (void)lout; (void)n; (void)sx; (void)idigit; (void)ifmt; (void)ifmt_len;
    return 0;
}

/* Single precision vector output */
int svout_(integer *lout, integer *n, real *sx, integer *idigit,
           char *ifmt, ftnlen ifmt_len) {
    (void)lout; (void)n; (void)sx; (void)idigit; (void)ifmt; (void)ifmt_len;
    return 0;
}

/* Complex single precision vector output */
int cvout_(integer *lout, integer *n, complex *cx, integer *idigit,
           char *ifmt, ftnlen ifmt_len) {
    (void)lout; (void)n; (void)cx; (void)idigit; (void)ifmt; (void)ifmt_len;
    return 0;
}

/* Complex double precision vector output */
int zvout_(integer *lout, integer *n, doublecomplex *cx, integer *idigit,
           char *ifmt, ftnlen ifmt_len) {
    (void)lout; (void)n; (void)cx; (void)idigit; (void)ifmt; (void)ifmt_len;
    return 0;
}

/* Double precision matrix output */
int dmout_(integer *lout, integer *m, integer *n, doublereal *a,
           integer *lda, integer *idigit, char *ifmt, ftnlen ifmt_len) {
    (void)lout; (void)m; (void)n; (void)a; (void)lda; (void)idigit;
    (void)ifmt; (void)ifmt_len;
    return 0;
}

/* Single precision matrix output */
int smout_(integer *lout, integer *m, integer *n, real *a,
           integer *lda, integer *idigit, char *ifmt, ftnlen ifmt_len) {
    (void)lout; (void)m; (void)n; (void)a; (void)lda; (void)idigit;
    (void)ifmt; (void)ifmt_len;
    return 0;
}

/* Complex single precision matrix output */
int cmout_(integer *lout, integer *m, integer *n, complex *a,
           integer *lda, integer *idigit, char *ifmt, ftnlen ifmt_len) {
    (void)lout; (void)m; (void)n; (void)a; (void)lda; (void)idigit;
    (void)ifmt; (void)ifmt_len;
    return 0;
}

/* Complex double precision matrix output */
int zmout_(integer *lout, integer *m, integer *n, doublecomplex *a,
           integer *lda, integer *idigit, char *ifmt, ftnlen ifmt_len) {
    (void)lout; (void)m; (void)n; (void)a; (void)lda; (void)idigit;
    (void)ifmt; (void)ifmt_len;
    return 0;
}
EOF

# Create f2c runtime library
echo "Creating f2c runtime library..."
cat > "$ARPACK_OUT/f2c_runtime.c" << 'EOF'
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
EOF

# Summary
echo ""
echo "========================================"
echo "Conversion Complete!"
echo "========================================"
echo ""
echo "Output directory: $ARPACK_OUT"
echo ""
echo "Files created:"
echo "  BLAS/:    $(ls "$ARPACK_OUT/BLAS"/*.c 2>/dev/null | wc -l) C files"
echo "  LAPACK/:  $(ls "$ARPACK_OUT/LAPACK"/*.c 2>/dev/null | wc -l) C files"
echo "  UTIL/:    $(ls "$ARPACK_OUT/UTIL"/*.c 2>/dev/null | wc -l) C files"
echo "  SRC/:     $(ls "$ARPACK_OUT/SRC"/*.c 2>/dev/null | wc -l) C files"
echo "  include/: f2c.h, arpack_config.h"
echo "  arpack_stubs.c, f2c_runtime.c"
echo ""
echo "Next steps:"
echo "  1. Ensure build-wasm.sh includes ARPACK compilation"
echo "  2. Run: npm run build:wasm"
