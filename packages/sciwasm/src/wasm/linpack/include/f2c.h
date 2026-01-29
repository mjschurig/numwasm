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
