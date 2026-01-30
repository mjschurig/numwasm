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
