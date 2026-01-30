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
