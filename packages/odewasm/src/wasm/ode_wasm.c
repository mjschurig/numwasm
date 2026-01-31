/*
 * ODE WASM Wrapper
 *
 * Provides clean C entry points for ODE solvers compiled from Fortran via f2c:
 * - Hairer solvers: DOPRI5, DOP853, RADAU5
 * - Netlib solvers: RKF45, DVERK, ODE, VODE, ZVODE, EPSODE, VODPK, RKSUITE, RKC
 *
 * The key challenge is bridging JavaScript callbacks to Fortran subroutine conventions.
 * We use global function pointers that are set before calling the solver.
 */

/* Include system headers BEFORE f2c.h to avoid macro conflicts */
#include <stdlib.h>
#include <string.h>
#include "include/f2c.h"

#ifdef __EMSCRIPTEN__
#include <emscripten.h>
#define EXPORT EMSCRIPTEN_KEEPALIVE
#else
#define EXPORT
#endif

/* ============================================
 * External declarations for f2c-generated solvers
 * ============================================ */

/* DOPRI5: Dormand-Prince 5(4) explicit Runge-Kutta */
extern int dopri5_(integer *n, U_fp fcn, doublereal *x, doublereal *y,
                   doublereal *xend, doublereal *rtol, doublereal *atol,
                   integer *itol, U_fp solout, integer *iout,
                   doublereal *work, integer *lwork, integer *iwork,
                   integer *liwork, doublereal *rpar, integer *ipar,
                   integer *idid);

/* DOP853: Dormand-Prince 8(5,3) explicit Runge-Kutta */
extern int dop853_(integer *n, U_fp fcn, doublereal *x, doublereal *y,
                   doublereal *xend, doublereal *rtol, doublereal *atol,
                   integer *itol, U_fp solout, integer *iout,
                   doublereal *work, integer *lwork, integer *iwork,
                   integer *liwork, doublereal *rpar, integer *ipar,
                   integer *idid);

/* RADAU5: Implicit Runge-Kutta for stiff systems */
extern int radau5_(integer *n, U_fp fcn, doublereal *x, doublereal *y,
                   doublereal *xend, doublereal *h__, doublereal *rtol,
                   doublereal *atol, integer *itol, U_fp jac, integer *ijac,
                   integer *mljac, integer *mujac, U_fp mas, integer *imas,
                   integer *mlmas, integer *mumas, U_fp solout, integer *iout,
                   doublereal *work, integer *lwork, integer *iwork,
                   integer *liwork, doublereal *rpar, integer *ipar,
                   integer *idid);

/* Dense output interpolation functions */
extern doublereal contd5_(integer *ii, doublereal *x, doublereal *con,
                          integer *icomp, integer *nd);
extern doublereal contd8_(integer *ii, doublereal *x, doublereal *con,
                          integer *icomp, integer *nd);
extern doublereal contr5_(integer *ii, doublereal *s, doublereal *cont,
                          integer *lrc);

/* ============================================
 * External declarations for Netlib solvers
 * ============================================ */

/* RKF45: Runge-Kutta-Fehlberg 4(5) */
extern int rkf45_(U_fp f, integer *neqn, doublereal *y, doublereal *t,
                  doublereal *tout, doublereal *relerr, doublereal *abserr,
                  integer *iflag, doublereal *work, integer *iwork);

/* DVERK: Verner 6(5) Runge-Kutta */
extern int dverk_(integer *n, S_fp fcn, doublereal *x, doublereal *y,
                  doublereal *xend, doublereal *tol, integer *ind,
                  doublereal *c__, integer *nw, doublereal *w);

/* ODE: Adams-Bashforth-Moulton */
extern int ode_(U_fp f, integer *neqn, doublereal *y, doublereal *t,
                doublereal *tout, doublereal *relerr, doublereal *abserr,
                integer *iflag, doublereal *work, integer *iwork);

/* VODE: Variable-coefficient ODE solver */
extern int dvode_(S_fp f, integer *neq, doublereal *y, doublereal *t,
                  doublereal *tout, integer *itol, doublereal *rtol,
                  doublereal *atol, integer *itask, integer *istate,
                  integer *iopt, doublereal *rwork, integer *lrw,
                  integer *iwork, integer *liw, U_fp jac, integer *mf,
                  doublereal *rpar, integer *ipar);

/* ZVODE: VODE for complex-valued ODEs */
extern int zvode_(S_fp f, integer *neq, doublecomplex *y, doublereal *t,
                  doublereal *tout, integer *itol, doublereal *rtol,
                  doublereal *atol, integer *itask, integer *istate,
                  integer *iopt, doublecomplex *zwork, integer *lzw,
                  doublereal *rwork, integer *lrw, integer *iwork,
                  integer *liw, U_fp jac, integer *mf,
                  doublereal *rpar, integer *ipar);

/* EPSODE: Episode package for stiff ODEs */
extern int epsode_(integer *n, doublereal *t0, doublereal *hmax,
                   doublereal *y, doublereal *tout, doublereal *eps,
                   integer *mf, integer *index__, doublereal *work);

/* VODPK: VODE with Krylov methods */
extern int dvodpk_(S_fp f, integer *neq, doublereal *y, doublereal *t,
                   doublereal *tout, integer *itol, doublereal *rtol,
                   doublereal *atol, integer *itask, integer *istate,
                   integer *iopt, doublereal *rwork, integer *lrw,
                   integer *iwork, integer *liw, U_fp jac, U_fp psol,
                   integer *mf, doublereal *rpar, integer *ipar);

/* RKSUITE: RK suite functions */
extern int setup_(integer *neq, doublereal *tstart, doublereal *ystart,
                  doublereal *tend, doublereal *tol, doublereal *thres,
                  integer *method, char *task, logical *errass,
                  doublereal *hstart, doublereal *work, integer *lenwrk,
                  logical *mesage, ftnlen task_len);
extern int ut_(U_fp f, doublereal *twant, doublereal *tgot, doublereal *ygot,
               doublereal *ypgot, doublereal *ymax, doublereal *work,
               integer *uflag);
extern int ct_(U_fp f, doublereal *tnow, doublereal *ynow, doublereal *ypnow,
               doublereal *work, integer *cflag);
extern int intrp_(doublereal *twant, char *reqest, integer *nwant,
                  doublereal *ywant, doublereal *ypwant, U_fp f,
                  doublereal *work, doublereal *wrkint, integer *lenint,
                  ftnlen reqest_len);
extern int stat_(integer *totfcn, integer *stpcst, doublereal *waste,
                 integer *stpsok, doublereal *hnext);

/* RKC: Runge-Kutta-Chebyshev order 2 */
extern int rkc_(integer *neqn, U_fp f, doublereal *y, doublereal *t,
                doublereal *tend, doublereal *rtol, doublereal *atol,
                integer *info, doublereal *work, integer *idid);
extern int rkcint_(doublereal *t, doublereal *yint, integer *neqn,
                   doublereal *work);

/* ============================================
 * JavaScript callback function pointers
 * ============================================ */

/* Callback type for ODE function: f(n, t, y, f_out) */
typedef void (*js_fcn_callback_t)(integer n, doublereal t, doublereal *y, doublereal *f);

/* Callback type for solution output: solout(nr, told, t, y, n, irtrn) */
typedef void (*js_solout_callback_t)(integer nr, doublereal told, doublereal t,
                                     doublereal *y, integer n, integer *irtrn);

/* Callback type for Jacobian: jac(n, t, y, dfy, ldfy) */
typedef void (*js_jac_callback_t)(integer n, doublereal t, doublereal *y,
                                  doublereal *dfy, integer ldfy);

/* Callback type for VODE/VODPK Jacobian: jac(n, t, y, ml, mu, pd, nrowpd) */
typedef void (*js_jac_vode_callback_t)(integer n, doublereal t, doublereal *y,
                                        integer ml, integer mu, doublereal *pd,
                                        integer nrowpd);

/* Callback type for VODPK preconditioner solve */
typedef int (*js_psol_callback_t)(integer n, doublereal t, doublereal *y,
                                   doublereal *savf, doublereal *wk,
                                   doublereal hrl1, doublereal *wp,
                                   integer *iwp, doublereal *b, integer lr);

/* Callback type for spectral radius estimation (RKC) */
typedef doublereal (*js_spcrad_callback_t)(integer n, doublereal t, doublereal *y);

/* Global callback pointers */
static js_fcn_callback_t g_fcn_callback = NULL;
static js_solout_callback_t g_solout_callback = NULL;
static js_jac_callback_t g_jac_callback = NULL;
static js_jac_vode_callback_t g_jac_vode_callback = NULL;
static js_psol_callback_t g_psol_callback = NULL;
static js_spcrad_callback_t g_spcrad_callback = NULL;

/* ============================================
 * Callback setters (called from JavaScript)
 * ============================================ */

EXPORT void wasm_set_fcn_callback(js_fcn_callback_t cb) {
    g_fcn_callback = cb;
}

EXPORT void wasm_set_solout_callback(js_solout_callback_t cb) {
    g_solout_callback = cb;
}

EXPORT void wasm_set_jac_callback(js_jac_callback_t cb) {
    g_jac_callback = cb;
}

EXPORT void wasm_set_jac_vode_callback(js_jac_vode_callback_t cb) {
    g_jac_vode_callback = cb;
}

EXPORT void wasm_set_psol_callback(js_psol_callback_t cb) {
    g_psol_callback = cb;
}

EXPORT void wasm_set_spcrad_callback(js_spcrad_callback_t cb) {
    g_spcrad_callback = cb;
}

/* ============================================
 * Fortran-compatible thunk functions
 * These adapt our JS callbacks to the Fortran calling convention
 * ============================================ */

/*
 * FCN thunk: Fortran calls FCN(N, X, Y, F, RPAR, IPAR)
 * We translate to our JS callback: fcn(n, x, y, f)
 */
static int fcn_thunk_(integer *n, doublereal *x, doublereal *y,
                      doublereal *f, doublereal *rpar, integer *ipar) {
    (void)rpar; (void)ipar;
    if (g_fcn_callback) {
        g_fcn_callback(*n, *x, y, f);
    }
    return 0;
}

/*
 * SOLOUT thunk for DOPRI5/DOP853: SOLOUT(NR, XOLD, X, Y, N, CON, ICOMP, ND, RPAR, IPAR, IRTRN)
 * We translate to: solout(nr, xold, x, y, n, irtrn)
 */
static int solout_dopri_thunk_(integer *nr, doublereal *xold, doublereal *x,
                               doublereal *y, integer *n, doublereal *con,
                               integer *icomp, integer *nd, doublereal *rpar,
                               integer *ipar, integer *irtrn) {
    (void)con; (void)icomp; (void)nd; (void)rpar; (void)ipar;
    if (g_solout_callback) {
        g_solout_callback(*nr, *xold, *x, y, *n, irtrn);
    }
    return 0;
}

/*
 * SOLOUT thunk for RADAU5: SOLOUT(NR, XOLD, X, Y, CONT, LRC, N, RPAR, IPAR, IRTRN)
 * We translate to: solout(nr, xold, x, y, n, irtrn)
 */
static int solout_radau_thunk_(integer *nr, doublereal *xold, doublereal *x,
                               doublereal *y, doublereal *cont, integer *lrc,
                               integer *n, doublereal *rpar, integer *ipar,
                               integer *irtrn) {
    (void)cont; (void)lrc; (void)rpar; (void)ipar;
    if (g_solout_callback) {
        g_solout_callback(*nr, *xold, *x, y, *n, irtrn);
    }
    return 0;
}

/*
 * JAC thunk for RADAU5: JAC(N, X, Y, DFY, LDFY, RPAR, IPAR)
 * We translate to: jac(n, x, y, dfy, ldfy)
 */
static int jac_thunk_(integer *n, doublereal *x, doublereal *y,
                      doublereal *dfy, integer *ldfy, doublereal *rpar,
                      integer *ipar) {
    (void)rpar; (void)ipar;
    if (g_jac_callback) {
        g_jac_callback(*n, *x, y, dfy, *ldfy);
    }
    return 0;
}

/*
 * MAS thunk for RADAU5 mass matrix (we don't support mass matrices yet)
 * Returns identity mass matrix behavior
 */
static int mas_thunk_(integer *n, doublereal *am, integer *lmas,
                      doublereal *rpar, integer *ipar) {
    (void)n; (void)am; (void)lmas; (void)rpar; (void)ipar;
    /* Identity mass matrix - nothing to do */
    return 0;
}

/*
 * FCN thunk for RKF45/ODE: F(T, Y, YP) - different signature from Hairer solvers
 * Fortran calls F(T, Y, YP) with T as scalar
 */
static int fcn_rkf45_thunk_(doublereal *t, doublereal *y, doublereal *yp) {
    if (g_fcn_callback) {
        /* We need to know n, but it's not passed - use a global or extract from work */
        /* For now, assume n is stored somewhere accessible */
        g_fcn_callback(0, *t, y, yp);  /* n=0 as placeholder, TS wrapper knows n */
    }
    return 0;
}

/*
 * FCN thunk for DVERK: FCN(N, X, Y, YP) - needs S_fp signature
 */
static int fcn_dverk_thunk_(integer *n, doublereal *x, doublereal *y, doublereal *yp) {
    if (g_fcn_callback) {
        g_fcn_callback(*n, *x, y, yp);
    }
    return 0;
}

/*
 * FCN thunk for VODE: F(NEQ, T, Y, YDOT, RPAR, IPAR)
 */
static int fcn_vode_thunk_(integer *neq, doublereal *t, doublereal *y,
                           doublereal *ydot, doublereal *rpar, integer *ipar) {
    (void)rpar; (void)ipar;
    if (g_fcn_callback) {
        g_fcn_callback(*neq, *t, y, ydot);
    }
    return 0;
}

/*
 * JAC thunk for VODE: JAC(NEQ, T, Y, ML, MU, PD, NROWPD, RPAR, IPAR)
 */
static int jac_vode_thunk_(integer *neq, doublereal *t, doublereal *y,
                           integer *ml, integer *mu, doublereal *pd,
                           integer *nrowpd, doublereal *rpar, integer *ipar) {
    (void)rpar; (void)ipar;
    if (g_jac_vode_callback) {
        g_jac_vode_callback(*neq, *t, y, *ml, *mu, pd, *nrowpd);
    }
    return 0;
}

/*
 * PSOL thunk for VODPK preconditioner solve
 */
static int psol_vodpk_thunk_(integer *neq, doublereal *t, doublereal *y,
                              doublereal *savf, doublereal *wk,
                              doublereal *hrl1, doublereal *wp, integer *iwp,
                              doublereal *b, integer *lr, integer *ier) {
    if (g_psol_callback) {
        *ier = g_psol_callback(*neq, *t, y, savf, wk, *hrl1, wp, iwp, b, *lr);
    } else {
        *ier = 0;
    }
    return 0;
}

/*
 * Dummy JAC for VODPK (preconditioner setup - often not needed)
 */
static int jac_vodpk_thunk_(S_fp f, integer *neq, doublereal *t, doublereal *y,
                            doublereal *ysv, doublereal *rewt, doublereal *fty,
                            doublereal *v, doublereal *hrl1, doublereal *wp,
                            integer *iwp, integer *ier, doublereal *rpar,
                            integer *ipar) {
    (void)f; (void)neq; (void)t; (void)y; (void)ysv; (void)rewt;
    (void)fty; (void)v; (void)hrl1; (void)wp; (void)iwp;
    (void)rpar; (void)ipar;
    *ier = 0;  /* Success - use identity preconditioner by default */
    return 0;
}

/*
 * SPCRAD thunk for RKC: Spectral radius estimation
 * RKC calls SPCRAD(N, T, Y) to get an estimate of the spectral radius.
 * If user provides callback, use it. Otherwise return 0 (RKC will estimate internally).
 */
doublereal spcrad_(integer *n, doublereal *t, doublereal *y) {
    if (g_spcrad_callback) {
        return g_spcrad_callback(*n, *t, y);
    }
    /* Return 0 to let RKC estimate spectral radius internally */
    return 0.0;
}

/* ============================================
 * WASM Entry Points
 * ============================================ */

/*
 * wasm_dopri5 - Solve ODE using DOPRI5 (Dormand-Prince 5(4))
 *
 * Parameters:
 *   n      - System dimension (number of equations)
 *   x      - Initial time (modified on output to final time)
 *   y      - State vector (n elements, modified on output)
 *   xend   - Target end time
 *   rtol   - Relative tolerance (scalar or array based on itol)
 *   atol   - Absolute tolerance (scalar or array based on itol)
 *   itol   - 0: scalar tolerances, 1: vector tolerances
 *   iout   - Output mode: 0=none, 1=after steps, 2=dense output
 *   work   - Work array (at least 8*n + 5*nrdens + 21 doubles)
 *   lwork  - Length of work array
 *   iwork  - Integer work array (at least nrdens + 21)
 *   liwork - Length of iwork array
 *   idid   - Output: status code (1=success, negative=error)
 */
EXPORT void wasm_dopri5(
    integer n, doublereal *x, doublereal *y, doublereal xend,
    doublereal *rtol, doublereal *atol, integer itol, integer iout,
    doublereal *work, integer lwork, integer *iwork, integer liwork,
    integer *idid)
{
    doublereal rpar = 0;
    integer ipar = 0;

    dopri5_(&n, (U_fp)fcn_thunk_, x, y, &xend, rtol, atol, &itol,
            (U_fp)solout_dopri_thunk_, &iout, work, &lwork, iwork, &liwork,
            &rpar, &ipar, idid);
}

/*
 * wasm_dop853 - Solve ODE using DOP853 (Dormand-Prince 8(5,3))
 *
 * Same interface as wasm_dopri5 but higher order method.
 * Work array needs at least 11*n + 8*nrdens + 21 doubles.
 */
EXPORT void wasm_dop853(
    integer n, doublereal *x, doublereal *y, doublereal xend,
    doublereal *rtol, doublereal *atol, integer itol, integer iout,
    doublereal *work, integer lwork, integer *iwork, integer liwork,
    integer *idid)
{
    doublereal rpar = 0;
    integer ipar = 0;

    dop853_(&n, (U_fp)fcn_thunk_, x, y, &xend, rtol, atol, &itol,
            (U_fp)solout_dopri_thunk_, &iout, work, &lwork, iwork, &liwork,
            &rpar, &ipar, idid);
}

/*
 * wasm_radau5 - Solve ODE using RADAU5 (implicit Runge-Kutta for stiff systems)
 *
 * Additional parameters vs DOPRI5:
 *   h      - Initial step size (can be 0 for automatic)
 *   ijac   - 0: internal numerical Jacobian, 1: user-supplied Jacobian
 *   mljac  - Lower bandwidth of Jacobian (n for full)
 *   mujac  - Upper bandwidth of Jacobian (n for full)
 *   imas   - 0: identity mass matrix, 1: user-supplied
 *   mlmas  - Lower bandwidth of mass matrix
 *   mumas  - Upper bandwidth of mass matrix
 *
 * Work array needs at least:
 *   Full Jacobian: n*(ljac+lmas+3*le+12)+20
 *   Banded: n*(ljac+lmas+3*le+12)+20 where ljac, lmas, le depend on bandwidths
 */
EXPORT void wasm_radau5(
    integer n, doublereal *x, doublereal *y, doublereal xend,
    doublereal h, doublereal *rtol, doublereal *atol, integer itol,
    integer ijac, integer mljac, integer mujac,
    integer imas, integer mlmas, integer mumas, integer iout,
    doublereal *work, integer lwork, integer *iwork, integer liwork,
    integer *idid)
{
    doublereal rpar = 0;
    integer ipar = 0;

    radau5_(&n, (U_fp)fcn_thunk_, x, y, &xend, &h, rtol, atol, &itol,
            (U_fp)jac_thunk_, &ijac, &mljac, &mujac,
            (U_fp)mas_thunk_, &imas, &mlmas, &mumas,
            (U_fp)solout_radau_thunk_, &iout, work, &lwork, iwork, &liwork,
            &rpar, &ipar, idid);
}

/* ============================================
 * Netlib Solver Entry Points
 * ============================================ */

/*
 * wasm_rkf45 - Solve ODE using RKF45 (Runge-Kutta-Fehlberg 4(5))
 *
 * Parameters:
 *   neqn   - Number of equations
 *   y      - State vector (modified on output)
 *   t      - Current time (modified on output)
 *   tout   - Target end time
 *   relerr - Relative error tolerance
 *   abserr - Absolute error tolerance
 *   iflag  - Status flag (input: 1=first call, output: 2=success)
 *   work   - Work array (at least 3+6*neqn doubles)
 *   iwork  - Integer work array (at least 5 integers)
 */
EXPORT void wasm_rkf45(
    integer neqn, doublereal *y, doublereal *t, doublereal tout,
    doublereal *relerr, doublereal *abserr, integer *iflag,
    doublereal *work, integer *iwork)
{
    rkf45_((U_fp)fcn_thunk_, &neqn, y, t, &tout, relerr, abserr,
           iflag, work, iwork);
}

/*
 * wasm_dverk - Solve ODE using DVERK (Verner 6(5) Runge-Kutta)
 *
 * Parameters:
 *   n      - Number of equations
 *   x      - Current time (modified on output)
 *   y      - State vector (modified on output)
 *   xend   - Target end time
 *   tol    - Tolerance
 *   ind    - Indicator (input: 1=first call, output: 3=success)
 *   c      - Options array (24 or n+30 doubles depending on c[0])
 *   nw     - Length of w array
 *   w      - Work array
 */
EXPORT void wasm_dverk(
    integer n, doublereal *x, doublereal *y, doublereal xend,
    doublereal tol, integer *ind, doublereal *c__, integer nw, doublereal *w)
{
    dverk_(&n, (S_fp)fcn_dverk_thunk_, x, y, &xend, &tol, ind, c__, &nw, w);
}

/*
 * wasm_ode - Solve ODE using Adams-Bashforth-Moulton
 *
 * Same interface as RKF45.
 */
EXPORT void wasm_ode(
    integer neqn, doublereal *y, doublereal *t, doublereal tout,
    doublereal *relerr, doublereal *abserr, integer *iflag,
    doublereal *work, integer *iwork)
{
    ode_((U_fp)fcn_thunk_, &neqn, y, t, &tout, relerr, abserr,
         iflag, work, iwork);
}

/*
 * wasm_vode - Solve ODE using VODE (Variable-coefficient ODE solver)
 *
 * Parameters:
 *   neq    - Number of equations
 *   y      - State vector (modified on output)
 *   t      - Current time (modified on output)
 *   tout   - Target end time
 *   itol   - Tolerance type (1=scalar, 2=array)
 *   rtol   - Relative tolerance
 *   atol   - Absolute tolerance
 *   itask  - Task indicator (1=normal, 2=one step)
 *   istate - State indicator (input: 1=first call, output: 2=success)
 *   iopt   - Options flag (0=no options, 1=use options in rwork/iwork)
 *   rwork  - Real work array
 *   lrw    - Length of rwork
 *   iwork  - Integer work array
 *   liw    - Length of iwork
 *   mf     - Method flag (10=Adams, 21/22=BDF with full/banded Jacobian)
 */
EXPORT void wasm_vode(
    integer neq, doublereal *y, doublereal *t, doublereal tout,
    integer itol, doublereal *rtol, doublereal *atol,
    integer itask, integer *istate, integer iopt,
    doublereal *rwork, integer lrw, integer *iwork, integer liw,
    integer mf)
{
    doublereal rpar = 0;
    integer ipar = 0;

    dvode_((S_fp)fcn_vode_thunk_, &neq, y, t, &tout, &itol, rtol, atol,
           &itask, istate, &iopt, rwork, &lrw, iwork, &liw,
           (U_fp)jac_vode_thunk_, &mf, &rpar, &ipar);
}

/*
 * wasm_zvode - Solve complex-valued ODE using ZVODE
 *
 * Same as VODE but with complex arrays.
 * y and zwork are interleaved real/imag pairs.
 */
EXPORT void wasm_zvode(
    integer neq, doublereal *y, doublereal *t, doublereal tout,
    integer itol, doublereal *rtol, doublereal *atol,
    integer itask, integer *istate, integer iopt,
    doublereal *zwork, integer lzw,
    doublereal *rwork, integer lrw, integer *iwork, integer liw,
    integer mf)
{
    doublereal rpar = 0;
    integer ipar = 0;

    /* y and zwork are treated as doublecomplex (pairs of doubles) */
    zvode_((S_fp)fcn_vode_thunk_, &neq, (doublecomplex*)y, t, &tout,
           &itol, rtol, atol, &itask, istate, &iopt,
           (doublecomplex*)zwork, &lzw, rwork, &lrw, iwork, &liw,
           (U_fp)jac_vode_thunk_, &mf, &rpar, &ipar);
}

/*
 * wasm_vodpk - Solve ODE using VODPK (VODE with Krylov methods)
 *
 * Same interface as VODE but uses preconditioned Krylov methods.
 */
EXPORT void wasm_vodpk(
    integer neq, doublereal *y, doublereal *t, doublereal tout,
    integer itol, doublereal *rtol, doublereal *atol,
    integer itask, integer *istate, integer iopt,
    doublereal *rwork, integer lrw, integer *iwork, integer liw,
    integer mf)
{
    doublereal rpar = 0;
    integer ipar = 0;

    dvodpk_((S_fp)fcn_vode_thunk_, &neq, y, t, &tout, &itol, rtol, atol,
            &itask, istate, &iopt, rwork, &lrw, iwork, &liw,
            (U_fp)jac_vodpk_thunk_, (U_fp)psol_vodpk_thunk_, &mf,
            &rpar, &ipar);
}

/*
 * wasm_rksuite_setup - Initialize RKSUITE solver
 *
 * Parameters:
 *   neq     - Number of equations
 *   tstart  - Initial time
 *   ystart  - Initial state vector
 *   tend    - Target end time
 *   tol     - Tolerance
 *   thres   - Threshold array (or scalar pointer)
 *   method  - Method (1=(2,3), 2=(4,5), 3=(7,8))
 *   task    - 'U' for UT mode, 'C' for CT mode
 *   errass  - Enable error assessment (0=no, 1=yes)
 *   hstart  - Initial step size (0 for automatic)
 *   work    - Work array
 *   lenwrk  - Length of work array
 *   mesage  - Enable messages (0=no, 1=yes)
 */
EXPORT void wasm_rksuite_setup(
    integer neq, doublereal tstart, doublereal *ystart, doublereal tend,
    doublereal tol, doublereal *thres, integer method,
    integer task_code, integer errass, doublereal hstart,
    doublereal *work, integer lenwrk, integer mesage)
{
    char task[2];
    logical errass_log = errass ? TRUE_ : FALSE_;
    logical mesage_log = mesage ? TRUE_ : FALSE_;

    task[0] = (task_code == 0) ? 'U' : 'C';
    task[1] = '\0';

    setup_(&neq, &tstart, ystart, &tend, &tol, thres, &method,
           task, &errass_log, &hstart, work, &lenwrk, &mesage_log, (ftnlen)1);
}

/*
 * wasm_rksuite_ut - Advance RKSUITE integration (UT mode)
 *
 * Parameters:
 *   twant  - Desired time (may not be reached exactly)
 *   tgot   - Time actually reached (output)
 *   ygot   - Solution at tgot (output)
 *   ypgot  - Derivative at tgot (output)
 *   ymax   - Maximum y values seen (output)
 *   work   - Work array from setup
 *   uflag  - Status flag (output: 1=success, others=error)
 */
EXPORT void wasm_rksuite_ut(
    doublereal twant, doublereal *tgot, doublereal *ygot,
    doublereal *ypgot, doublereal *ymax, doublereal *work, integer *uflag)
{
    ut_((U_fp)fcn_thunk_, &twant, tgot, ygot, ypgot, ymax, work, uflag);
}

/*
 * wasm_rksuite_ct - Advance RKSUITE one step (CT mode)
 */
EXPORT void wasm_rksuite_ct(
    doublereal *tnow, doublereal *ynow, doublereal *ypnow,
    doublereal *work, integer *cflag)
{
    ct_((U_fp)fcn_thunk_, tnow, ynow, ypnow, work, cflag);
}

/*
 * wasm_rksuite_stat - Get RKSUITE statistics
 */
EXPORT void wasm_rksuite_stat(
    integer *totfcn, integer *stpcst, doublereal *waste,
    integer *stpsok, doublereal *hnext)
{
    stat_(totfcn, stpcst, waste, stpsok, hnext);
}

/*
 * wasm_rkc - Solve ODE using RKC (Runge-Kutta-Chebyshev order 2)
 *
 * Parameters:
 *   neqn  - Number of equations
 *   y     - State vector (modified on output)
 *   t     - Current time (modified on output)
 *   tend  - Target end time
 *   rtol  - Relative tolerance
 *   atol  - Absolute tolerance (scalar or array based on info[1])
 *   info  - Options array (4 integers)
 *   work  - Work array
 *   idid  - Status flag (output: 1=success)
 */
EXPORT void wasm_rkc(
    integer neqn, doublereal *y, doublereal *t, doublereal tend,
    doublereal rtol, doublereal *atol, integer *info,
    doublereal *work, integer *idid)
{
    rkc_(&neqn, (U_fp)fcn_thunk_, y, t, &tend, &rtol, atol, info, work, idid);
}

/*
 * wasm_rkc_int - Dense output interpolation for RKC
 */
EXPORT void wasm_rkc_int(
    doublereal t, doublereal *yint, integer neqn, doublereal *work)
{
    rkcint_(&t, yint, &neqn, work);
}

/* ============================================
 * Dense Output Functions
 * ============================================ */

/*
 * wasm_contd5 - Dense output for DOPRI5
 *
 * Parameters:
 *   ii   - Component index (1-based Fortran indexing)
 *   x    - Evaluation point (must be in [xold, x] from last step)
 *   con  - Continuous output coefficients (from work array)
 *   icomp - Component indices for dense output
 *   nd   - Number of components for dense output
 *
 * Returns: Interpolated value of component ii at time x
 */
EXPORT doublereal wasm_contd5(integer ii, doublereal x, doublereal *con,
                               integer *icomp, integer nd) {
    return contd5_(&ii, &x, con, icomp, &nd);
}

/*
 * wasm_contd8 - Dense output for DOP853
 */
EXPORT doublereal wasm_contd8(integer ii, doublereal x, doublereal *con,
                               integer *icomp, integer nd) {
    return contd8_(&ii, &x, con, icomp, &nd);
}

/*
 * wasm_contr5 - Dense output for RADAU5
 */
EXPORT doublereal wasm_contr5(integer ii, doublereal s, doublereal *cont,
                               integer lrc) {
    return contr5_(&ii, &s, cont, &lrc);
}

/* ============================================
 * Utility Functions
 * ============================================ */

/*
 * Get minimum work array size for DOPRI5
 * n: system dimension
 * nrdens: number of components for dense output (0 if not used)
 */
EXPORT integer wasm_dopri5_work_size(integer n, integer nrdens) {
    return 8 * n + 5 * nrdens + 21;
}

EXPORT integer wasm_dopri5_iwork_size(integer nrdens) {
    return nrdens + 21;
}

/*
 * Get minimum work array size for DOP853
 */
EXPORT integer wasm_dop853_work_size(integer n, integer nrdens) {
    return 11 * n + 8 * nrdens + 21;
}

EXPORT integer wasm_dop853_iwork_size(integer nrdens) {
    return nrdens + 21;
}

/*
 * Get minimum work array size for RADAU5 with full Jacobian
 */
EXPORT integer wasm_radau5_work_size(integer n) {
    /* For full Jacobian: ljac=n, le=n, lmas=0 (identity mass) */
    /* Formula: n*(ljac+lmas+3*le+12)+20 = n*(n+0+3*n+12)+20 = n*(4*n+12)+20 */
    return n * (4 * n + 12) + 20;
}

EXPORT integer wasm_radau5_iwork_size(integer n) {
    /* For full Jacobian: 3*n + 20 */
    return 3 * n + 20;
}

/*
 * Get minimum work array size for RKF45
 */
EXPORT integer wasm_rkf45_work_size(integer n) {
    return 3 + 6 * n;
}

EXPORT integer wasm_rkf45_iwork_size(void) {
    return 5;
}

/*
 * Get minimum work array size for DVERK
 * error_control: c[0] value (0-5)
 */
EXPORT integer wasm_dverk_c_size(integer n, integer error_control) {
    /* c array is 24 normally, or n+30 for error_control=4 or 5 */
    if (error_control == 4 || error_control == 5) {
        return n + 30;
    }
    return 24;
}

EXPORT integer wasm_dverk_work_size(integer n) {
    return 9 * n;
}

/*
 * Get minimum work array size for ODE (Adams-Bashforth-Moulton)
 */
EXPORT integer wasm_ode_work_size(integer n) {
    return 100 + 21 * n;
}

EXPORT integer wasm_ode_iwork_size(void) {
    return 5;
}

/*
 * Get minimum work array size for VODE
 * mf: method flag
 *   10 = Adams (non-stiff), no Jacobian
 *   21 = BDF, full Jacobian
 *   22 = BDF, full Jacobian (internal)
 *   24 = BDF, banded Jacobian
 *   25 = BDF, banded Jacobian (internal)
 */
EXPORT integer wasm_vode_rwork_size(integer n, integer mf) {
    if (mf == 10) {
        return 20 + 16 * n;
    } else if (mf == 21 || mf == 22) {
        return 22 + 9 * n + 2 * n * n;
    } else {
        /* Banded: 22 + 10*n + (2*ml + mu)*n, use worst case */
        return 22 + 9 * n + 2 * n * n;
    }
}

EXPORT integer wasm_vode_iwork_size(integer n, integer mf) {
    if (mf == 10) {
        return 30;
    } else {
        return 30 + n;
    }
}

/*
 * Get minimum work array size for ZVODE (complex VODE)
 */
EXPORT integer wasm_zvode_zwork_size(integer n, integer mf) {
    if (mf == 10) {
        return 15 * n;
    } else if (mf == 21 || mf == 22) {
        return 8 * n + 2 * n * n;
    } else {
        return 8 * n + 2 * n * n;
    }
}

EXPORT integer wasm_zvode_rwork_size(integer n, integer mf) {
    return 20 + n;
}

/*
 * Get minimum work array size for VODPK
 */
EXPORT integer wasm_vodpk_rwork_size(integer n, integer maxl, integer maxp, integer lwp) {
    /* Base: 20 + 16*n + n for SAVF */
    /* Krylov: maxl*n + (maxp+3)*n + (maxl+3)*maxl + 1 */
    /* Preconditioner: lwp */
    return 20 + 17 * n + maxl * n + (maxp + 3) * n + (maxl + 3) * maxl + 1 + lwp;
}

EXPORT integer wasm_vodpk_iwork_size(integer n, integer liwp) {
    return 30 + liwp;
}

/*
 * Get minimum work array size for RKSUITE
 * method: 1=(2,3), 2=(4,5), 3=(7,8)
 */
EXPORT integer wasm_rksuite_work_size(integer n, integer method, integer errass) {
    /* Conservative estimate: 32*n covers all cases */
    (void)method; (void)errass;
    return 32 * n;
}

/*
 * Get minimum work array size for RKC
 * use_spcrad: 0 if user provides spectral radius function, 1 otherwise
 */
EXPORT integer wasm_rkc_work_size(integer n, integer use_spcrad) {
    if (use_spcrad) {
        return 8 + 5 * n;  /* RKC estimates spectral radius internally */
    } else {
        return 8 + 4 * n;  /* User provides spectral radius */
    }
}

/* ============================================
 * Memory helpers for JavaScript
 * ============================================ */

EXPORT void* wasm_malloc(size_t size) {
    return malloc(size);
}

EXPORT void wasm_free(void *ptr) {
    free(ptr);
}
