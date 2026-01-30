/*
 * ODE WASM Wrapper
 *
 * Provides clean C entry points for the Hairer ODE solvers (DOPRI5, DOP853, RADAU5)
 * compiled from Fortran via f2c.
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

/* Global callback pointers */
static js_fcn_callback_t g_fcn_callback = NULL;
static js_solout_callback_t g_solout_callback = NULL;
static js_jac_callback_t g_jac_callback = NULL;

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

/* ============================================
 * Memory helpers for JavaScript
 * ============================================ */

EXPORT void* wasm_malloc(size_t size) {
    return malloc(size);
}

EXPORT void wasm_free(void *ptr) {
    free(ptr);
}
