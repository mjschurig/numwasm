/*
 * vode_common.h - Shared declarations for VODE family of ODE solvers
 *
 * This header declares the common blocks and shared utility functions
 * used by vode.c, vodpk.c, and zvode.c to avoid duplicate symbol errors.
 */

#ifndef VODE_COMMON_H
#define VODE_COMMON_H

#include "f2c.h"

/* ========================================================================
 * Common Block Declarations
 * ======================================================================== */

/* DVOD01 - Primary VODE state common block (real: 48 doubles, int: 33 ints) */
extern union {
    struct {
        doublereal acnrm, ccmxj, conp, crate, drc, el[13], eta, etamax, h__,
                hmin, hmxi, hnew, hscal, prl1, rc, rl1, tau[13], tq[5], tn,
                uround;
        integer icf, init, ipup, jcur, jstart, jsv, kflag, kuth, l, lmax, lyh,
                lewt, lacor, lsavf, lwm, liwm, locjs, maxord, meth, miter,
                msbj, mxhnil, mxstep, n, newh, newq, nhnil, nq, nqnyh, nqwait,
                nslj, nslp, nyh;
    } _1;
    struct {
        doublereal rvod1[48];
        integer ivod1[33];
    } _2;
} dvod01_;

#define dvod01_1 (dvod01_._1)
#define dvod01_2 (dvod01_._2)

/* DVOD02 - VODE statistics common block (real: 1 double, int: 8 ints)
   Note: The 4th integer is 'nje' for VODE/ZVODE but 'npe' for VODPK.
   Both refer to the same memory location, so we provide both names. */
extern union {
    struct {
        doublereal hu;
        integer ncfn, netf, nfe, nje, nlu, nni, nqu, nst;
    } _1;
    struct {
        doublereal rvod2[1];
        integer ivod2[8];
    } _2;
    struct {
        doublereal hu;
        integer ncfn, netf, nfe, npe, nlu, nni, nqu, nst;
    } _3;
} dvod02_;

#define dvod02_1 (dvod02_._1)
#define dvod02_2 (dvod02_._2)
#define dvod02_3 (dvod02_._3)

/* ========================================================================
 * Shared Utility Function Prototypes
 * ======================================================================== */

/* Error message handler */
extern int xerrwd_(char *msg, integer *nmes, integer *nerr, integer *level,
                   integer *ni, integer *i1, integer *i2, integer *nr,
                   doublereal *r1, doublereal *r2, ftnlen msg_len);

/* Machine epsilon computation */
extern doublereal dumach_(void);

/* Weighted root-mean-square vector norm */
extern doublereal dvnorm_(integer *n, doublereal *v, doublereal *w);

/* Error weight vector setup */
extern int dewset_(integer *n, integer *itol, doublereal *rtol,
                   doublereal *atol, doublereal *ycur, doublereal *ewt);

/* Initial step size computation */
extern int dvhin_(integer *n, doublereal *t0, doublereal *y0,
                  doublereal *ydot, S_fp f, doublereal *rpar, integer *ipar,
                  doublereal *tout, doublereal *uround, doublereal *ewt,
                  integer *itol, doublereal *atol, doublereal *y,
                  doublereal *temp, doublereal *h0, integer *niter, integer *ier);

/* Interpolation routine for solution output */
extern int dvindy_(doublereal *t, integer *k, doublereal *yh,
                   integer *ldyh, doublereal *dky, integer *iflag);

/* Core integration step */
extern int dvstep_(doublereal *y, doublereal *yh, integer *ldyh,
                   doublereal *yh1, doublereal *ewt, doublereal *savf,
                   doublereal *vsav, doublereal *acor, doublereal *wm,
                   integer *iwm, S_fp f, U_fp jac, S_fp psol, S_fp vnls,
                   doublereal *rpar, integer *ipar);

/* Method parameter setup */
extern int dvset_(void);

/* Array adjustment routine */
extern int dvjust_(doublereal *yh, integer *ldyh, integer *iord);

#endif /* VODE_COMMON_H */
