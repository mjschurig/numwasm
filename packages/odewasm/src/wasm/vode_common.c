/*
 * vode_common.c - Shared routines for VODE family of ODE solvers
 *
 * This file contains common utility functions and data structures shared
 * between vode.c, vodpk.c, and zvode.c. Extracting these avoids duplicate
 * symbol errors when linking all three solvers together.
 *
 * Functions included:
 *   - xerrwd_  : Error message handler
 *   - dumach_  : Machine epsilon computation
 *   - dumsum_  : Helper for dumach_
 *   - dvnorm_  : Weighted RMS vector norm
 *   - dewset_  : Error weight vector setup
 *   - dvhin_   : Initial step size computation
 *   - dvindy_  : Interpolation for solution output
 *   - dvstep_  : Core integration step
 *   - dvset_   : Method parameter setup
 *   - dvjust_  : Array adjustment routine
 *   - ixsav_   : Save/restore integer values
 *   - xsetf_   : Set error print flag
 *   - xsetun_  : Set error unit number
 *   - iumach_  : Get machine unit number
 *
 * Common blocks defined here:
 *   - dvod01_  : Primary VODE state
 *   - dvod02_  : VODE statistics
 */

#include "f2c.h"

/* ========================================================================
 * Common Block Definitions
 * ======================================================================== */

/* DVOD01 - Primary VODE state common block */
union {
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

/* DVOD02 - VODE statistics common block */
union {
    struct {
        doublereal hu;
        integer ncfn, netf, nfe, nje, nlu, nni, nqu, nst;
    } _1;
    struct {
        doublereal rvod2[1];
        integer ivod2[8];
    } _2;
} dvod02_;

/* ========================================================================
 * Static Constants (used by multiple functions)
 * ======================================================================== */

static integer c__0 = 0;
static integer c__1 = 1;
static integer c__2 = 2;
static integer c__30 = 30;
static integer c__51 = 51;
static integer c__52 = 52;
static integer c__60 = 60;
static integer c_n1 = -1;
static logical c_false = FALSE_;
static logical c_true = TRUE_;
static doublereal c_b_one = 1.;

/* ========================================================================
 * IXSAV - Save and restore integer values
 * ======================================================================== */

integer ixsav_(integer *ipar, integer *ivalue, logical *iset)
{
    static integer lunit = 6;
    static integer mesflg = 1;
    integer ret_val;

    if (*ipar == 1) {
        if (*iset) {
            lunit = *ivalue;
        }
        ret_val = lunit;
    } else if (*ipar == 2) {
        if (*iset) {
            mesflg = *ivalue;
        }
        ret_val = mesflg;
    }
    return ret_val;
}

/* ========================================================================
 * XERRWD - Write error message with values
 * ======================================================================== */

int xerrwd_(char *msg, integer *nmes, integer *nerr, integer *level,
            integer *ni, integer *i1, integer *i2, integer *nr,
            doublereal *r1, doublereal *r2, ftnlen msg_len)
{
    static char fmt_10[] = "(1x,a)";
    static char fmt_20[] = "(6x,\"In above message,  I1 =\",i10)";
    static char fmt_30[] = "(6x,\"In above message,  I1 =\",i10,3x,\"I2 =\",i10)";
    static char fmt_40[] = "(6x,\"In above message,  R1 =\",d21.13)";
    static char fmt_50[] = "(6x,\"In above,  R1 =\",d21.13,3x,\"R2 =\",d21.13)";

    integer lunit, mesflg;
    static cilist io___1 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___2 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___3 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___4 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___5 = { 0, 0, 0, fmt_50, 0 };

    lunit = ixsav_(&c__1, &c__0, &c_false);
    mesflg = ixsav_(&c__2, &c__0, &c_false);
    if (mesflg == 0) {
        goto L100;
    }

    io___1.ciunit = lunit;
    s_wsfe(&io___1);
    do_fio(&c__1, msg, msg_len);
    e_wsfe();
    if (*ni == 1) {
        io___2.ciunit = lunit;
        s_wsfe(&io___2);
        do_fio(&c__1, (char *)&(*i1), (ftnlen)sizeof(integer));
        e_wsfe();
    }
    if (*ni == 2) {
        io___3.ciunit = lunit;
        s_wsfe(&io___3);
        do_fio(&c__1, (char *)&(*i1), (ftnlen)sizeof(integer));
        do_fio(&c__1, (char *)&(*i2), (ftnlen)sizeof(integer));
        e_wsfe();
    }
    if (*nr == 1) {
        io___4.ciunit = lunit;
        s_wsfe(&io___4);
        do_fio(&c__1, (char *)&(*r1), (ftnlen)sizeof(doublereal));
        e_wsfe();
    }
    if (*nr == 2) {
        io___5.ciunit = lunit;
        s_wsfe(&io___5);
        do_fio(&c__1, (char *)&(*r1), (ftnlen)sizeof(doublereal));
        do_fio(&c__1, (char *)&(*r2), (ftnlen)sizeof(doublereal));
        e_wsfe();
    }

L100:
    if (*level != 2) {
        return 0;
    }
    s_stop("", (ftnlen)0);
    return 0;
}

/* ========================================================================
 * XSETF - Reset the error print control flag
 * ======================================================================== */

int xsetf_(integer *mflag)
{
    integer junk;
    if (*mflag == 0 || *mflag == 1) {
        junk = ixsav_(&c__2, mflag, &c_true);
    }
    return 0;
}

/* ========================================================================
 * XSETUN - Reset the logical unit number for error messages
 * ======================================================================== */

int xsetun_(integer *lun)
{
    integer junk;
    if (*lun > 0) {
        junk = ixsav_(&c__1, lun, &c_true);
    }
    return 0;
}

/* ========================================================================
 * IUMACH - Get standard output unit number
 * ======================================================================== */

integer iumach_(void)
{
    return 6;
}

/* ========================================================================
 * DUMSUM - Helper routine for DUMACH
 * ======================================================================== */

int dumsum_(doublereal *a, doublereal *b, doublereal *c__)
{
    *c__ = *a + *b;
    return 0;
}

/* ========================================================================
 * DUMACH - Compute the unit roundoff of the machine
 * ======================================================================== */

doublereal dumach_(void)
{
    doublereal ret_val;
    doublereal u, comp;

    u = 1.;
L10:
    u *= .5;
    dumsum_(&c_b_one, &u, &comp);
    if (comp != 1.) {
        goto L10;
    }
    ret_val = u * 2.;
    return ret_val;
}

/* ========================================================================
 * DVNORM - Weighted root-mean-square vector norm
 * ======================================================================== */

doublereal dvnorm_(integer *n, doublereal *v, doublereal *w)
{
    integer i__1;
    doublereal ret_val, d__1;
    integer i__;
    doublereal sum;

    --w;
    --v;

    sum = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        d__1 = v[i__] * w[i__];
        sum += d__1 * d__1;
    }
    ret_val = sqrt(sum / *n);
    return ret_val;
}

/* ========================================================================
 * DEWSET - Set error weight vector
 * ======================================================================== */

int dewset_(integer *n, integer *itol, doublereal *rtol,
            doublereal *atol, doublereal *ycur, doublereal *ewt)
{
    integer i__1;
    doublereal d__1;
    integer i__;

    --ewt;
    --ycur;
    --rtol;
    --atol;

    switch (*itol) {
        case 1:  goto L10;
        case 2:  goto L20;
        case 3:  goto L30;
        case 4:  goto L40;
    }
L10:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        ewt[i__] = rtol[1] * (d__1 = ycur[i__], abs(d__1)) + atol[1];
    }
    return 0;
L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        ewt[i__] = rtol[1] * (d__1 = ycur[i__], abs(d__1)) + atol[i__];
    }
    return 0;
L30:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        ewt[i__] = rtol[i__] * (d__1 = ycur[i__], abs(d__1)) + atol[1];
    }
    return 0;
L40:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        ewt[i__] = rtol[i__] * (d__1 = ycur[i__], abs(d__1)) + atol[i__];
    }
    return 0;
}

/* ========================================================================
 * DVHIN - Compute initial step size
 * ======================================================================== */

int dvhin_(integer *n, doublereal *t0, doublereal *y0,
           doublereal *ydot, S_fp f, doublereal *rpar, integer *ipar,
           doublereal *tout, doublereal *uround, doublereal *ewt,
           integer *itol, doublereal *atol, doublereal *y,
           doublereal *temp, doublereal *h0, integer *niter, integer *ier)
{
    static doublereal half = .5;
    static doublereal hun = 100.;
    static doublereal pt1 = .1;
    static doublereal two = 2.;

    integer i__1;
    doublereal d__1, d__2;

    doublereal h__;
    integer i__;
    doublereal t1, hg, afi, hlb, hub, hrat, hnew;
    integer iter;
    doublereal delyi, atoli, tdist, yddnrm;
    doublereal tround;

    --temp;
    --y;
    --atol;
    --ewt;
    --ipar;
    --rpar;
    --ydot;
    --y0;

    *niter = 0;
    tdist = (d__1 = *tout - *t0, abs(d__1));
    d__1 = abs(*t0), d__2 = abs(*tout);
    tround = *uround * max(d__1,d__2);
    if (tdist < two * tround) {
        goto L100;
    }

    hlb = hun * tround;
    hub = pt1 * tdist;
    atoli = atol[1];
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        if (*itol == 2 || *itol == 4) {
            atoli = atol[i__];
        }
        delyi = pt1 * (d__1 = y0[i__], abs(d__1)) + atoli;
        afi = (d__1 = ydot[i__], abs(d__1));
        if (afi * hub > delyi) {
            hub = delyi / afi;
        }
    }

    iter = 0;
    hg = sqrt(hlb * hub);
    if (hub < hlb) {
        *h0 = hg;
        goto L90;
    }

L50:
    d__1 = *tout - *t0;
    h__ = d_sign(&hg, &d__1);
    t1 = *t0 + h__;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        y[i__] = y0[i__] + h__ * ydot[i__];
    }
    (*f)(n, &t1, &y[1], &temp[1], &rpar[1], &ipar[1]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        temp[i__] = (temp[i__] - ydot[i__]) / h__;
    }
    yddnrm = dvnorm_(n, &temp[1], &ewt[1]);
    if (yddnrm * hub * hub > two) {
        hnew = sqrt(two / yddnrm);
    } else {
        hnew = sqrt(hg * hub);
    }
    ++iter;
    if (iter >= 4) {
        goto L80;
    }
    hrat = hnew / hg;
    if (hrat > half && hrat < two) {
        goto L80;
    }
    if (iter >= 2 && hnew > two * hg) {
        hnew = hg;
        goto L80;
    }
    hg = hnew;
    goto L50;

L80:
    *h0 = hnew * half;
    if (*h0 < hlb) {
        *h0 = hlb;
    }
    if (*h0 > hub) {
        *h0 = hub;
    }
L90:
    d__1 = *tout - *t0;
    *h0 = d_sign(h0, &d__1);
    *niter = iter;
    *ier = 0;
    return 0;
L100:
    *ier = -1;
    return 0;
}

/* ========================================================================
 * DVINDY - Interpolation routine
 * ======================================================================== */

int dvindy_(doublereal *t, integer *k, doublereal *yh,
            integer *ldyh, doublereal *dky, integer *iflag)
{
    static doublereal hun = 100.;
    static doublereal zero = 0.;

    integer yh_dim1, yh_offset, i__1, i__2;
    doublereal d__1;

    doublereal c__;
    integer i__, j;
    doublereal r__, s;
    integer ic, jb, jj;
    doublereal tp;
    integer jb2, jj1, jp1;
    doublereal tn1;
    char msg[80];
    doublereal tfuzz;

    yh_dim1 = *ldyh;
    yh_offset = 1 + yh_dim1;
    yh -= yh_offset;
    --dky;

    *iflag = 0;
    if (*k < 0 || *k > dvod01_._1.nq) {
        goto L80;
    }
    d__1 = abs(dvod01_._1.tn) + abs(dvod02_._1.hu);
    tfuzz = hun * dvod01_._1.uround * d_sign(&d__1, &dvod02_._1.hu);
    tp = dvod01_._1.tn - dvod02_._1.hu - tfuzz;
    tn1 = dvod01_._1.tn + tfuzz;
    if ((*t - tp) * (*t - tn1) > zero) {
        goto L90;
    }

    s = (*t - dvod01_._1.tn) / dvod01_._1.h__;
    ic = 1;
    if (*k == 0) {
        goto L15;
    }
    jj1 = dvod01_._1.l - *k;
    i__1 = dvod01_._1.nq;
    for (jj = jj1; jj <= i__1; ++jj) {
        ic *= jj;
    }
L15:
    c__ = (real) ic;
    i__1 = dvod01_._1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
        dky[i__] = c__ * yh[i__ + dvod01_._1.l * yh_dim1];
    }
    if (*k == dvod01_._1.nq) {
        goto L55;
    }
    jb2 = dvod01_._1.nq - *k;
    i__1 = jb2;
    for (jb = 1; jb <= i__1; ++jb) {
        j = dvod01_._1.nq - jb;
        jp1 = j + 1;
        ic = 1;
        if (*k == 0) {
            goto L35;
        }
        jj1 = jp1 - *k;
        i__2 = j;
        for (jj = jj1; jj <= i__2; ++jj) {
            ic *= jj;
        }
L35:
        c__ = (real) ic;
        i__2 = dvod01_._1.n;
        for (i__ = 1; i__ <= i__2; ++i__) {
            dky[i__] = c__ * yh[i__ + jp1 * yh_dim1] + s * dky[i__];
        }
    }
    if (*k == 0) {
        return 0;
    }
L55:
    i__1 = -(*k);
    r__ = pow_di(&dvod01_._1.h__, &i__1);
    dscal_(&dvod01_._1.n, &r__, &dky[1], &c__1);
    return 0;

L80:
    s_copy(msg, "DVINDY-- K (=I1) illegal      ", (ftnlen)80, (ftnlen)30);
    xerrwd_(msg, &c__30, &c__51, &c__1, &c__1, k, &c__0, &c__0, &zero, &zero,
            (ftnlen)80);
    *iflag = -1;
    return 0;
L90:
    s_copy(msg, "DVINDY-- T (=R1) illegal      ", (ftnlen)80, (ftnlen)30);
    xerrwd_(msg, &c__30, &c__52, &c__1, &c__0, &c__0, &c__0, &c__1, t, &zero,
            (ftnlen)80);
    s_copy(msg, "      T not in interval TCUR - HU (= R1) to TCUR (=R2)      ",
            (ftnlen)80, (ftnlen)60);
    xerrwd_(msg, &c__60, &c__52, &c__1, &c__0, &c__0, &c__0, &c__2, &tp,
            &dvod01_._1.tn, (ftnlen)80);
    *iflag = -2;
    return 0;
}
