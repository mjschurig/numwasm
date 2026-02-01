/* ../reference/netlib/epsode.f -- translated by f2c (version 20240504).
   You must link the resulting object file with libf2c:
	on Microsoft Windows system, link with libf2c.lib;
	on Linux or Unix systems, link with .../path/to/libf2c.a -lm
	or, if you install libf2c.a in a standard place, with -lf2c -lm
	-- in that order, at the end of the command line, as in
		cc *.o -lf2c -lm
	Source for libf2c is in /netlib/f2c/libf2c.zip, e.g.,

		http://netlib.org/f2c/libf2c.zip
*/

#include "f2c.h"

/* Common Block Declarations */

union {
    struct {
	doublereal t, h__, hmin, hmax, epsc, ss, uround;
	integer nc, mfc, kflag, jstart;
    } _1;
    struct {
	doublereal t, h__, hmin, hmax, eps, ss, uround;
	integer n, mf, kflag, jstart;
    } _2;
} epcom1_;

#define epcom1_1 (epcom1_._1)
#define epcom1_2 (epcom1_._2)

union {
    struct {
	doublereal ymax[20];
    } _1;
    struct {
	doublereal ymax[1];
    } _2;
} epcom2_;

#define epcom2_1 (epcom2_._1)
#define epcom2_2 (epcom2_._2)

union {
    struct {
	doublereal error[20];
    } _1;
    struct {
	doublereal error[1];
    } _2;
} epcom3_;

#define epcom3_1 (epcom3_._1)
#define epcom3_2 (epcom3_._2)

union {
    struct {
	doublereal save1[20];
    } _1;
    struct {
	doublereal save1[1];
    } _2;
} epcom4_;

#define epcom4_1 (epcom4_._1)
#define epcom4_2 (epcom4_._2)

union {
    struct {
	doublereal save2[20];
    } _1;
    struct {
	doublereal save2[1];
    } _2;
} epcom5_;

#define epcom5_1 (epcom5_._1)
#define epcom5_2 (epcom5_._2)

union {
    struct {
	doublereal pw[400];
    } _1;
    struct {
	doublereal pw[1];
    } _2;
} epcom6_;

#define epcom6_1 (epcom6_._1)
#define epcom6_2 (epcom6_._2)

union {
    struct {
	integer ipiv[20];
    } _1;
    struct {
	integer ipiv[1];
    } _2;
} epcom7_;

#define epcom7_1 (epcom7_._1)
#define epcom7_2 (epcom7_._2)

struct {
    doublereal epsj;
    integer nsq;
} epcom8_;

#define epcom8_1 epcom8_

struct {
    doublereal hused;
    integer nqused, nstep, nfe, nje;
} epcom9_;

#define epcom9_1 epcom9_

struct {
    doublereal tau[13], el[13], tq[5];
    integer lmax, meth, nq, l, nqindx;
} epcm10_;

#define epcm10_1 epcm10_

/* Table of constant values */

static integer c__4 = 4;
static integer c__1 = 1;

/* Subroutine */ int epsode_(U_fp diffun, U_fp pederv, integer *n, doublereal 
	*t0, doublereal *h0, doublereal *y0, doublereal *tout, doublereal *
	eps, integer *ierror, integer *mf, integer *index)
{
    /* Initialized data */

    static integer lout = 6;
    static doublereal hcut = .1;
    static doublereal four = 4.;
    static doublereal hundrd = 100.;
    static doublereal one = 1.;
    static doublereal ten = 10.;
    static doublereal zero = 0.;

    /* Format strings */
    static char fmt_15[] = "(/\002---  message from subroutine epsode in epi"
	    "sode,\002,\002 the o.d.e. solver.  ---\002/\002 warning.. t + h "
	    "= t =\002,e18.8,\002 in the next step.\002/)";
    static char fmt_101[] = "(/\002---  message from subroutine epsode in ep"
	    "isode,\002,\002 the o.d.e. solver.  ---\002/)";
    static char fmt_105[] = "(//\002 kflag = -1 from integrator at t = \002,"
	    "e18.8/\002  error test failed with abs(h) = hmin =\002,e18.8/)";
    static char fmt_115[] = "(\002  h has been reduced to \002,e18.8,\002  a"
	    "nd step will be retried\002//)";
    static char fmt_155[] = "(//\002 problem appears unsolvable with given i"
	    "nput\002//)";
    static char fmt_205[] = "(//\002 kflag = -2  t =\002,e18.8,\002 h =\002,"
	    "e18.8,\002 eps =\002,e18.8/\002  the requested error is too smal"
	    "l for integrator.\002//)";
    static char fmt_255[] = "(//\002 integration halted by subroutine epsode"
	    " at t =\002,e18.8/\002 eps is too small for machine precision and"
	    "\002/\002 problem being solved.  eps =\002,e18.8//)";
    static char fmt_305[] = "(//\002 kflag = -3 from integrator at t =\002,e"
	    "18.8/\002  corrector convergence could not be achieved\002/)";
    static char fmt_405[] = "(//\002 illegal input.. eps .le. 0. eps = \002,"
	    "e18.8//)";
    static char fmt_415[] = "(//\002 illegal input.. n .le. 0. n = \002,i8//)"
	    ;
    static char fmt_425[] = "(//\002 illegal input.. (t0 - tout)*h0 .ge. 0"
	    ".\002/\002 t0 =\002,e18.8,\002 tout =\002,e18.8,\002 h0 =\002,e1"
	    "8.8//)";
    static char fmt_435[] = "(//\002 illegal input.. index =\002,i8//)";
    static char fmt_445[] = "(//\002 illegal input.  the number of ordinar"
	    "y\002/\002 differential equations being solved is n =\002,i6/"
	    "\002 storage allocation in subroutine epsode is\002/\002 too sma"
	    "ll.  see writeups#epsode public........\002/)";
    static char fmt_455[] = "(//\002 index = -1 on input with (t - tout)*h ."
	    "ge. 0.\002/\002 interpolation was done as on normal return.\002"
	    "/\002 desired parameter changes were not made.\002/\002 t =\002,"
	    "e18.8,\002 tout =\002,e18.8,\002 h =\002,e18.8//)";
    static char fmt_465[] = "(//\002 index = 2 on input with (t - tout)*h .g"
	    "e. 0.\002/\002 t =\002,e18.8,\002 tout =\002,e18.8,\002 h =\002,"
	    "e18.8//)";

    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    double sqrt(doublereal);

    /* Local variables */
    doublereal d__;
    integer i__;
    doublereal y[260]	/* was [20][13] */;
    integer n0;
    doublereal t0p;
    integer kgo;
    doublereal ayi;
    integer nhcut;
    extern doublereal d1mach_(integer *);
    extern /* Subroutine */ int tstep_(U_fp, U_fp, doublereal *, integer *), 
	    interp_(doublereal *, doublereal *, integer *, doublereal *);

    /* Fortran I/O blocks */
    static cilist io___10 = { 0, 0, 0, fmt_15, 0 };
    static cilist io___14 = { 0, 0, 0, fmt_15, 0 };
    static cilist io___18 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___19 = { 0, 0, 0, fmt_105, 0 };
    static cilist io___20 = { 0, 0, 0, fmt_115, 0 };
    static cilist io___21 = { 0, 0, 0, fmt_155, 0 };
    static cilist io___22 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___23 = { 0, 0, 0, fmt_205, 0 };
    static cilist io___24 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___25 = { 0, 0, 0, fmt_255, 0 };
    static cilist io___26 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___27 = { 0, 0, 0, fmt_305, 0 };
    static cilist io___28 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___29 = { 0, 0, 0, fmt_405, 0 };
    static cilist io___30 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___31 = { 0, 0, 0, fmt_415, 0 };
    static cilist io___32 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___33 = { 0, 0, 0, fmt_425, 0 };
    static cilist io___34 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___35 = { 0, 0, 0, fmt_435, 0 };
    static cilist io___36 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___37 = { 0, 0, 0, fmt_445, 0 };
    static cilist io___38 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___39 = { 0, 0, 0, fmt_455, 0 };
    static cilist io___40 = { 0, 0, 0, fmt_101, 0 };
    static cilist io___41 = { 0, 0, 0, fmt_465, 0 };


/* this is the june 24, 1975  version of */
/* episode..  experimental package for integration of */
/* systems of ordinary differential equations, */
/*    dy/dt = f(y,t),  y = (y(1),y(2),...,y(n)) transpose, */
/* given the initial value of y. */
/* this code is for the ibm 370/195 at argonne national laboratory */
/* and is a modification of earlier versions by g.d.byrne */
/* and a.c.hindmarsh. */
/*                 references */
/* 1.  g. d. byrne and a. c. hindmarsh, a polyalgorithm for the */
/*       numerical solution of ordinary differential equations, */
/*       ucrl-75652, lawrence livermore laboratory, p. o. box 808, */
/*       livermore, ca 94550, april 1974. also in acm transactions */
/*       on mathematical software, 1 (1975), pp. 71-96. */

/* 2.  a. c. hindmarsh and g. d. byrne, episode.. an experimental */
/*       package for the integration of systems of ordinary */
/*       differential equations, ucid-30112, l.l.l., may, 1975. */

/* 3.  a. c. hindmarsh, gear.. ordinary differential equation */
/*       system solver, ucid-30001, rev. 3, l.l.l., december, 1974. */

/* ----------------------------------------------------------------------- */
/* epsode is a driver subroutine for the episode package. */
/* epsode is to be called once for each output value of t. */
/* it then makes repeated calls to the core integrator */
/* subroutine, tstep. */

/* the input parameters are as follows. */
/*   diffun=  name of double precision subroutine diffun(n,t,y,ydot) */
/*            which the user must supply and declare external. see */
/*            note below. */
/*   pederv=  name of double precision subroutine pederv(n,t,y,pd,n0) */
/*            which the user must supply and declare external. see */
/*            note below. */
/*   n     =  the number of differential equations (used only on */
/*              first call, unless index = -1).  n must never be */
/*              increased during a given problem. */
/*   t0    =  the initial value of t, the independent variable */
/*              (used for input only on first call). */
/*   h0    =  the step size in t (used for input only on the */
/*              first call, unless index = 3 on input).  when */
/*              index = 3, h0 is the maximum absolute value of */
/*              the step size to be used. */
/*   y0    =  a vector of length n containing the initial values of */
/*              y (used for input only on first call). */
/*   tout  =  the value of t at which output is desired next. */
/*              integration will normally go beyond tout and */
/*              interpolate to t = tout.  (used only for input.) */
/*   eps   =  the relative error bound (used only on first call, */
/*              unless index = -1).  this bound is used as follows. */
/*              let r(i) denote the estimated relative local error */
/*              in y(i), i.e. the error relative to ymax(i), as */
/*              measured per step (of size h) or per ss units of t. */
/*              then eps is a bound on the root-mean-square norm */
/*              of the vector r, i.e. */
/*                      n */
/*              sqrt ( sum ( r(i)**2 )/n ) .lt. eps. */
/*                     i=1 */
/*              the vector ymax is computed in epsode as described */
/*              under ierror below. */
/*              if error control per ss units of t is desired, set ss */
/*              to a positive number after statement 10 (where it is */
/*              now set to zero) and update it after statement 60. */
/*              see also the comments on ss and ymax below. */
/*  ierror =  the error flag with values and meanings as follow. */
/*         1    absolute error is controlled.  ymax(i) = 1.0. */
/*         2    error relative to abs(y) is controlled.  if y(i) = 0.0 */
/*              a divide error will occur.  ymax(i) = abs(y(i)). */
/*         3    error relative to the largest value of abs(y(i)) seen */
/*              so far is controlled.  if the initial value of y(i) is */
/*              0.0, then ymax(i) is set to 1.0 initially and remains */
/*              at least 1.0. */
/*   mf    =  the method flag (used only on first call, unless */
/*              index = -1).  allowed values are 10, 11, 12, 13, */
/*              20, 21, 22, 23.  mf is an integer with two decimal */
/*              digits, meth and miter (mf = 10*meth + miter).  (mf */
/*              can be thought of as the ordered pair (meth,miter).) */
/*              meth is the basic method indicator. */
/*                meth = 1 indicates variable-step size, variable- */
/*                         order adams method, suitable for non- */
/*                         stiff problems. */
/*                meth = 2 indicates variable-step size, variable- */
/*                         order backward differentiation method, */
/*                         suitable for stiff problems. */
/*              miter indicates the method of iterative correction */
/*                (nonlinear system solution). */
/*                miter = 0 indicates functional iteration (no */
/*                          partial derivatives needed). */
/*                miter = 1 indicates a chord or semi-stationary */
/*                          newton method with closed form (exact) */
/*                          jacobian, which is computed in the */
/*                          user supplied subroutine */
/*                          pederv(n,t,y,pd,n0) described below. */
/*                miter = 2 indicates a chord  or semi-stationary */
/*                          newton method with an internally */
/*                          computed finite difference approximation */
/*                          to the jacobian. */
/*                miter = 3 indicates a chord or semi-stationary */
/*                          newton method with an internally */
/*                          computed diagonal matrix approximation */
/*                          to the jacobian, based on a directional */
/*                          derivative. */
/*   index =  integer used on input to indicate type of call, */
/*              with the following values and meanings.. */
/*         1    this is the first call for this problem. */
/*         0    this is not the first call for this problem, */
/*              and integration is to continue. */
/*        -1    this is not the first call for the problem, */
/*              and the user has reset n, eps, and/or mf. */
/*         2    same as 0 except that tout is to be hit */
/*              exactly (no interpolation is done). */
/*              assumes tout .ge. the current t. */
/*         3    same as 0 except control returns to calling */
/*              program after one step.  tout is ignored. */
/*            since the normal output value of index is 0, */
/*            it need not be reset for normal continuation. */

/* after the initial call, if a normal return occurred and a normal */
/* continuation is desired, simply reset tout and call again. */
/* all other parameters will be ready for the next call. */
/* a change of parameters with index = -1 can be made after */
/* either a successful or an unsuccessful return. */

/* the output parameters are as follows. */
/*   t0    =  the output value of t.  if integration was successful, */
/*            t0 = tout.  otherwise, t0 is the last value of t */
/*            reached successfully. */
/*   h0    =  the step size h used last, whether successfully or not. */
/*   y0    =  the computed values of y at t = t0. */
/*   index =  integer used on output to indicate results, */
/*              with the following values and meanings.. */
/*         0    integration was completed to tout or beyond. */
/*        -1    the integration was halted after failing to pass the */
/*              error test even after reducing h by a factor of */
/*              1.e10 from its initial value. */
/*        -2    after some initial success, the integration was */
/*              halted either by repeated error test failures or */
/*              by a test on eps.  possibly too much accuracy has */
/*              been requested, or a bad choice of mf was made. */
/*        -3    the integration was halted after failing to achieve */
/*              corrector convergence even after reducing h by a */
/*              factor of 1.e10 from its initial value. */
/*        -4    immediate halt because of illegal values of input */
/*              parameters.  see printed message. */
/*        -5    index was -1 on input, but the desired changes of */
/*              parameters were not implemented because tout */
/*              was not beyond t.  interpolation to t = tout was */
/*              performed as on a normal return.  to continue, */
/*              simply call again with index = -1 and a new tout. */
/*        -6    index was 2 on input, but tout was not beyond t. */
/*              no action was taken. */

/* in addition to epsode, the following subroutines are used by and */
/* provided in this package: */
/*   interp(tout,y,n0,y0)  interpolates to give output values at */
/*                         t = tout  by using data in the y array. */
/*   tstep(diffun,pederv,y,n0)  is the core integration subroutine, which */
/*                integrates over a single step and does associated error */
/*                control. */
/*   coset  sets coefficients for use in tstep. */
/*   adjust(y,n0)  adjusts the history array y on reduction of order. */
/*   pset(diffun,pederv,y,n0,con,miter,ier)  computes and processes the */
/*                             jacobian matrix, j = df/dy. */
/*   dec(n,n0,a,ip,ier)  performs the lu decomposition of a matrix. */
/*   sol(n,n0,a,b,ip)  solves a linear system a*x = b, after dec */
/*                       has been called for the matrix a. */
/* note:  pset, dec, and sol are called if and only if miter = 1 */
/*        or miter = 2. */

/* the user must furnish the following double precision subroutines */
/*      and declare them external in his calling program: */
/*   diffun(n,t,y,ydot)  computes the function  ydot = f(y,t), */
/*                       the right hand side of the ordinary */
/*                       differential equation system, where y */
/*                       and ydot are vectors of length n. */
/*   pederv(n,t,y,pd,n0)  computes the n by n jacobian matrix of */
/*                        partial derivatives and stores it in pd as */
/*                        an n0 by n0 array.   pd(i,j) is to be set */
/*                        to the partial derivative of ydot(i) with */
/*                        respect to y(j).  pederv is called if and */
/*                        only if miter = 1. for other values of */
/*                        miter, pederv can be a dummy subroutine. */

/* caution:  at the present time the maximum number of differential */
/*           equations, which can be solved by episode, is 20.  to */
/*           change this number to a new value, say nmax, change */
/*           y(20,13) to y(nmax,13), ymax(20) to ymax(nmax), */
/*           error(20) to error(nmax), save1(20) to save1(nmax), */
/*           save2(20) to save2(nmax), pw(400) to pw(nmax*nmax), */
/*           and ipiv(20) to ipiv(nmax) in the common and dimension */
/*           statements below.  also change the argument in the */
/*           if...go to 440 statement (after the common statements) */
/*           from 20 to nmax.  no other changes need to be made to */
/*           any other subroutine in this package when the maximum */
/*           number of equations is changed.  elsewhere, the column */
/*           length of the y array is n0 instead of 20.  the row */
/*           length of y can be reduced from 13 to 6 if meth = 2. */
/*           the array ipiv is used if and only if miter = 1 or */
/*           miter = 2.  the size of the pw array can be reduced */
/*           to 1 if miter = 0 or to n if miter = 3. */

/* the common block epcom2 can be accessed externally by the user, */
/* if he desires to modify the error control behaviour of the package. */
/* the block contains ymax(20) used to control the error. see above */
/* and reference (2). */

/* the common block epcom9 can be accessed externally by the user, */
/* if he desires.  it contains the step size last used successfully */
/* (hused), the order last used successfully (nqused), the */
/* number of steps taken so far (nstep), the number of function */
/* evaluations (diffun calls) so far (nfe), and the number of */
/* jacobian evaluations so far (nje). */

/* in a data statement below, lout is set to the logical unit number */
/* for the output of messages during integration.  currently, lout */
/* = 6. */
/* ----------------------------------------------------------------------- */



    /* Parameter adjustments */
    --y0;

    /* Function Body */
    if (*index == 0) {
	goto L20;
    }
    if (*index == 2) {
	goto L25;
    }
    if (*index == -1) {
	goto L30;
    }
    if (*index == 3) {
	goto L40;
    }
    if (*index != 1) {
	goto L430;
    }
    if (*eps <= zero) {
	goto L400;
    }
    if (*n <= 0) {
	goto L410;
    }
    if (*n > 20) {
	goto L440;
    }
    if ((*t0 - *tout) * *h0 >= zero) {
	goto L420;
    }
/* ----------------------------------------------------------------------- */
/* if initial values for ymax other than those below are desired, */
/* they should be set here.  all ymax(i) must be positive.  if */
/* values for hmin or hmax, the bounds on the absolute value of h, */
/* other than those below, are desired, they also should be set here. */
/* if error per ss units of t is to be controlled, ss should be set */
/* to a positive value below.  error per unit step is controlled */
/* when ss = 1.  the default value for ss is 0 and yields control */
/* of error per step. */
/* ----------------------------------------------------------------------- */
/* set uround, the machine roundoff constant, here. */
/* use statement below for short precision on ibm 360 or 370. */
/*     uround = 9.53674e-7 */
/* use statement below for single precision on cdc 7600 or 6600. */
/*     uround = 7.105427406e-15 */
/* use statement below for long precision on ibm 360 or 370. */
/*     uround = 2.220446047d-16 */
    epcom1_1.uround = d1mach_(&c__4);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	switch (*ierror) {
	    case 1:  goto L5;
	    case 2:  goto L6;
	    case 3:  goto L7;
	}
/* ierror   =   1, 2, 3  ------------------------------------------------ */
L5:
	epcom2_1.ymax[i__ - 1] = one;
	goto L10;
L6:
	epcom2_1.ymax[i__ - 1] = (d__1 = y0[i__], abs(d__1));
	goto L10;
L7:
	epcom2_1.ymax[i__ - 1] = (d__1 = y0[i__], abs(d__1));
	if (epcom2_1.ymax[i__ - 1] == zero) {
	    epcom2_1.ymax[i__ - 1] = one;
	}
L10:
	y[i__ - 1] = y0[i__];
    }
    epcom1_1.nc = *n;
    epcom1_1.t = *t0;
    epcom1_1.h__ = *h0;
    if (epcom1_1.t + epcom1_1.h__ == epcom1_1.t) {
	io___10.ciunit = lout;
	s_wsfe(&io___10);
	do_fio(&c__1, (char *)&epcom1_1.t, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    epcom1_1.hmin = abs(*h0);
    epcom1_1.hmax = (d__1 = *t0 - *tout, abs(d__1)) * ten;
    epcom1_1.epsc = *eps;
    epcom1_1.mfc = *mf;
    epcom1_1.jstart = 0;
    epcom1_1.ss = zero;
    n0 = *n;
    epcom8_1.nsq = n0 * n0;
    epcom8_1.epsj = sqrt(epcom1_1.uround);
    nhcut = 0;
    goto L50;
/* t0p is the previous output value of t0 for use in hmax. -------------- */
L20:
    epcom1_1.hmax = (d__1 = *tout - t0p, abs(d__1)) * ten;
    goto L80;
L25:
    epcom1_1.hmax = (d__1 = *tout - t0p, abs(d__1)) * ten;

    if ((epcom1_1.t - *tout) * epcom1_1.h__ >= zero) {
	goto L460;
    }
    goto L85;

L30:
    if ((epcom1_1.t - *tout) * epcom1_1.h__ >= zero) {
	goto L450;
    }
    if (*mf != epcom1_1.mfc) {
	epcom1_1.jstart = -1;
    }
    epcom1_1.nc = *n;
    epcom1_1.epsc = *eps;
    epcom1_1.mfc = *mf;
    goto L45;

L40:
    epcom1_1.hmax = *h0;

L45:
    if (epcom1_1.t + epcom1_1.h__ == epcom1_1.t) {
	io___14.ciunit = lout;
	s_wsfe(&io___14);
	do_fio(&c__1, (char *)&epcom1_1.t, (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

L50:
    tstep_((U_fp)diffun, (U_fp)pederv, y, &n0);

    kgo = 1 - epcom1_1.kflag;
    switch (kgo) {
	case 1:  goto L60;
	case 2:  goto L100;
	case 3:  goto L200;
	case 4:  goto L300;
    }
/* kflag  =   0,  -1,  -2,  -3  ----------------------------------------- */

L60:
/* ----------------------------------------------------------------------- */
/* normal return from tstep. */

/* the weights ymax(i) are updated.  if different values are desired, */
/* they should be set here.  if ss is to be updated for control of */
/* error per ss units of t, it should also be done here.  a test is */
/* made to determine if eps is too small for machine precision. */

/* any other tests or calculations that are required after each step */
/* should be inserted here. */

/* if index = 3, y0 is set to the current y values on return. */
/* if index = 2, h is controlled to hit tout (within roundoff */
/* error), and then the current y values are put in y0 on */
/* return.  for any other value of index, control returns to */
/* the integrator unless tout has been reached.  then */
/* interpolated values of y are computed and stored in y0 on */
/* return. */
/* if interpolation is not desired, the call to interp should */
/* be deleted and control transferred to statement 500 instead */
/* of 520. */
/* ----------------------------------------------------------------------- */
    d__ = zero;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ayi = (d__1 = y[i__ - 1], abs(d__1));
	switch (*ierror) {
	    case 1:  goto L70;
	    case 2:  goto L66;
	    case 3:  goto L67;
	}
/* ierror  =     1,  2,  3  --------------------------------------------- */
L66:
	epcom2_1.ymax[i__ - 1] = ayi;
	goto L70;
L67:
/* Computing MAX */
	d__1 = epcom2_1.ymax[i__ - 1];
	epcom2_1.ymax[i__ - 1] = max(d__1,ayi);
L70:
/* Computing 2nd power */
	d__1 = ayi / epcom2_1.ymax[i__ - 1];
	d__ += d__1 * d__1;
    }
/* Computing 2nd power */
    d__1 = epcom1_1.uround / *eps;
    d__ *= d__1 * d__1;
    if (d__ > (doublereal) (*n)) {
	goto L250;
    }
    if (*index == 3) {
	goto L500;
    }
    if (*index == 2) {
	goto L85;
    }
L80:
    if ((epcom1_1.t - *tout) * epcom1_1.h__ < zero) {
	goto L45;
    }
    interp_(tout, y, &n0, &y0[1]);
    *t0 = *tout;
    goto L520;
L85:
    if ((epcom1_1.t + epcom1_1.h__ - *tout) * epcom1_1.h__ <= zero) {
	goto L45;
    }
    if ((d__1 = epcom1_1.t - *tout, abs(d__1)) <= hundrd * epcom1_1.uround * 
	    epcom1_1.hmax) {
	goto L500;
    }
    if ((epcom1_1.t - *tout) * epcom1_1.h__ >= zero) {
	goto L500;
    }
    epcom1_1.h__ = (*tout - epcom1_1.t) * (one - four * epcom1_1.uround);
    epcom1_1.jstart = -1;
    goto L45;
/* ----------------------------------------------------------------------- */
/* on an error return from tstep, an immediate return occurs if */
/* kflag = -2, and recovery attempts are made otherwise. */
/* to recover, h and hmin are reduced by a factor of .1 up to 10 */
/* times before giving up. */
/* ----------------------------------------------------------------------- */
L100:
    io___18.ciunit = lout;
    s_wsfe(&io___18);
    e_wsfe();
    io___19.ciunit = lout;
    s_wsfe(&io___19);
    do_fio(&c__1, (char *)&epcom1_1.t, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&epcom1_1.hmin, (ftnlen)sizeof(doublereal));
    e_wsfe();
L110:
    if (nhcut == 10) {
	goto L150;
    }
    ++nhcut;
    epcom1_1.hmin = hcut * epcom1_1.hmin;
    epcom1_1.h__ = hcut * epcom1_1.h__;
    io___20.ciunit = lout;
    s_wsfe(&io___20);
    do_fio(&c__1, (char *)&epcom1_1.h__, (ftnlen)sizeof(doublereal));
    e_wsfe();
    epcom1_1.jstart = -1;
    goto L45;

L150:
    io___21.ciunit = lout;
    s_wsfe(&io___21);
    e_wsfe();
    goto L500;

L200:
    io___22.ciunit = lout;
    s_wsfe(&io___22);
    e_wsfe();
    io___23.ciunit = lout;
    s_wsfe(&io___23);
    do_fio(&c__1, (char *)&epcom1_1.t, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&epcom1_1.h__, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&(*eps), (ftnlen)sizeof(doublereal));
    e_wsfe();
    goto L500;

L250:
    io___24.ciunit = lout;
    s_wsfe(&io___24);
    e_wsfe();
    io___25.ciunit = lout;
    s_wsfe(&io___25);
    do_fio(&c__1, (char *)&epcom1_1.t, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&(*eps), (ftnlen)sizeof(doublereal));
    e_wsfe();
    epcom1_1.kflag = -2;
    goto L500;

L300:
    io___26.ciunit = lout;
    s_wsfe(&io___26);
    e_wsfe();
    io___27.ciunit = lout;
    s_wsfe(&io___27);
    do_fio(&c__1, (char *)&epcom1_1.t, (ftnlen)sizeof(doublereal));
    e_wsfe();
    goto L110;

L400:
    io___28.ciunit = lout;
    s_wsfe(&io___28);
    e_wsfe();
    io___29.ciunit = lout;
    s_wsfe(&io___29);
    do_fio(&c__1, (char *)&(*eps), (ftnlen)sizeof(doublereal));
    e_wsfe();
    *index = -4;
    return 0;

L410:
    io___30.ciunit = lout;
    s_wsfe(&io___30);
    e_wsfe();
    io___31.ciunit = lout;
    s_wsfe(&io___31);
    do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
    e_wsfe();
    *index = -4;
    return 0;

L420:
    io___32.ciunit = lout;
    s_wsfe(&io___32);
    e_wsfe();
    io___33.ciunit = lout;
    s_wsfe(&io___33);
    do_fio(&c__1, (char *)&(*t0), (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&(*tout), (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&(*h0), (ftnlen)sizeof(doublereal));
    e_wsfe();
    *index = -4;
    return 0;

L430:
    io___34.ciunit = lout;
    s_wsfe(&io___34);
    e_wsfe();
    io___35.ciunit = lout;
    s_wsfe(&io___35);
    do_fio(&c__1, (char *)&(*index), (ftnlen)sizeof(integer));
    e_wsfe();
    *index = -4;
    return 0;

L440:
    io___36.ciunit = lout;
    s_wsfe(&io___36);
    e_wsfe();
    io___37.ciunit = lout;
    s_wsfe(&io___37);
    do_fio(&c__1, (char *)&(*n), (ftnlen)sizeof(integer));
    e_wsfe();
    *index = -4;
    return 0;

L450:
    io___38.ciunit = lout;
    s_wsfe(&io___38);
    e_wsfe();
    io___39.ciunit = lout;
    s_wsfe(&io___39);
    do_fio(&c__1, (char *)&epcom1_1.t, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&(*tout), (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&epcom1_1.h__, (ftnlen)sizeof(doublereal));
    e_wsfe();
    interp_(tout, y, &n0, &y0[1]);
    *t0 = *tout;
    *index = -5;
    return 0;

L460:
    io___40.ciunit = lout;
    s_wsfe(&io___40);
    e_wsfe();
    io___41.ciunit = lout;
    s_wsfe(&io___41);
    do_fio(&c__1, (char *)&epcom1_1.t, (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&(*tout), (ftnlen)sizeof(doublereal));
    do_fio(&c__1, (char *)&epcom1_1.h__, (ftnlen)sizeof(doublereal));
    e_wsfe();
    *index = -6;
    return 0;

L500:
    *t0 = epcom1_1.t;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L510: */
	y0[i__] = y[i__ - 1];
    }
L520:
    *index = epcom1_1.kflag;
    t0p = *t0;
    *h0 = epcom9_1.hused;
    if (epcom1_1.kflag != 0) {
	*h0 = epcom1_1.h__;
    }
    return 0;
/* -----------------------  end of subroutine epsode --------------------- */
} /* epsode_ */

/* Subroutine */ int interp_(doublereal *tout, doublereal *y, integer *n0, 
	doublereal *y0)
{
    /* Initialized data */

    static doublereal one = 1.;

    /* System generated locals */
    integer y_dim1, y_offset, i__1, i__2;

    /* Local variables */
    integer i__, j, l;
    doublereal s, s1;

/* ----------------------------------------------------------------------- */
/* subroutine interp computes interpolated values of the dependent */
/* variable y and stores them in y0.  the interpolation is to the */
/* point t = tout and uses the nordsieck history array y as follows.. */
/*                             nq */
/*                  y0(i)  =  sum  y(i,j+1)*s**j , */
/*                            j=0 */
/* where s = -(t-tout)/h. */
/* ----------------------------------------------------------------------- */
/* caution:  not all members of epcom1 are used in this subroutine. */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --y0;
    y_dim1 = *n0;
    y_offset = 1 + y_dim1;
    y -= y_offset;

    /* Function Body */
    i__1 = epcom1_2.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	y0[i__] = y[i__ + y_dim1];
    }
    l = epcom1_2.jstart + 1;
    s = (*tout - epcom1_2.t) / epcom1_2.h__;
    s1 = one;
    i__1 = l;
    for (j = 2; j <= i__1; ++j) {
	s1 *= s;
	i__2 = epcom1_2.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L20: */
	    y0[i__] += s1 * y[i__ + j * y_dim1];
	}
/* L30: */
    }
    return 0;
/* -----------------------  end of subroutine interp  -------------------- */
} /* interp_ */

/* Subroutine */ int tstep_(S_fp diffun, U_fp pederv, doublereal *y, integer *
	n0)
{
    /* Initialized data */

    static integer istepj = 20;
    static integer kfc = -3;
    static integer kfh = -7;
    static integer maxcor = 3;
    static doublereal addon = 1e-6;
    static doublereal bias1 = 25.;
    static doublereal bias2 = 25.;
    static doublereal bias3 = 100.;
    static doublereal crdown = .1;
    static doublereal delrc = .3;
    static doublereal etacf = .25;
    static doublereal etamin = .1;
    static doublereal etamxf = .2;
    static doublereal etamx1 = 1e4;
    static doublereal etamx2 = 10.;
    static doublereal etamx3 = 1.5;
    static doublereal onepsm = 1.00001;
    static doublereal short__ = .1;
    static doublereal thresh = 1.3;
    static doublereal one = 1.;
    static doublereal pt5 = .5;
    static doublereal zero = 0.;

    /* System generated locals */
    integer y_dim1, y_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), pow_di(doublereal *, integer *)
	    ;

    /* Local variables */
    doublereal d__, e;
    integer i__, j, m;
    doublereal r__, d1;
    integer j1, j2;
    doublereal r0, r1, rc, rl1, bnd, edn, drc, eta;
    integer ier;
    doublereal con;
    integer mio;
    doublereal eup;
    extern /* Subroutine */ int sol_(integer *, integer *, doublereal *, 
	    doublereal *, integer *);
    doublereal hrl1, prl1, hold, etaq, conp, told;
    integer newj;
    extern /* Subroutine */ int pset_(S_fp, U_fp, doublereal *, integer *, 
	    doublereal *, integer *, integer *);
    doublereal phrl1;
    integer iback;
    doublereal crate;
    integer mfold, iredo;
    extern /* Subroutine */ int coset_(void);
    integer miter;
    doublereal flotl, flotn, etaqm1, etaqp1;
    integer miter1;
    doublereal etamax;
    extern /* Subroutine */ int adjust_(doublereal *, integer *);
    integer nstepj;
    doublereal cnquot;

/* ----------------------------------------------------------------------- */
/* tstep performs one step of the integration of an initial value */
/* problem for a system of ordinary differential equations. */
/* communication with tstep is via the following variables.. */

/*   y       an n0 by lmax array containing the dependent variables */
/*             and their scaled derivatives.  lmax is currently 6 for */
/*             the variable step backward differentiation formulas, */
/*             and 13 for the variable step adams formulas. */
/*             (lmax -1) = maxder, the maximum order used. */
/*             see subroutine coset.  y(i,j+1) contains the */
/*             j-th derivative of y(i), scaled by h**j/factorial(j) */
/*             for j = 0,1,...,nq, where nq is the current order. */
/*   n0      a constant integer .ge. n, used for dimensioning */
/*             purposes. */
/*   t       the independent variable, updated on each step taken. */
/*   h       the step size to be attempted on the next step. */
/*             h is altered by the error control algorithm during */
/*             the solution of the problem. h can be either positive */
/*             or negative, but its sign must remain constant */
/*             throughout the problem run. */
/*   hmin,   the minimum and maximum absolute values of the step */
/*    hmax     size to be used for the step. these may be changed at */
/*             any time, but the change will not take effect until the */
/*             next change in h is made. */
/*   eps     the relative error bound.  see description in */
/*             subroutine epsode. */
/*   ss      the size of the time interval to be used for error */
/*             control.  a default value of 0 is used to produce */
/*             control of error per step.  see subroutine epsode. */
/*   uround  the unit of roundoff for the computer being used. */
/*   n       the number of first order ordinary differential */
/*             equations being solved. */
/*   mf      the method flag.  see description in subroutine epsode. */
/*   kflag   a completion code with the following meanings.. */
/*                     0  the step was succesful. */
/*                    -1  the requested error could not be achieved */
/*                          with abs(h) = hmin. */
/*                    -2  the requested error is smaller than can */
/*                          be handled for this problem. */
/*                    -3  corrector convergence could not be */
/*                          achieved for abs(h) = hmin. */
/*             on a return with kflag negative, the values of t and */
/*             the y array are as of the beginning of the last */
/*             step and h is the last step size attempted. */
/*   jstart  an integer used on input and output. */
/*             on input, it has the following values and meanings.. */
/*                     0  perform the first step. */
/*                 .gt.0  take a new step continuing from the last. */
/*                 .lt.0  take the next step with a new value of */
/*                          h and/or mf. */
/*             on exit, jstart is set to nq, the current order of the */
/*             method. */
/*   ymax      an array of n elements with which the estimated local */
/*               errors in y are compared. */
/*   error     an array of n elements.  error(i)/tq(2) is the */
/*               estimated local error in y(i) per ss units of */
/*               t or per step (of size h). */
/*   save1,    two arrays for working storage, */
/*     save2     each of length n. */
/*   pw        a block of locations used for the partial derivatives */
/*               of f with respect to y, if miter is not o.  see */
/*               description in subroutine epsode. */
/*   ipiv      an integer array of length n, which is used for pivot */
/*               information for the linear algebraic system in the */
/*               correction process, when miter = 1 or 2. */

/* the common block epcm10, declared below, is primarily intended */
/* for internal use, but it can be accessed externally. */
/* ----------------------------------------------------------------------- */


    /* Parameter adjustments */
    y_dim1 = *n0;
    y_offset = 1 + y_dim1;
    y -= y_offset;

    /* Function Body */
    epcom1_2.kflag = 0;
    told = epcom1_2.t;
    flotn = (doublereal) epcom1_2.n;
    if (epcom1_2.jstart > 0) {
	goto L200;
    }
    if (epcom1_2.jstart != 0) {
	goto L150;
    }
/* ----------------------------------------------------------------------- */
/* on the first call, the order is set to 1 and the initial */
/* derivatives are calculated.  etamax is the maximum ratio by */
/* which h can be increased in a single step.  it is 1.e04 for the */
/* first step to compensate for the small initial h, then 10 for */
/* the next 10 steps, and then 1.5 thereafter.  if a failure */
/* occurs (in corrector convergence or error test), etamax is set at 1 */
/* for the next increase.  etamin = .1 is the minimum ratio by which */
/* h can be reduced on any retry of a step. */
/* ----------------------------------------------------------------------- */
    (*diffun)(&epcom1_2.n, &epcom1_2.t, &y[y_offset], epcom4_2.save1);
    i__1 = epcom1_2.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L110: */
	y[i__ + (y_dim1 << 1)] = epcom1_2.h__ * epcom4_2.save1[i__ - 1];
    }
    epcm10_1.meth = epcom1_2.mf / 10;
    miter = epcom1_2.mf - epcm10_1.meth * 10;
    miter1 = miter + 1;
    mfold = epcom1_2.mf;
    epcm10_1.nq = 1;
    epcm10_1.l = 2;
    epcm10_1.tau[0] = epcom1_2.h__;
    prl1 = one;
    rc = zero;
    etamax = etamx1;
    epcm10_1.nqindx = 2;
    epcom9_1.nstep = 0;
    nstepj = 0;
    epcom9_1.nfe = 1;
    epcom9_1.nje = 0;
    goto L200;
/* ----------------------------------------------------------------------- */
/* if the user has changed h, then y must be rescaled.  if the */
/* user has changed miter, then newj is set to miter to force */
/* the partial derivativees to be updated, if they are being used. */
/* ----------------------------------------------------------------------- */
L150:
    if (epcom1_2.mf == mfold) {
	goto L170;
    }
    mio = miter;
    epcm10_1.meth = epcom1_2.mf / 10;
    miter = epcom1_2.mf - epcm10_1.meth * 10;
    mfold = epcom1_2.mf;
    if (miter == mio) {
	goto L170;
    }
    newj = miter;
    miter1 = miter + 1;
L170:
    if (epcom1_2.h__ == hold) {
	goto L200;
    }
    eta = epcom1_2.h__ / hold;
    epcom1_2.h__ = hold;
    iredo = 3;
    goto L185;
L180:
/* Computing MAX */
    d__1 = eta, d__2 = epcom1_2.hmin / abs(epcom1_2.h__), d__1 = max(d__1,
	    d__2);
    eta = max(d__1,etamin);
L185:
/* Computing MIN */
    d__1 = eta, d__2 = epcom1_2.hmax / abs(epcom1_2.h__), d__1 = min(d__1,
	    d__2);
    eta = min(d__1,etamax);
    r1 = one;
    i__1 = epcm10_1.l;
    for (j = 2; j <= i__1; ++j) {
	r1 *= eta;
	i__2 = epcom1_2.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L190: */
	    y[i__ + j * y_dim1] *= r1;
	}
    }
    epcom1_2.h__ *= eta;
    rc *= eta;
    if (iredo == 0) {
	goto L690;
    }
/* ----------------------------------------------------------------------- */
/* this section computes the predicted values by effectively */
/* multiplying the y array by the pascal triangle matrix.  then */
/* coset is called to obtain el, the vector of coefficients of */
/* length nq + 1.  rc is the ratio of new to old values of the */
/* coefficient h/el(2).  when rc differs from 1 by more than */
/* delrc, newj is set to miter to force the partial derivatives */
/* to be updated, if used.  delrc is 0.3.  in any case, the partial */
/* derivatives are updated at least every 20-th step. */
/* ----------------------------------------------------------------------- */
L200:
    epcom1_2.t += epcom1_2.h__;
    i__2 = epcm10_1.nq;
    for (j1 = 1; j1 <= i__2; ++j1) {
	i__1 = epcm10_1.nq;
	for (j2 = j1; j2 <= i__1; ++j2) {
	    j = epcm10_1.nq + j1 - j2;
	    i__3 = epcom1_2.n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/* L210: */
		y[i__ + j * y_dim1] += y[i__ + (j + 1) * y_dim1];
	    }
	}
    }
    coset_();
/* Computing 2nd power */
    d__1 = epcm10_1.tq[3] * epcom1_2.eps;
    bnd = flotn * (d__1 * d__1);
    rl1 = one / epcm10_1.el[1];
    rc *= rl1 / prl1;
    prl1 = rl1;
    if (epcom9_1.nstep >= nstepj + istepj) {
	newj = miter;
    }
    drc = (d__1 = rc - one, abs(d__1));
    if (drc <= delrc) {
	goto L215;
    }
    newj = miter;
    crate = one;
    rc = one;
    goto L220;
L215:
    if (miter != 0 && drc != zero) {
	crate = one;
    }
/* ----------------------------------------------------------------------- */
/* up to 3 corrector iterations are taken.  a convergence test is made */
/* on the root mean square norm of each correction, using bnd, which */
/* is dependent on eps.  the sum of the corrections is accumulated in */
/* the vector error.  the y array is not altered in the corrector */
/* loop.  the updated y vector is stored temporarily in save1. */
/* ----------------------------------------------------------------------- */
L220:
    i__3 = epcom1_2.n;
    for (i__ = 1; i__ <= i__3; ++i__) {
/* L230: */
	epcom3_2.error[i__ - 1] = zero;
    }
    m = 0;
    (*diffun)(&epcom1_2.n, &epcom1_2.t, &y[y_offset], epcom5_2.save2);
    ++epcom9_1.nfe;
    if (newj <= 0) {
	goto L290;
    }
/* ----------------------------------------------------------------------- */
/* if indicated, the matrix p = i - h*rl1*j is reevaluated before */
/* starting the corrector iteration.  newj is set to 0 as an */
/* indicator that this has been done.  if miter = 1 or 2, p is */
/* computed and processed in pset.  if miter = 3, the matrix  is */
/* p = i - h*rl1*d, where d is a diagonal matrix.  rl1 is 1/el(2). */
/* ----------------------------------------------------------------------- */
    newj = 0;
    rc = one;
    ++epcom9_1.nje;
    nstepj = epcom9_1.nstep;
    switch (miter) {
	case 1:  goto L250;
	case 2:  goto L240;
	case 3:  goto L260;
    }
L240:
    epcom9_1.nfe += epcom1_2.n;
L250:
    con = -epcom1_2.h__ * rl1;
    pset_((S_fp)diffun, (U_fp)pederv, &y[y_offset], n0, &con, &miter, &ier);
    if (ier != 0) {
	goto L420;
    }
    goto L350;
L260:
    r__ = rl1 * short__;
    i__3 = epcom1_2.n;
    for (i__ = 1; i__ <= i__3; ++i__) {
/* L270: */
	epcom6_2.pw[i__ - 1] = y[i__ + y_dim1] + r__ * (epcom1_2.h__ * 
		epcom5_2.save2[i__ - 1] - y[i__ + (y_dim1 << 1)]);
    }
    (*diffun)(&epcom1_2.n, &epcom1_2.t, epcom6_2.pw, epcom4_2.save1);
    ++epcom9_1.nfe;
    hrl1 = epcom1_2.h__ * rl1;
    i__3 = epcom1_2.n;
    for (i__ = 1; i__ <= i__3; ++i__) {
	r0 = epcom1_2.h__ * epcom5_2.save2[i__ - 1] - y[i__ + (y_dim1 << 1)];
	epcom6_2.pw[i__ - 1] = one;
	d__ = short__ * r0 - epcom1_2.h__ * (epcom4_2.save1[i__ - 1] - 
		epcom5_2.save2[i__ - 1]);
	epcom4_2.save1[i__ - 1] = zero;
	if (abs(r0) < epcom1_2.uround * epcom2_2.ymax[i__ - 1]) {
	    goto L280;
	}
	if (abs(d__) == zero) {
	    goto L420;
	}
	epcom6_2.pw[i__ - 1] = short__ * r0 / d__;
	epcom4_2.save1[i__ - 1] = epcom6_2.pw[i__ - 1] * rl1 * r0;
L280:
	;
    }
    goto L370;
L290:
    switch (miter1) {
	case 1:  goto L295;
	case 2:  goto L350;
	case 3:  goto L350;
	case 4:  goto L310;
    }
/* ----------------------------------------------------------------------- */
/* in the case of functional iteration, y is updated directly from */
/* the result of the last diffun call. */
/* ----------------------------------------------------------------------- */
L295:
    d__ = zero;
    i__3 = epcom1_2.n;
    for (i__ = 1; i__ <= i__3; ++i__) {
	r__ = rl1 * (epcom1_2.h__ * epcom5_2.save2[i__ - 1] - y[i__ + (y_dim1 
		<< 1)]);
/* Computing 2nd power */
	d__1 = (r__ - epcom3_2.error[i__ - 1]) / epcom2_2.ymax[i__ - 1];
	d__ += d__1 * d__1;
	epcom4_2.save1[i__ - 1] = y[i__ + y_dim1] + r__;
/* L300: */
	epcom3_2.error[i__ - 1] = r__;
    }
    goto L400;
/* ----------------------------------------------------------------------- */
/* in the case of a chord method, the residual -g(y sub n(m)) */
/* is computed and the linear system with that as right-hand side */
/* and p as coefficient matrix is solved, using the lu decomposition */
/* of p if miter = 1 or 2.  if miter = 3 the scalar h*rl1 is updated. */
/* ----------------------------------------------------------------------- */
L310:
    phrl1 = hrl1;
    hrl1 = epcom1_2.h__ * rl1;
    if (hrl1 == phrl1) {
	goto L330;
    }
    r__ = hrl1 / phrl1;
    i__3 = epcom1_2.n;
    for (i__ = 1; i__ <= i__3; ++i__) {
	d__ = one - r__ * (one - one / epcom6_2.pw[i__ - 1]);
	if (abs(d__) == zero) {
	    goto L440;
	}
/* L320: */
	epcom6_2.pw[i__ - 1] = one / d__;
    }
L330:
    i__3 = epcom1_2.n;
    for (i__ = 1; i__ <= i__3; ++i__) {
/* L340: */
	epcom4_2.save1[i__ - 1] = epcom6_2.pw[i__ - 1] * (rl1 * epcom1_2.h__ *
		 epcom5_2.save2[i__ - 1] - (rl1 * y[i__ + (y_dim1 << 1)] + 
		epcom3_2.error[i__ - 1]));
    }
    goto L370;
L350:
    i__3 = epcom1_2.n;
    for (i__ = 1; i__ <= i__3; ++i__) {
/* L360: */
	epcom4_2.save1[i__ - 1] = rl1 * epcom1_2.h__ * epcom5_2.save2[i__ - 1]
		 - (rl1 * y[i__ + (y_dim1 << 1)] + epcom3_2.error[i__ - 1]);
    }
    sol_(&epcom1_2.n, n0, epcom6_2.pw, epcom4_2.save1, epcom7_2.ipiv);
L370:
    d__ = zero;
    i__3 = epcom1_2.n;
    for (i__ = 1; i__ <= i__3; ++i__) {
	epcom3_2.error[i__ - 1] += epcom4_2.save1[i__ - 1];
/* Computing 2nd power */
	d__1 = epcom4_2.save1[i__ - 1] / epcom2_2.ymax[i__ - 1];
	d__ += d__1 * d__1;
/* L380: */
	epcom4_2.save1[i__ - 1] = y[i__ + y_dim1] + epcom3_2.error[i__ - 1];
    }
/* ----------------------------------------------------------------------- */
/* test for convergence.  if m .gt. 0, an estimate of the square of */
/* the convergence rate constant is stored in crate, and this is used */
/* in the test. */
/* ----------------------------------------------------------------------- */
L400:
    if (m != 0) {
/* Computing MAX */
	d__1 = crdown * crate, d__2 = d__ / d1;
	crate = max(d__1,d__2);
    }
    if (d__ * min(one,crate) <= bnd) {
	goto L450;
    }
    d1 = d__;
    ++m;
    if (m == maxcor) {
	goto L410;
    }
    (*diffun)(&epcom1_2.n, &epcom1_2.t, epcom4_2.save1, epcom5_2.save2);
    switch (miter1) {
	case 1:  goto L295;
	case 2:  goto L350;
	case 3:  goto L350;
	case 4:  goto L310;
    }
/* ----------------------------------------------------------------------- */
/* the corrector iteration failed to converge in 3 tries. if partial */
/* derivatives are involved but are not up to date, they are */
/* reevaluated for the next try.  otherwise the y array is restored */
/* to its values before prediction, and h is reduced, */
/* if possible.  if not, a no-convergence exit is taken. */
/* ----------------------------------------------------------------------- */
L410:
    epcom9_1.nfe = epcom9_1.nfe + maxcor - 1;
    if (newj == -1) {
	goto L440;
    }
L420:
    epcom1_2.t = told;
    etamax = one;
    i__3 = epcm10_1.nq;
    for (j1 = 1; j1 <= i__3; ++j1) {
	i__1 = epcm10_1.nq;
	for (j2 = j1; j2 <= i__1; ++j2) {
	    j = epcm10_1.nq + j1 - j2;
	    i__2 = epcom1_2.n;
	    for (i__ = 1; i__ <= i__2; ++i__) {
/* L430: */
		y[i__ + j * y_dim1] -= y[i__ + (j + 1) * y_dim1];
	    }
	}
    }
    if (abs(epcom1_2.h__) <= epcom1_2.hmin * onepsm) {
	goto L680;
    }
    eta = etacf;
    iredo = 1;
    goto L180;
L440:
    newj = miter;
    goto L220;
/* ----------------------------------------------------------------------- */
/* the corrector has converged.  newj is set to -1 if partial */
/* derivatives were used, to signal that they may need updating on */
/* subsequent steps.  the error test is made and control passes to */
/* statement 500 if it fails. */
/* ----------------------------------------------------------------------- */
L450:
    if (miter != 0) {
	newj = -1;
    }
    epcom9_1.nfe += m;
    d__ = zero;
    i__2 = epcom1_2.n;
    for (i__ = 1; i__ <= i__2; ++i__) {
/* L460: */
/* Computing 2nd power */
	d__1 = epcom3_2.error[i__ - 1] / epcom2_2.ymax[i__ - 1];
	d__ += d__1 * d__1;
    }
/* Computing 2nd power */
    d__1 = epcm10_1.tq[1] * epcom1_2.eps;
    e = flotn * (d__1 * d__1);
    if (d__ > e) {
	goto L500;
    }
/* ----------------------------------------------------------------------- */
/* after a successful step, the y array, tau, nstep, and nqindx are */
/* updated, and a new value of h at order nq is computed. */
/* the vector tau contains the nq + 1 most recent values of h. */
/* a change in nq up or down by 1 is considered if nqindx = 0. */
/* if nqindx = 1 and nq .lt. maxder, then error is saved */
/* for use in a possible order increase on the next step. */
/* a change in h or nq is made only of the increase in h */
/* is by a factor of at least 1.3. */
/* if not, nqindx is set to 2 to prevent testing for that many */
/* steps.  if nq is changed, nqindx is set to nq + 1 (new value). */
/* ----------------------------------------------------------------------- */
    epcom1_2.kflag = 0;
    iredo = 0;
    ++epcom9_1.nstep;
    epcom9_1.hused = epcom1_2.h__;
    epcom9_1.nqused = epcm10_1.nq;
    i__2 = epcm10_1.nq;
    for (iback = 1; iback <= i__2; ++iback) {
	i__ = epcm10_1.l - iback;
/* L470: */
	epcm10_1.tau[i__] = epcm10_1.tau[i__ - 1];
    }
    epcm10_1.tau[0] = epcom1_2.h__;
    i__2 = epcm10_1.l;
    for (j = 1; j <= i__2; ++j) {
	i__1 = epcom1_2.n;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* L480: */
	    y[i__ + j * y_dim1] += epcom3_2.error[i__ - 1] * epcm10_1.el[j - 
		    1];
	}
    }
    --epcm10_1.nqindx;
    if (epcm10_1.l == epcm10_1.lmax || epcm10_1.nqindx != 1) {
	goto L495;
    }
    i__1 = epcom1_2.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L490: */
	y[i__ + epcm10_1.lmax * y_dim1] = epcom3_2.error[i__ - 1];
    }
    conp = epcm10_1.tq[4];
L495:
    if (etamax != one) {
	goto L520;
    }
    if (epcm10_1.nqindx < 2) {
	epcm10_1.nqindx = 2;
    }
    goto L690;
/* ----------------------------------------------------------------------- */
/* the error test failed.  kflag keeps track of multiple failures. */
/* t and the y array are restored to their previous values.  a new */
/* h for a retry of the step is computed.  the order is kept fixed. */
/* ----------------------------------------------------------------------- */
L500:
    --epcom1_2.kflag;
    epcom1_2.t = told;
    i__1 = epcm10_1.nq;
    for (j1 = 1; j1 <= i__1; ++j1) {
	i__2 = epcm10_1.nq;
	for (j2 = j1; j2 <= i__2; ++j2) {
	    j = epcm10_1.nq + j1 - j2;
	    i__3 = epcom1_2.n;
	    for (i__ = 1; i__ <= i__3; ++i__) {
/* L510: */
		y[i__ + j * y_dim1] -= y[i__ + (j + 1) * y_dim1];
	    }
	}
    }
    newj = miter;
    etamax = one;
    if (abs(epcom1_2.h__) <= epcom1_2.hmin * onepsm) {
	goto L660;
    }
    if (epcom1_2.kflag <= kfc) {
	goto L630;
    }
    iredo = 2;
/* compute ratio of new h to current h at the current order. ------------ */
L520:
    flotl = (doublereal) epcm10_1.l;
    d__1 = bias2 * d__ / e;
    d__2 = pt5 / flotl;
    etaq = one / (pow_dd(&d__1, &d__2) + addon);
    if (epcm10_1.nqindx != 0 || epcom1_2.kflag != 0) {
	goto L580;
    }
    etaqm1 = zero;
    if (epcm10_1.nq == 1) {
	goto L540;
    }
/* compute ratio of new h to current h at the current order less one. --- */
    d__ = zero;
    i__3 = epcom1_2.n;
    for (i__ = 1; i__ <= i__3; ++i__) {
/* L530: */
/* Computing 2nd power */
	d__1 = y[i__ + epcm10_1.l * y_dim1] / epcom2_2.ymax[i__ - 1];
	d__ += d__1 * d__1;
    }
/* Computing 2nd power */
    d__1 = epcm10_1.tq[0] * epcom1_2.eps;
    edn = flotn * (d__1 * d__1);
    d__1 = bias1 * d__ / edn;
    d__2 = pt5 / (flotl - one);
    etaqm1 = one / (pow_dd(&d__1, &d__2) + addon);
L540:
    etaqp1 = zero;
    if (epcm10_1.l == epcm10_1.lmax) {
	goto L560;
    }
/* compute ratio of new h to current h at current order plus one. ------- */
    d__1 = epcom1_2.h__ / epcm10_1.tau[1];
    cnquot = epcm10_1.tq[4] / conp * pow_di(&d__1, &epcm10_1.l);
    d__ = zero;
    i__3 = epcom1_2.n;
    for (i__ = 1; i__ <= i__3; ++i__) {
/* L550: */
/* Computing 2nd power */
	d__1 = (epcom3_2.error[i__ - 1] - cnquot * y[i__ + epcm10_1.lmax * 
		y_dim1]) / epcom2_2.ymax[i__ - 1];
	d__ += d__1 * d__1;
    }
/* Computing 2nd power */
    d__1 = epcm10_1.tq[2] * epcom1_2.eps;
    eup = flotn * (d__1 * d__1);
    d__1 = bias3 * d__ / eup;
    d__2 = pt5 / (flotl + one);
    etaqp1 = one / (pow_dd(&d__1, &d__2) + addon);
L560:
    epcm10_1.nqindx = 2;
    if (etaq >= etaqp1) {
	goto L570;
    }
    if (etaqp1 > etaqm1) {
	goto L600;
    }
    goto L590;
L570:
    if (etaq < etaqm1) {
	goto L590;
    }
L580:
    if (etaq < thresh && epcom1_2.kflag == 0) {
	goto L690;
    }
    eta = etaq;
    if (epcom1_2.kflag <= -2 && eta > etamxf) {
	eta = etamxf;
    }
    goto L180;
L590:
    if (etaqm1 < thresh) {
	goto L690;
    }
    adjust_(&y[y_offset], n0);
    epcm10_1.l = epcm10_1.nq;
    --epcm10_1.nq;
    eta = etaqm1;
    epcm10_1.nqindx = epcm10_1.l;
    goto L180;
L600:
    if (etaqp1 < thresh) {
	goto L690;
    }
    epcm10_1.nq = epcm10_1.l;
    eta = etaqp1;
    ++epcm10_1.l;
    i__3 = epcom1_2.n;
    for (i__ = 1; i__ <= i__3; ++i__) {
/* L610: */
	y[i__ + epcm10_1.l * y_dim1] = zero;
    }
    epcm10_1.nqindx = epcm10_1.l;
    goto L180;
/* ----------------------------------------------------------------------- */
/* control reaches this section if 3 or more consecutive failures */
/* have occurred.  it is assumed that the elements of the y array */
/* have accumulated errors of the wrong order.  the order is reduced */
/* by one, if possible.  then h is reduced by a factor of 0.1 and */
/* the step is retried.  after a total of 7 consecutive failures, */
/* an exit is taken with kflag = -2. */
/* ----------------------------------------------------------------------- */
L630:
    if (epcom1_2.kflag == kfh) {
	goto L670;
    }
    if (epcm10_1.nq == 1) {
	goto L640;
    }
    eta = etamin;
    adjust_(&y[y_offset], n0);
    epcm10_1.l = epcm10_1.nq;
    --epcm10_1.nq;
    epcm10_1.nqindx = epcm10_1.l;
    goto L180;
L640:
/* Computing MAX */
    d__1 = etamin, d__2 = epcom1_2.hmin / abs(epcom1_2.h__);
    eta = max(d__1,d__2);
    epcom1_2.h__ *= eta;
    (*diffun)(&epcom1_2.n, &epcom1_2.t, &y[y_offset], epcom4_2.save1);
    ++epcom9_1.nfe;
    i__3 = epcom1_2.n;
    for (i__ = 1; i__ <= i__3; ++i__) {
/* L650: */
	y[i__ + (y_dim1 << 1)] = epcom1_2.h__ * epcom4_2.save1[i__ - 1];
    }
    epcm10_1.nqindx = 10;
    goto L200;
/* ----------------------------------------------------------------------- */
/* all returns are made through this section.  h is saved in hold */
/* to allow the caller to change h on the next step. */
/* ----------------------------------------------------------------------- */
L660:
    epcom1_2.kflag = -1;
    goto L700;
L670:
    epcom1_2.kflag = -2;
    goto L700;
L680:
    epcom1_2.kflag = -3;
    goto L700;
L690:
    etamax = etamx3;
    if (epcom9_1.nstep <= 10) {
	etamax = etamx2;
    }
L700:
    hold = epcom1_2.h__;
    epcom1_2.jstart = epcm10_1.nq;
    return 0;
/* -----------------------  end of subroutine tstep  --------------------- */
} /* tstep_ */

/* Subroutine */ int coset_(void)
{
    /* Initialized data */

    static doublereal cortes = .1;
    static doublereal one = 1.;
    static doublereal six = 6.;
    static doublereal two = 2.;
    static doublereal zero = 0.;

    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer i__, j;
    doublereal s, em[13], xi, em0;
    integer jp1;
    doublereal elp, rxi;
    integer nqm1;
    doublereal prod, csum, hsum, cnqm1, hsum1;
    integer iback;
    doublereal ahdss, floti, flotl;
    integer maxder;
    doublereal flotnq;

/* ----------------------------------------------------------------------- */
/* coset is called by tstep and sets coefficients for use there. */

/* for each order nq, the coefficients in el are calculated by use of */
/*  the generating polynomial lambda(x), with coefficients el(i): */
/*      lambda(x) = el(1) + el(2)*x + ... + el(nq+1)*(x**nq). */
/* for the backward differentiation formulas, */
/*                    nq */
/*      lambda(x) = product (1 + x/xi(i) ) . */
/*                   i = 1 */
/* for the adams formulas, */
/*                              nq-1 */
/*      (d/dx) lambda(x) = c * product (1 + x/xi(i) ) , */
/*                              i = 1 */
/*      lambda(-1) = 0,    lambda(0) = 1, */
/* where c is a normalization constant. */
/* in both cases, xi(i) is defined by */
/*      h*xi(i) = t sub n  -  t sub (n-i) */
/*              = h + tau(1) + tau(2) + ... tau(i-1). */

/* coset also sets maxder, the maximum order of the formulas */
/* available. currently this is 5 for the backward differentiation */
/* formulas, and 12 for the adams formulas.  to use different */
/* values (.le. 13),  change the numbers in statements 1 and 2 below. */

/* in addition to variables described previously, communication */
/* with coset uses the following.. */
/*   tau    = a vector of length 13 containing the past nq values */
/*            of h. */
/*   el     = a vector of length 13 in which coset stores the */
/*            coefficients for the corrector formula. */
/*   tq     = a vector of length 5 in which coset stores constants */
/*            used for the convergence test, the error test, and */
/*            selection of h at a new order. */
/*   lmax   = maxder + 1, where maxder is the maximum order */
/*            available.  lmax is the maximum number of columns */
/*            of the y array to be used. */
/*   meth   = the basic method indicator. */
/*   nq     = the current order. */
/*   l      = nq + 1, the length of the vector stored in el, and */
/*            the number of columns of the y array being used. */
/*   nqindx = a counter controlling the frequency of order changes. */
/*            an order change is about to be considered if */
/*            nqindx = 1. */
/* ----------------------------------------------------------------------- */
/* caution:  not all members of epcom1 are used in this subroutine. */
/* ----------------------------------------------------------------------- */

    ahdss = one;
    if (epcom1_2.ss != zero) {
	ahdss = abs(epcom1_2.h__) / epcom1_2.ss;
    }
    flotl = (doublereal) epcm10_1.l;
    nqm1 = epcm10_1.nq - 1;
    switch (epcm10_1.meth) {
	case 1:  goto L1;
	case 2:  goto L2;
    }
L1:
    maxder = 12;
    goto L100;

L2:
    maxder = 5;
    goto L200;

L100:
    if (epcm10_1.nq != 1) {
	goto L110;
    }
    epcm10_1.el[0] = one;
    epcm10_1.el[1] = one;
    epcm10_1.tq[0] = one;
    epcm10_1.tq[1] = two * ahdss;
    epcm10_1.tq[2] = six * epcm10_1.tq[1];
    epcm10_1.tq[4] = one;
    goto L300;
L110:
    hsum = epcom1_2.h__;
    em[0] = one;
    flotnq = flotl - one;
    i__1 = epcm10_1.l;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* L115: */
	em[i__ - 1] = zero;
    }
    i__1 = nqm1;
    for (j = 1; j <= i__1; ++j) {
	if (j != nqm1 || epcm10_1.nqindx != 1) {
	    goto L130;
	}
	s = one;
	csum = zero;
	i__2 = nqm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    csum += s * em[i__ - 1] / (doublereal) (i__ + 1);
/* L120: */
	    s = -s;
	}
	epcm10_1.tq[0] = ahdss * em[nqm1 - 1] / (flotnq * csum);
L130:
	rxi = epcom1_2.h__ / hsum;
	i__2 = j;
	for (iback = 1; iback <= i__2; ++iback) {
	    i__ = j + 2 - iback;
/* L140: */
	    em[i__ - 1] += em[i__ - 2] * rxi;
	}
/* L150: */
	hsum += epcm10_1.tau[j - 1];
    }
/* compute integral from -1 to 0 of polynomial and of x times it. ------- */
    s = one;
    em0 = zero;
    csum = zero;
    i__1 = epcm10_1.nq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	floti = (doublereal) i__;
	em0 += s * em[i__ - 1] / floti;
	csum += s * em[i__ - 1] / (floti + 1);
/* L160: */
	s = -s;
    }
/* in el, form coefficients of normalized integrated polynomial. -------- */
    s = one / em0;
    epcm10_1.el[0] = one;
    i__1 = epcm10_1.nq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L170: */
	epcm10_1.el[i__] = s * em[i__ - 1] / (doublereal) i__;
    }
    xi = hsum / epcom1_2.h__;
    epcm10_1.tq[1] = ahdss * xi * em0 / csum;
    epcm10_1.tq[4] = xi / epcm10_1.el[epcm10_1.l - 1];
    if (epcm10_1.nqindx != 1) {
	goto L300;
    }
/* for higher order control constant, multiply polynomial by 1+x/xi(q). - */
    rxi = one / xi;
    i__1 = epcm10_1.nq;
    for (iback = 1; iback <= i__1; ++iback) {
	i__ = epcm10_1.l + 1 - iback;
/* L180: */
	em[i__ - 1] += em[i__ - 2] * rxi;
    }
/* compute integral of polynomial. -------------------------------------- */
    s = one;
    csum = zero;
    i__1 = epcm10_1.l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	csum += s * em[i__ - 1] / (doublereal) (i__ + 1);
/* L190: */
	s = -s;
    }
    epcm10_1.tq[2] = ahdss * flotl * em0 / csum;
    goto L300;

L200:
    i__1 = epcm10_1.l;
    for (i__ = 3; i__ <= i__1; ++i__) {
/* L210: */
	epcm10_1.el[i__ - 1] = zero;
    }
    epcm10_1.el[0] = one;
    epcm10_1.el[1] = one;
    hsum = epcom1_2.h__;
    hsum1 = zero;
    prod = one;
    rxi = one;
    if (epcm10_1.nq == 1) {
	goto L240;
    }
    i__1 = nqm1;
    for (j = 1; j <= i__1; ++j) {
/* in el, construct coefficients of (1+x/xi(1))*...*(1+x/xi(j+1)). ------ */
	hsum += epcm10_1.tau[j - 1];
	hsum1 += epcm10_1.tau[j - 1];
	prod *= hsum / hsum1;
	rxi = epcom1_2.h__ / hsum;
	jp1 = j + 1;
	i__2 = jp1;
	for (iback = 1; iback <= i__2; ++iback) {
	    i__ = j + 3 - iback;
/* L220: */
	    epcm10_1.el[i__ - 1] += epcm10_1.el[i__ - 2] * rxi;
	}
/* L230: */
    }
L240:
    epcm10_1.tq[1] = ahdss * epcm10_1.el[1] * (one + prod);
    epcm10_1.tq[4] = (one + prod) / epcm10_1.el[epcm10_1.l - 1];
    if (epcm10_1.nqindx != 1) {
	goto L300;
    }
    cnqm1 = rxi / epcm10_1.el[epcm10_1.l - 1];
    elp = epcm10_1.el[1] - rxi;
    epcm10_1.tq[0] = ahdss * elp / cnqm1;
    hsum += epcm10_1.tau[epcm10_1.nq - 1];
    rxi = epcom1_2.h__ / hsum;
    elp = epcm10_1.el[1] + rxi;
    epcm10_1.tq[2] = ahdss * elp * rxi * (one + prod) * (flotl + one);
L300:
    epcm10_1.tq[3] = cortes * epcm10_1.tq[1];
    epcm10_1.lmax = maxder + 1;
    return 0;
/* -----------------------  end of subroutine coset  --------------------- */
} /* coset_ */

/* Subroutine */ int adjust_(doublereal *y, integer *n0)
{
    /* Initialized data */

    static doublereal one = 1.;
    static doublereal zero = 0.;

    /* System generated locals */
    integer y_dim1, y_offset, i__1, i__2;

    /* Local variables */
    integer i__, j;
    doublereal xi;
    integer jp1, nqm1, nqm2;
    doublereal hsum;
    integer iback;

/* ----------------------------------------------------------------------- */
/* this subroutine adjusts the y array on reduction of order. */
/* see reference 1 for details. */
/* ----------------------------------------------------------------------- */
/* caution:  not all members of epcom1 are used in this subroutine. */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    y_dim1 = *n0;
    y_offset = 1 + y_dim1;
    y -= y_offset;

    /* Function Body */
    if (epcm10_1.nq == 2) {
	return 0;
    }
    nqm1 = epcm10_1.nq - 1;
    nqm2 = epcm10_1.nq - 2;
    switch (epcm10_1.meth) {
	case 1:  goto L100;
	case 2:  goto L200;
    }

L100:
    i__1 = epcm10_1.lmax;
    for (j = 1; j <= i__1; ++j) {
/* L110: */
	epcm10_1.el[j - 1] = zero;
    }
    epcm10_1.el[1] = one;
    hsum = zero;
    i__1 = nqm2;
    for (j = 1; j <= i__1; ++j) {
/* construct coefficients of x*(x+xi(1))*...*(x+xi(j)). ----------------- */
	hsum += epcm10_1.tau[j - 1];
	xi = hsum / epcom1_2.h__;
	jp1 = j + 1;
	i__2 = jp1;
	for (iback = 1; iback <= i__2; ++iback) {
	    i__ = j + 3 - iback;
/* L120: */
	    epcm10_1.el[i__ - 1] = epcm10_1.el[i__ - 1] * xi + epcm10_1.el[
		    i__ - 2];
	}
/* L130: */
    }
/* construct coefficients of integrated polynomial. --------------------- */
    i__1 = nqm1;
    for (j = 2; j <= i__1; ++j) {
/* L140: */
	epcm10_1.el[j] = (doublereal) epcm10_1.nq * epcm10_1.el[j - 1] / (
		doublereal) j;
    }
    goto L300;

L200:
    i__1 = epcm10_1.lmax;
    for (j = 1; j <= i__1; ++j) {
/* L210: */
	epcm10_1.el[j - 1] = zero;
    }
    epcm10_1.el[2] = one;
    hsum = zero;
    i__1 = nqm2;
    for (j = 1; j <= i__1; ++j) {
/* construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). --------------- */
	hsum += epcm10_1.tau[j - 1];
	xi = hsum / epcom1_2.h__;
	jp1 = j + 1;
	i__2 = jp1;
	for (iback = 1; iback <= i__2; ++iback) {
	    i__ = j + 4 - iback;
/* L220: */
	    epcm10_1.el[i__ - 1] = epcm10_1.el[i__ - 1] * xi + epcm10_1.el[
		    i__ - 2];
	}
/* L230: */
    }

/* subtract correction terms from y array. ------------------------------ */
L300:
    i__1 = epcm10_1.nq;
    for (j = 3; j <= i__1; ++j) {
	i__2 = epcom1_2.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L310: */
	    y[i__ + j * y_dim1] -= y[i__ + epcm10_1.l * y_dim1] * epcm10_1.el[
		    j - 1];
	}
/* L320: */
    }
    return 0;
/* -----------------------  end of subroutine adjust  -------------------- */
} /* adjust_ */

/* Subroutine */ int pset_(S_fp diffun, S_fp pederv, doublereal *y, integer *
	n0, doublereal *con, integer *miter, integer *ier)
{
    /* Initialized data */

    static doublereal one = 1.;
    static doublereal rep = .001;
    static doublereal zero = 0.;

    /* System generated locals */
    integer y_dim1, y_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    doublereal d__;
    integer i__, j;
    doublereal r__;
    integer j1;
    doublereal r0, yj;
    extern /* Subroutine */ int dec_(integer *, integer *, doublereal *, 
	    integer *, integer *);

/* ----------------------------------------------------------------------- */
/* pset is called by tstep to compute and to process the matrix */
/* p = i - (h/el(2))*j, where j is an approximation to the */
/* jacobian.  j is computed by either the user supplied */
/* subroutine pederv, when miter = 1, or by finite differences, */
/* when miter = 2.  j is stored in pw and replaced by p, using */
/* con = -h/el(2).  then p is subjected to an lu decomposition */
/* for later solution of linear algebraic systems with p as the */
/* coefficient matrix. */

/* in addition to variables described previously, communication */
/* with pset uses the following.. */
/*   epsj = sqrt(uround), used in the numerical jacobian increments. */
/*   nsq  = n0**2. */
/* ----------------------------------------------------------------------- */
/* caution:  not all epcom1 variables are used inthis subroutine. */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    y_dim1 = *n0;
    y_offset = 1 + y_dim1;
    y -= y_offset;

    /* Function Body */
    if (*miter == 2) {
	goto L20;
    }
/* if miter = 1, call pederv and multiply by a scalar.  ----------------- */
    (*pederv)(&epcom1_2.n, &epcom1_2.t, &y[y_offset], epcom6_2.pw, n0);
    i__1 = epcom8_1.nsq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	epcom6_2.pw[i__ - 1] *= *con;
    }
    goto L60;
/* if miter = 2, make n calls to diffun to approximate j. --------------- */
L20:
    d__ = zero;
    i__1 = epcom1_2.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L30: */
/* Computing 2nd power */
	d__1 = epcom5_2.save2[i__ - 1];
	d__ += d__1 * d__1;
    }
    r0 = abs(epcom1_2.h__) * sqrt(d__) * epcom1_2.uround / rep;
    j1 = 0;
    i__1 = epcom1_2.n;
    for (j = 1; j <= i__1; ++j) {
	yj = y[j + y_dim1];
	r__ = epcom8_1.epsj * epcom2_2.ymax[j - 1];
	r__ = max(r__,r0);
	y[j + y_dim1] += r__;
	d__ = *con / r__;
	(*diffun)(&epcom1_2.n, &epcom1_2.t, &y[y_offset], epcom4_2.save1);
	i__2 = epcom1_2.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L40: */
	    epcom6_2.pw[i__ + j1 - 1] = (epcom4_2.save1[i__ - 1] - 
		    epcom5_2.save2[i__ - 1]) * d__;
	}
	y[j + y_dim1] = yj;
	j1 += *n0;
/* L50: */
    }
/* add on the identity matrix.  ----------------------------------------- */
L60:
    j = 1;
    i__1 = epcom1_2.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	epcom6_2.pw[j - 1] += one;
/* L70: */
	j += *n0 + 1;
    }
/* get lu decomposition of p. ------------------------------------------- */
    dec_(&epcom1_2.n, n0, epcom6_2.pw, epcom7_2.ipiv, ier);
    return 0;
/* -----------------------  end of subroutine pset ----------------------- */
} /* pset_ */

#ifndef DECSOL_DEFINED
/* dec_ and sol_ are also defined in decsol.c. When linking together,
   we use the versions from decsol.c and skip these. */

/* Subroutine */ int dec_(integer *n, integer *ndim, doublereal *a, integer *
	ip, integer *ier)
{
    /* Initialized data */

    static doublereal one = 1.;
    static doublereal zero = 0.;

    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Local variables */
    integer i__, j, k, m;
    doublereal t;
    integer kp1, nm1;

/* ----------------------------------------------------------------------- */
/* matrix triangularization by gaussian elimination. */
/* input.. */
/*    n = order of matrix. */
/*    ndim = declared dimension of array  a . */
/*    a = matrix to be triangularized. */
/* output.. */
/*    a(i,j), i.le.j = upper triangular factor, u . */
/*    a(i,j), i.gt.j = multipliers = lower triangular factor, i - l. */
/*    ip(k), k.lt.n = index of k-th pivot row. */
/*    ip(n) = (-1)**(number of interchanges) or o . */
/*    ier = 0 if a nonsingular, or k if a found to be */
/*          singular at stage k. */
/* use  sol  to obtain solution of linear system. */
/* determ(a) = ip(n)*a(1,1)*a(2,2)*...*a(n,n). */
/* if ip(n)=0, a is singular, sol will divide by zero. */
/* interchanges finished in u , only partly in l . */

/* reference.. */
/* c. b. moler, algorithm 423, linear equation solver, */
/* comm. assoc. comput. mach., 15 (1972), p. 274. */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --ip;
    a_dim1 = *ndim;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    *ier = 0;
    ip[*n] = 1;
    if (*n == 1) {
	goto L70;
    }
    nm1 = *n - 1;
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;
	m = k;
	i__2 = *n;
	for (i__ = kp1; i__ <= i__2; ++i__) {
/* L10: */
	    if ((d__1 = a[i__ + k * a_dim1], abs(d__1)) > (d__2 = a[m + k * 
		    a_dim1], abs(d__2))) {
		m = i__;
	    }
	}
	ip[k] = m;
	t = a[m + k * a_dim1];
	if (m == k) {
	    goto L20;
	}
	ip[*n] = -ip[*n];
	a[m + k * a_dim1] = a[k + k * a_dim1];
	a[k + k * a_dim1] = t;
L20:
	if (t == zero) {
	    goto L80;
	}
	t = one / t;
	i__2 = *n;
	for (i__ = kp1; i__ <= i__2; ++i__) {
/* L30: */
	    a[i__ + k * a_dim1] = -a[i__ + k * a_dim1] * t;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    t = a[m + j * a_dim1];
	    a[m + j * a_dim1] = a[k + j * a_dim1];
	    a[k + j * a_dim1] = t;
	    if (t == zero) {
		goto L50;
	    }
	    i__3 = *n;
	    for (i__ = kp1; i__ <= i__3; ++i__) {
/* L40: */
		a[i__ + j * a_dim1] += a[i__ + k * a_dim1] * t;
	    }
L50:
	    ;
	}
/* L60: */
    }
L70:
    k = *n;
    if (a[*n + *n * a_dim1] == zero) {
	goto L80;
    }
    return 0;
L80:
    *ier = k;
    ip[*n] = 0;
    return 0;
/* -----------------------  end of subroutine dec  ----------------------- */
} /* dec_ */

/* Subroutine */ int sol_(integer *n, integer *ndim, doublereal *a, 
	doublereal *b, integer *ip)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    integer i__, k, m;
    doublereal t;
    integer kb, km1, kp1, nm1;

/* ----------------------------------------------------------------------- */
/* solution of linear system, a*x = b . */
/* input.. */
/*   n = order of matrix. */
/*   ndim = declared dimension of array  a . */
/*   a = triangularized matrix obtained from dec. */
/*   b = right hand side vector. */
/*   ip = pivot vector obtained from dec. */
/* do not use if dec has set ier .ne. 0. */
/* output.. */
/*   b = solution vector, x . */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --ip;
    --b;
    a_dim1 = *ndim;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    if (*n == 1) {
	goto L50;
    }
    nm1 = *n - 1;
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;
	m = ip[k];
	t = b[m];
	b[m] = b[k];
	b[k] = t;
	i__2 = *n;
	for (i__ = kp1; i__ <= i__2; ++i__) {
/* L10: */
	    b[i__] += a[i__ + k * a_dim1] * t;
	}
/* L20: */
    }
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
	km1 = *n - kb;
	k = km1 + 1;
	b[k] /= a[k + k * a_dim1];
	t = -b[k];
	i__2 = km1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L30: */
	    b[i__] += a[i__ + k * a_dim1] * t;
	}
/* L40: */
    }
L50:
    b[1] /= a[a_dim1 + 1];
    return 0;
/* -----------------------  end of subroutine sol  ----------------------- */
} /* sol_ */

#endif /* DECSOL_DEFINED */

