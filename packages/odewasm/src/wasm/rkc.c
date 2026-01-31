/* ../reference/netlib/rkc.f -- translated by f2c (version 20240504).
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

struct {
    integer nfe, nsteps, naccpt, nrejct, nfesig, maxm;
} rkcdid_;

#define rkcdid_1 rkcdid_

/* Table of constant values */

static integer c__9 = 9;
static integer c__1 = 1;
static doublereal c_b9 = 1.;
static doublereal c_b14 = .33333333333333331;
static doublereal c_b18 = .66666666666666663;

/* Subroutine */ int rkc_(integer *neqn, U_fp f, doublereal *y, doublereal *t,
	 doublereal *tend, doublereal *rtol, doublereal *atol, integer *info, 
	doublereal *work, integer *idid)
{
    /* System generated locals */
    integer i__1;

    /* Builtin functions */
    double sqrt(doublereal);
    integer s_wsle(cilist *), do_lio(integer *, integer *, char *, ftnlen), 
	    e_wsle(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    static integer i__, ptr1, ptr2, ptr3, ptr4;
    static logical valid, array;
    extern /* Subroutine */ int rkclow_(integer *, doublereal *, doublereal *,
	     doublereal *, U_fp, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, integer *);

    /* Fortran I/O blocks */
    static cilist io___8 = { 0, 6, 0, 0, 0 };


/* -------------------------------------------------------------------------- */

/*  ABSTRACT:  RKC integrates initial value problems for systems of first */
/*  order ordinary differential equations.  It is based on a family of */
/*  explicit Runge-Kutta-Chebyshev formulas of order two.  The stability */
/*  of members of the family increases quadratically in the number of */
/*  stages m. An estimate of the spectral radius is used at each step to */
/*  select the smallest m resulting in a stable integration. RKC is */
/*  appropriate for the solution to modest accuracy of mildly stiff problems */
/*  with eigenvalues of Jacobians that are close to the negative real axis. */
/*  For such problems it has the advantages of explicit one-step methods and */
/*  very low storage. If it should turn out that RKC is using m far beyond */
/*  100, the problem is not mildly stiff and alternative methods should be */
/*  considered.  Answers can be obtained cheaply anywhere in the interval */
/*  of integration by means of a continuous extension evaluated in the */
/*  subroutine RKCINT. */

/*  The initial value problems arising from semi-discretization of */
/*  diffusion-dominated parabolic partial differential equations and of */
/*  reaction-diffusion equations, especially in two and three spatial */
/*  variables, exemplify the problems for which RKC was designed.  Two */
/*  example programs, ExA and ExB, are provided that show how to use RKC. */

/* --------------------------------------------------------------------------- */
/*  USAGE:  RKC integrates a system of NEQN first order ordinary differential */
/*  equations specified by a subroutine F from T to TEND.  The initial values */
/*  at T are input in Y(*).  On all returns from RKC, Y(*) is an approximate */
/*  solution at T.  In the computation of Y(*), the local error has been */
/*  controlled at each step to satisfy a relative error tolerance RTOL and */
/*  absolute error tolerances ATOL(*).  The array INFO(*) specifies the way */
/*  the problem is to be solved.  WORK(*) is a work array. IDID reports */
/*  success or the reason the computation has been terminated. */

/*  FIRST CALL TO RKC */

/*  You must provide storage in your calling program for the arrays in the */
/*  call list -- Y(NEQN), INFO(4), WORK(8+5*NEQN).  If INFO(2) = 1, you can */
/*      [Note: originally this said INFO(2) = 0, but that is reportedly */
/*       an error. See readme. Edited 2025-01-16.] */
/*  reduce the storage for the work array to WORK(8+4*NEQN).  ATOL may be */
/*  a scalar or an array.  If it is an array, you must provide storage for */
/*  ATOL(NEQN).  You must declare F in an external statement, supply the */
/*  subroutine F and the function SPCRAD, and initialize the following */
/*  quantities: */

/*    NEQN:  The number of differential equations. Integer. */

/*    T:     The initial point of the integration. Double precision. */
/*           Must be a variable. */

/*    TEND:  The end of the interval of integration.  Double precision. */
/*           TEND may be less than T. */

/*    Y(*):  The initial value of the solution.  Double precision array */
/*           of length NEQN. */

/*    F:     The name of a subroutine for evaluating the differential */
/*           equation.  It must have the form */

/*             subroutine f(neqn,t,y,dy) */
/*             integer          neqn */
/*             double precision t,y(neqn),dy(neqn) */
/*             dy(1)    = ... */
/*             ... */
/*             dy(neqn) = ... */
/*             return */
/*             end */

/*  RTOL, */
/*  ATOL(*):  At each step of the integration the local error is controlled */
/*            so that its RMS norm is no larger than tolerances RTOL, ATOL(*). */
/*            RTOL is a double precision scalar. ATOL(*) is either a double */
/*            precision scalar or a double precision array of length NEQN. */
/*            RKC is designed for the solution of problems to modest accuracy. */
/*            Because it is based on a method of order 2, it is relatively */
/*            expensive to achieve high accuracy. */

/*            RTOL is a relative error tolerance.  You must ask for some */
/*            relative accuracy, but you cannot ask for too much for the */
/*            precision available.  Accordingly, it is required that */
/*            0.1 >= RTOL >= 10*uround. (See below for the machine and */
/*            precision dependent quantity uround.) */

/*            ATOL is an absolute error tolerance that can be either a */
/*            scalar or an array.  When it is an array, the tolerances are */
/*            applied to corresponding components of the solution and when */
/*            it is a scalar, it is applied to all components.  A scalar */
/*            tolerance is reasonable only when all solution components are */
/*            scaled to be of comparable size.  A scalar tolerance saves a */
/*            useful amount of storage and is convenient.  Use INFO(*) to */
/*            tell RKC whether ATOL is a scalar or an array. */

/*            The absolute error tolerances ATOL(*) must satisfy ATOL(i) >= 0 */
/*            for i = 1,...,NEQN.  ATOL(j)= 0 specifies a pure relative error */
/*            test on component j of the solution, so it is an error if this */
/*            component vanishes in the course of the integration. */

/*            If all is going well, reducing the tolerances by a factor of */
/*            0.1 will reduce the error in the computed solution by a factor */
/*            of roughly 0.2. */

/*  INFO(*)   Integer array of length 4 that specifies how the problem */
/*            is to be solved. */

/*  INFO(1):  RKC integrates the initial value problem from T to TEND. */
/*            This is done by computing approximate solutions at points */
/*            chosen automatically throughout [T, TEND].  Ordinarily RKC */
/*            returns at each step with an approximate solution. These */
/*            approximations show how y behaves throughout the interval. */
/*            The subroutine RKCINT can be used to obtain answers anywhere */
/*            in the span of a step very inexpensively. This makes it */
/*            possible to obtain answers at specific points in [T, TEND] */
/*            and to obtain many answers very cheaply when attempting to */
/*            locating where some function of the solution has a zero */
/*            (event location).  Sometimes you will be interested only in */
/*            a solution at TEND, so you can suppress the returns at each */
/*            step along the way if you wish. */

/*  INFO(1)  = 0 Return after each step on the way to TEND with a */
/*               solution Y(*) at the output value of T. */

/*           = 1 Compute a solution Y(*) at TEND only. */

/*  INFO(2):  RKC needs an estimate of the spectral radius of the Jacobian. */
/*            You must provide a function that must be called SPCRAD and */
/*            have the form */

/*              double precision function spcrad(neqn,t,y) */
/*              integer          neqn */
/*              double precision t,y(neqn) */

/*              spcrad = < expression depending on info(2) > */

/*              return */
/*              end */

/*            You can provide a dummy function and let RKC compute the */
/*            estimate. Sometimes it is convenient for you to compute in */
/*            SPCRAD a reasonably close upper bound on the spectral radius, */
/*            using, e.g., Gershgorin's theorem.  This may be faster and/or */
/*            more reliable than having RKC compute one. */

/*  INFO(2)  = 0  RKC is to compute the estimate internally. */
/*                Assign any value to SPCRAD. */

/*           = 1  SPCRAD returns an upper bound on the spectral */
/*                radius of the Jacobian of f at (t,y). */

/*  INFO(3):  If you know that the Jacobian is constant, you should say so. */

/*  INFO(3)  = 0  The Jacobian may not be constant. */

/*           = 1  The Jacobian is constant. */

/*  INFO(4):  You must tell RKC whether ATOL is a scalar or an array. */

/*  INFO(4)  = 0  ATOL is a double precision scalar. */

/*           = 1  ATOL is a double precision array of length NEQN. */

/*  WORK(*):  Work array.  Double precision array of length at least */
/*            8 + 5*NEQN if INFO(2) = 0 and otherwise, 8 + 4*NEQN. */

/*  IDID:     Set IDID = 0 to initialize the integration. */



/*  RETURNS FROM RKC */

/*  T:        The integration has advanced to T. */

/*  Y(*):     The solution at T. */

/*  IDID:     The value of IDID reports what happened. */

/*                          SUCCESS */

/*    IDID     = 1 T = TEND, so the integration is complete. */

/*             = 2 Took a step to the output value of T.  To continue on */
/*                 towards TEND, just call RKC again.   WARNING:  Do not */
/*                 alter any argument between calls. */

/*                 The last step, HLAST, is returned as WORK(1). RKCINT */
/*                 can be used to approximate the solution anywhere in */
/*                 [T-HLAST, T] very inexpensively using data in WORK(*). */

/*                 The work can be monitored by inspecting data in RKCDID. */

/*                          FAILURE */

/*             = 3 Improper error control: For some j, ATOL(j) = 0 */
/*                 and Y(j) = 0. */

/*             = 4 Unable to achieve the desired accuracy with the */
/*                 precision available.  A severe lack of smoothness in */
/*                 the solution y(t) or the function f(t,y) is likely. */

/*             = 5 Invalid input parameters:  NEQN <= 0, RTOL > 0.1, */
/*                 RTOL < 10*UROUND, or ATOL(i) < 0 for some i. */

/*             = 6 The method used by RKC to estimate the spectral */
/*                 radius of the Jacobian failed to converge. */

/*  RKCDID is a labelled common block that communicates statistics */
/*         about the integration process: */
/*         common /rkcdid/   nfe,nsteps,naccpt,nrejct,nfesig,maxm */

/*         The integer counters are: */

/*        NFE      number of evaluations of F used */
/*                   to integrate the initial value problem */
/*        NSTEPS   number of integration steps */
/*        NACCPT   number of accepted steps */
/*        NREJCT   number of rejected steps */
/*        NFESIG   number of evaluations of F used */
/*                   to estimate the spectral radius */
/*        MAXM     maximum number of stages used */

/*        This data can be used to monitor the work and terminate a run */
/*        that proves to be unacceptably expensive.  Also, if MAXM should */
/*        be far beyond 100, the problem is too expensive for RKC and */
/*        alternative methods should be considered. */

/* -------------------------------------------------------------------------- */

/*   CAUTION: MACHINE/PRECISION ISSUES */

/*     UROUND (the machine precision) is the smallest number such that */
/*     1 + UROUND > 1, where 1 is a floating point number in the working */
/*     precision. UROUND is set in a parameter statement in RKC. Its */
/*     value depends on both the precision and the machine used, so it */
/*     must be set appropriately.  UROUND is the only constant in RKC */
/*     that depends on the precision. */

/*     This version of RKC is written in double precision. It can be changed */
/*     to single precision by replacing DOUBLE PRECISION in the declarations */
/*     by REAL and changing the type of the floating point constants set in */
/*     PARAMETER statements from double precision to real. */

/* -------------------------------------------------------------------------- */

/*  Authors: B.P. Sommeijer and J.G. Verwer */
/*           Centre for Mathematics and Computer Science (CWI) */
/*           Kruislaan 413 */
/*           1098 SJ  Amsterdam */
/*           The Netherlands */
/*           e-mail: bsom@cwi.nl */

/*           L.F. Shampine */
/*           Mathematics Department */
/*           Southern Methodist University */
/*           Dallas, Texas 75275-0156 */
/*           USA */
/*           e-mail: lshampin@mail.smu.edu */

/*  Details of the methods used and the performance of RKC can be */
/*  found in */

/*         B.P. Sommeijer, L.F. Shampine and J.G. Verwer */
/*         RKC: an Explicit Solver for Parabolic PDEs. */
/*              Technical Report MAS-R9715, CWI, Amsterdam, 1997 */

/*  This source code for RKC and some examples, as well as the */
/*  reference solution to the second example can also be obtained */
/*  by anonymous ftp from the address ftp://ftp.cwi.nl/pub/bsom/rkc */
/* ------------------------------------------------------------------ */

/* ********************************************************************* */
/*  uround is set here for IEEE double precision arithmetic. */
/* ********************************************************************* */


    /* Parameter adjustments */
    --y;
    --atol;
    --info;
    --work;

    /* Function Body */
    if (*idid == 0) {
/* ---------------------- */
/*  Test the input data. */
/* ---------------------- */
	array = info[4] == 1;
	valid = *neqn > 0;
	if (*rtol > .1 || *rtol < 2.2200000000000002e-15) {
	    valid = FALSE_;
	}
	if (atol[1] < 0.) {
	    valid = FALSE_;
	}
	if (array) {
	    i__1 = *neqn;
	    for (i__ = 2; i__ <= i__1; ++i__) {
		if (atol[i__] < 0.) {
		    valid = FALSE_;
		}
/* L10: */
	    }
	}
	if (! valid) {
	    *idid = 5;
	    return 0;
	}
/* ----------------------------------- */
/*  Initialize counters and pointers. */
/* ----------------------------------- */
	rkcdid_1.nfe = 0;
	rkcdid_1.nsteps = 0;
	rkcdid_1.naccpt = 0;
	rkcdid_1.nrejct = 0;
	rkcdid_1.nfesig = 0;
	rkcdid_1.maxm = 0;
/* ----------------------------------------------------------- */
/*  work(*) contains information needed for interpolation, */
/*  continuation after a return, and working storage. Items */
/*  relevant here are: */

/*  The last step taken, hlast, is work(1). */
/*  The current t is work(2). */
/*  The number of equations, neqn, is work(3). */
/*  The unit roundoff, uround, is work(4). */
/*  The square root of uround, sqrtu, is work(5). */
/*  The maximum step size, hmax, is work(6). */
/*  The base address for the solution is ptr1 = nint(work(7)). */
/*  The solution at t starts at ptr1. */
/*  The derivative of the solution at t starts at ptr2. */
/*  The solution at t-hlast starts at ptr3. */
/*  The derivative of the solution at t-hlast starts at ptr4. */
/*  The estimated dominant eigenvector starts at ptr4 + neqn. */
/* ------------------------------------------------------------ */
	work[2] = *t;
	work[3] = (doublereal) (*neqn);
	work[4] = 2.22e-16;
	work[5] = sqrt(2.22e-16);
	ptr1 = 8;
	work[7] = (doublereal) ptr1;
	ptr2 = ptr1 + *neqn;
	ptr3 = ptr2 + *neqn;
	ptr4 = ptr3 + *neqn;
    } else if (*idid != 2) {
	s_wsle(&io___8);
	do_lio(&c__9, &c__1, " RKC was called with an illegal value of IDID.",
		 (ftnlen)46);
	e_wsle();
	s_stop("", (ftnlen)0);
    }

    rkclow_(neqn, t, tend, &y[1], (U_fp)f, &info[1], rtol, &atol[1], &work[1],
	     &work[ptr1], &work[ptr2], &work[ptr3], &work[ptr4], idid);
    return 0;
} /* rkc_ */

/* Subroutine */ int rkclow_(integer *neqn, doublereal *t, doublereal *tend, 
	doublereal *y, S_fp f, integer *info, doublereal *rtol, doublereal *
	atol, doublereal *work, doublereal *yn, doublereal *fn, doublereal *
	vtemp1, doublereal *vtemp2, integer *idid)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double sqrt(doublereal);
    integer i_dnnt(doublereal *);
    double d_sign(doublereal *, doublereal *), pow_dd(doublereal *, 
	    doublereal *);

    /* Local variables */
    static doublereal h__;
    static integer i__, m;
    static doublereal at, wt, fac, err, est, absh, hold, hmin, hmax;
    static integer mmax;
    static doublereal tdir;
    static logical last;
    extern /* Subroutine */ int step_(integer *, S_fp, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, doublereal *,
	     doublereal *, doublereal *);
    static doublereal temp1, temp2, sprad;
    static logical array;
    static doublereal ylast;
    static logical jacatt;
    extern doublereal spcrad_(integer *, doublereal *, doublereal *);
    static doublereal errold;
    extern /* Subroutine */ int rkcrho_(integer *, doublereal *, S_fp, 
	    doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *);
    static logical newspc;
    static integer nstsig;
    static doublereal uround, yplast;

/* ---------------------------------------------------------------------- */
/*  RKC is an interface to RKCLOW where the actual solution takes place. */
/* ---------------------------------------------------------------------- */


/* --------------------------------- */
/*    Initialize on the first call. */
/* --------------------------------- */
    /* Parameter adjustments */
    --vtemp2;
    --vtemp1;
    --fn;
    --yn;
    --work;
    --atol;
    --info;
    --y;

    /* Function Body */
    if (*idid == 0) {
	array = info[4] == 1;
	uround = work[4];
	d__1 = sqrt(*rtol / (uround * 10.));
	mmax = i_dnnt(&d__1);
	mmax = max(mmax,2);
	newspc = TRUE_;
	jacatt = FALSE_;
	nstsig = 0;
	i__1 = *neqn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    yn[i__] = y[i__];
/* L10: */
	}
	(*f)(neqn, t, &yn[1], &fn[1]);
	++rkcdid_1.nfe;
	d__1 = *tend - *t;
	tdir = d_sign(&c_b9, &d__1);
	hmax = (d__1 = *tend - *t, abs(d__1));
	work[6] = hmax;
/* Computing MAX */
	d__1 = abs(*t);
	hmin = uround * 10. * max(d__1,hmax);
    }
/* ------------------------------------ */
/*  Start of loop for taking one step. */
/* ------------------------------------ */
L20:
/* ---------------------------------------------- */
/*  Estimate the spectral radius of the Jacobian */
/*  when newspc = .true..  A convergence failure */
/*  in rkcrho is reported by idid = 6. */
/* ---------------------------------------------- */
    if (newspc) {
	if (info[2] == 1) {
	    sprad = spcrad_(neqn, t, &yn[1]);
	} else {
	    rkcrho_(neqn, t, (S_fp)f, &yn[1], &fn[1], &vtemp1[1], &vtemp2[1], 
		    &work[1], &sprad, idid);
	    if (*idid == 6) {
		return 0;
	    }
	}
	jacatt = TRUE_;
    }
/* ------------------------------- */
/*  Compute an initial step size. */
/* ------------------------------- */
    if (rkcdid_1.nsteps == 0) {
	absh = hmax;
	if (sprad * absh > 1.) {
	    absh = 1. / sprad;
	}
	absh = max(absh,hmin);
	i__1 = *neqn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    vtemp1[i__] = yn[i__] + absh * fn[i__];
/* L30: */
	}
	d__1 = *t + absh;
	(*f)(neqn, &d__1, &vtemp1[1], &vtemp2[1]);
	++rkcdid_1.nfe;
	est = 0.;
	at = atol[1];
	i__1 = *neqn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    if (array) {
		at = atol[i__];
	    }
	    wt = at + *rtol * (d__1 = yn[i__], abs(d__1));
	    if (wt == 0.) {
		*idid = 3;
		return 0;
	    }
/* Computing 2nd power */
	    d__1 = (vtemp2[i__] - fn[i__]) / wt;
	    est += d__1 * d__1;
/* L40: */
	}
	est = absh * sqrt(est / *neqn);
	if (absh * .1 < hmax * sqrt(est)) {
/* Computing MAX */
	    d__1 = absh * .1 / sqrt(est);
	    absh = max(d__1,hmin);
	} else {
	    absh = hmax;
	}
    }
/* ------------------------------------------------------------ */
/*  Adjust the step size and determine the number of stages m. */
/* ------------------------------------------------------------ */
    last = FALSE_;
    if (absh * 1.1 >= (d__1 = *tend - *t, abs(d__1))) {
	absh = (d__1 = *tend - *t, abs(d__1));
	last = TRUE_;
    }
    m = (integer) sqrt(absh * 1.54 * sprad + 1.) + 1;
/* ---------------------------------------------------------- */
/*  Limit m to mmax to control the growth of roundoff error. */
/* ---------------------------------------------------------- */
    if (m > mmax) {
	m = mmax;
/* Computing 2nd power */
	i__1 = m;
	absh = (i__1 * i__1 - 1) / (sprad * 1.54);
	last = FALSE_;
    }
    rkcdid_1.maxm = max(m,rkcdid_1.maxm);
/* -------------------------------------------- */
/*  A tentative solution at t+h is returned in */
/*  y and its slope is evaluated in vtemp1(*). */
/* -------------------------------------------- */
    h__ = tdir * absh;
/* Computing MAX */
    d__2 = abs(*t), d__3 = (d__1 = *t + h__, abs(d__1));
    hmin = uround * 10. * max(d__2,d__3);
    step_(neqn, (S_fp)f, t, &yn[1], &fn[1], &h__, &m, &y[1], &vtemp1[1], &
	    vtemp2[1]);
    d__1 = *t + h__;
    (*f)(neqn, &d__1, &y[1], &vtemp1[1]);
    rkcdid_1.nfe += m;
    ++rkcdid_1.nsteps;
/* ------------------------------------------------------------- */
/*  Estimate the local error and compute its weighted RMS norm. */
/* ------------------------------------------------------------- */
    err = 0.;
    at = atol[1];
    i__1 = *neqn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (array) {
	    at = atol[i__];
	}
/* Computing MAX */
	d__3 = (d__1 = y[i__], abs(d__1)), d__4 = (d__2 = yn[i__], abs(d__2));
	wt = at + *rtol * max(d__3,d__4);
	if (wt == 0.) {
	    *idid = 3;
	    return 0;
	}
	est = (yn[i__] - y[i__]) * .8 + h__ * .4 * (fn[i__] + vtemp1[i__]);
/* Computing 2nd power */
	d__1 = est / wt;
	err += d__1 * d__1;
/* L50: */
    }
    err = sqrt(err / *neqn);

    if (err > 1.) {
/* ------------------- */
/*  Step is rejected. */
/* ------------------- */
	++rkcdid_1.nrejct;
	absh = absh * .8 / pow_dd(&err, &c_b14);
	if (absh < hmin) {
	    *idid = 4;
	    return 0;
	} else {
	    newspc = ! jacatt;
	    goto L20;
	}
    }
/* ------------------- */
/*  Step is accepted. */
/* ------------------- */
    ++rkcdid_1.naccpt;
    *t += h__;
    jacatt = info[3] == 1;
    nstsig = (nstsig + 1) % 25;
    newspc = FALSE_;
    if (info[2] == 1 || nstsig == 0) {
	newspc = ! jacatt;
    }
/* ------------------------------------------------------ */
/*  Update the data for interpolation stored in work(*). */
/* ------------------------------------------------------ */
    work[1] = h__;
    work[2] = *t;
    i__1 = *neqn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ylast = yn[i__];
	yplast = fn[i__];
	yn[i__] = y[i__];
	fn[i__] = vtemp1[i__];
	vtemp1[i__] = ylast;
	vtemp2[i__] = yplast;
/* L60: */
    }
    fac = 10.;
    if (rkcdid_1.naccpt == 1) {
	temp2 = pow_dd(&err, &c_b14);
	if (.8 < fac * temp2) {
	    fac = .8 / temp2;
	}
    } else {
	temp1 = absh * .8 * pow_dd(&errold, &c_b14);
	temp2 = abs(hold) * pow_dd(&err, &c_b18);
	if (temp1 < fac * temp2) {
	    fac = temp1 / temp2;
	}
    }
    absh = max(.1,fac) * absh;
/* Computing MAX */
    d__1 = hmin, d__2 = min(hmax,absh);
    absh = max(d__1,d__2);
    errold = err;
    hold = h__;
    h__ = tdir * absh;
    if (last) {
	*idid = 1;
	return 0;
    } else if (info[1] == 0) {
	*idid = 2;
	return 0;
    } else {
	goto L20;
    }
    return 0;
} /* rkclow_ */

/* Subroutine */ int step_(integer *neqn, S_fp f, doublereal *t, doublereal *
	yn, doublereal *fn, doublereal *h__, integer *m, doublereal *y, 
	doublereal *yjm1, doublereal *yjm2)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal), log(doublereal), sinh(doublereal), cosh(
	    doublereal);

    /* Local variables */
    integer i__, j;
    doublereal w0, w1, bj, mu, nu, zj, arg, thj, dzj, mus, ajm1, bjm1, bjm2, 
	    d2zj, zjm1, zjm2, thjm1, thjm2, dzjm1, dzjm2, temp1, temp2, 
	    d2zjm1, d2zjm2;

/* -------------------------------------------------- */
/*  Take a step of size H from T to T+H to get Y(*). */
/* -------------------------------------------------- */


    /* Parameter adjustments */
    --yjm2;
    --yjm1;
    --y;
    --fn;
    --yn;

    /* Function Body */
/* Computing 2nd power */
    i__1 = *m;
    w0 = 2. / (i__1 * i__1 * 13.) + 1.;
/* Computing 2nd power */
    d__1 = w0;
    temp1 = d__1 * d__1 - 1.;
    temp2 = sqrt(temp1);
    arg = *m * log(w0 + temp2);
    w1 = sinh(arg) * temp1 / (cosh(arg) * *m * temp2 - w0 * sinh(arg));
/* Computing 2nd power */
    d__1 = w0 * 2.;
    bjm1 = 1. / (d__1 * d__1);
    bjm2 = bjm1;
/* --------------------------- */
/*  Evaluate the first stage. */
/* --------------------------- */
    i__1 = *neqn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	yjm2[i__] = yn[i__];
/* L10: */
    }
    mus = w1 * bjm1;
    i__1 = *neqn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	yjm1[i__] = yn[i__] + *h__ * mus * fn[i__];
/* L20: */
    }
    thjm2 = 0.;
    thjm1 = mus;
    zjm1 = w0;
    zjm2 = 1.;
    dzjm1 = 1.;
    dzjm2 = 0.;
    d2zjm1 = 0.;
    d2zjm2 = 0.;
/* ------------------------------ */
/*  Evaluate stages j = 2,...,m. */
/* ------------------------------ */
    i__1 = *m;
    for (j = 2; j <= i__1; ++j) {
	zj = w0 * 2. * zjm1 - zjm2;
	dzj = w0 * 2. * dzjm1 - dzjm2 + zjm1 * 2.;
	d2zj = w0 * 2. * d2zjm1 - d2zjm2 + dzjm1 * 4.;
/* Computing 2nd power */
	d__1 = dzj;
	bj = d2zj / (d__1 * d__1);
	ajm1 = 1. - zjm1 * bjm1;
	mu = w0 * 2. * bj / bjm1;
	nu = -bj / bjm2;
	mus = mu * w1 / w0;
/* --------------------------------------------- */
/*  Use the y array for temporary storage here. */
/* --------------------------------------------- */
	d__1 = *t + *h__ * thjm1;
	(*f)(neqn, &d__1, &yjm1[1], &y[1]);
	i__2 = *neqn;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    y[i__] = mu * yjm1[i__] + nu * yjm2[i__] + (1. - mu - nu) * yn[
		    i__] + *h__ * mus * (y[i__] - ajm1 * fn[i__]);
/* L30: */
	}
	thj = mu * thjm1 + nu * thjm2 + mus * (1. - ajm1);
/* ------------------------------------ */
/*  Shift the data for the next stage. */
/* ------------------------------------ */
	if (j < *m) {
	    i__2 = *neqn;
	    for (i__ = 1; i__ <= i__2; ++i__) {
		yjm2[i__] = yjm1[i__];
		yjm1[i__] = y[i__];
/* L40: */
	    }
	    thjm2 = thjm1;
	    thjm1 = thj;
	    bjm2 = bjm1;
	    bjm1 = bj;
	    zjm2 = zjm1;
	    zjm1 = zj;
	    dzjm2 = dzjm1;
	    dzjm1 = dzj;
	    d2zjm2 = d2zjm1;
	    d2zjm1 = d2zj;
	}
/* L50: */
    }
    return 0;
} /* step_ */

/* Subroutine */ int rkcint_(doublereal *work, doublereal *arg, doublereal *
	yarg)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer i_dnnt(doublereal *);

    /* Local variables */
    integer i__;
    doublereal s, t, a1, a2, b1, b2;
    integer ptr1, ptr2, ptr3, ptr4, neqn;
    doublereal hlast, tlast;

/* ------------------------------------------------------------------------- */
/*  RKCINT is used to compute approximate solutions at specific t and to */
/*  compute cheaply the large number of approximations that may be needed */
/*  for plotting or locating when events occur. */

/*  After a step to T, RKC provides HLAST, the step just taken, in WORK(1). */
/*  In other entries of WORK(*) it provides the data needed to interpolate */
/*  anywhere in [T-HLAST, T]. YARG(*), the approximate solution at t = ARG */
/*  computed by interpolation in RKCINT has the same order of accuracy as */
/*  the Y(*) computed directly by RKC. */

/*  INPUT: */

/*    WORK(*)   Double precision array returned by RKC. */

/*    ARG       The point at which a solution is desired. Double precision. */

/*  OUTPUT: */

/*    YARG(*)   The approximate solution at t = ARG.  Double precision */
/*              array of length neqn. */
/* -------------------------------------------------------------------------- */


/* --------------------------------------------------------------------- */
/*  The data needed for interpolation are stored in work(*) as follows: */

/*  The last step taken, hlast, is work(1). */
/*  The current t is work(2). */
/*  The number of equations, neqn, is work(3). */
/*  The base address for the solution is ptr1 = nint(work(7)) */
/*  The solution at t starts at ptr1. */
/*  The derivative of the solution at t starts at ptr2. */
/*  The solution at t-hlast starts at ptr3. */
/*  The derivative of the solution at t-hlast starts at ptr4. */
/* --------------------------------------------------------------------- */
    /* Parameter adjustments */
    --yarg;
    --work;

    /* Function Body */
    hlast = work[1];
    t = work[2];
    tlast = t - hlast;
    neqn = i_dnnt(&work[3]);
    ptr1 = i_dnnt(&work[7]);
    ptr2 = ptr1 + neqn;
    ptr3 = ptr2 + neqn;
    ptr4 = ptr3 + neqn;

    s = (*arg - tlast) / hlast;
/* Computing 2nd power */
    d__1 = s - 1.;
    a1 = (s * 2. + 1.) * (d__1 * d__1);
/* Computing 2nd power */
    d__1 = s;
    a2 = (3. - s * 2.) * (d__1 * d__1);
/* Computing 2nd power */
    d__1 = s - 1.;
    b1 = hlast * s * (d__1 * d__1);
/* Computing 2nd power */
    d__1 = s;
    b2 = hlast * (s - 1.) * (d__1 * d__1);

    i__1 = neqn;
    for (i__ = 1; i__ <= i__1; ++i__) {
	yarg[i__] = a1 * work[ptr3 + i__ - 1] + a2 * work[ptr1 + i__ - 1] + 
		b1 * work[ptr4 + i__ - 1] + b2 * work[ptr2 + i__ - 1];
/* L10: */
    }
    return 0;
} /* rkcint_ */

/* Subroutine */ int rkcrho_(integer *neqn, doublereal *t, S_fp f, doublereal 
	*yn, doublereal *fn, doublereal *v, doublereal *fv, doublereal *work, 
	doublereal *sprad, integer *idid)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Builtin functions */
    integer i_dnnt(doublereal *);
    double sqrt(doublereal);

    /* Local variables */
    integer i__, ptr5, iter;
    doublereal vnrm, ynrm, sigma, dfnrm;
    integer index;
    doublereal small, dynrm, sqrtu, sigmal, uround;

/* --------------------------------------------------------------- */
/*  RKCRHO attempts to compute a close upper bound, SPRAD, on */
/*  the spectral radius of the Jacobian matrix using a nonlinear */
/*  power method.  A convergence failure is reported by IDID = 6. */
/* --------------------------------------------------------------- */


    /* Parameter adjustments */
    --fv;
    --v;
    --fn;
    --yn;
    --work;

    /* Function Body */
    uround = work[4];
    sqrtu = work[5];
/* ------------------------------------------------------------ */
/*  hmax = work(6).  sprad smaller than small = 1/hmax are not */
/*  interesting because they do not constrain the step size. */
/* ------------------------------------------------------------ */
    small = 1. / work[6];
/* --------------------------------------------------------- */
/*  The initial slope is used as guess when nsteps = 0 and */
/*  thereafter the last computed eigenvector.  Some care */
/*  is needed to deal with special cases. Approximations to */
/*  the eigenvector are normalized so that their Euclidean */
/*  norm has the constant value dynrm. */
/* --------------------------------------------------------- */
    ptr5 = i_dnnt(&work[7]) + (*neqn << 2);
    if (rkcdid_1.nsteps == 0) {
	i__1 = *neqn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    v[i__] = fn[i__];
/* L10: */
	}
    } else {
	i__1 = *neqn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    v[i__] = work[ptr5 + i__ - 1];
/* L20: */
	}
    }
    ynrm = 0.;
    vnrm = 0.;
    i__1 = *neqn;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = yn[i__];
	ynrm += d__1 * d__1;
/* Computing 2nd power */
	d__1 = v[i__];
	vnrm += d__1 * d__1;
/* L30: */
    }
    ynrm = sqrt(ynrm);
    vnrm = sqrt(vnrm);
    if (ynrm != 0. && vnrm != 0.) {
	dynrm = ynrm * sqrtu;
	i__1 = *neqn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    v[i__] = yn[i__] + v[i__] * (dynrm / vnrm);
/* L40: */
	}
    } else if (ynrm != 0.) {
	dynrm = ynrm * sqrtu;
	i__1 = *neqn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    v[i__] = yn[i__] + yn[i__] * sqrtu;
/* L50: */
	}
    } else if (vnrm != 0.) {
	dynrm = uround;
	i__1 = *neqn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    v[i__] *= dynrm / vnrm;
/* L60: */
	}
    } else {
	dynrm = uround;
	i__1 = *neqn;
	for (i__ = 1; i__ <= i__1; ++i__) {
	    v[i__] = dynrm;
/* L70: */
	}
    }
/* -------------------------------------------- */
/*  Now iterate with a nonlinear power method. */
/* -------------------------------------------- */
    sigma = 0.;
    for (iter = 1; iter <= 50; ++iter) {
	(*f)(neqn, t, &v[1], &fv[1]);
	++rkcdid_1.nfesig;
	dfnrm = 0.;
	i__1 = *neqn;
	for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	    d__1 = fv[i__] - fn[i__];
	    dfnrm += d__1 * d__1;
/* L80: */
	}
	dfnrm = sqrt(dfnrm);
	sigmal = sigma;
	sigma = dfnrm / dynrm;
/* ---------------------------------------------------------- */
/*  sprad is a little bigger than the estimate sigma of the */
/*  spectral radius, so is more likely to be an upper bound. */
/* ---------------------------------------------------------- */
	*sprad = sigma * 1.2;
	if (iter >= 2 && (d__1 = sigma - sigmal, abs(d__1)) <= max(sigma,
		small) * .01) {
	    i__1 = *neqn;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		work[ptr5 + i__ - 1] = v[i__] - yn[i__];
/* L90: */
	    }
	    return 0;
	}
/* -------------------------------------- */
/*  The next v(*) is the change in f */
/*  scaled so that norm(v - yn) = dynrm. */
/* -------------------------------------- */
	if (dfnrm != 0.) {
	    i__1 = *neqn;
	    for (i__ = 1; i__ <= i__1; ++i__) {
		v[i__] = yn[i__] + (fv[i__] - fn[i__]) * (dynrm / dfnrm);
/* L100: */
	    }
	} else {
/* ------------------------------------------------------- */
/*  The new v(*) degenerated to yn(*)--"randomly" perturb */
/*  current approximation to the eigenvector by changing */
/*  the sign of one component. */
/* ------------------------------------------------------- */
	    index = iter % *neqn + 1;
	    v[index] = yn[index] - (v[index] - yn[index]);
	}
/* L110: */
    }
/* ------------------------------------------- */
/*  Set flag to report a convergence failure. */
/* ------------------------------------------- */
    *idid = 6;
    return 0;
} /* rkcrho_ */

