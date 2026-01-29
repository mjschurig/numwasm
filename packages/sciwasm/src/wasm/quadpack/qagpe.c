/* qagpe.f -- translated by f2c (version 20240504).
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

/* Table of constant values */

static integer c__4 = 4;
static integer c__1 = 1;
static integer c__2 = 2;

/* Subroutine */ int qagpe_(E_fp f, real *a, real *b, integer *npts2, real *
	points, real *epsabs, real *epsrel, integer *limit, real *result, 
	real *abserr, integer *neval, integer *ier, real *alist__, real *
	blist, real *rlist, real *elist, real *pts, integer *iord, integer *
	level, integer *ndin, integer *last)
{
    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2;

    /* Local variables */
    integer i__, j, k;
    real a1, a2, b1, b2;
    integer id, ip1;
    extern /* Subroutine */ int qk21_(E_fp, real *, real *, real *, real *, 
	    real *, real *);
    integer ind1, ind2;
    real area;
    extern /* Subroutine */ int qelg_(integer *, real *, real *, real *, real 
	    *, integer *);
    real resa, dres, sign;
    integer ksgn;
    real temp;
    integer nres, nint, jlow, npts;
    real area1, area2, area12, erro12;
    integer ierro;
    real defab1, defab2;
    integer ktmin, nrmax;
    real oflow, uflow;
    logical noext;
    extern /* Subroutine */ int qpsrt_(integer *, integer *, integer *, real *
	    , real *, integer *, integer *);
    extern doublereal r1mach_(integer *);
    integer iroff1, iroff2, iroff3;
    real res3la[3];
    integer nintp1;
    real error1, error2, rlist2[52];
    integer numrl2;
    real defabs, epmach, erlarg, abseps, correc, errbnd, resabs;
    integer jupbnd;
    real erlast;
    integer levmax;
    real errmax;
    integer maxerr, levcur;
    real reseps;
    logical extrap;
    real ertest, errsum;

/* ***begin prologue  qagpe */
/* ***date written   800101   (yymmdd) */
/* ***revision date  830518   (yymmdd) */
/* ***category no.  h2a2a1 */
/* ***keywords  automatic integrator, general-purpose, */
/*             singularities at user specified points, */
/*             extrapolation, globally adaptive. */
/* ***author  piessens,robert ,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
/* ***purpose  the routine calculates an approximation result to a given */
/*            definite integral i = integral of f over (a,b),hopefully */
/*            satisfying following claim for accuracy abs(i-result).le. */
/*            max(epsabs,epsrel*abs(i)). break points of the integration */
/*            interval, where local difficulties of the integrand may */
/*            occur(e.g. singularities,discontinuities),provided by user. */
/* ***description */

/*        computation of a definite integral */
/*        standard fortran subroutine */
/*        real version */

/*        parameters */
/*         on entry */
/*            f      - real */
/*                     function subprogram defining the integrand */
/*                     function f(x). the actual name for f needs to be */
/*                     declared e x t e r n a l in the driver program. */

/*            a      - real */
/*                     lower limit of integration */

/*            b      - real */
/*                     upper limit of integration */

/*            npts2  - integer */
/*                     number equal to two more than the number of */
/*                     user-supplied break points within the integration */
/*                     range, npts2.ge.2. */
/*                     if npts2.lt.2, the routine will end with ier = 6. */

/*            points - real */
/*                     vector of dimension npts2, the first (npts2-2) */
/*                     elements of which are the user provided break */
/*                     points. if these points do not constitute an */
/*                     ascending sequence there will be an automatic */
/*                     sorting. */

/*            epsabs - real */
/*                     absolute accuracy requested */
/*            epsrel - real */
/*                     relative accuracy requested */
/*                     if  epsabs.le.0 */
/*                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
/*                     the routine will end with ier = 6. */

/*            limit  - integer */
/*                     gives an upper bound on the number of subintervals */
/*                     in the partition of (a,b), limit.ge.npts2 */
/*                     if limit.lt.npts2, the routine will end with */
/*                     ier = 6. */

/*         on return */
/*            result - real */
/*                     approximation to the integral */

/*            abserr - real */
/*                     estimate of the modulus of the absolute error, */
/*                     which should equal or exceed abs(i-result) */

/*            neval  - integer */
/*                     number of integrand evaluations */

/*            ier    - integer */
/*                     ier = 0 normal and reliable termination of the */
/*                             routine. it is assumed that the requested */
/*                             accuracy has been achieved. */
/*                     ier.gt.0 abnormal termination of the routine. */
/*                             the estimates for integral and error are */
/*                             less reliable. it is assumed that the */
/*                             requested accuracy has not been achieved. */
/*            error messages */
/*                     ier = 1 maximum number of subdivisions allowed */
/*                             has been achieved. one can allow more */
/*                             subdivisions by increasing the value of */
/*                             limit (and taking the according dimension */
/*                             adjustments into account). however, if */
/*                             this yields no improvement it is advised */
/*                             to analyze the integrand in order to */
/*                             determine the integration difficulties. if */
/*                             the position of a local difficulty can be */
/*                             determined (i.e. singularity, */
/*                             discontinuity within the interval), it */
/*                             should be supplied to the routine as an */
/*                             element of the vector points. if necessary */
/*                             an appropriate special-purpose integrator */
/*                             must be used, which is designed for */
/*                             handling the type of difficulty involved. */
/*                         = 2 the occurrence of roundoff error is */
/*                             detected, which prevents the requested */
/*                             tolerance from being achieved. */
/*                             the error may be under-estimated. */
/*                         = 3 extremely bad integrand behaviour occurs */
/*                             at some points of the integration */
/*                             interval. */
/*                         = 4 the algorithm does not converge. */
/*                             roundoff error is detected in the */
/*                             extrapolation table. it is presumed that */
/*                             the requested tolerance cannot be */
/*                             achieved, and that the returned result is */
/*                             the best which can be obtained. */
/*                         = 5 the integral is probably divergent, or */
/*                             slowly convergent. it must be noted that */
/*                             divergence can occur with any other value */
/*                             of ier.gt.0. */
/*                         = 6 the input is invalid because */
/*                             npts2.lt.2 or */
/*                             break points are specified outside */
/*                             the integration range or */
/*                             (epsabs.le.0 and */
/*                              epsrel.lt.max(50*rel.mach.acc.,0.5d-28)) */
/*                             or limit.lt.npts2. */
/*                             result, abserr, neval, last, rlist(1), */
/*                             and elist(1) are set to zero. alist(1) and */
/*                             blist(1) are set to a and b respectively. */

/*            alist  - real */
/*                     vector of dimension at least limit, the first */
/*                      last  elements of which are the left end points */
/*                     of the subintervals in the partition of the given */
/*                     integration range (a,b) */

/*            blist  - real */
/*                     vector of dimension at least limit, the first */
/*                      last  elements of which are the right end points */
/*                     of the subintervals in the partition of the given */
/*                     integration range (a,b) */

/*            rlist  - real */
/*                     vector of dimension at least limit, the first */
/*                      last  elements of which are the integral */
/*                     approximations on the subintervals */

/*            elist  - real */
/*                     vector of dimension at least limit, the first */
/*                      last  elements of which are the moduli of the */
/*                     absolute error estimates on the subintervals */

/*            pts    - real */
/*                     vector of dimension at least npts2, containing the */
/*                     integration limits and the break points of the */
/*                     interval in ascending sequence. */

/*            level  - integer */
/*                     vector of dimension at least limit, containing the */
/*                     subdivision levels of the subinterval, i.e. if */
/*                     (aa,bb) is a subinterval of (p1,p2) where p1 as */
/*                     well as p2 is a user-provided break point or */
/*                     integration limit, then (aa,bb) has level l if */
/*                     abs(bb-aa) = abs(p2-p1)*2**(-l). */

/*            ndin   - integer */
/*                     vector of dimension at least npts2, after first */
/*                     integration over the intervals (pts(i)),pts(i+1), */
/*                     i = 0,1, ..., npts2-2, the error estimates over */
/*                     some of the intervals may have been increased */
/*                     artificially, in order to put their subdivision */
/*                     forward. if this happens for the subinterval */
/*                     numbered k, ndin(k) is put to 1, otherwise */
/*                     ndin(k) = 0. */

/*            iord   - integer */
/*                     vector of dimension at least limit, the first k */
/*                     elements of which are pointers to the */
/*                     error estimates over the subintervals, */
/*                     such that elist(iord(1)), ..., elist(iord(k)) */
/*                     form a decreasing sequence, with k = last */
/*                     if last.le.(limit/2+2), and k = limit+1-last */
/*                     otherwise */

/*            last   - integer */
/*                     number of subintervals actually produced in the */
/*                     subdivisions process */

/* ***references  (none) */
/* ***routines called  qelg,qk21,qpsrt,r1mach */
/* ***end prologue  qagpe */




/*            the dimension of rlist2 is determined by the value of */
/*            limexp in subroutine epsalg (rlist2 should be of dimension */
/*            (limexp+2) at least). */


/*            list of major variables */
/*            ----------------------- */

/*           alist     - list of left end points of all subintervals */
/*                       considered up to now */
/*           blist     - list of right end points of all subintervals */
/*                       considered up to now */
/*           rlist(i)  - approximation to the integral over */
/*                       (alist(i),blist(i)) */
/*           rlist2    - array of dimension at least limexp+2 */
/*                       containing the part of the epsilon table which */
/*                       is still needed for further computations */
/*           elist(i)  - error estimate applying to rlist(i) */
/*           maxerr    - pointer to the interval with largest error */
/*                       estimate */
/*           errmax    - elist(maxerr) */
/*           erlast    - error on the interval currently subdivided */
/*                       (before that subdivision has taken place) */
/*           area      - sum of the integrals over the subintervals */
/*           errsum    - sum of the errors over the subintervals */
/*           errbnd    - requested accuracy max(epsabs,epsrel* */
/*                       abs(result)) */
/*           *****1    - variable for the left subinterval */
/*           *****2    - variable for the right subinterval */
/*           last      - index for subdivision */
/*           nres      - number of calls to the extrapolation routine */
/*           numrl2    - number of elements in rlist2. if an */
/*                       appropriate approximation to the compounded */
/*                       integral has been obtained, it is put in */
/*                       rlist2(numrl2) after numrl2 has been increased */
/*                       by one. */
/*           erlarg    - sum of the errors over the intervals larger */
/*                       than the smallest interval considered up to now */
/*           extrap    - logical variable denoting that the routine */
/*                       is attempting to perform extrapolation. i.e. */
/*                       before subdividing the smallest interval we */
/*                       try to decrease the value of erlarg. */
/*           noext     - logical variable denoting that extrapolation is */
/*                       no longer allowed (true-value) */

/*            machine dependent constants */
/*            --------------------------- */

/*           epmach is the largest relative spacing. */
/*           uflow is the smallest positive magnitude. */
/*           oflow is the largest positive magnitude. */

/* ***first executable statement  qagpe */
    /* Parameter adjustments */
    --ndin;
    --pts;
    --points;
    --level;
    --iord;
    --elist;
    --rlist;
    --blist;
    --alist__;

    /* Function Body */
    epmach = r1mach_(&c__4);

/*            test on validity of parameters */
/*            ----------------------------- */

    *ier = 0;
    *neval = 0;
    *last = 0;
    *result = 0.f;
    *abserr = 0.f;
    alist__[1] = *a;
    blist[1] = *b;
    rlist[1] = 0.f;
    elist[1] = 0.f;
    iord[1] = 0;
    level[1] = 0;
    npts = *npts2 - 2;
/* Computing MAX */
    r__1 = epmach * 50.f;
    if (*npts2 < 2 || *limit <= npts || *epsabs <= 0.f && *epsrel < dmax(r__1,
	    5e-15f)) {
	*ier = 6;
    }
    if (*ier == 6) {
	goto L210;
    }

/*            if any break points are provided, sort them into an */
/*            ascending sequence. */

    sign = 1.f;
    if (*a > *b) {
	sign = -1.f;
    }
    pts[1] = dmin(*a,*b);
    if (npts == 0) {
	goto L15;
    }
    i__1 = npts;
    for (i__ = 1; i__ <= i__1; ++i__) {
	pts[i__ + 1] = points[i__];
/* L10: */
    }
L15:
    pts[npts + 2] = dmax(*a,*b);
    nint = npts + 1;
    a1 = pts[1];
    if (npts == 0) {
	goto L40;
    }
    nintp1 = nint + 1;
    i__1 = nint;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ip1 = i__ + 1;
	i__2 = nintp1;
	for (j = ip1; j <= i__2; ++j) {
	    if (pts[i__] <= pts[j]) {
		goto L20;
	    }
	    temp = pts[i__];
	    pts[i__] = pts[j];
	    pts[j] = temp;
L20:
	    ;
	}
    }
    if (pts[1] != dmin(*a,*b) || pts[nintp1] != dmax(*a,*b)) {
	*ier = 6;
    }
    if (*ier == 6) {
	goto L999;
    }

/*            compute first integral and error approximations. */
/*            ------------------------------------------------ */

L40:
    resabs = 0.f;
    i__2 = nint;
    for (i__ = 1; i__ <= i__2; ++i__) {
	b1 = pts[i__ + 1];
	qk21_((E_fp)f, &a1, &b1, &area1, &error1, &defabs, &resa);
	*abserr += error1;
	*result += area1;
	ndin[i__] = 0;
	if (error1 == resa && error1 != 0.f) {
	    ndin[i__] = 1;
	}
	resabs += defabs;
	level[i__] = 0;
	elist[i__] = error1;
	alist__[i__] = a1;
	blist[i__] = b1;
	rlist[i__] = area1;
	iord[i__] = i__;
	a1 = b1;
/* L50: */
    }
    errsum = 0.f;
    i__2 = nint;
    for (i__ = 1; i__ <= i__2; ++i__) {
	if (ndin[i__] == 1) {
	    elist[i__] = *abserr;
	}
	errsum += elist[i__];
/* L55: */
    }

/*           test on accuracy. */

    *last = nint;
    *neval = nint * 21;
    dres = dabs(*result);
/* Computing MAX */
    r__1 = *epsabs, r__2 = *epsrel * dres;
    errbnd = dmax(r__1,r__2);
    if (*abserr <= epmach * 100.f * resabs && *abserr > errbnd) {
	*ier = 2;
    }
    if (nint == 1) {
	goto L80;
    }
    i__2 = npts;
    for (i__ = 1; i__ <= i__2; ++i__) {
	jlow = i__ + 1;
	ind1 = iord[i__];
	i__1 = nint;
	for (j = jlow; j <= i__1; ++j) {
	    ind2 = iord[j];
	    if (elist[ind1] > elist[ind2]) {
		goto L60;
	    }
	    ind1 = ind2;
	    k = j;
L60:
	    ;
	}
	if (ind1 == iord[i__]) {
	    goto L70;
	}
	iord[k] = iord[i__];
	iord[i__] = ind1;
L70:
	;
    }
    if (*limit < *npts2) {
	*ier = 1;
    }
L80:
    if (*ier != 0 || *abserr <= errbnd) {
	goto L999;
    }

/*           initialization */
/*           -------------- */

    rlist2[0] = *result;
    maxerr = iord[1];
    errmax = elist[maxerr];
    area = *result;
    nrmax = 1;
    nres = 0;
    numrl2 = 1;
    ktmin = 0;
    extrap = FALSE_;
    noext = FALSE_;
    erlarg = errsum;
    ertest = errbnd;
    levmax = 1;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ierro = 0;
    uflow = r1mach_(&c__1);
    oflow = r1mach_(&c__2);
    *abserr = oflow;
    ksgn = -1;
    if (dres >= (1.f - epmach * 50.f) * resabs) {
	ksgn = 1;
    }

/*           main do-loop */
/*           ------------ */

    i__2 = *limit;
    for (*last = *npts2; *last <= i__2; ++(*last)) {

/*           bisect the subinterval with the nrmax-th largest */
/*           error estimate. */

	levcur = level[maxerr] + 1;
	a1 = alist__[maxerr];
	b1 = (alist__[maxerr] + blist[maxerr]) * .5f;
	a2 = b1;
	b2 = blist[maxerr];
	erlast = errmax;
	qk21_((E_fp)f, &a1, &b1, &area1, &error1, &resa, &defab1);
	qk21_((E_fp)f, &a2, &b2, &area2, &error2, &resa, &defab2);

/*           improve previous approximations to integral */
/*           and error and test for accuracy. */

	*neval += 42;
	area12 = area1 + area2;
	erro12 = error1 + error2;
	errsum = errsum + erro12 - errmax;
	area = area + area12 - rlist[maxerr];
	if (defab1 == error1 || defab2 == error2) {
	    goto L95;
	}
	if ((r__1 = rlist[maxerr] - area12, dabs(r__1)) > dabs(area12) * 
		1e-5f || erro12 < errmax * .99f) {
	    goto L90;
	}
	if (extrap) {
	    ++iroff2;
	}
	if (! extrap) {
	    ++iroff1;
	}
L90:
	if (*last > 10 && erro12 > errmax) {
	    ++iroff3;
	}
L95:
	level[maxerr] = levcur;
	level[*last] = levcur;
	rlist[maxerr] = area1;
	rlist[*last] = area2;
/* Computing MAX */
	r__1 = *epsabs, r__2 = *epsrel * dabs(area);
	errbnd = dmax(r__1,r__2);

/*           test for roundoff error and eventually */
/*           set error flag. */

	if (iroff1 + iroff2 >= 10 || iroff3 >= 20) {
	    *ier = 2;
	}
	if (iroff2 >= 5) {
	    ierro = 3;
	}

/*           set error flag in the case that the number of */
/*           subintervals equals limit. */

	if (*last == *limit) {
	    *ier = 1;
	}

/*           set error flag in the case of bad integrand behaviour */
/*           at a point of the integration range */

/* Computing MAX */
	r__1 = dabs(a1), r__2 = dabs(b2);
	if (dmax(r__1,r__2) <= (epmach * 100.f + 1.f) * (dabs(a2) + uflow * 
		1e3f)) {
	    *ier = 4;
	}

/*           append the newly-created intervals to the list. */

	if (error2 > error1) {
	    goto L100;
	}
	alist__[*last] = a2;
	blist[maxerr] = b1;
	blist[*last] = b2;
	elist[maxerr] = error1;
	elist[*last] = error2;
	goto L110;
L100:
	alist__[maxerr] = a2;
	alist__[*last] = a1;
	blist[*last] = b1;
	rlist[maxerr] = area2;
	rlist[*last] = area1;
	elist[maxerr] = error2;
	elist[*last] = error1;

/*           call subroutine qpsrt to maintain the descending ordering */
/*           in the list of error estimates and select the */
/*           subinterval with nrmax-th largest error estimate (to be */
/*           bisected next). */

L110:
	qpsrt_(limit, last, &maxerr, &errmax, &elist[1], &iord[1], &nrmax);
/* ***jump out of do-loop */
	if (errsum <= errbnd) {
	    goto L190;
	}
/* ***jump out of do-loop */
	if (*ier != 0) {
	    goto L170;
	}
	if (noext) {
	    goto L160;
	}
	erlarg -= erlast;
	if (levcur + 1 <= levmax) {
	    erlarg += erro12;
	}
	if (extrap) {
	    goto L120;
	}

/*           test whether the interval to be bisected next is the */
/*           smallest interval. */

	if (level[maxerr] + 1 <= levmax) {
	    goto L160;
	}
	extrap = TRUE_;
	nrmax = 2;
L120:
	if (ierro == 3 || erlarg <= ertest) {
	    goto L140;
	}

/*           the smallest interval has the largest error. */
/*           before bisecting decrease the sum of the errors */
/*           over the larger intervals (erlarg) and perform */
/*           extrapolation. */

	id = nrmax;
	jupbnd = *last;
	if (*last > *limit / 2 + 2) {
	    jupbnd = *limit + 3 - *last;
	}
	i__1 = jupbnd;
	for (k = id; k <= i__1; ++k) {
	    maxerr = iord[nrmax];
	    errmax = elist[maxerr];
/* ***jump out of do-loop */
	    if (level[maxerr] + 1 <= levmax) {
		goto L160;
	    }
	    ++nrmax;
/* L130: */
	}

/*           perform extrapolation. */

L140:
	++numrl2;
	rlist2[numrl2 - 1] = area;
	if (numrl2 <= 2) {
	    goto L155;
	}
	qelg_(&numrl2, rlist2, &reseps, &abseps, res3la, &nres);
	++ktmin;
	if (ktmin > 5 && *abserr < errsum * .001f) {
	    *ier = 5;
	}
	if (abseps >= *abserr) {
	    goto L150;
	}
	ktmin = 0;
	*abserr = abseps;
	*result = reseps;
	correc = erlarg;
/* Computing MAX */
	r__1 = *epsabs, r__2 = *epsrel * dabs(reseps);
	ertest = dmax(r__1,r__2);
/* ***jump out of do-loop */
	if (*abserr < ertest) {
	    goto L170;
	}

/*           prepare bisection of the smallest interval. */

L150:
	if (numrl2 == 1) {
	    noext = TRUE_;
	}
	if (*ier >= 5) {
	    goto L170;
	}
L155:
	maxerr = iord[1];
	errmax = elist[maxerr];
	nrmax = 1;
	extrap = FALSE_;
	++levmax;
	erlarg = errsum;
L160:
	;
    }

/*           set the final result. */
/*           --------------------- */


L170:
    if (*abserr == oflow) {
	goto L190;
    }
    if (*ier + ierro == 0) {
	goto L180;
    }
    if (ierro == 3) {
	*abserr += correc;
    }
    if (*ier == 0) {
	*ier = 3;
    }
    if (*result != 0.f && area != 0.f) {
	goto L175;
    }
    if (*abserr > errsum) {
	goto L190;
    }
    if (area == 0.f) {
	goto L210;
    }
    goto L180;
L175:
    if (*abserr / dabs(*result) > errsum / dabs(area)) {
	goto L190;
    }

/*           test on divergence. */

L180:
/* Computing MAX */
    r__1 = dabs(*result), r__2 = dabs(area);
    if (ksgn == -1 && dmax(r__1,r__2) <= resabs * .01f) {
	goto L210;
    }
    if (.01f > *result / area || *result / area > 100.f || errsum > dabs(area)
	    ) {
	*ier = 6;
    }
    goto L210;

/*           compute global integral sum. */

L190:
    *result = 0.f;
    i__2 = *last;
    for (k = 1; k <= i__2; ++k) {
	*result += rlist[k];
/* L200: */
    }
    *abserr = errsum;
L210:
    if (*ier > 2) {
	--(*ier);
    }
    *result *= sign;
L999:
    return 0;
} /* qagpe_ */

