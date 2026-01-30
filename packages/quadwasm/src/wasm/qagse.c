/* qagse.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int qagse_(E_fp f, real *a, real *b, real *epsabs, real *
	epsrel, integer *limit, real *result, real *abserr, integer *neval, 
	integer *ier, real *alist__, real *blist, real *rlist, real *elist, 
	integer *iord, integer *last)
{
    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2;

    /* Local variables */
    integer k;
    real a1, a2, b1, b2;
    integer id;
    extern /* Subroutine */ int qk21_(E_fp, real *, real *, real *, real *, 
	    real *, real *);
    real area;
    extern /* Subroutine */ int qelg_(integer *, real *, real *, real *, real 
	    *, integer *);
    real dres;
    integer ksgn, nres;
    real area1, area2, area12, small, erro12;
    integer ierro;
    real defab1, defab2;
    integer ktmin, nrmax;
    real oflow, uflow;
    logical noext;
    extern /* Subroutine */ int qpsrt_(integer *, integer *, integer *, real *
	    , real *, integer *, integer *);
    extern doublereal r1mach_(integer *);
    integer iroff1, iroff2, iroff3;
    real res3la[3], error1, error2, rlist2[52];
    integer numrl2;
    real defabs, epmach, erlarg, abseps, correc, errbnd, resabs;
    integer jupbnd;
    real erlast, errmax;
    integer maxerr;
    real reseps;
    logical extrap;
    real ertest, errsum;

/* ***begin prologue  qagse */
/* ***date written   800101   (yymmdd) */
/* ***revision date  830518   (yymmdd) */
/* ***category no.  h2a1a1 */
/* ***keywords  automatic integrator, general-purpose, */
/*             (end point) singularities, extrapolation, */
/*             globally adaptive */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
/* ***purpose  the routine calculates an approximation result to a given */
/*            definite integral i = integral of f over (a,b), */
/*            hopefully satisfying following claim for accuracy */
/*            abs(i-result).le.max(epsabs,epsrel*abs(i)). */
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

/*            epsabs - real */
/*                     absolute accuracy requested */
/*            epsrel - real */
/*                     relative accuracy requested */
/*                     if  epsabs.le.0 */
/*                     and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
/*                     the routine will end with ier = 6. */

/*            limit  - integer */
/*                     gives an upperbound on the number of subintervals */
/*                     in the partition of (a,b) */

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
/*                     ier.gt.0 abnormal termination of the routine */
/*                             the estimates for integral and error are */
/*                             less reliable. it is assumed that the */
/*                             requested accuracy has not been achieved. */
/*            error messages */
/*                         = 1 maximum number of subdivisions allowed */
/*                             has been achieved. one can allow more sub- */
/*                             divisions by increasing the value of limit */
/*                             (and taking the according dimension */
/*                             adjustments into account). however, if */
/*                             this yields no improvement it is advised */
/*                             to analyze the integrand in order to */
/*                             determine the integration difficulties. if */
/*                             the position of a local difficulty can be */
/*                             determined (e.g. singularity, */
/*                             discontinuity within the interval) one */
/*                             will probably gain from splitting up the */
/*                             interval at this point and calling the */
/*                             integrator on the subranges. if possible, */
/*                             an appropriate special-purpose integrator */
/*                             should be used, which is designed for */
/*                             handling the type of difficulty involved. */
/*                         = 2 the occurrence of roundoff error is detec- */
/*                             ted, which prevents the requested */
/*                             tolerance from being achieved. */
/*                             the error may be under-estimated. */
/*                         = 3 extremely bad integrand behaviour */
/*                             occurs at some points of the integration */
/*                             interval. */
/*                         = 4 the algorithm does not converge. */
/*                             roundoff error is detected in the */
/*                             extrapolation table. */
/*                             it is presumed that the requested */
/*                             tolerance cannot be achieved, and that the */
/*                             returned result is the best which can be */
/*                             obtained. */
/*                         = 5 the integral is probably divergent, or */
/*                             slowly convergent. it must be noted that */
/*                             divergence can occur with any other value */
/*                             of ier. */
/*                         = 6 the input is invalid, because */
/*                             epsabs.le.0 and */
/*                             epsrel.lt.max(50*rel.mach.acc.,0.5d-28). */
/*                             result, abserr, neval, last, rlist(1), */
/*                             iord(1) and elist(1) are set to zero. */
/*                             alist(1) and blist(1) are set to a and b */
/*                             respectively. */

/*            alist  - real */
/*                     vector of dimension at least limit, the first */
/*                      last  elements of which are the left end points */
/*                     of the subintervals in the partition of the */
/*                     given integration range (a,b) */

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
/*                     subdivision process */

/* ***references  (none) */
/* ***routines called  qelg,qk21,qpsrt,r1mach */
/* ***end prologue  qagse */




/*            the dimension of rlist2 is determined by the value of */
/*            limexp in subroutine qelg (rlist2 should be of dimension */
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
/*                       containing the part of the epsilon table */
/*                       which is still needed for further computations */
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
/*           *****1    - variable for the left interval */
/*           *****2    - variable for the right interval */
/*           last      - index for subdivision */
/*           nres      - number of calls to the extrapolation routine */
/*           numrl2    - number of elements currently in rlist2. if an */
/*                       appropriate approximation to the compounded */
/*                       integral has been obtained it is put in */
/*                       rlist2(numrl2) after numrl2 has been increased */
/*                       by one. */
/*           small     - length of the smallest interval considered */
/*                       up to now, multiplied by 1.5 */
/*           erlarg    - sum of the errors over the intervals larger */
/*                       than the smallest interval considered up to now */
/*           extrap    - logical variable denoting that the routine */
/*                       is attempting to perform extrapolation */
/*                       i.e. before subdividing the smallest interval */
/*                       we try to decrease the value of erlarg. */
/*           noext     - logical variable denoting that extrapolation */
/*                       is no longer allowed (true value) */

/*            machine dependent constants */
/*            --------------------------- */

/*           epmach is the largest relative spacing. */
/*           uflow is the smallest positive magnitude. */
/*           oflow is the largest positive magnitude. */

/* ***first executable statement  qagse */
    /* Parameter adjustments */
    --iord;
    --elist;
    --rlist;
    --blist;
    --alist__;

    /* Function Body */
    epmach = r1mach_(&c__4);

/*            test on validity of parameters */
/*            ------------------------------ */
    *ier = 0;
    *neval = 0;
    *last = 0;
    *result = 0.f;
    *abserr = 0.f;
    alist__[1] = *a;
    blist[1] = *b;
    rlist[1] = 0.f;
    elist[1] = 0.f;
/* Computing MAX */
    r__1 = epmach * 50.f;
    if (*epsabs <= 0.f && *epsrel < dmax(r__1,5e-15f)) {
	*ier = 6;
    }
    if (*ier == 6) {
	goto L999;
    }

/*           first approximation to the integral */
/*           ----------------------------------- */

    uflow = r1mach_(&c__1);
    oflow = r1mach_(&c__2);
    ierro = 0;
    qk21_((E_fp)f, a, b, result, abserr, &defabs, &resabs);

/*           test on accuracy. */

    dres = dabs(*result);
/* Computing MAX */
    r__1 = *epsabs, r__2 = *epsrel * dres;
    errbnd = dmax(r__1,r__2);
    *last = 1;
    rlist[1] = *result;
    elist[1] = *abserr;
    iord[1] = 1;
    if (*abserr <= epmach * 100.f * defabs && *abserr > errbnd) {
	*ier = 2;
    }
    if (*limit == 1) {
	*ier = 1;
    }
    if (*ier != 0 || *abserr <= errbnd && *abserr != resabs || *abserr == 0.f)
	     {
	goto L140;
    }

/*           initialization */
/*           -------------- */

    rlist2[0] = *result;
    errmax = *abserr;
    maxerr = 1;
    area = *result;
    errsum = *abserr;
    *abserr = oflow;
    nrmax = 1;
    nres = 0;
    numrl2 = 2;
    ktmin = 0;
    extrap = FALSE_;
    noext = FALSE_;
    iroff1 = 0;
    iroff2 = 0;
    iroff3 = 0;
    ksgn = -1;
    if (dres >= (1.f - epmach * 50.f) * defabs) {
	ksgn = 1;
    }

/*           main do-loop */
/*           ------------ */

    i__1 = *limit;
    for (*last = 2; *last <= i__1; ++(*last)) {

/*           bisect the subinterval with the nrmax-th largest */
/*           error estimate. */

	a1 = alist__[maxerr];
	b1 = (alist__[maxerr] + blist[maxerr]) * .5f;
	a2 = b1;
	b2 = blist[maxerr];
	erlast = errmax;
	qk21_((E_fp)f, &a1, &b1, &area1, &error1, &resabs, &defab1);
	qk21_((E_fp)f, &a2, &b2, &area2, &error2, &resabs, &defab2);

/*           improve previous approximations to integral */
/*           and error and test for accuracy. */

	area12 = area1 + area2;
	erro12 = error1 + error2;
	errsum = errsum + erro12 - errmax;
	area = area + area12 - rlist[maxerr];
	if (defab1 == error1 || defab2 == error2) {
	    goto L15;
	}
	if ((r__1 = rlist[maxerr] - area12, dabs(r__1)) > dabs(area12) * 
		1e-5f || erro12 < errmax * .99f) {
	    goto L10;
	}
	if (extrap) {
	    ++iroff2;
	}
	if (! extrap) {
	    ++iroff1;
	}
L10:
	if (*last > 10 && erro12 > errmax) {
	    ++iroff3;
	}
L15:
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
/*           at a point of the integration range. */

/* Computing MAX */
	r__1 = dabs(a1), r__2 = dabs(b2);
	if (dmax(r__1,r__2) <= (epmach * 100.f + 1.f) * (dabs(a2) + uflow * 
		1e3f)) {
	    *ier = 4;
	}

/*           append the newly-created intervals to the list. */

	if (error2 > error1) {
	    goto L20;
	}
	alist__[*last] = a2;
	blist[maxerr] = b1;
	blist[*last] = b2;
	elist[maxerr] = error1;
	elist[*last] = error2;
	goto L30;
L20:
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

L30:
	qpsrt_(limit, last, &maxerr, &errmax, &elist[1], &iord[1], &nrmax);
/* ***jump out of do-loop */
	if (errsum <= errbnd) {
	    goto L115;
	}
/* ***jump out of do-loop */
	if (*ier != 0) {
	    goto L100;
	}
	if (*last == 2) {
	    goto L80;
	}
	if (noext) {
	    goto L90;
	}
	erlarg -= erlast;
	if ((r__1 = b1 - a1, dabs(r__1)) > small) {
	    erlarg += erro12;
	}
	if (extrap) {
	    goto L40;
	}

/*           test whether the interval to be bisected next is the */
/*           smallest interval. */

	if ((r__1 = blist[maxerr] - alist__[maxerr], dabs(r__1)) > small) {
	    goto L90;
	}
	extrap = TRUE_;
	nrmax = 2;
L40:
	if (ierro == 3 || erlarg <= ertest) {
	    goto L60;
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
	i__2 = jupbnd;
	for (k = id; k <= i__2; ++k) {
	    maxerr = iord[nrmax];
	    errmax = elist[maxerr];
/* ***jump out of do-loop */
	    if ((r__1 = blist[maxerr] - alist__[maxerr], dabs(r__1)) > small) 
		    {
		goto L90;
	    }
	    ++nrmax;
/* L50: */
	}

/*           perform extrapolation. */

L60:
	++numrl2;
	rlist2[numrl2 - 1] = area;
	qelg_(&numrl2, rlist2, &reseps, &abseps, res3la, &nres);
	++ktmin;
	if (ktmin > 5 && *abserr < errsum * .001f) {
	    *ier = 5;
	}
	if (abseps >= *abserr) {
	    goto L70;
	}
	ktmin = 0;
	*abserr = abseps;
	*result = reseps;
	correc = erlarg;
/* Computing MAX */
	r__1 = *epsabs, r__2 = *epsrel * dabs(reseps);
	ertest = dmax(r__1,r__2);
/* ***jump out of do-loop */
	if (*abserr <= ertest) {
	    goto L100;
	}

/*           prepare bisection of the smallest interval. */

L70:
	if (numrl2 == 1) {
	    noext = TRUE_;
	}
	if (*ier == 5) {
	    goto L100;
	}
	maxerr = iord[1];
	errmax = elist[maxerr];
	nrmax = 1;
	extrap = FALSE_;
	small *= .5f;
	erlarg = errsum;
	goto L90;
L80:
	small = (r__1 = *b - *a, dabs(r__1)) * .375f;
	erlarg = errsum;
	ertest = errbnd;
	rlist2[1] = area;
L90:
	;
    }

/*           set final result and error estimate. */
/*           ------------------------------------ */

L100:
    if (*abserr == oflow) {
	goto L115;
    }
    if (*ier + ierro == 0) {
	goto L110;
    }
    if (ierro == 3) {
	*abserr += correc;
    }
    if (*ier == 0) {
	*ier = 3;
    }
    if (*result != 0.f && area != 0.f) {
	goto L105;
    }
    if (*abserr > errsum) {
	goto L115;
    }
    if (area == 0.f) {
	goto L130;
    }
    goto L110;
L105:
    if (*abserr / dabs(*result) > errsum / dabs(area)) {
	goto L115;
    }

/*           test on divergence. */

L110:
/* Computing MAX */
    r__1 = dabs(*result), r__2 = dabs(area);
    if (ksgn == -1 && dmax(r__1,r__2) <= defabs * .01f) {
	goto L130;
    }
    if (.01f > *result / area || *result / area > 100.f || errsum > dabs(area)
	    ) {
	*ier = 6;
    }
    goto L130;

/*           compute global integral sum. */

L115:
    *result = 0.f;
    i__1 = *last;
    for (k = 1; k <= i__1; ++k) {
	*result += rlist[k];
/* L120: */
    }
    *abserr = errsum;
L130:
    if (*ier > 2) {
	--(*ier);
    }
L140:
    *neval = *last * 42 - 21;
L999:
    return 0;
} /* qagse_ */

