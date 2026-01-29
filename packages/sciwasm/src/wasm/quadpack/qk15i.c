/* qk15i.f -- translated by f2c (version 20240504).
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
static doublereal c_b6 = 1.5;

/* Subroutine */ int qk15i_(E_fp f, real *boun, integer *inf, real *a, real *
	b, real *result, real *abserr, real *resabs, real *resasc)
{
    /* Initialized data */

    static real xgk[8] = { .9914553711208126f,.9491079123427585f,
	    .8648644233597691f,.7415311855993944f,.5860872354676911f,
	    .4058451513773972f,.2077849550078985f,0.f };
    static real wgk[8] = { .02293532201052922f,.06309209262997855f,
	    .1047900103222502f,.1406532597155259f,.1690047266392679f,
	    .1903505780647854f,.2044329400752989f,.2094821410847278f };
    static real wg[8] = { 0.f,.1294849661688697f,0.f,.2797053914892767f,0.f,
	    .3818300505051189f,0.f,.4179591836734694f };

    /* System generated locals */
    real r__1, r__2;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    integer j;
    real fc, fv1[7], fv2[7], absc, dinf, resg, resk, fsum, absc1, absc2, 
	    fval1, fval2, hlgth, centr, reskh, uflow;
    extern doublereal r1mach_(integer *);
    real tabsc1, tabsc2, epmach;

/* ***begin prologue  qk15i */
/* ***date written   800101   (yymmdd) */
/* ***revision date  830518   (yymmdd) */
/* ***category no.  h2a3a2,h2a4a2 */
/* ***keywords  15-point transformed gauss-kronrod rules */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
/* ***purpose  the original (infinite integration range is mapped */
/*            onto the interval (0,1) and (a,b) is a part of (0,1). */
/*            it is the purpose to compute */
/*            i = integral of transformed integrand over (a,b), */
/*            j = integral of abs(transformed integrand) over (a,b). */
/* ***description */

/*           integration rule */
/*           standard fortran subroutine */
/*           real version */

/*           parameters */
/*            on entry */
/*              f      - real */
/*                       fuction subprogram defining the integrand */
/*                       function f(x). the actual name for f needs to be */
/*                       declared e x t e r n a l in the calling program. */

/*              boun   - real */
/*                       finite bound of original integration */
/*                       range (set to zero if inf = +2) */

/*              inf    - integer */
/*                       if inf = -1, the original interval is */
/*                                   (-infinity,bound), */
/*                       if inf = +1, the original interval is */
/*                                   (bound,+infinity), */
/*                       if inf = +2, the original interval is */
/*                                   (-infinity,+infinity) and */
/*                       the integral is computed as the sum of two */
/*                       integrals, one over (-infinity,0) and one over */
/*                       (0,+infinity). */

/*              a      - real */
/*                       lower limit for integration over subrange */
/*                       of (0,1) */

/*              b      - real */
/*                       upper limit for integration over subrange */
/*                       of (0,1) */

/*            on return */
/*              result - real */
/*                       approximation to the integral i */
/*                       result is computed by applying the 15-point */
/*                       kronrod rule(resk) obtained by optimal addition */
/*                       of abscissae to the 7-point gauss rule(resg). */

/*              abserr - real */
/*                       estimate of the modulus of the absolute error, */
/*                       which should equal or exceed abs(i-result) */

/*              resabs - real */
/*                       approximation to the integral j */

/*              resasc - real */
/*                       approximation to the integral of */
/*                       abs((transformed integrand)-i/(b-a)) over (a,b) */

/* ***references  (none) */
/* ***routines called  r1mach */
/* ***end prologue  qk15i */



/*           the abscissae and weights are supplied for the interval */
/*           (-1,1).  because of symmetry only the positive abscissae and */
/*           their corresponding weights are given. */

/*           xgk    - abscissae of the 15-point kronrod rule */
/*                    xgk(2), xgk(4), ... abscissae of the 7-point */
/*                    gauss rule */
/*                    xgk(1), xgk(3), ...  abscissae which are optimally */
/*                    added to the 7-point gauss rule */

/*           wgk    - weights of the 15-point kronrod rule */

/*           wg     - weights of the 7-point gauss rule, corresponding */
/*                    to the abscissae xgk(2), xgk(4), ... */
/*                    wg(1), wg(3), ... are set to zero. */





/*           list of major variables */
/*           ----------------------- */

/*           centr  - mid point of the interval */
/*           hlgth  - half-length of the interval */
/*           absc*  - abscissa */
/*           tabsc* - transformed abscissa */
/*           fval*  - function value */
/*           resg   - result of the 7-point gauss formula */
/*           resk   - result of the 15-point kronrod formula */
/*           reskh  - approximation to the mean value of the transformed */
/*                    integrand over (a,b), i.e. to i/(b-a) */

/*           machine dependent constants */
/*           --------------------------- */

/*           epmach is the largest relative spacing. */
/*           uflow is the smallest positive magnitude. */

/* ***first executable statement  qk15i */
    epmach = r1mach_(&c__4);
    uflow = r1mach_(&c__1);
    dinf = (real) min(1,*inf);

    centr = (*a + *b) * .5f;
    hlgth = (*b - *a) * .5f;
    tabsc1 = *boun + dinf * (1.f - centr) / centr;
    fval1 = (*f)(&tabsc1);
    if (*inf == 2) {
	r__1 = -tabsc1;
	fval1 += (*f)(&r__1);
    }
    fc = fval1 / centr / centr;

/*           compute the 15-point kronrod approximation to */
/*           the integral, and estimate the error. */

    resg = wg[7] * fc;
    resk = wgk[7] * fc;
    *resabs = dabs(resk);
    for (j = 1; j <= 7; ++j) {
	absc = hlgth * xgk[j - 1];
	absc1 = centr - absc;
	absc2 = centr + absc;
	tabsc1 = *boun + dinf * (1.f - absc1) / absc1;
	tabsc2 = *boun + dinf * (1.f - absc2) / absc2;
	fval1 = (*f)(&tabsc1);
	fval2 = (*f)(&tabsc2);
	if (*inf == 2) {
	    r__1 = -tabsc1;
	    fval1 += (*f)(&r__1);
	}
	if (*inf == 2) {
	    r__1 = -tabsc2;
	    fval2 += (*f)(&r__1);
	}
	fval1 = fval1 / absc1 / absc1;
	fval2 = fval2 / absc2 / absc2;
	fv1[j - 1] = fval1;
	fv2[j - 1] = fval2;
	fsum = fval1 + fval2;
	resg += wg[j - 1] * fsum;
	resk += wgk[j - 1] * fsum;
	*resabs += wgk[j - 1] * (dabs(fval1) + dabs(fval2));
/* L10: */
    }
    reskh = resk * .5f;
    *resasc = wgk[7] * (r__1 = fc - reskh, dabs(r__1));
    for (j = 1; j <= 7; ++j) {
	*resasc += wgk[j - 1] * ((r__1 = fv1[j - 1] - reskh, dabs(r__1)) + (
		r__2 = fv2[j - 1] - reskh, dabs(r__2)));
/* L20: */
    }
    *result = resk * hlgth;
    *resasc *= hlgth;
    *resabs *= hlgth;
    *abserr = (r__1 = (resk - resg) * hlgth, dabs(r__1));
    if (*resasc != 0.f && *abserr != 0.f) {
/* Computing MIN */
	d__1 = (doublereal) (*abserr * 200.f / *resasc);
	r__1 = 1.f, r__2 = pow_dd(&d__1, &c_b6);
	*abserr = *resasc * dmin(r__1,r__2);
    }
    if (*resabs > uflow / (epmach * 50.f)) {
/* Computing MAX */
	r__1 = epmach * 50.f * *resabs;
	*abserr = dmax(r__1,*abserr);
    }
    return 0;
} /* qk15i_ */

