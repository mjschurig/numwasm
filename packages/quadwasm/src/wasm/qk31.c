/* qk31.f -- translated by f2c (version 20240504).
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
static doublereal c_b7 = 1.5;

/* Subroutine */ int qk31_(E_fp f, real *a, real *b, real *result, real *
	abserr, real *resabs, real *resasc)
{
    /* Initialized data */

    static real xgk[16] = { .9980022986933971f,.9879925180204854f,
	    .9677390756791391f,.9372733924007059f,.8972645323440819f,
	    .8482065834104272f,.7904185014424659f,.72441773136017f,
	    .650996741297417f,.5709721726085388f,.4850818636402397f,
	    .3941513470775634f,.2991800071531688f,.2011940939974345f,
	    .1011420669187175f,0.f };
    static real wgk[16] = { .005377479872923349f,.01500794732931612f,
	    .02546084732671532f,.03534636079137585f,.04458975132476488f,
	    .05348152469092809f,.06200956780067064f,.06985412131872826f,
	    .07684968075772038f,.08308050282313302f,.08856444305621177f,
	    .09312659817082532f,.09664272698362368f,.09917359872179196f,
	    .1007698455238756f,.1013300070147915f };
    static real wg[8] = { .03075324199611727f,.07036604748810812f,
	    .1071592204671719f,.1395706779261543f,.1662692058169939f,
	    .1861610000155622f,.1984314853271116f,.2025782419255613f };

    /* System generated locals */
    real r__1, r__2;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    integer j;
    real fc, fv1[15], fv2[15];
    integer jtw;
    real absc, resg, resk, fsum, fval1, fval2;
    integer jtwm1;
    real hlgth, centr, reskh, uflow;
    extern doublereal r1mach_(integer *);
    real epmach, dhlgth;

/* ***begin prologue  qk31 */
/* ***date written   800101   (yymmdd) */
/* ***revision date  830518   (yymmdd) */
/* ***category no.  h2a1a2 */
/* ***keywords  31-point gauss-kronrod rules */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
/* ***purpose  to compute i = integral of f over (a,b) with error */
/*                           estimate */
/*                       j = integral of abs(f) over (a,b) */
/* ***description */

/*           integration rules */
/*           standard fortran subroutine */
/*           real version */

/*           parameters */
/*            on entry */
/*              f      - real */
/*                       function subprogram defining the integrand */
/*                       function f(x). the actual name for f needs to be */
/*                       declared e x t e r n a l in the calling program. */

/*              a      - real */
/*                       lower limit of integration */

/*              b      - real */
/*                       upper limit of integration */

/*            on return */
/*              result - real */
/*                       approximation to the integral i */
/*                       result is computed by applying the 31-point */
/*                       gauss-kronrod rule (resk), obtained by optimal */
/*                       addition of abscissae to the 15-point gauss */
/*                       rule (resg). */

/*              abserr - real */
/*                       estimate of the modulus of the modulus, */
/*                       which should not exceed abs(i-result) */

/*              resabs - real */
/*                       approximation to the integral j */

/*              resasc - real */
/*                       approximation to the integral of abs(f-i/(b-a)) */
/*                       over (a,b) */

/* ***references  (none) */
/* ***routines called  r1mach */
/* ***end prologue  qk31 */


/*           the abscissae and weights are given for the interval (-1,1). */
/*           because of symmetry only the positive abscissae and their */
/*           corresponding weights are given. */

/*           xgk    - abscissae of the 31-point kronrod rule */
/*                    xgk(2), xgk(4), ...  abscissae of the 15-point */
/*                    gauss rule */
/*                    xgk(1), xgk(3), ...  abscissae which are optimally */
/*                    added to the 15-point gauss rule */

/*           wgk    - weights of the 31-point kronrod rule */

/*           wg     - weights of the 15-point gauss rule */



/*           list of major variables */
/*           ----------------------- */
/*           centr  - mid point of the interval */
/*           hlgth  - half-length of the interval */
/*           absc   - abscissa */
/*           fval*  - function value */
/*           resg   - result of the 15-point gauss formula */
/*           resk   - result of the 31-point kronrod formula */
/*           reskh  - approximation to the mean value of f over (a,b), */
/*                    i.e. to i/(b-a) */

/*           machine dependent constants */
/*           --------------------------- */
/*           epmach is the largest relative spacing. */
/*           uflow is the smallest positive magnitude. */

/* ***first executable statement  qk31 */
    epmach = r1mach_(&c__4);
    uflow = r1mach_(&c__1);

    centr = (*a + *b) * .5f;
    hlgth = (*b - *a) * .5f;
    dhlgth = dabs(hlgth);

/*           compute the 31-point kronrod approximation to */
/*           the integral, and estimate the absolute error. */

    fc = (*f)(&centr);
    resg = wg[7] * fc;
    resk = wgk[15] * fc;
    *resabs = dabs(resk);
    for (j = 1; j <= 7; ++j) {
	jtw = j << 1;
	absc = hlgth * xgk[jtw - 1];
	r__1 = centr - absc;
	fval1 = (*f)(&r__1);
	r__1 = centr + absc;
	fval2 = (*f)(&r__1);
	fv1[jtw - 1] = fval1;
	fv2[jtw - 1] = fval2;
	fsum = fval1 + fval2;
	resg += wg[j - 1] * fsum;
	resk += wgk[jtw - 1] * fsum;
	*resabs += wgk[jtw - 1] * (dabs(fval1) + dabs(fval2));
/* L10: */
    }
    for (j = 1; j <= 8; ++j) {
	jtwm1 = (j << 1) - 1;
	absc = hlgth * xgk[jtwm1 - 1];
	r__1 = centr - absc;
	fval1 = (*f)(&r__1);
	r__1 = centr + absc;
	fval2 = (*f)(&r__1);
	fv1[jtwm1 - 1] = fval1;
	fv2[jtwm1 - 1] = fval2;
	fsum = fval1 + fval2;
	resk += wgk[jtwm1 - 1] * fsum;
	*resabs += wgk[jtwm1 - 1] * (dabs(fval1) + dabs(fval2));
/* L15: */
    }
    reskh = resk * .5f;
    *resasc = wgk[15] * (r__1 = fc - reskh, dabs(r__1));
    for (j = 1; j <= 15; ++j) {
	*resasc += wgk[j - 1] * ((r__1 = fv1[j - 1] - reskh, dabs(r__1)) + (
		r__2 = fv2[j - 1] - reskh, dabs(r__2)));
/* L20: */
    }
    *result = resk * hlgth;
    *resabs *= dhlgth;
    *resasc *= dhlgth;
    *abserr = (r__1 = (resk - resg) * hlgth, dabs(r__1));
    if (*resasc != 0.f && *abserr != 0.f) {
/* Computing MIN */
	d__1 = (doublereal) (*abserr * 200.f / *resasc);
	r__1 = 1.f, r__2 = pow_dd(&d__1, &c_b7);
	*abserr = *resasc * dmin(r__1,r__2);
    }
    if (*resabs > uflow / (epmach * 50.f)) {
/* Computing MAX */
	r__1 = epmach * 50.f * *resabs;
	*abserr = dmax(r__1,*abserr);
    }
    return 0;
} /* qk31_ */

