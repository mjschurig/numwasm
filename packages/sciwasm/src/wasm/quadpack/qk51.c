/* qk51.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int qk51_(E_fp f, real *a, real *b, real *result, real *
	abserr, real *resabs, real *resasc)
{
    /* Initialized data */

    static real xgk[26] = { .9992621049926098f,.9955569697904981f,
	    .9880357945340772f,.9766639214595175f,.9616149864258425f,
	    .9429745712289743f,.9207471152817016f,.8949919978782754f,
	    .8658470652932756f,.833442628760834f,.7978737979985001f,
	    .7592592630373576f,.7177664068130844f,.6735663684734684f,
	    .6268100990103174f,.577662930241223f,.5263252843347192f,
	    .473002731445715f,.4178853821930377f,.3611723058093878f,
	    .3030895389311078f,.2438668837209884f,.1837189394210489f,
	    .1228646926107104f,.06154448300568508f,0.f };
    static real wgk[26] = { .001987383892330316f,.005561932135356714f,
	    .009473973386174152f,.01323622919557167f,.0168478177091283f,
	    .02043537114588284f,.02400994560695322f,.02747531758785174f,
	    .03079230016738749f,.03400213027432934f,.03711627148341554f,
	    .04008382550403238f,.04287284502017005f,.04550291304992179f,
	    .04798253713883671f,.05027767908071567f,.05236288580640748f,
	    .05425112988854549f,.05595081122041232f,.05743711636156783f,
	    .05868968002239421f,.05972034032417406f,.06053945537604586f,
	    .06112850971705305f,.06147118987142532f,.06158081806783294f };
    static real wg[13] = { .01139379850102629f,.02635498661503214f,
	    .04093915670130631f,.05490469597583519f,.06803833381235692f,
	    .08014070033500102f,.09102826198296365f,.1005359490670506f,
	    .1085196244742637f,.1148582591457116f,.1194557635357848f,
	    .12224244299031f,.1231760537267155f };

    /* System generated locals */
    real r__1, r__2;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    integer j;
    real fc, fv1[25], fv2[25];
    integer jtw;
    real absc, resg, resk, fsum, fval1, fval2;
    integer jtwm1;
    real hlgth, centr, reskh, uflow;
    extern doublereal r1mach_(integer *);
    real epmach, dhlgth;

/* ***begin prologue  qk51 */
/* ***date written   800101   (yymmdd) */
/* ***revision date  830518   (yymmdd) */
/* ***category no.  h2a1a2 */
/* ***keywords  51-point gauss-kronrod rules */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math & progr. div. - k.u.leuven */
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
/*                       function subroutine defining the integrand */
/*                       function f(x). the actual name for f needs to be */
/*                       declared e x t e r n a l in the calling program. */

/*              a      - real */
/*                       lower limit of integration */

/*              b      - real */
/*                       upper limit of integration */

/*            on return */
/*              result - real */
/*                       approximation to the integral i */
/*                       result is computed by applying the 51-point */
/*                       kronrod rule (resk) obtained by optimal addition */
/*                       of abscissae to the 25-point gauss rule (resg). */

/*              abserr - real */
/*                       estimate of the modulus of the absolute error, */
/*                       which should not exceed abs(i-result) */

/*              resabs - real */
/*                       approximation to the integral j */

/*              resasc - real */
/*                       approximation to the integral of abs(f-i/(b-a)) */
/*                       over (a,b) */

/* ***references  (none) */
/* ***routines called  r1mach */
/* ***end prologue  qk51 */



/*           the abscissae and weights are given for the interval (-1,1). */
/*           because of symmetry only the positive abscissae and their */
/*           corresponding weights are given. */

/*           xgk    - abscissae of the 51-point kronrod rule */
/*                    xgk(2), xgk(4), ...  abscissae of the 25-point */
/*                    gauss rule */
/*                    xgk(1), xgk(3), ...  abscissae which are optimally */
/*                    added to the 25-point gauss rule */

/*           wgk    - weights of the 51-point kronrod rule */

/*           wg     - weights of the 25-point gauss rule */



/*           list of major variables */
/*           ----------------------- */

/*           centr  - mid point of the interval */
/*           hlgth  - half-length of the interval */
/*           absc   - abscissa */
/*           fval*  - function value */
/*           resg   - result of the 25-point gauss formula */
/*           resk   - result of the 51-point kronrod formula */
/*           reskh  - approximation to the mean value of f over (a,b), */
/*                    i.e. to i/(b-a) */

/*           machine dependent constants */
/*           --------------------------- */

/*           epmach is the largest relative spacing. */
/*           uflow is the smallest positive magnitude. */

/* ***first executable statement  qk51 */
    epmach = r1mach_(&c__4);
    uflow = r1mach_(&c__1);

    centr = (*a + *b) * .5f;
    hlgth = (*b - *a) * .5f;
    dhlgth = dabs(hlgth);

/*           compute the 51-point kronrod approximation to */
/*           the integral, and estimate the absolute error. */

    fc = (*f)(&centr);
    resg = wg[12] * fc;
    resk = wgk[25] * fc;
    *resabs = dabs(resk);
    for (j = 1; j <= 12; ++j) {
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
    for (j = 1; j <= 13; ++j) {
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
    *resasc = wgk[25] * (r__1 = fc - reskh, dabs(r__1));
    for (j = 1; j <= 25; ++j) {
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
} /* qk51_ */

