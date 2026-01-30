/* qk61.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int qk61_(E_fp f, real *a, real *b, real *result, real *
	abserr, real *resabs, real *resasc)
{
    /* Initialized data */

    static real xgk[31] = { .9994844100504906f,.9968934840746495f,
	    .9916309968704046f,.9836681232797472f,.9731163225011263f,
	    .9600218649683075f,.94437444474856f,.9262000474292743f,
	    .9055733076999078f,.8825605357920527f,.8572052335460611f,
	    .8295657623827684f,.7997278358218391f,.7677774321048262f,
	    .7337900624532268f,.6978504947933158f,.660061064126627f,
	    .6205261829892429f,.5793452358263617f,.5366241481420199f,
	    .4924804678617786f,.4470337695380892f,.4004012548303944f,
	    .3527047255308781f,.3040732022736251f,.2546369261678898f,
	    .2045251166823099f,.1538699136085835f,.102806937966737f,
	    .0514718425553177f,0.f };
    static real wgk[31] = { .001389013698677008f,.003890461127099884f,
	    .006630703915931292f,.009273279659517763f,.01182301525349634f,
	    .0143697295070458f,.01692088918905327f,.01941414119394238f,
	    .02182803582160919f,.0241911620780806f,.0265099548823331f,
	    .02875404876504129f,.03090725756238776f,.03298144705748373f,
	    .03497933802806002f,.03688236465182123f,.03867894562472759f,
	    .04037453895153596f,.04196981021516425f,.04345253970135607f,
	    .04481480013316266f,.04605923827100699f,.04718554656929915f,
	    .04818586175708713f,.04905543455502978f,.04979568342707421f,
	    .05040592140278235f,.05088179589874961f,.05122154784925877f,
	    .05142612853745903f,.05149472942945157f };
    static real wg[15] = { .007968192496166606f,.01846646831109096f,
	    .02878470788332337f,.03879919256962705f,.04840267283059405f,
	    .05749315621761907f,.0659742298821805f,.07375597473770521f,
	    .08075589522942022f,.08689978720108298f,.09212252223778613f,
	    .09636873717464426f,.09959342058679527f,.1017623897484055f,
	    .1028526528935588f };

    /* System generated locals */
    real r__1, r__2;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    integer j;
    real fc, fv1[30], fv2[30];
    integer jtw;
    real absc, resg, resk, fsum, fval1, fval2;
    integer jtwm1;
    real hlgth, centr, reskh, uflow;
    extern doublereal r1mach_(integer *);
    real epmach, dhlgth;

/* ***begin prologue  qk61 */
/* ***date written   800101   (yymmdd) */
/* ***revision date  830518   (yymmdd) */
/* ***category no.  h2a1a2 */
/* ***keywords  61-point gauss-kronrod rules */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
/* ***purpose  to compute i = integral of f over (a,b) with error */
/*                           estimate */
/*                       j = integral of dabs(f) over (a,b) */
/* ***description */

/*        integration rule */
/*        standard fortran subroutine */
/*        real version */


/*        parameters */
/*         on entry */
/*           f      - real */
/*                    function subprogram defining the integrand */
/*                    function f(x). the actual name for f needs to be */
/*                    declared e x t e r n a l in the calling program. */

/*           a      - real */
/*                    lower limit of integration */

/*           b      - real */
/*                    upper limit of integration */

/*         on return */
/*           result - real */
/*                    approximation to the integral i */
/*                    result is computed by applying the 61-point */
/*                    kronrod rule (resk) obtained by optimal addition of */
/*                    abscissae to the 30-point gauss rule (resg). */

/*           abserr - real */
/*                    estimate of the modulus of the absolute error, */
/*                    which should equal or exceed dabs(i-result) */

/*           resabs - real */
/*                    approximation to the integral j */

/*           resasc - real */
/*                    approximation to the integral of dabs(f-i/(b-a)) */


/* ***references  (none) */
/* ***routines called  r1mach */
/* ***end prologue  qk61 */



/*           the abscissae and weights are given for the */
/*           interval (-1,1). because of symmetry only the positive */
/*           abscissae and their corresponding weights are given. */

/*           xgk   - abscissae of the 61-point kronrod rule */
/*                   xgk(2), xgk(4)  ... abscissae of the 30-point */
/*                   gauss rule */
/*                   xgk(1), xgk(3)  ... optimally added abscissae */
/*                   to the 30-point gauss rule */

/*           wgk   - weights of the 61-point kronrod rule */

/*           wg    - weigths of the 30-point gauss rule */


/*           list of major variables */
/*           ----------------------- */

/*           centr  - mid point of the interval */
/*           hlgth  - half-length of the interval */
/*           absc   - abscissa */
/*           fval*  - function value */
/*           resg   - result of the 30-point gauss rule */
/*           resk   - result of the 61-point kronrod rule */
/*           reskh  - approximation to the mean value of f */
/*                    over (a,b), i.e. to i/(b-a) */

/*           machine dependent constants */
/*           --------------------------- */

/*           epmach is the largest relative spacing. */
/*           uflow is the smallest positive magnitude. */

/* ***first executable statement  qk61 */
    epmach = r1mach_(&c__4);
    uflow = r1mach_(&c__1);

    centr = (*b + *a) * .5f;
    hlgth = (*b - *a) * .5f;
    dhlgth = dabs(hlgth);

/*           compute the 61-point kronrod approximation to the */
/*           integral, and estimate the absolute error. */

    resg = 0.f;
    fc = (*f)(&centr);
    resk = wgk[30] * fc;
    *resabs = dabs(resk);
    for (j = 1; j <= 15; ++j) {
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
    for (j = 1; j <= 15; ++j) {
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
    *resasc = wgk[30] * (r__1 = fc - reskh, dabs(r__1));
    for (j = 1; j <= 30; ++j) {
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
} /* qk61_ */

