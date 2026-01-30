/* qk41.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int qk41_(E_fp f, real *a, real *b, real *result, real *
	abserr, real *resabs, real *resasc)
{
    /* Initialized data */

    static real xgk[21] = { .9988590315882777f,.9931285991850949f,
	    .9815078774502503f,.9639719272779138f,.9408226338317548f,
	    .9122344282513259f,.878276811252282f,.8391169718222188f,
	    .7950414288375512f,.7463319064601508f,.6932376563347514f,
	    .636053680726515f,.5751404468197103f,.5108670019508271f,
	    .4435931752387251f,.3737060887154196f,.301627868114913f,
	    .2277858511416451f,.1526054652409227f,.07652652113349733f,0.f };
    static real wgk[21] = { .003073583718520532f,.008600269855642942f,
	    .01462616925697125f,.02038837346126652f,.02588213360495116f,
	    .0312873067770328f,.0366001697582008f,.04166887332797369f,
	    .04643482186749767f,.05094457392372869f,.05519510534828599f,
	    .05911140088063957f,.06265323755478117f,.06583459713361842f,
	    .06864867292852162f,.07105442355344407f,.07303069033278667f,
	    .07458287540049919f,.07570449768455667f,.07637786767208074f,
	    .07660071191799966f };
    static real wg[10] = { .01761400713915212f,.04060142980038694f,
	    .06267204833410906f,.08327674157670475f,.1019301198172404f,
	    .1181945319615184f,.1316886384491766f,.1420961093183821f,
	    .1491729864726037f,.1527533871307259f };

    /* System generated locals */
    real r__1, r__2;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    integer j;
    real fc, fv1[20], fv2[20];
    integer jtw;
    real absc, resg, resk, fsum, fval1, fval2;
    integer jtwm1;
    real hlgth, centr, reskh, uflow;
    extern doublereal r1mach_(integer *);
    real epmach, dhlgth;

/* ***begin prologue  qk41 */
/* ***date written   800101   (yymmdd) */
/* ***revision date  830518   (yymmdd) */
/* ***category no.  h2a1a2 */
/* ***keywords  41-point gauss-kronrod rules */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
/* ***purpose  to compute i = integral of f over (a,b), with error */
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
/*                       result is computed by applying the 41-point */
/*                       gauss-kronrod rule (resk) obtained by optimal */
/*                       addition of abscissae to the 20-point gauss */
/*                       rule (resg). */

/*              abserr - real */
/*                       estimate of the modulus of the absolute error, */
/*                       which should not exceed abs(i-result) */

/*              resabs - real */
/*                       approximation to the integral j */

/*              resasc - real */
/*                       approximation to the integal of abs(f-i/(b-a)) */
/*                       over (a,b) */

/* ***references  (none) */
/* ***routines called  r1mach */
/* ***end prologue  qk41 */



/*           the abscissae and weights are given for the interval (-1,1). */
/*           because of symmetry only the positive abscissae and their */
/*           corresponding weights are given. */

/*           xgk    - abscissae of the 41-point gauss-kronrod rule */
/*                    xgk(2), xgk(4), ...  abscissae of the 20-point */
/*                    gauss rule */
/*                    xgk(1), xgk(3), ...  abscissae which are optimally */
/*                    added to the 20-point gauss rule */

/*           wgk    - weights of the 41-point gauss-kronrod rule */

/*           wg     - weights of the 20-point gauss rule */



/*           list of major variables */
/*           ----------------------- */

/*           centr  - mid point of the interval */
/*           hlgth  - half-length of the interval */
/*           absc   - abscissa */
/*           fval*  - function value */
/*           resg   - result of the 20-point gauss formula */
/*           resk   - result of the 41-point kronrod formula */
/*           reskh  - approximation to mean value of f over (a,b), i.e. */
/*                    to i/(b-a) */

/*           machine dependent constants */
/*           --------------------------- */

/*           epmach is the largest relative spacing. */
/*           uflow is the smallest positive magnitude. */

/* ***first executable statement  qk41 */
    epmach = r1mach_(&c__4);
    uflow = r1mach_(&c__1);

    centr = (*a + *b) * .5f;
    hlgth = (*b - *a) * .5f;
    dhlgth = dabs(hlgth);

/*           compute the 41-point gauss-kronrod approximation to */
/*           the integral, and estimate the absolute error. */

    resg = 0.f;
    fc = (*f)(&centr);
    resk = wgk[20] * fc;
    *resabs = dabs(resk);
    for (j = 1; j <= 10; ++j) {
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
    for (j = 1; j <= 10; ++j) {
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
    *resasc = wgk[20] * (r__1 = fc - reskh, dabs(r__1));
    for (j = 1; j <= 20; ++j) {
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
} /* qk41_ */

