/* qng.f -- translated by f2c (version 20240504).
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
static doublereal c_b17 = 1.5;
static integer c__26 = 26;
static integer c__0 = 0;

/* Subroutine */ int qng_(E_fp f, real *a, real *b, real *epsabs, real *
	epsrel, real *result, real *abserr, integer *neval, integer *ier)
{
    /* Initialized data */

    static real x1[5] = { .9739065285171717f,.8650633666889845f,
	    .6794095682990244f,.4333953941292472f,.1488743389816312f };
    static real w87a[21] = { .008148377384149173f,.01876143820156282f,
	    .02734745105005229f,.03367770731163793f,.03693509982042791f,
	    .002884872430211531f,.0136859460227127f,.02328041350288831f,
	    .03087249761171336f,.03569363363941877f,9.152833452022414e-4f,
	    .005399280219300471f,.01094767960111893f,.01629873169678734f,
	    .02108156888920384f,.02537096976925383f,.02918969775647575f,
	    .03237320246720279f,.03478309895036514f,.03641222073135179f,
	    .03725387550304771f };
    static real w87b[23] = { 2.741455637620724e-4f,.001807124155057943f,
	    .004096869282759165f,.006758290051847379f,.009549957672201647f,
	    .01232944765224485f,.01501044734638895f,.01754896798624319f,
	    .01993803778644089f,.02219493596101229f,.02433914712600081f,
	    .02637450541483921f,.0282869107887712f,.0300525811280927f,
	    .03164675137143993f,.0330504134199785f,.03425509970422606f,
	    .03526241266015668f,.0360769896228887f,.03669860449845609f,
	    .03712054926983258f,.03733422875193504f,.03736107376267902f };
    static real x2[5] = { .9956571630258081f,.9301574913557082f,
	    .7808177265864169f,.5627571346686047f,.2943928627014602f };
    static real x3[11] = { .9993333609019321f,.9874334029080889f,
	    .9548079348142663f,.9001486957483283f,.8251983149831142f,
	    .732148388989305f,.6228479705377252f,.4994795740710565f,
	    .3649016613465808f,.2222549197766013f,.07465061746138332f };
    static real x4[22] = { .9999029772627292f,.9979898959866787f,
	    .9921754978606872f,.9813581635727128f,.9650576238583846f,
	    .9431676131336706f,.9158064146855072f,.8832216577713165f,
	    .8457107484624157f,.803557658035231f,.7570057306854956f,
	    .7062732097873218f,.6515894665011779f,.5932233740579611f,
	    .5314936059708319f,.4667636230420228f,.3994248478592188f,
	    .3298748771061883f,.2585035592021616f,.1856953965683467f,
	    .1118422131799075f,.03735212339461987f };
    static real w10[5] = { .06667134430868814f,.1494513491505806f,
	    .219086362515982f,.2692667193099964f,.2955242247147529f };
    static real w21a[5] = { .03255816230796473f,.07503967481091995f,
	    .1093871588022976f,.1347092173114733f,.1477391049013385f };
    static real w21b[6] = { .01169463886737187f,.054755896574352f,
	    .09312545458369761f,.1234919762620659f,.1427759385770601f,
	    .1494455540029169f };
    static real w43a[10] = { .01629673428966656f,.0375228761208695f,
	    .05469490205825544f,.06735541460947809f,.07387019963239395f,
	    .005768556059769796f,.02737189059324884f,.04656082691042883f,
	    .06174499520144256f,.0713872672686934f };
    static real w43b[12] = { .001844477640212414f,.01079868958589165f,
	    .02189536386779543f,.03259746397534569f,.04216313793519181f,
	    .05074193960018458f,.05837939554261925f,.06474640495144589f,
	    .06956619791235648f,.07282444147183321f,.07450775101417512f,
	    .07472214751740301f };

    /* System generated locals */
    real r__1, r__2, r__3, r__4;
    doublereal d__1;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *);

    /* Local variables */
    integer k, l;
    real fv1[5], fv2[5], fv3[5], fv4[5];
    integer ipx;
    real absc, fval, res10, res21, res43, res87, fval1, fval2, hlgth, centr, 
	    reskh, uflow;
    extern doublereal r1mach_(integer *);
    real epmach, dhlgth, resabs, resasc, fcentr, savfun[21];
    extern /* Subroutine */ int xerror_(char *, integer *, integer *, integer 
	    *, ftnlen);

/* ***begin prologue  qng */
/* ***date written   800101   (yymmdd) */
/* ***revision date  830518   (yymmdd) */
/* ***category no.  h2a1a1 */
/* ***keywords  automatic integrator, smooth integrand, */
/*             non-adaptive, gauss-kronrod(patterson) */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl math & progr. div. - k.u.leuven */
/*           kahaner,david,nbs - modified (2/82) */
/* ***purpose  the routine calculates an approximation result to a */
/*            given definite integral i = integral of f over (a,b), */
/*            hopefully satisfying following claim for accuracy */
/*            abs(i-result).le.max(epsabs,epsrel*abs(i)). */
/* ***description */

/* non-adaptive integration */
/* standard fortran subroutine */
/* real version */

/*           f      - real version */
/*                    function subprogram defining the integrand function */
/*                    f(x). the actual name for f needs to be declared */
/*                    e x t e r n a l in the driver program. */

/*           a      - real version */
/*                    lower limit of integration */

/*           b      - real version */
/*                    upper limit of integration */

/*           epsabs - real */
/*                    absolute accuracy requested */
/*           epsrel - real */
/*                    relative accuracy requested */
/*                    if  epsabs.le.0 */
/*                    and epsrel.lt.max(50*rel.mach.acc.,0.5d-28), */
/*                    the routine will end with ier = 6. */

/*         on return */
/*           result - real */
/*                    approximation to the integral i */
/*                    result is obtained by applying the 21-point */
/*                    gauss-kronrod rule (res21) obtained by optimal */
/*                    addition of abscissae to the 10-point gauss rule */
/*                    (res10), or by applying the 43-point rule (res43) */
/*                    obtained by optimal addition of abscissae to the */
/*                    21-point gauss-kronrod rule, or by applying the */
/*                    87-point rule (res87) obtained by optimal addition */
/*                    of abscissae to the 43-point rule. */

/*           abserr - real */
/*                    estimate of the modulus of the absolute error, */
/*                    which should equal or exceed abs(i-result) */

/*           neval  - integer */
/*                    number of integrand evaluations */

/*           ier    - ier = 0 normal and reliable termination of the */
/*                            routine. it is assumed that the requested */
/*                            accuracy has been achieved. */
/*                    ier.gt.0 abnormal termination of the routine. it is */
/*                            assumed that the requested accuracy has */
/*                            not been achieved. */
/*           error messages */
/*                    ier = 1 the maximum number of steps has been */
/*                            executed. the integral is probably too */
/*                            difficult to be calculated by dqng. */
/*                        = 6 the input is invalid, because */
/*                            epsabs.le.0 and */
/*                            epsrel.lt.max(50*rel.mach.acc.,0.5d-28). */
/*                            result, abserr and neval are set to zero. */

/* ***references  (none) */
/* ***routines called  r1mach,xerror */
/* ***end prologue  qng */



/*           the following data statements contain the */
/*           abscissae and weights of the integration rules used. */

/*           x1      abscissae common to the 10-, 21-, 43- */
/*                   and 87-point rule */
/*           x2      abscissae common to the 21-, 43- and */
/*                   87-point rule */
/*           x3      abscissae common to the 43- and 87-point */
/*                   rule */
/*           x4      abscissae of the 87-point rule */
/*           w10     weights of the 10-point formula */
/*           w21a    weights of the 21-point formula for */
/*                   abscissae x1 */
/*           w21b    weights of the 21-point formula for */
/*                   abscissae x2 */
/*           w43a    weights of the 43-point formula for */
/*                   abscissae x1, x3 */
/*           w43b    weights of the 43-point formula for */
/*                   abscissae x3 */
/*           w87a    weights of the 87-point formula for */
/*                   abscissae x1, x2, x3 */
/*           w87b    weights of the 87-point formula for */
/*                   abscissae x4 */


/*           list of major variables */
/*           ----------------------- */

/*           centr  - mid point of the integration interval */
/*           hlgth  - half-length of the integration interval */
/*           fcentr - function value at mid point */
/*           absc   - abscissa */
/*           fval   - function value */
/*           savfun - array of function values which */
/*                    have already been computed */
/*           res10  - 10-point gauss result */
/*           res21  - 21-point kronrod result */
/*           res43  - 43-point result */
/*           res87  - 87-point result */
/*           resabs - approximation to the integral of abs(f) */
/*           resasc - approximation to the integral of abs(f-i/(b-a)) */

/*           machine dependent constants */
/*           --------------------------- */

/*           epmach is the largest relative spacing. */
/*           uflow is the smallest positive magnitude. */

/* ***first executable statement  qng */
    epmach = r1mach_(&c__4);
    uflow = r1mach_(&c__1);

/*           test on validity of parameters */
/*           ------------------------------ */

    *result = 0.f;
    *abserr = 0.f;
    *neval = 0;
    *ier = 6;
/* Computing MAX */
    r__1 = 5e-15f, r__2 = epmach * 50.f;
    if (*epsabs <= 0.f && *epsrel < dmax(r__1,r__2)) {
	goto L80;
    }
    hlgth = (*b - *a) * .5f;
    dhlgth = dabs(hlgth);
    centr = (*b + *a) * .5f;
    fcentr = (*f)(&centr);
    *neval = 21;
    *ier = 1;

/*          compute the integral using the 10- and 21-point formula. */

    for (l = 1; l <= 3; ++l) {
	switch (l) {
	    case 1:  goto L5;
	    case 2:  goto L25;
	    case 3:  goto L45;
	}
L5:
	res10 = 0.f;
	res21 = w21b[5] * fcentr;
	resabs = w21b[5] * dabs(fcentr);
	for (k = 1; k <= 5; ++k) {
	    absc = hlgth * x1[k - 1];
	    r__1 = centr + absc;
	    fval1 = (*f)(&r__1);
	    r__1 = centr - absc;
	    fval2 = (*f)(&r__1);
	    fval = fval1 + fval2;
	    res10 += w10[k - 1] * fval;
	    res21 += w21a[k - 1] * fval;
	    resabs += w21a[k - 1] * (dabs(fval1) + dabs(fval2));
	    savfun[k - 1] = fval;
	    fv1[k - 1] = fval1;
	    fv2[k - 1] = fval2;
/* L10: */
	}
	ipx = 5;
	for (k = 1; k <= 5; ++k) {
	    ++ipx;
	    absc = hlgth * x2[k - 1];
	    r__1 = centr + absc;
	    fval1 = (*f)(&r__1);
	    r__1 = centr - absc;
	    fval2 = (*f)(&r__1);
	    fval = fval1 + fval2;
	    res21 += w21b[k - 1] * fval;
	    resabs += w21b[k - 1] * (dabs(fval1) + dabs(fval2));
	    savfun[ipx - 1] = fval;
	    fv3[k - 1] = fval1;
	    fv4[k - 1] = fval2;
/* L15: */
	}

/*          test for convergence. */

	*result = res21 * hlgth;
	resabs *= dhlgth;
	reskh = res21 * .5f;
	resasc = w21b[5] * (r__1 = fcentr - reskh, dabs(r__1));
	for (k = 1; k <= 5; ++k) {
	    resasc = resasc + w21a[k - 1] * ((r__1 = fv1[k - 1] - reskh, dabs(
		    r__1)) + (r__2 = fv2[k - 1] - reskh, dabs(r__2))) + w21b[
		    k - 1] * ((r__3 = fv3[k - 1] - reskh, dabs(r__3)) + (r__4 
		    = fv4[k - 1] - reskh, dabs(r__4)));
/* L20: */
	}
	*abserr = (r__1 = (res21 - res10) * hlgth, dabs(r__1));
	resasc *= dhlgth;
	goto L65;

/*          compute the integral using the 43-point formula. */

L25:
	res43 = w43b[11] * fcentr;
	*neval = 43;
	for (k = 1; k <= 10; ++k) {
	    res43 += savfun[k - 1] * w43a[k - 1];
/* L30: */
	}
	for (k = 1; k <= 11; ++k) {
	    ++ipx;
	    absc = hlgth * x3[k - 1];
	    r__1 = absc + centr;
	    r__2 = centr - absc;
	    fval = (*f)(&r__1) + (*f)(&r__2);
	    res43 += fval * w43b[k - 1];
	    savfun[ipx - 1] = fval;
/* L40: */
	}

/*          test for convergence. */

	*result = res43 * hlgth;
	*abserr = (r__1 = (res43 - res21) * hlgth, dabs(r__1));
	goto L65;

/*          compute the integral using the 87-point formula. */

L45:
	res87 = w87b[22] * fcentr;
	*neval = 87;
	for (k = 1; k <= 21; ++k) {
	    res87 += savfun[k - 1] * w87a[k - 1];
/* L50: */
	}
	for (k = 1; k <= 22; ++k) {
	    absc = hlgth * x4[k - 1];
	    r__1 = absc + centr;
	    r__2 = centr - absc;
	    res87 += w87b[k - 1] * ((*f)(&r__1) + (*f)(&r__2));
/* L60: */
	}
	*result = res87 * hlgth;
	*abserr = (r__1 = (res87 - res43) * hlgth, dabs(r__1));
L65:
	if (resasc != 0.f && *abserr != 0.f) {
/* Computing MIN */
	    d__1 = (doublereal) (*abserr * 200.f / resasc);
	    r__1 = 1.f, r__2 = pow_dd(&d__1, &c_b17);
	    *abserr = resasc * dmin(r__1,r__2);
	}
	if (resabs > uflow / (epmach * 50.f)) {
/* Computing MAX */
	    r__1 = epmach * 50.f * resabs;
	    *abserr = dmax(r__1,*abserr);
	}
/* Computing MAX */
	r__1 = *epsabs, r__2 = *epsrel * dabs(*result);
	if (*abserr <= dmax(r__1,r__2)) {
	    *ier = 0;
	}
/* ***jump out of do-loop */
	if (*ier == 0) {
	    goto L999;
	}
/* L70: */
    }
L80:
    xerror_("abnormal return from  qng ", &c__26, ier, &c__0, (ftnlen)26);
L999:
    return 0;
} /* qng_ */

