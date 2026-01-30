/* qwgts.f -- translated by f2c (version 20240504).
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

doublereal qwgts_(real *x, real *a, real *b, real *alfa, real *beta, integer *
	integr)
{
    /* System generated locals */
    real ret_val;
    doublereal d__1, d__2, d__3, d__4;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), log(doublereal);

    /* Local variables */
    real xma, bmx;

/* ***begin prologue  qwgts */
/* ***refer to qk15w */
/* ***routines called  (none) */
/* ***revision date  810101   (yymmdd) */
/* ***keywords  weight function, algebraico-logarithmic */
/*             end-point singularities */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
/* ***purpose  this function subprogram is used together with the */
/*            routine qaws and defines the weight function. */
/* ***end prologue  qwgts */

/* ***first executable statement */
    xma = *x - *a;
    bmx = *b - *x;
    d__1 = (doublereal) xma;
    d__2 = (doublereal) (*alfa);
    d__3 = (doublereal) bmx;
    d__4 = (doublereal) (*beta);
    ret_val = pow_dd(&d__1, &d__2) * pow_dd(&d__3, &d__4);
    switch (*integr) {
	case 1:  goto L40;
	case 2:  goto L10;
	case 3:  goto L20;
	case 4:  goto L30;
    }
L10:
    ret_val *= log(xma);
    goto L40;
L20:
    ret_val *= log(bmx);
    goto L40;
L30:
    ret_val = ret_val * log(xma) * log(bmx);
L40:
    return ret_val;
} /* qwgts_ */

