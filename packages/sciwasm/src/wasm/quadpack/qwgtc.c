/* qwgtc.f -- translated by f2c (version 20240504).
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

doublereal qwgtc_(real *x, real *c__, real *p2, real *p3, real *p4, integer *
	kp)
{
    /* System generated locals */
    real ret_val;

/* ***begin prologue  qwgtc */
/* ***refer to qk15w */
/* ***routines called  (none) */
/* ***revision date  810101   (yymmdd) */
/* ***keywords  weight function, cauchy principal value */
/* ***author  piessens,robert,appl. math. & progr. div. - k.u.leuven */
/*           de doncker,elise,appl. math. & progr. div. - k.u.leuven */
/* ***purpose  this function subprogram is used together with the */
/*            routine qawc and defines the weight function. */
/* ***end prologue  qwgtc */

/* ***first executable statement */
    ret_val = 1.f / (*x - *c__);
    return ret_val;
} /* qwgtc_ */

