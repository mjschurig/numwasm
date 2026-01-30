/* dptsl.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dptsl_(integer *n, doublereal *d__, doublereal *e, 
	doublereal *b)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer k;
    doublereal t1, t2;
    integer ke, kf, kp1, nm1, kbm1, nm1d2;


/*     dptsl given a positive definite tridiagonal matrix and a right */
/*     hand side will find the solution. */

/*     on entry */

/*        n        integer */
/*                 is the order of the tridiagonal matrix. */

/*        d        double precision(n) */
/*                 is the diagonal of the tridiagonal matrix. */
/*                 on output d is destroyed. */

/*        e        double precision(n) */
/*                 is the offdiagonal of the tridiagonal matrix. */
/*                 e(1) through e(n-1) should contain the */
/*                 offdiagonal. */

/*        b        double precision(n) */
/*                 is the right hand side vector. */

/*     on return */

/*        b        contains the soultion. */

/*     linpack. this version dated 08/14/78 . */
/*     jack dongarra, argonne national laboratory. */

/*     no externals */
/*     fortran mod */

/*     internal variables */


/*     check for 1 x 1 case */

    /* Parameter adjustments */
    --b;
    --e;
    --d__;

    /* Function Body */
    if (*n != 1) {
	goto L10;
    }
    b[1] /= d__[1];
    goto L70;
L10:
    nm1 = *n - 1;
    nm1d2 = nm1 / 2;
    if (*n == 2) {
	goto L30;
    }
    kbm1 = *n - 1;

/*           zero top half of subdiagonal and bottom half of */
/*           superdiagonal */

    i__1 = nm1d2;
    for (k = 1; k <= i__1; ++k) {
	t1 = e[k] / d__[k];
	d__[k + 1] -= t1 * e[k];
	b[k + 1] -= t1 * b[k];
	t2 = e[kbm1] / d__[kbm1 + 1];
	d__[kbm1] -= t2 * e[kbm1];
	b[kbm1] -= t2 * b[kbm1 + 1];
	--kbm1;
/* L20: */
    }
L30:
    kp1 = nm1d2 + 1;

/*        clean up for possible 2 x 2 block at center */

    if (*n % 2 != 0) {
	goto L40;
    }
    t1 = e[kp1] / d__[kp1];
    d__[kp1 + 1] -= t1 * e[kp1];
    b[kp1 + 1] -= t1 * b[kp1];
    ++kp1;
L40:

/*        back solve starting at the center, going towards the top */
/*        and bottom */

    b[kp1] /= d__[kp1];
    if (*n == 2) {
	goto L60;
    }
    k = kp1 - 1;
    ke = kp1 + nm1d2 - 1;
    i__1 = ke;
    for (kf = kp1; kf <= i__1; ++kf) {
	b[k] = (b[k] - e[k] * b[k + 1]) / d__[k];
	b[kf + 1] = (b[kf + 1] - e[kf] * b[kf]) / d__[kf + 1];
	--k;
/* L50: */
    }
L60:
    if (*n % 2 == 0) {
	b[1] = (b[1] - e[1] * b[2]) / d__[1];
    }
L70:
    return 0;
} /* dptsl_ */

