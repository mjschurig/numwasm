/* dppdi.f -- translated by f2c (version 20240504).
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

static integer c__1 = 1;

/* Subroutine */ int dppdi_(doublereal *ap, integer *n, doublereal *det, 
	integer *job)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    integer i__, j, k;
    doublereal s, t;
    integer j1, k1, ii, jj, kj, kk, jm1, kp1;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), daxpy_(integer *, doublereal *, doublereal *, integer 
	    *, doublereal *, integer *);


/*     dppdi computes the determinant and inverse */
/*     of a double precision symmetric positive definite matrix */
/*     using the factors computed by dppco or dppfa . */

/*     on entry */

/*        ap      double precision (n*(n+1)/2) */
/*                the output from dppco or dppfa. */

/*        n       integer */
/*                the order of the matrix  a . */

/*        job     integer */
/*                = 11   both determinant and inverse. */
/*                = 01   inverse only. */
/*                = 10   determinant only. */

/*     on return */

/*        ap      the upper triangular half of the inverse . */
/*                the strict lower triangle is unaltered. */

/*        det     double precision(2) */
/*                determinant of original matrix if requested. */
/*                otherwise not referenced. */
/*                determinant = det(1) * 10.0**det(2) */
/*                with  1.0 .le. det(1) .lt. 10.0 */
/*                or  det(1) .eq. 0.0 . */

/*     error condition */

/*        a division by zero will occur if the input factor contains */
/*        a zero on the diagonal and the inverse is requested. */
/*        it will not occur if the subroutines are called correctly */
/*        and if dpoco or dpofa has set info .eq. 0 . */

/*     linpack.  this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas daxpy,dscal */
/*     fortran mod */

/*     internal variables */


/*     compute determinant */

    /* Parameter adjustments */
    --det;
    --ap;

    /* Function Body */
    if (*job / 10 == 0) {
	goto L70;
    }
    det[1] = 1.;
    det[2] = 0.;
    s = 10.;
    ii = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii += i__;
/* Computing 2nd power */
	d__1 = ap[ii];
	det[1] = d__1 * d__1 * det[1];
/*        ...exit */
	if (det[1] == 0.) {
	    goto L60;
	}
L10:
	if (det[1] >= 1.) {
	    goto L20;
	}
	det[1] = s * det[1];
	det[2] += -1.;
	goto L10;
L20:
L30:
	if (det[1] < s) {
	    goto L40;
	}
	det[1] /= s;
	det[2] += 1.;
	goto L30;
L40:
/* L50: */
	;
    }
L60:
L70:

/*     compute inverse(r) */

    if (*job % 10 == 0) {
	goto L140;
    }
    kk = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	k1 = kk + 1;
	kk += k;
	ap[kk] = 1. / ap[kk];
	t = -ap[kk];
	i__2 = k - 1;
	dscal_(&i__2, &t, &ap[k1], &c__1);
	kp1 = k + 1;
	j1 = kk + 1;
	kj = kk + k;
	if (*n < kp1) {
	    goto L90;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    t = ap[kj];
	    ap[kj] = 0.;
	    daxpy_(&k, &t, &ap[k1], &c__1, &ap[j1], &c__1);
	    j1 += j;
	    kj += j;
/* L80: */
	}
L90:
/* L100: */
	;
    }

/*        form  inverse(r) * trans(inverse(r)) */

    jj = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	j1 = jj + 1;
	jj += j;
	jm1 = j - 1;
	k1 = 1;
	kj = j1;
	if (jm1 < 1) {
	    goto L120;
	}
	i__2 = jm1;
	for (k = 1; k <= i__2; ++k) {
	    t = ap[kj];
	    daxpy_(&k, &t, &ap[j1], &c__1, &ap[k1], &c__1);
	    k1 += k;
	    ++kj;
/* L110: */
	}
L120:
	t = ap[jj];
	dscal_(&j, &t, &ap[j1], &c__1);
/* L130: */
    }
L140:
    return 0;
} /* dppdi_ */

