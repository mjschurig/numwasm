/* sppfa.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int sppfa_(real *ap, integer *n, integer *info)
{
    /* System generated locals */
    integer i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer j, k;
    real s, t;
    integer jj, kj, kk, jm1;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);


/*     sppfa factors a real symmetric positive definite matrix */
/*     stored in packed form. */

/*     sppfa is usually called by sppco, but it can be called */
/*     directly with a saving in time if  rcond  is not needed. */
/*     (time for sppco) = (1 + 18/n)*(time for sppfa) . */

/*     on entry */

/*        ap      real (n*(n+1)/2) */
/*                the packed form of a symmetric matrix  a .  the */
/*                columns of the upper triangle are stored sequentially */
/*                in a one-dimensional array of length  n*(n+1)/2 . */
/*                see comments below for details. */

/*        n       integer */
/*                the order of the matrix  a . */

/*     on return */

/*        ap      an upper triangular matrix  r , stored in packed */
/*                form, so that  a = trans(r)*r . */

/*        info    integer */
/*                = 0  for normal return. */
/*                = k  if the leading minor of order  k  is not */
/*                     positive definite. */


/*     packed storage */

/*          the following program segment will pack the upper */
/*          triangle of a symmetric matrix. */

/*                k = 0 */
/*                do 20 j = 1, n */
/*                   do 10 i = 1, j */
/*                      k = k + 1 */
/*                      ap(k) = a(i,j) */
/*             10    continue */
/*             20 continue */

/*     linpack.  this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas sdot */
/*     fortran sqrt */

/*     internal variables */

/*     begin block with ...exits to 40 */


    /* Parameter adjustments */
    --ap;

    /* Function Body */
    jj = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	*info = j;
	s = 0.f;
	jm1 = j - 1;
	kj = jj;
	kk = 0;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (k = 1; k <= i__2; ++k) {
	    ++kj;
	    i__3 = k - 1;
	    t = ap[kj] - sdot_(&i__3, &ap[kk + 1], &c__1, &ap[jj + 1], &c__1);
	    kk += k;
	    t /= ap[kk];
	    ap[kj] = t;
	    s += t * t;
/* L10: */
	}
L20:
	jj += j;
	s = ap[jj] - s;
/*     ......exit */
	if (s <= 0.f) {
	    goto L40;
	}
	ap[jj] = sqrt(s);
/* L30: */
    }
    *info = 0;
L40:
    return 0;
} /* sppfa_ */

