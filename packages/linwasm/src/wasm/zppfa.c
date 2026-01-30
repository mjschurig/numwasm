/* zppfa.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int zppfa_(doublecomplex *ap, integer *n, integer *info)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);
    double sqrt(doublereal);

    /* Local variables */
    integer j, k;
    doublereal s;
    doublecomplex t;
    integer jj, kj, kk, jm1;
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);


/*     zppfa factors a complex*16 hermitian positive definite matrix */
/*     stored in packed form. */

/*     zppfa is usually called by zppco, but it can be called */
/*     directly with a saving in time if  rcond  is not needed. */
/*     (time for zppco) = (1 + 18/n)*(time for zppfa) . */

/*     on entry */

/*        ap      complex*16 (n*(n+1)/2) */
/*                the packed form of a hermitian matrix  a .  the */
/*                columns of the upper triangle are stored sequentially */
/*                in a one-dimensional array of length  n*(n+1)/2 . */
/*                see comments below for details. */

/*        n       integer */
/*                the order of the matrix  a . */

/*     on return */

/*        ap      an upper triangular matrix  r , stored in packed */
/*                form, so that  a = ctrans(r)*r . */

/*        info    integer */
/*                = 0  for normal return. */
/*                = k  if the leading minor of order  k  is not */
/*                     positive definite. */


/*     packed storage */

/*          the following program segment will pack the upper */
/*          triangle of a hermitian matrix. */

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

/*     blas zdotc */
/*     fortran dcmplx,dconjg,dsqrt */

/*     internal variables */

/*     begin block with ...exits to 40 */


    /* Parameter adjustments */
    --ap;

    /* Function Body */
    jj = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	*info = j;
	s = 0.;
	jm1 = j - 1;
	kj = jj;
	kk = 0;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (k = 1; k <= i__2; ++k) {
	    ++kj;
	    i__3 = kj;
	    i__4 = k - 1;
	    zdotc_(&z__2, &i__4, &ap[kk + 1], &c__1, &ap[jj + 1], &c__1);
	    z__1.r = ap[i__3].r - z__2.r, z__1.i = ap[i__3].i - z__2.i;
	    t.r = z__1.r, t.i = z__1.i;
	    kk += k;
	    z_div(&z__1, &t, &ap[kk]);
	    t.r = z__1.r, t.i = z__1.i;
	    i__3 = kj;
	    ap[i__3].r = t.r, ap[i__3].i = t.i;
	    d_cnjg(&z__3, &t);
	    z__2.r = t.r * z__3.r - t.i * z__3.i, z__2.i = t.r * z__3.i + t.i 
		    * z__3.r;
	    z__1.r = z__2.r, z__1.i = z__2.i;
	    s += z__1.r;
/* L10: */
	}
L20:
	jj += j;
	i__2 = jj;
	s = ap[i__2].r - s;
/*     ......exit */
	i__2 = jj;
	z__1.r = ap[i__2].r * 0. - ap[i__2].i * -1., z__1.i = ap[i__2].i * 0. 
		+ ap[i__2].r * -1.;
	if (s <= 0. || z__1.r != 0.) {
	    goto L40;
	}
	i__2 = jj;
	d__1 = sqrt(s);
	z__1.r = d__1, z__1.i = 0.;
	ap[i__2].r = z__1.r, ap[i__2].i = z__1.i;
/* L30: */
    }
    *info = 0;
L40:
    return 0;
} /* zppfa_ */

