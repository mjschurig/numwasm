/* zpofa.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int zpofa_(doublecomplex *a, integer *lda, integer *n, 
	integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
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
    integer jm1;
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);


/*     zpofa factors a complex*16 hermitian positive definite matrix. */

/*     zpofa is usually called by zpoco, but it can be called */
/*     directly with a saving in time if  rcond  is not needed. */
/*     (time for zpoco) = (1 + 18/n)*(time for zpofa) . */

/*     on entry */

/*        a       complex*16(lda, n) */
/*                the hermitian matrix to be factored.  only the */
/*                diagonal and upper triangle are used. */

/*        lda     integer */
/*                the leading dimension of the array  a . */

/*        n       integer */
/*                the order of the matrix  a . */

/*     on return */

/*        a       an upper triangular matrix  r  so that  a = */
/*                ctrans(r)*r where  ctrans(r)  is the conjugate */
/*                transpose.  the strict lower triangle is unaltered. */
/*                if  info .ne. 0 , the factorization is not complete. */

/*        info    integer */
/*                = 0  for normal return. */
/*                = k  signals an error condition.  the leading minor */
/*                     of order  k  is not positive definite. */

/*     linpack.  this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas zdotc */
/*     fortran dcmplx,dconjg,dsqrt */

/*     internal variables */

/*     begin block with ...exits to 40 */


    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	*info = j;
	s = 0.;
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (k = 1; k <= i__2; ++k) {
	    i__3 = k + j * a_dim1;
	    i__4 = k - 1;
	    zdotc_(&z__2, &i__4, &a[k * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1]
		    , &c__1);
	    z__1.r = a[i__3].r - z__2.r, z__1.i = a[i__3].i - z__2.i;
	    t.r = z__1.r, t.i = z__1.i;
	    z_div(&z__1, &t, &a[k + k * a_dim1]);
	    t.r = z__1.r, t.i = z__1.i;
	    i__3 = k + j * a_dim1;
	    a[i__3].r = t.r, a[i__3].i = t.i;
	    d_cnjg(&z__3, &t);
	    z__2.r = t.r * z__3.r - t.i * z__3.i, z__2.i = t.r * z__3.i + t.i 
		    * z__3.r;
	    z__1.r = z__2.r, z__1.i = z__2.i;
	    s += z__1.r;
/* L10: */
	}
L20:
	i__2 = j + j * a_dim1;
	s = a[i__2].r - s;
/*     ......exit */
	i__2 = j + j * a_dim1;
	z__1.r = a[i__2].r * 0. - a[i__2].i * -1., z__1.i = a[i__2].i * 0. + 
		a[i__2].r * -1.;
	if (s <= 0. || z__1.r != 0.) {
	    goto L40;
	}
	i__2 = j + j * a_dim1;
	d__1 = sqrt(s);
	z__1.r = d__1, z__1.i = 0.;
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
/* L30: */
    }
    *info = 0;
L40:
    return 0;
} /* zpofa_ */

