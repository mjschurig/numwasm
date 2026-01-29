/* zpbfa.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int zpbfa_(doublecomplex *abd, integer *lda, integer *n, 
	integer *m, integer *info)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3, i__4;
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
    integer ik, jk, mu;
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);


/*     zpbfa factors a complex*16 hermitian positive definite matrix */
/*     stored in band form. */

/*     zpbfa is usually called by zpbco, but it can be called */
/*     directly with a saving in time if  rcond  is not needed. */

/*     on entry */

/*        abd     complex*16(lda, n) */
/*                the matrix to be factored.  the columns of the upper */
/*                triangle are stored in the columns of abd and the */
/*                diagonals of the upper triangle are stored in the */
/*                rows of abd .  see the comments below for details. */

/*        lda     integer */
/*                the leading dimension of the array  abd . */
/*                lda must be .ge. m + 1 . */

/*        n       integer */
/*                the order of the matrix  a . */

/*        m       integer */
/*                the number of diagonals above the main diagonal. */
/*                0 .le. m .lt. n . */

/*     on return */

/*        abd     an upper triangular matrix  r , stored in band */
/*                form, so that  a = ctrans(r)*r . */

/*        info    integer */
/*                = 0  for normal return. */
/*                = k  if the leading minor of order  k  is not */
/*                     positive definite. */

/*     band storage */

/*           if  a  is a hermitian positive definite band matrix, */
/*           the following program segment will set up the input. */

/*                   m = (band width above diagonal) */
/*                   do 20 j = 1, n */
/*                      i1 = max0(1, j-m) */
/*                      do 10 i = i1, j */
/*                         k = i-j+m+1 */
/*                         abd(k,j) = a(i,j) */
/*                10    continue */
/*                20 continue */

/*     linpack.  this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas zdotc */
/*     fortran dcmplx,dconjg,max0,dsqrt */

/*     internal variables */

/*     begin block with ...exits to 40 */


    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	*info = j;
	s = 0.;
	ik = *m + 1;
/* Computing MAX */
	i__2 = j - *m;
	jk = max(i__2,1);
/* Computing MAX */
	i__2 = *m + 2 - j;
	mu = max(i__2,1);
	if (*m < mu) {
	    goto L20;
	}
	i__2 = *m;
	for (k = mu; k <= i__2; ++k) {
	    i__3 = k + j * abd_dim1;
	    i__4 = k - mu;
	    zdotc_(&z__2, &i__4, &abd[ik + jk * abd_dim1], &c__1, &abd[mu + j 
		    * abd_dim1], &c__1);
	    z__1.r = abd[i__3].r - z__2.r, z__1.i = abd[i__3].i - z__2.i;
	    t.r = z__1.r, t.i = z__1.i;
	    z_div(&z__1, &t, &abd[*m + 1 + jk * abd_dim1]);
	    t.r = z__1.r, t.i = z__1.i;
	    i__3 = k + j * abd_dim1;
	    abd[i__3].r = t.r, abd[i__3].i = t.i;
	    d_cnjg(&z__3, &t);
	    z__2.r = t.r * z__3.r - t.i * z__3.i, z__2.i = t.r * z__3.i + t.i 
		    * z__3.r;
	    z__1.r = z__2.r, z__1.i = z__2.i;
	    s += z__1.r;
	    --ik;
	    ++jk;
/* L10: */
	}
L20:
	i__2 = *m + 1 + j * abd_dim1;
	s = abd[i__2].r - s;
/*     ......exit */
	i__2 = *m + 1 + j * abd_dim1;
	z__1.r = abd[i__2].r * 0. - abd[i__2].i * -1., z__1.i = abd[i__2].i * 
		0. + abd[i__2].r * -1.;
	if (s <= 0. || z__1.r != 0.) {
	    goto L40;
	}
	i__2 = *m + 1 + j * abd_dim1;
	d__1 = sqrt(s);
	z__1.r = d__1, z__1.i = 0.;
	abd[i__2].r = z__1.r, abd[i__2].i = z__1.i;
/* L30: */
    }
    *info = 0;
L40:
    return 0;
} /* zpbfa_ */

