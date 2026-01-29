/* cgbdi.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int cgbdi_(complex *abd, integer *lda, integer *n, integer *
	ml, integer *mu, integer *ipvt, complex *det)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2;
    real r__1, r__2;
    complex q__1, q__2;

    /* Builtin functions */
    double r_imag(complex *);
    void c_div(complex *, complex *, complex *);

    /* Local variables */
    integer i__, m;
    real ten;


/*     cgbdi computes the determinant of a band matrix */
/*     using the factors computed by cgbco or cgbfa. */
/*     if the inverse is needed, use cgbsl  n  times. */

/*     on entry */

/*        abd     complex(lda, n) */
/*                the output from cgbco or cgbfa. */

/*        lda     integer */
/*                the leading dimension of the array  abd . */

/*        n       integer */
/*                the order of the original matrix. */

/*        ml      integer */
/*                number of diagonals below the main diagonal. */

/*        mu      integer */
/*                number of diagonals above the main diagonal. */

/*        ipvt    integer(n) */
/*                the pivot vector from cgbco or cgbfa. */

/*     on return */

/*        det     complex(2) */
/*                determinant of original matrix. */
/*                determinant = det(1) * 10.0**det(2) */
/*                with  1.0 .le. cabs1(det(1)) .lt. 10.0 */
/*                or  det(1) = 0.0 . */

/*     linpack. this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     fortran abs,aimag,cmplx,real */

/*     internal variables */



    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --ipvt;
    --det;

    /* Function Body */
    m = *ml + *mu + 1;
    det[1].r = 1.f, det[1].i = 0.f;
    det[2].r = 0.f, det[2].i = 0.f;
    ten = 10.f;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (ipvt[i__] != i__) {
	    q__1.r = -det[1].r, q__1.i = -det[1].i;
	    det[1].r = q__1.r, det[1].i = q__1.i;
	}
	i__2 = m + i__ * abd_dim1;
	q__1.r = abd[i__2].r * det[1].r - abd[i__2].i * det[1].i, q__1.i = 
		abd[i__2].r * det[1].i + abd[i__2].i * det[1].r;
	det[1].r = q__1.r, det[1].i = q__1.i;
/*     ...exit */
	if ((r__1 = det[1].r, dabs(r__1)) + (r__2 = r_imag(&det[1]), dabs(
		r__2)) == 0.f) {
	    goto L60;
	}
L10:
	if ((r__1 = det[1].r, dabs(r__1)) + (r__2 = r_imag(&det[1]), dabs(
		r__2)) >= 1.f) {
	    goto L20;
	}
	q__2.r = ten, q__2.i = 0.f;
	q__1.r = q__2.r * det[1].r - q__2.i * det[1].i, q__1.i = q__2.r * det[
		1].i + q__2.i * det[1].r;
	det[1].r = q__1.r, det[1].i = q__1.i;
	q__1.r = det[2].r - 1.f, q__1.i = det[2].i - 0.f;
	det[2].r = q__1.r, det[2].i = q__1.i;
	goto L10;
L20:
L30:
	if ((r__1 = det[1].r, dabs(r__1)) + (r__2 = r_imag(&det[1]), dabs(
		r__2)) < ten) {
	    goto L40;
	}
	q__2.r = ten, q__2.i = 0.f;
	c_div(&q__1, &det[1], &q__2);
	det[1].r = q__1.r, det[1].i = q__1.i;
	q__1.r = det[2].r + 1.f, q__1.i = det[2].i + 0.f;
	det[2].r = q__1.r, det[2].i = q__1.i;
	goto L30;
L40:
/* L50: */
	;
    }
L60:
    return 0;
} /* cgbdi_ */

