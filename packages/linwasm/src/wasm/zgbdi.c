/* zgbdi.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int zgbdi_(doublecomplex *abd, integer *lda, integer *n, 
	integer *ml, integer *mu, integer *ipvt, doublecomplex *det)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    integer i__, m;
    doublereal ten;


/*     zgbdi computes the determinant of a band matrix */
/*     using the factors computed by zgbco or zgbfa. */
/*     if the inverse is needed, use zgbsl  n  times. */

/*     on entry */

/*        abd     complex*16(lda, n) */
/*                the output from zgbco or zgbfa. */

/*        lda     integer */
/*                the leading dimension of the array  abd . */

/*        n       integer */
/*                the order of the original matrix. */

/*        ml      integer */
/*                number of diagonals below the main diagonal. */

/*        mu      integer */
/*                number of diagonals above the main diagonal. */

/*        ipvt    integer(n) */
/*                the pivot vector from zgbco or zgbfa. */

/*     on return */

/*        det     complex*16(2) */
/*                determinant of original matrix. */
/*                determinant = det(1) * 10.0**det(2) */
/*                with  1.0 .le. cabs1(det(1)) .lt. 10.0 */
/*                or  det(1) = 0.0 . */

/*     linpack. this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     fortran dabs,dcmplx */

/*     internal variables */



    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --ipvt;
    --det;

    /* Function Body */
    m = *ml + *mu + 1;
    det[1].r = 1., det[1].i = 0.;
    det[2].r = 0., det[2].i = 0.;
    ten = 10.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (ipvt[i__] != i__) {
	    z__1.r = -det[1].r, z__1.i = -det[1].i;
	    det[1].r = z__1.r, det[1].i = z__1.i;
	}
	i__2 = m + i__ * abd_dim1;
	z__1.r = abd[i__2].r * det[1].r - abd[i__2].i * det[1].i, z__1.i = 
		abd[i__2].r * det[1].i + abd[i__2].i * det[1].r;
	det[1].r = z__1.r, det[1].i = z__1.i;
/*     ...exit */
	z__1.r = det[1].r * 0. - det[1].i * -1., z__1.i = det[1].i * 0. + det[
		1].r * -1.;
	if ((d__1 = det[1].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) == 0.) {
	    goto L60;
	}
L10:
	z__1.r = det[1].r * 0. - det[1].i * -1., z__1.i = det[1].i * 0. + det[
		1].r * -1.;
	if ((d__1 = det[1].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) >= 1.) {
	    goto L20;
	}
	z__2.r = ten, z__2.i = 0.;
	z__1.r = z__2.r * det[1].r - z__2.i * det[1].i, z__1.i = z__2.r * det[
		1].i + z__2.i * det[1].r;
	det[1].r = z__1.r, det[1].i = z__1.i;
	z__1.r = det[2].r - 1., z__1.i = det[2].i - 0.;
	det[2].r = z__1.r, det[2].i = z__1.i;
	goto L10;
L20:
L30:
	z__1.r = det[1].r * 0. - det[1].i * -1., z__1.i = det[1].i * 0. + det[
		1].r * -1.;
	if ((d__1 = det[1].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) < ten) {
	    goto L40;
	}
	z__2.r = ten, z__2.i = 0.;
	z_div(&z__1, &det[1], &z__2);
	det[1].r = z__1.r, det[1].i = z__1.i;
	z__1.r = det[2].r + 1., z__1.i = det[2].i + 0.;
	det[2].r = z__1.r, det[2].i = z__1.i;
	goto L30;
L40:
/* L50: */
	;
    }
L60:
    return 0;
} /* zgbdi_ */

