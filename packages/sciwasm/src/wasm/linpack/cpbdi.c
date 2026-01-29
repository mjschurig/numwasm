/* cpbdi.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int cpbdi_(complex *abd, integer *lda, integer *n, integer *
	m, real *det)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2;
    real r__1;

    /* Local variables */
    integer i__;
    real s;


/*     cpbdi computes the determinant */
/*     of a complex hermitian positive definite band matrix */
/*     using the factors computed by cpbco or cpbfa. */
/*     if the inverse is needed, use cpbsl  n  times. */

/*     on entry */

/*        abd     complex(lda, n) */
/*                the output from cpbco or cpbfa. */

/*        lda     integer */
/*                the leading dimension of the array  abd . */

/*        n       integer */
/*                the order of the matrix  a . */

/*        m       integer */
/*                the number of diagonals above the main diagonal. */

/*     on return */

/*        det     real(2) */
/*                determinant of original matrix in the form */
/*                determinant = det(1) * 10.0**det(2) */
/*                with  1.0 .le. det(1) .lt. 10.0 */
/*                or  det(1) .eq. 0.0 . */

/*     linpack.  this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     fortran real */

/*     internal variables */


/*     compute determinant */

    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --det;

    /* Function Body */
    det[1] = 1.f;
    det[2] = 0.f;
    s = 10.f;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = *m + 1 + i__ * abd_dim1;
/* Computing 2nd power */
	r__1 = abd[i__2].r;
	det[1] = r__1 * r__1 * det[1];
/*     ...exit */
	if (det[1] == 0.f) {
	    goto L60;
	}
L10:
	if (det[1] >= 1.f) {
	    goto L20;
	}
	det[1] = s * det[1];
	det[2] += -1.f;
	goto L10;
L20:
L30:
	if (det[1] < s) {
	    goto L40;
	}
	det[1] /= s;
	det[2] += 1.f;
	goto L30;
L40:
/* L50: */
	;
    }
L60:
    return 0;
} /* cpbdi_ */

