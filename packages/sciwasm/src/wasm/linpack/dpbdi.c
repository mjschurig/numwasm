/* dpbdi.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dpbdi_(doublereal *abd, integer *lda, integer *n, 
	integer *m, doublereal *det)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1;
    doublereal d__1;

    /* Local variables */
    integer i__;
    doublereal s;


/*     dpbdi computes the determinant */
/*     of a double precision symmetric positive definite band matrix */
/*     using the factors computed by dpbco or dpbfa. */
/*     if the inverse is needed, use dpbsl  n  times. */

/*     on entry */

/*        abd     double precision(lda, n) */
/*                the output from dpbco or dpbfa. */

/*        lda     integer */
/*                the leading dimension of the array  abd . */

/*        n       integer */
/*                the order of the matrix  a . */

/*        m       integer */
/*                the number of diagonals above the main diagonal. */

/*     on return */

/*        det     double precision(2) */
/*                determinant of original matrix in the form */
/*                determinant = det(1) * 10.0**det(2) */
/*                with  1.0 .le. det(1) .lt. 10.0 */
/*                or  det(1) .eq. 0.0 . */

/*     linpack.  this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */


/*     internal variables */


/*     compute determinant */

    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --det;

    /* Function Body */
    det[1] = 1.;
    det[2] = 0.;
    s = 10.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* Computing 2nd power */
	d__1 = abd[*m + 1 + i__ * abd_dim1];
	det[1] = d__1 * d__1 * det[1];
/*     ...exit */
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
    return 0;
} /* dpbdi_ */

