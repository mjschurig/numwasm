/* sgbdi.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int sgbdi_(real *abd, integer *lda, integer *n, integer *ml, 
	integer *mu, integer *ipvt, real *det)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1;

    /* Local variables */
    integer i__, m;
    real ten;


/*     sgbdi computes the determinant of a band matrix */
/*     using the factors computed by sgbco or sgbfa. */
/*     if the inverse is needed, use sgbsl  n  times. */

/*     on entry */

/*        abd     real(lda, n) */
/*                the output from sgbco or sgbfa. */

/*        lda     integer */
/*                the leading dimension of the array  abd . */

/*        n       integer */
/*                the order of the original matrix. */

/*        ml      integer */
/*                number of diagonals below the main diagonal. */

/*        mu      integer */
/*                number of diagonals above the main diagonal. */

/*        ipvt    integer(n) */
/*                the pivot vector from sgbco or sgbfa. */

/*     on return */

/*        det     real(2) */
/*                determinant of original matrix. */
/*                determinant = det(1) * 10.0**det(2) */
/*                with  1.0 .le. abs(det(1)) .lt. 10.0 */
/*                or  det(1) = 0.0 . */

/*     linpack. this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     fortran abs */

/*     internal variables */



    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --ipvt;
    --det;

    /* Function Body */
    m = *ml + *mu + 1;
    det[1] = 1.f;
    det[2] = 0.f;
    ten = 10.f;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (ipvt[i__] != i__) {
	    det[1] = -det[1];
	}
	det[1] = abd[m + i__ * abd_dim1] * det[1];
/*     ...exit */
	if (det[1] == 0.f) {
	    goto L60;
	}
L10:
	if (dabs(det[1]) >= 1.f) {
	    goto L20;
	}
	det[1] = ten * det[1];
	det[2] += -1.f;
	goto L10;
L20:
L30:
	if (dabs(det[1]) < ten) {
	    goto L40;
	}
	det[1] /= ten;
	det[2] += 1.f;
	goto L30;
L40:
/* L50: */
	;
    }
L60:
    return 0;
} /* sgbdi_ */

