/* sgedi.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int sgedi_(real *a, integer *lda, integer *n, integer *ipvt, 
	real *det, real *work, integer *job)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    integer i__, j, k, l;
    real t;
    integer kb, kp1, nm1;
    real ten;
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *), 
	    sswap_(integer *, real *, integer *, real *, integer *), saxpy_(
	    integer *, real *, real *, integer *, real *, integer *);


/*     sgedi computes the determinant and inverse of a matrix */
/*     using the factors computed by sgeco or sgefa. */

/*     on entry */

/*        a       real(lda, n) */
/*                the output from sgeco or sgefa. */

/*        lda     integer */
/*                the leading dimension of the array  a . */

/*        n       integer */
/*                the order of the matrix  a . */

/*        ipvt    integer(n) */
/*                the pivot vector from sgeco or sgefa. */

/*        work    real(n) */
/*                work vector.  contents destroyed. */

/*        job     integer */
/*                = 11   both determinant and inverse. */
/*                = 01   inverse only. */
/*                = 10   determinant only. */

/*     on return */

/*        a       inverse of original matrix if requested. */
/*                otherwise unchanged. */

/*        det     real(2) */
/*                determinant of original matrix if requested. */
/*                otherwise not referenced. */
/*                determinant = det(1) * 10.0**det(2) */
/*                with  1.0 .le. abs(det(1)) .lt. 10.0 */
/*                or  det(1) .eq. 0.0 . */

/*     error condition */

/*        a division by zero will occur if the input factor contains */
/*        a zero on the diagonal and the inverse is requested. */
/*        it will not occur if the subroutines are called correctly */
/*        and if sgeco has set rcond .gt. 0.0 or sgefa has set */
/*        info .eq. 0 . */

/*     linpack. this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas saxpy,sscal,sswap */
/*     fortran abs,mod */

/*     internal variables */



/*     compute determinant */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipvt;
    --det;
    --work;

    /* Function Body */
    if (*job / 10 == 0) {
	goto L70;
    }
    det[1] = 1.f;
    det[2] = 0.f;
    ten = 10.f;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (ipvt[i__] != i__) {
	    det[1] = -det[1];
	}
	det[1] = a[i__ + i__ * a_dim1] * det[1];
/*        ...exit */
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
L70:

/*     compute inverse(u) */

    if (*job % 10 == 0) {
	goto L150;
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	a[k + k * a_dim1] = 1.f / a[k + k * a_dim1];
	t = -a[k + k * a_dim1];
	i__2 = k - 1;
	sscal_(&i__2, &t, &a[k * a_dim1 + 1], &c__1);
	kp1 = k + 1;
	if (*n < kp1) {
	    goto L90;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    t = a[k + j * a_dim1];
	    a[k + j * a_dim1] = 0.f;
	    saxpy_(&k, &t, &a[k * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &
		    c__1);
/* L80: */
	}
L90:
/* L100: */
	;
    }

/*        form inverse(u)*inverse(l) */

    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L140;
    }
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n - kb;
	kp1 = k + 1;
	i__2 = *n;
	for (i__ = kp1; i__ <= i__2; ++i__) {
	    work[i__] = a[i__ + k * a_dim1];
	    a[i__ + k * a_dim1] = 0.f;
/* L110: */
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    t = work[j];
	    saxpy_(n, &t, &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &
		    c__1);
/* L120: */
	}
	l = ipvt[k];
	if (l != k) {
	    sswap_(n, &a[k * a_dim1 + 1], &c__1, &a[l * a_dim1 + 1], &c__1);
	}
/* L130: */
    }
L140:
L150:
    return 0;
} /* sgedi_ */

