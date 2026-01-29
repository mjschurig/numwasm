/* cgedi.f -- translated by f2c (version 20240504).
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

static complex c_b3 = {1.f,0.f};
static integer c__1 = 1;

/* Subroutine */ int cgedi_(complex *a, integer *lda, integer *n, integer *
	ipvt, complex *det, complex *work, integer *job)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2;
    complex q__1, q__2;

    /* Builtin functions */
    double r_imag(complex *);
    void c_div(complex *, complex *, complex *);

    /* Local variables */
    integer i__, j, k, l;
    complex t;
    integer kb, kp1, nm1;
    real ten;
    extern /* Subroutine */ int cscal_(integer *, complex *, complex *, 
	    integer *), cswap_(integer *, complex *, integer *, complex *, 
	    integer *), caxpy_(integer *, complex *, complex *, integer *, 
	    complex *, integer *);


/*     cgedi computes the determinant and inverse of a matrix */
/*     using the factors computed by cgeco or cgefa. */

/*     on entry */

/*        a       complex(lda, n) */
/*                the output from cgeco or cgefa. */

/*        lda     integer */
/*                the leading dimension of the array  a . */

/*        n       integer */
/*                the order of the matrix  a . */

/*        ipvt    integer(n) */
/*                the pivot vector from cgeco or cgefa. */

/*        work    complex(n) */
/*                work vector.  contents destroyed. */

/*        job     integer */
/*                = 11   both determinant and inverse. */
/*                = 01   inverse only. */
/*                = 10   determinant only. */

/*     on return */

/*        a       inverse of original matrix if requested. */
/*                otherwise unchanged. */

/*        det     complex(2) */
/*                determinant of original matrix if requested. */
/*                otherwise not referenced. */
/*                determinant = det(1) * 10.0**det(2) */
/*                with  1.0 .le. cabs1(det(1)) .lt. 10.0 */
/*                or  det(1) .eq. 0.0 . */

/*     error condition */

/*        a division by zero will occur if the input factor contains */
/*        a zero on the diagonal and the inverse is requested. */
/*        it will not occur if the subroutines are called correctly */
/*        and if cgeco has set rcond .gt. 0.0 or cgefa has set */
/*        info .eq. 0 . */

/*     linpack. this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas caxpy,cscal,cswap */
/*     fortran abs,aimag,cmplx,mod,real */

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
    det[1].r = 1.f, det[1].i = 0.f;
    det[2].r = 0.f, det[2].i = 0.f;
    ten = 10.f;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (ipvt[i__] != i__) {
	    q__1.r = -det[1].r, q__1.i = -det[1].i;
	    det[1].r = q__1.r, det[1].i = q__1.i;
	}
	i__2 = i__ + i__ * a_dim1;
	q__1.r = a[i__2].r * det[1].r - a[i__2].i * det[1].i, q__1.i = a[i__2]
		.r * det[1].i + a[i__2].i * det[1].r;
	det[1].r = q__1.r, det[1].i = q__1.i;
/*        ...exit */
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
L70:

/*     compute inverse(u) */

    if (*job % 10 == 0) {
	goto L150;
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k + k * a_dim1;
	c_div(&q__1, &c_b3, &a[k + k * a_dim1]);
	a[i__2].r = q__1.r, a[i__2].i = q__1.i;
	i__2 = k + k * a_dim1;
	q__1.r = -a[i__2].r, q__1.i = -a[i__2].i;
	t.r = q__1.r, t.i = q__1.i;
	i__2 = k - 1;
	cscal_(&i__2, &t, &a[k * a_dim1 + 1], &c__1);
	kp1 = k + 1;
	if (*n < kp1) {
	    goto L90;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    i__3 = k + j * a_dim1;
	    t.r = a[i__3].r, t.i = a[i__3].i;
	    i__3 = k + j * a_dim1;
	    a[i__3].r = 0.f, a[i__3].i = 0.f;
	    caxpy_(&k, &t, &a[k * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &
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
	    i__3 = i__;
	    i__4 = i__ + k * a_dim1;
	    work[i__3].r = a[i__4].r, work[i__3].i = a[i__4].i;
	    i__3 = i__ + k * a_dim1;
	    a[i__3].r = 0.f, a[i__3].i = 0.f;
/* L110: */
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    i__3 = j;
	    t.r = work[i__3].r, t.i = work[i__3].i;
	    caxpy_(n, &t, &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &
		    c__1);
/* L120: */
	}
	l = ipvt[k];
	if (l != k) {
	    cswap_(n, &a[k * a_dim1 + 1], &c__1, &a[l * a_dim1 + 1], &c__1);
	}
/* L130: */
    }
L140:
L150:
    return 0;
} /* cgedi_ */

