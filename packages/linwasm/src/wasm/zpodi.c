/* zpodi.f -- translated by f2c (version 20240504).
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

static doublecomplex c_b11 = {1.,0.};
static integer c__1 = 1;

/* Subroutine */ int zpodi_(doublecomplex *a, integer *lda, integer *n, 
	doublereal *det, integer *job)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1;
    doublecomplex z__1;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    integer i__, j, k;
    doublereal s;
    doublecomplex t;
    integer jm1, kp1;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);


/*     zpodi computes the determinant and inverse of a certain */
/*     complex*16 hermitian positive definite matrix (see below) */
/*     using the factors computed by zpoco, zpofa or zqrdc. */

/*     on entry */

/*        a       complex*16(lda, n) */
/*                the output  a  from zpoco or zpofa */
/*                or the output  x  from zqrdc. */

/*        lda     integer */
/*                the leading dimension of the array  a . */

/*        n       integer */
/*                the order of the matrix  a . */

/*        job     integer */
/*                = 11   both determinant and inverse. */
/*                = 01   inverse only. */
/*                = 10   determinant only. */

/*     on return */

/*        a       if zpoco or zpofa was used to factor  a  then */
/*                zpodi produces the upper half of inverse(a) . */
/*                if zqrdc was used to decompose  x  then */
/*                zpodi produces the upper half of inverse(ctrans(x)*x) */
/*                where ctrans(x) is the conjugate transpose. */
/*                elements of  a  below the diagonal are unchanged. */
/*                if the units digit of job is zero,  a  is unchanged. */

/*        det     double precision(2) */
/*                determinant of  a  or of  ctrans(x)*x  if requested. */
/*                otherwise not referenced. */
/*                determinant = det(1) * 10.0**det(2) */
/*                with  1.0 .le. det(1) .lt. 10.0 */
/*                or  det(1) .eq. 0.0 . */

/*     error condition */

/*        a division by zero will occur if the input factor contains */
/*        a zero on the diagonal and the inverse is requested. */
/*        it will not occur if the subroutines are called correctly */
/*        and if zpoco or zpofa has set info .eq. 0 . */

/*     linpack.  this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas zaxpy,zscal */
/*     fortran dconjg,mod */

/*     internal variables */


/*     compute determinant */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --det;

    /* Function Body */
    if (*job / 10 == 0) {
	goto L70;
    }
    det[1] = 1.;
    det[2] = 0.;
    s = 10.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__ + i__ * a_dim1;
/* Computing 2nd power */
	d__1 = a[i__2].r;
	det[1] = d__1 * d__1 * det[1];
/*        ...exit */
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
L70:

/*     compute inverse(r) */

    if (*job % 10 == 0) {
	goto L140;
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k + k * a_dim1;
	z_div(&z__1, &c_b11, &a[k + k * a_dim1]);
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
	i__2 = k + k * a_dim1;
	z__1.r = -a[i__2].r, z__1.i = -a[i__2].i;
	t.r = z__1.r, t.i = z__1.i;
	i__2 = k - 1;
	zscal_(&i__2, &t, &a[k * a_dim1 + 1], &c__1);
	kp1 = k + 1;
	if (*n < kp1) {
	    goto L90;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    i__3 = k + j * a_dim1;
	    t.r = a[i__3].r, t.i = a[i__3].i;
	    i__3 = k + j * a_dim1;
	    a[i__3].r = 0., a[i__3].i = 0.;
	    zaxpy_(&k, &t, &a[k * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &
		    c__1);
/* L80: */
	}
L90:
/* L100: */
	;
    }

/*        form  inverse(r) * ctrans(inverse(r)) */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L120;
	}
	i__2 = jm1;
	for (k = 1; k <= i__2; ++k) {
	    d_cnjg(&z__1, &a[k + j * a_dim1]);
	    t.r = z__1.r, t.i = z__1.i;
	    zaxpy_(&k, &t, &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &
		    c__1);
/* L110: */
	}
L120:
	d_cnjg(&z__1, &a[j + j * a_dim1]);
	t.r = z__1.r, t.i = z__1.i;
	zscal_(&j, &t, &a[j * a_dim1 + 1], &c__1);
/* L130: */
    }
L140:
    return 0;
} /* zpodi_ */

