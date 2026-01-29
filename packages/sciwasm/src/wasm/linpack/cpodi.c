/* cpodi.f -- translated by f2c (version 20240504).
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

static complex c_b11 = {1.f,0.f};
static integer c__1 = 1;

/* Subroutine */ int cpodi_(complex *a, integer *lda, integer *n, real *det, 
	integer *job)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    real r__1;
    complex q__1;

    /* Builtin functions */
    void c_div(complex *, complex *, complex *), r_cnjg(complex *, complex *);

    /* Local variables */
    integer i__, j, k;
    real s;
    complex t;
    integer jm1, kp1;
    extern /* Subroutine */ int cscal_(integer *, complex *, complex *, 
	    integer *), caxpy_(integer *, complex *, complex *, integer *, 
	    complex *, integer *);


/*     cpodi computes the determinant and inverse of a certain */
/*     complex hermitian positive definite matrix (see below) */
/*     using the factors computed by cpoco, cpofa or cqrdc. */

/*     on entry */

/*        a       complex(lda, n) */
/*                the output  a  from cpoco or cpofa */
/*                or the output  x  from cqrdc. */

/*        lda     integer */
/*                the leading dimension of the array  a . */

/*        n       integer */
/*                the order of the matrix  a . */

/*        job     integer */
/*                = 11   both determinant and inverse. */
/*                = 01   inverse only. */
/*                = 10   determinant only. */

/*     on return */

/*        a       if cpoco or cpofa was used to factor  a  then */
/*                cpodi produces the upper half of inverse(a) . */
/*                if cqrdc was used to decompose  x  then */
/*                cpodi produces the upper half of inverse(ctrans(x)*x) */
/*                where ctrans(x) is the conjugate transpose. */
/*                elements of  a  below the diagonal are unchanged. */
/*                if the units digit of job is zero,  a  is unchanged. */

/*        det     real(2) */
/*                determinant of  a  or of  ctrans(x)*x  if requested. */
/*                otherwise not referenced. */
/*                determinant = det(1) * 10.0**det(2) */
/*                with  1.0 .le. det(1) .lt. 10.0 */
/*                or  det(1) .eq. 0.0 . */

/*     error condition */

/*        a division by zero will occur if the input factor contains */
/*        a zero on the diagonal and the inverse is requested. */
/*        it will not occur if the subroutines are called correctly */
/*        and if cpoco or cpofa has set info .eq. 0 . */

/*     linpack.  this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas caxpy,cscal */
/*     fortran conjg,mod,real */

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
    det[1] = 1.f;
    det[2] = 0.f;
    s = 10.f;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__ + i__ * a_dim1;
/* Computing 2nd power */
	r__1 = a[i__2].r;
	det[1] = r__1 * r__1 * det[1];
/*        ...exit */
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
L70:

/*     compute inverse(r) */

    if (*job % 10 == 0) {
	goto L140;
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k + k * a_dim1;
	c_div(&q__1, &c_b11, &a[k + k * a_dim1]);
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

/*        form  inverse(r) * ctrans(inverse(r)) */

    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L120;
	}
	i__2 = jm1;
	for (k = 1; k <= i__2; ++k) {
	    r_cnjg(&q__1, &a[k + j * a_dim1]);
	    t.r = q__1.r, t.i = q__1.i;
	    caxpy_(&k, &t, &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &
		    c__1);
/* L110: */
	}
L120:
	r_cnjg(&q__1, &a[j + j * a_dim1]);
	t.r = q__1.r, t.i = q__1.i;
	cscal_(&j, &t, &a[j * a_dim1 + 1], &c__1);
/* L130: */
    }
L140:
    return 0;
} /* cpodi_ */

