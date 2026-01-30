/* cppdi.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int cppdi_(complex *ap, integer *n, real *det, integer *job)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1;
    complex q__1;

    /* Builtin functions */
    void c_div(complex *, complex *, complex *), r_cnjg(complex *, complex *);

    /* Local variables */
    integer i__, j, k;
    real s;
    complex t;
    integer j1, k1, ii, jj, kj, kk, jm1, kp1;
    extern /* Subroutine */ int cscal_(integer *, complex *, complex *, 
	    integer *), caxpy_(integer *, complex *, complex *, integer *, 
	    complex *, integer *);


/*     cppdi computes the determinant and inverse */
/*     of a complex hermitian positive definite matrix */
/*     using the factors computed by cppco or cppfa . */

/*     on entry */

/*        ap      complex (n*(n+1)/2) */
/*                the output from cppco or cppfa. */

/*        n       integer */
/*                the order of the matrix  a . */

/*        job     integer */
/*                = 11   both determinant and inverse. */
/*                = 01   inverse only. */
/*                = 10   determinant only. */

/*     on return */

/*        ap      the upper triangular half of the inverse . */
/*                the strict lower triangle is unaltered. */

/*        det     real(2) */
/*                determinant of original matrix if requested. */
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
    --det;
    --ap;

    /* Function Body */
    if (*job / 10 == 0) {
	goto L70;
    }
    det[1] = 1.f;
    det[2] = 0.f;
    s = 10.f;
    ii = 0;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii += i__;
	i__2 = ii;
/* Computing 2nd power */
	r__1 = ap[i__2].r;
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
    kk = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	k1 = kk + 1;
	kk += k;
	i__2 = kk;
	c_div(&q__1, &c_b11, &ap[kk]);
	ap[i__2].r = q__1.r, ap[i__2].i = q__1.i;
	i__2 = kk;
	q__1.r = -ap[i__2].r, q__1.i = -ap[i__2].i;
	t.r = q__1.r, t.i = q__1.i;
	i__2 = k - 1;
	cscal_(&i__2, &t, &ap[k1], &c__1);
	kp1 = k + 1;
	j1 = kk + 1;
	kj = kk + k;
	if (*n < kp1) {
	    goto L90;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    i__3 = kj;
	    t.r = ap[i__3].r, t.i = ap[i__3].i;
	    i__3 = kj;
	    ap[i__3].r = 0.f, ap[i__3].i = 0.f;
	    caxpy_(&k, &t, &ap[k1], &c__1, &ap[j1], &c__1);
	    j1 += j;
	    kj += j;
/* L80: */
	}
L90:
/* L100: */
	;
    }

/*        form  inverse(r) * ctrans(inverse(r)) */

    jj = 0;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	j1 = jj + 1;
	jj += j;
	jm1 = j - 1;
	k1 = 1;
	kj = j1;
	if (jm1 < 1) {
	    goto L120;
	}
	i__2 = jm1;
	for (k = 1; k <= i__2; ++k) {
	    r_cnjg(&q__1, &ap[kj]);
	    t.r = q__1.r, t.i = q__1.i;
	    caxpy_(&k, &t, &ap[j1], &c__1, &ap[k1], &c__1);
	    k1 += k;
	    ++kj;
/* L110: */
	}
L120:
	r_cnjg(&q__1, &ap[jj]);
	t.r = q__1.r, t.i = q__1.i;
	cscal_(&j, &t, &ap[j1], &c__1);
/* L130: */
    }
L140:
    return 0;
} /* cppdi_ */

