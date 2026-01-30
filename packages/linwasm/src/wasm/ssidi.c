/* ssidi.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int ssidi_(real *a, integer *lda, integer *n, integer *kpvt, 
	real *det, integer *inert, real *work, integer *job)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    real r__1;

    /* Local variables */
    real d__;
    integer j, k;
    real t, ak;
    integer jb, ks, km1;
    real ten, akp1, temp;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    real akkp1;
    logical nodet;
    integer kstep;
    logical noert, noinv;
    extern /* Subroutine */ int scopy_(integer *, real *, integer *, real *, 
	    integer *), sswap_(integer *, real *, integer *, real *, integer *
	    ), saxpy_(integer *, real *, real *, integer *, real *, integer *)
	    ;


/*     ssidi computes the determinant, inertia and inverse */
/*     of a real symmetric matrix using the factors from ssifa. */

/*     on entry */

/*        a       real(lda,n) */
/*                the output from ssifa. */

/*        lda     integer */
/*                the leading dimension of the array a. */

/*        n       integer */
/*                the order of the matrix a. */

/*        kpvt    integer(n) */
/*                the pivot vector from ssifa. */

/*        work    real(n) */
/*                work vector.  contents destroyed. */

/*        job     integer */
/*                job has the decimal expansion  abc  where */
/*                   if  c .ne. 0, the inverse is computed, */
/*                   if  b .ne. 0, the determinant is computed, */
/*                   if  a .ne. 0, the inertia is computed. */

/*                for example, job = 111  gives all three. */

/*     on return */

/*        variables not requested by job are not used. */

/*        a      contains the upper triangle of the inverse of */
/*               the original matrix.  the strict lower triangle */
/*               is never referenced. */

/*        det    real(2) */
/*               determinant of original matrix. */
/*               determinant = det(1) * 10.0**det(2) */
/*               with 1.0 .le. abs(det(1)) .lt. 10.0 */
/*               or det(1) = 0.0. */

/*        inert  integer(3) */
/*               the inertia of the original matrix. */
/*               inert(1)  =  number of positive eigenvalues. */
/*               inert(2)  =  number of negative eigenvalues. */
/*               inert(3)  =  number of zero eigenvalues. */

/*     error condition */

/*        a division by zero may occur if the inverse is requested */
/*        and  ssico  has set rcond .eq. 0.0 */
/*        or  ssifa  has set  info .ne. 0 . */

/*     linpack. this version dated 08/14/78 . */
/*     james bunch, univ. calif. san diego, argonne nat. lab */

/*     subroutines and functions */

/*     blas saxpy,scopy,sdot,sswap */
/*     fortran abs,iabs,mod */

/*     internal variables. */


    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --kpvt;
    --det;
    --inert;
    --work;

    /* Function Body */
    noinv = *job % 10 == 0;
    nodet = *job % 100 / 10 == 0;
    noert = *job % 1000 / 100 == 0;

    if (nodet && noert) {
	goto L140;
    }
    if (noert) {
	goto L10;
    }
    inert[1] = 0;
    inert[2] = 0;
    inert[3] = 0;
L10:
    if (nodet) {
	goto L20;
    }
    det[1] = 1.f;
    det[2] = 0.f;
    ten = 10.f;
L20:
    t = 0.f;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	d__ = a[k + k * a_dim1];

/*           check if 1 by 1 */

	if (kpvt[k] > 0) {
	    goto L50;
	}

/*              2 by 2 block */
/*              use det (d  s)  =  (d/t * c - t) * t  ,  t = abs(s) */
/*                      (s  c) */
/*              to avoid underflow/overflow troubles. */
/*              take two passes through scaling.  use  t  for flag. */

	if (t != 0.f) {
	    goto L30;
	}
	t = (r__1 = a[k + (k + 1) * a_dim1], dabs(r__1));
	d__ = d__ / t * a[k + 1 + (k + 1) * a_dim1] - t;
	goto L40;
L30:
	d__ = t;
	t = 0.f;
L40:
L50:

	if (noert) {
	    goto L60;
	}
	if (d__ > 0.f) {
	    ++inert[1];
	}
	if (d__ < 0.f) {
	    ++inert[2];
	}
	if (d__ == 0.f) {
	    ++inert[3];
	}
L60:

	if (nodet) {
	    goto L120;
	}
	det[1] = d__ * det[1];
	if (det[1] == 0.f) {
	    goto L110;
	}
L70:
	if (dabs(det[1]) >= 1.f) {
	    goto L80;
	}
	det[1] = ten * det[1];
	det[2] += -1.f;
	goto L70;
L80:
L90:
	if (dabs(det[1]) < ten) {
	    goto L100;
	}
	det[1] /= ten;
	det[2] += 1.f;
	goto L90;
L100:
L110:
L120:
/* L130: */
	;
    }
L140:

/*     compute inverse(a) */

    if (noinv) {
	goto L270;
    }
    k = 1;
L150:
    if (k > *n) {
	goto L260;
    }
    km1 = k - 1;
    if (kpvt[k] < 0) {
	goto L180;
    }

/*              1 by 1 */

    a[k + k * a_dim1] = 1.f / a[k + k * a_dim1];
    if (km1 < 1) {
	goto L170;
    }
    scopy_(&km1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	a[j + k * a_dim1] = sdot_(&j, &a[j * a_dim1 + 1], &c__1, &work[1], &
		c__1);
	i__2 = j - 1;
	saxpy_(&i__2, &work[j], &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1],
		 &c__1);
/* L160: */
    }
    a[k + k * a_dim1] += sdot_(&km1, &work[1], &c__1, &a[k * a_dim1 + 1], &
	    c__1);
L170:
    kstep = 1;
    goto L220;
L180:

/*              2 by 2 */

    t = (r__1 = a[k + (k + 1) * a_dim1], dabs(r__1));
    ak = a[k + k * a_dim1] / t;
    akp1 = a[k + 1 + (k + 1) * a_dim1] / t;
    akkp1 = a[k + (k + 1) * a_dim1] / t;
    d__ = t * (ak * akp1 - 1.f);
    a[k + k * a_dim1] = akp1 / d__;
    a[k + 1 + (k + 1) * a_dim1] = ak / d__;
    a[k + (k + 1) * a_dim1] = -akkp1 / d__;
    if (km1 < 1) {
	goto L210;
    }
    scopy_(&km1, &a[(k + 1) * a_dim1 + 1], &c__1, &work[1], &c__1);
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	a[j + (k + 1) * a_dim1] = sdot_(&j, &a[j * a_dim1 + 1], &c__1, &work[
		1], &c__1);
	i__2 = j - 1;
	saxpy_(&i__2, &work[j], &a[j * a_dim1 + 1], &c__1, &a[(k + 1) * 
		a_dim1 + 1], &c__1);
/* L190: */
    }
    a[k + 1 + (k + 1) * a_dim1] += sdot_(&km1, &work[1], &c__1, &a[(k + 1) * 
	    a_dim1 + 1], &c__1);
    a[k + (k + 1) * a_dim1] += sdot_(&km1, &a[k * a_dim1 + 1], &c__1, &a[(k + 
	    1) * a_dim1 + 1], &c__1);
    scopy_(&km1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	a[j + k * a_dim1] = sdot_(&j, &a[j * a_dim1 + 1], &c__1, &work[1], &
		c__1);
	i__2 = j - 1;
	saxpy_(&i__2, &work[j], &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1],
		 &c__1);
/* L200: */
    }
    a[k + k * a_dim1] += sdot_(&km1, &work[1], &c__1, &a[k * a_dim1 + 1], &
	    c__1);
L210:
    kstep = 2;
L220:

/*           swap */

    ks = (i__1 = kpvt[k], abs(i__1));
    if (ks == k) {
	goto L250;
    }
    sswap_(&ks, &a[ks * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
    i__1 = k;
    for (jb = ks; jb <= i__1; ++jb) {
	j = k + ks - jb;
	temp = a[j + k * a_dim1];
	a[j + k * a_dim1] = a[ks + j * a_dim1];
	a[ks + j * a_dim1] = temp;
/* L230: */
    }
    if (kstep == 1) {
	goto L240;
    }
    temp = a[ks + (k + 1) * a_dim1];
    a[ks + (k + 1) * a_dim1] = a[k + (k + 1) * a_dim1];
    a[k + (k + 1) * a_dim1] = temp;
L240:
L250:
    k += kstep;
    goto L150;
L260:
L270:
    return 0;
} /* ssidi_ */

