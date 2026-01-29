/* zhidi.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int zhidi_(doublecomplex *a, integer *lda, integer *n, 
	integer *kpvt, doublereal *det, integer *inert, doublecomplex *work, 
	integer *job)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    doublereal d__1;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    double z_abs(doublecomplex *);
    void d_cnjg(doublecomplex *, doublecomplex *);

    /* Local variables */
    doublereal d__;
    integer j, k;
    doublereal t, ak;
    integer jb, ks, km1;
    doublereal ten, akp1;
    doublecomplex temp, akkp1;
    logical nodet;
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    integer kstep;
    logical noert, noinv;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zswap_(integer *, doublecomplex *, 
	    integer *, doublecomplex *, integer *), zaxpy_(integer *, 
	    doublecomplex *, doublecomplex *, integer *, doublecomplex *, 
	    integer *);


/*     zhidi computes the determinant, inertia and inverse */
/*     of a complex*16 hermitian matrix using the factors from zhifa. */

/*     on entry */

/*        a       complex*16(lda,n) */
/*                the output from zhifa. */

/*        lda     integer */
/*                the leading dimension of the array a. */

/*        n       integer */
/*                the order of the matrix a. */

/*        kpvt    integer(n) */
/*                the pivot vector from zhifa. */

/*        work    complex*16(n) */
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

/*        det    double precision(2) */
/*               determinant of original matrix. */
/*               determinant = det(1) * 10.0**det(2) */
/*               with 1.0 .le. dabs(det(1)) .lt. 10.0 */
/*               or det(1) = 0.0. */

/*        inert  integer(3) */
/*               the inertia of the original matrix. */
/*               inert(1)  =  number of positive eigenvalues. */
/*               inert(2)  =  number of negative eigenvalues. */
/*               inert(3)  =  number of zero eigenvalues. */

/*     error condition */

/*        a division by zero may occur if the inverse is requested */
/*        and  zhico  has set rcond .eq. 0.0 */
/*        or  zhifa  has set  info .ne. 0 . */

/*     linpack. this version dated 08/14/78 . */
/*     james bunch, univ. calif. san diego, argonne nat. lab */

/*     subroutines and functions */

/*     blas zaxpy,zcopy,zdotc,zswap */
/*     fortran dabs,cdabs,dcmplx,dconjg,iabs,mod */

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
    det[1] = 1.;
    det[2] = 0.;
    ten = 10.;
L20:
    t = 0.;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k + k * a_dim1;
	d__ = a[i__2].r;

/*           check if 1 by 1 */

	if (kpvt[k] > 0) {
	    goto L50;
	}

/*              2 by 2 block */
/*              use det (d  s)  =  (d/t * c - t) * t  ,  t = cdabs(s) */
/*                      (s  c) */
/*              to avoid underflow/overflow troubles. */
/*              take two passes through scaling.  use  t  for flag. */

	if (t != 0.) {
	    goto L30;
	}
	t = z_abs(&a[k + (k + 1) * a_dim1]);
	i__2 = k + 1 + (k + 1) * a_dim1;
	d__ = d__ / t * a[i__2].r - t;
	goto L40;
L30:
	d__ = t;
	t = 0.;
L40:
L50:

	if (noert) {
	    goto L60;
	}
	if (d__ > 0.) {
	    ++inert[1];
	}
	if (d__ < 0.) {
	    ++inert[2];
	}
	if (d__ == 0.) {
	    ++inert[3];
	}
L60:

	if (nodet) {
	    goto L120;
	}
	det[1] = d__ * det[1];
	if (det[1] == 0.) {
	    goto L110;
	}
L70:
	if (abs(det[1]) >= 1.) {
	    goto L80;
	}
	det[1] = ten * det[1];
	det[2] += -1.;
	goto L70;
L80:
L90:
	if (abs(det[1]) < ten) {
	    goto L100;
	}
	det[1] /= ten;
	det[2] += 1.;
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

    i__1 = k + k * a_dim1;
    i__2 = k + k * a_dim1;
    d__1 = 1. / a[i__2].r;
    z__1.r = d__1, z__1.i = 0.;
    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
    if (km1 < 1) {
	goto L170;
    }
    zcopy_(&km1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j + k * a_dim1;
	zdotc_(&z__1, &j, &a[j * a_dim1 + 1], &c__1, &work[1], &c__1);
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
	i__2 = j - 1;
	zaxpy_(&i__2, &work[j], &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1],
		 &c__1);
/* L160: */
    }
    zdotc_(&z__2, &km1, &work[1], &c__1, &a[k * a_dim1 + 1], &c__1);
    z__1.r = z__2.r, z__1.i = z__2.i;
    i__1 = k + k * a_dim1;
    i__2 = k + k * a_dim1;
    d__1 = z__1.r;
    z__4.r = d__1, z__4.i = 0.;
    z__3.r = a[i__2].r + z__4.r, z__3.i = a[i__2].i + z__4.i;
    a[i__1].r = z__3.r, a[i__1].i = z__3.i;
L170:
    kstep = 1;
    goto L220;
L180:

/*              2 by 2 */

    t = z_abs(&a[k + (k + 1) * a_dim1]);
    i__1 = k + k * a_dim1;
    ak = a[i__1].r / t;
    i__1 = k + 1 + (k + 1) * a_dim1;
    akp1 = a[i__1].r / t;
    i__1 = k + (k + 1) * a_dim1;
    z__1.r = a[i__1].r / t, z__1.i = a[i__1].i / t;
    akkp1.r = z__1.r, akkp1.i = z__1.i;
    d__ = t * (ak * akp1 - 1.);
    i__1 = k + k * a_dim1;
    d__1 = akp1 / d__;
    z__1.r = d__1, z__1.i = 0.;
    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
    i__1 = k + 1 + (k + 1) * a_dim1;
    d__1 = ak / d__;
    z__1.r = d__1, z__1.i = 0.;
    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
    i__1 = k + (k + 1) * a_dim1;
    z__2.r = -akkp1.r, z__2.i = -akkp1.i;
    z__1.r = z__2.r / d__, z__1.i = z__2.i / d__;
    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
    if (km1 < 1) {
	goto L210;
    }
    zcopy_(&km1, &a[(k + 1) * a_dim1 + 1], &c__1, &work[1], &c__1);
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j + (k + 1) * a_dim1;
	zdotc_(&z__1, &j, &a[j * a_dim1 + 1], &c__1, &work[1], &c__1);
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
	i__2 = j - 1;
	zaxpy_(&i__2, &work[j], &a[j * a_dim1 + 1], &c__1, &a[(k + 1) * 
		a_dim1 + 1], &c__1);
/* L190: */
    }
    zdotc_(&z__2, &km1, &work[1], &c__1, &a[(k + 1) * a_dim1 + 1], &c__1);
    z__1.r = z__2.r, z__1.i = z__2.i;
    i__1 = k + 1 + (k + 1) * a_dim1;
    i__2 = k + 1 + (k + 1) * a_dim1;
    d__1 = z__1.r;
    z__4.r = d__1, z__4.i = 0.;
    z__3.r = a[i__2].r + z__4.r, z__3.i = a[i__2].i + z__4.i;
    a[i__1].r = z__3.r, a[i__1].i = z__3.i;
    i__1 = k + (k + 1) * a_dim1;
    i__2 = k + (k + 1) * a_dim1;
    zdotc_(&z__2, &km1, &a[k * a_dim1 + 1], &c__1, &a[(k + 1) * a_dim1 + 1], &
	    c__1);
    z__1.r = a[i__2].r + z__2.r, z__1.i = a[i__2].i + z__2.i;
    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
    zcopy_(&km1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j + k * a_dim1;
	zdotc_(&z__1, &j, &a[j * a_dim1 + 1], &c__1, &work[1], &c__1);
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
	i__2 = j - 1;
	zaxpy_(&i__2, &work[j], &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1],
		 &c__1);
/* L200: */
    }
    zdotc_(&z__2, &km1, &work[1], &c__1, &a[k * a_dim1 + 1], &c__1);
    z__1.r = z__2.r, z__1.i = z__2.i;
    i__1 = k + k * a_dim1;
    i__2 = k + k * a_dim1;
    d__1 = z__1.r;
    z__4.r = d__1, z__4.i = 0.;
    z__3.r = a[i__2].r + z__4.r, z__3.i = a[i__2].i + z__4.i;
    a[i__1].r = z__3.r, a[i__1].i = z__3.i;
L210:
    kstep = 2;
L220:

/*           swap */

    ks = (i__1 = kpvt[k], abs(i__1));
    if (ks == k) {
	goto L250;
    }
    zswap_(&ks, &a[ks * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
    i__1 = k;
    for (jb = ks; jb <= i__1; ++jb) {
	j = k + ks - jb;
	d_cnjg(&z__1, &a[j + k * a_dim1]);
	temp.r = z__1.r, temp.i = z__1.i;
	i__2 = j + k * a_dim1;
	d_cnjg(&z__1, &a[ks + j * a_dim1]);
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
	i__2 = ks + j * a_dim1;
	a[i__2].r = temp.r, a[i__2].i = temp.i;
/* L230: */
    }
    if (kstep == 1) {
	goto L240;
    }
    i__1 = ks + (k + 1) * a_dim1;
    temp.r = a[i__1].r, temp.i = a[i__1].i;
    i__1 = ks + (k + 1) * a_dim1;
    i__2 = k + (k + 1) * a_dim1;
    a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
    i__1 = k + (k + 1) * a_dim1;
    a[i__1].r = temp.r, a[i__1].i = temp.i;
L240:
L250:
    k += kstep;
    goto L150;
L260:
L270:
    return 0;
} /* zhidi_ */

