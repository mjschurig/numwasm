/* zsidi.f -- translated by f2c (version 20240504).
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

static doublecomplex c_b3 = {1.,0.};
static integer c__1 = 1;

/* Subroutine */ int zsidi_(doublecomplex *a, integer *lda, integer *n, 
	integer *kpvt, doublecomplex *det, doublecomplex *work, integer *job)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    doublecomplex d__;
    integer j, k;
    doublecomplex t, ak;
    integer jb, ks, km1;
    doublereal ten;
    doublecomplex akp1, temp, akkp1;
    logical nodet;
    integer kstep;
    logical noinv;
    extern /* Subroutine */ int zcopy_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *);
    extern /* Double Complex */ VOID zdotu_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);


/*     zsidi computes the determinant and inverse */
/*     of a complex*16 symmetric matrix using the factors from zsifa. */

/*     on entry */

/*        a       complex*16(lda,n) */
/*                the output from zsifa. */

/*        lda     integer */
/*                the leading dimension of the array a. */

/*        n       integer */
/*                the order of the matrix a. */

/*        kpvt    integer(n) */
/*                the pivot vector from zsifa. */

/*        work    complex*16(n) */
/*                work vector.  contents destroyed. */

/*        job     integer */
/*                job has the decimal expansion  ab  where */
/*                   if  b .ne. 0, the inverse is computed, */
/*                   if  a .ne. 0, the determinant is computed, */

/*                for example, job = 11  gives both. */

/*     on return */

/*        variables not requested by job are not used. */

/*        a      contains the upper triangle of the inverse of */
/*               the original matrix.  the strict lower triangle */
/*               is never referenced. */

/*        det    complex*16(2) */
/*               determinant of original matrix. */
/*               determinant = det(1) * 10.0**det(2) */
/*               with 1.0 .le. dabs(det(1)) .lt. 10.0 */
/*               or det(1) = 0.0. */

/*     error condition */

/*        a division by zero may occur if the inverse is requested */
/*        and  zsico  has set rcond .eq. 0.0 */
/*        or  zsifa  has set  info .ne. 0 . */

/*     linpack. this version dated 08/14/78 . */
/*     james bunch, univ. calif. san diego, argonne nat. lab */

/*     subroutines and functions */

/*     blas zaxpy,zcopy,zdotu,zswap */
/*     fortran dabs,dcmplx,iabs,mod */

/*     internal variables. */



    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --kpvt;
    --det;
    --work;

    /* Function Body */
    noinv = *job % 10 == 0;
    nodet = *job % 100 / 10 == 0;

    if (nodet) {
	goto L100;
    }
    det[1].r = 1., det[1].i = 0.;
    det[2].r = 0., det[2].i = 0.;
    ten = 10.;
    t.r = 0., t.i = 0.;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k + k * a_dim1;
	d__.r = a[i__2].r, d__.i = a[i__2].i;

/*           check if 1 by 1 */

	if (kpvt[k] > 0) {
	    goto L30;
	}

/*              2 by 2 block */
/*              use det (d  t)  =  (d/t * c - t) * t */
/*                      (t  c) */
/*              to avoid underflow/overflow troubles. */
/*              take two passes through scaling.  use  t  for flag. */

	z__1.r = t.r * 0. - t.i * -1., z__1.i = t.i * 0. + t.r * -1.;
	if ((d__1 = t.r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) != 0.) {
	    goto L10;
	}
	i__2 = k + (k + 1) * a_dim1;
	t.r = a[i__2].r, t.i = a[i__2].i;
	z_div(&z__3, &d__, &t);
	i__2 = k + 1 + (k + 1) * a_dim1;
	z__2.r = z__3.r * a[i__2].r - z__3.i * a[i__2].i, z__2.i = z__3.r * a[
		i__2].i + z__3.i * a[i__2].r;
	z__1.r = z__2.r - t.r, z__1.i = z__2.i - t.i;
	d__.r = z__1.r, d__.i = z__1.i;
	goto L20;
L10:
	d__.r = t.r, d__.i = t.i;
	t.r = 0., t.i = 0.;
L20:
L30:

	z__1.r = d__.r * det[1].r - d__.i * det[1].i, z__1.i = d__.r * det[1]
		.i + d__.i * det[1].r;
	det[1].r = z__1.r, det[1].i = z__1.i;
	z__1.r = det[1].r * 0. - det[1].i * -1., z__1.i = det[1].i * 0. + det[
		1].r * -1.;
	if ((d__1 = det[1].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) == 0.) {
	    goto L80;
	}
L40:
	z__1.r = det[1].r * 0. - det[1].i * -1., z__1.i = det[1].i * 0. + det[
		1].r * -1.;
	if ((d__1 = det[1].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) >= 1.) {
	    goto L50;
	}
	z__2.r = ten, z__2.i = 0.;
	z__1.r = z__2.r * det[1].r - z__2.i * det[1].i, z__1.i = z__2.r * det[
		1].i + z__2.i * det[1].r;
	det[1].r = z__1.r, det[1].i = z__1.i;
	z__1.r = det[2].r - 1., z__1.i = det[2].i - 0.;
	det[2].r = z__1.r, det[2].i = z__1.i;
	goto L40;
L50:
L60:
	z__1.r = det[1].r * 0. - det[1].i * -1., z__1.i = det[1].i * 0. + det[
		1].r * -1.;
	if ((d__1 = det[1].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) < ten) {
	    goto L70;
	}
	z__2.r = ten, z__2.i = 0.;
	z_div(&z__1, &det[1], &z__2);
	det[1].r = z__1.r, det[1].i = z__1.i;
	z__1.r = det[2].r + 1., z__1.i = det[2].i + 0.;
	det[2].r = z__1.r, det[2].i = z__1.i;
	goto L60;
L70:
L80:
/* L90: */
	;
    }
L100:

/*     compute inverse(a) */

    if (noinv) {
	goto L230;
    }
    k = 1;
L110:
    if (k > *n) {
	goto L220;
    }
    km1 = k - 1;
    if (kpvt[k] < 0) {
	goto L140;
    }

/*              1 by 1 */

    i__1 = k + k * a_dim1;
    z_div(&z__1, &c_b3, &a[k + k * a_dim1]);
    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
    if (km1 < 1) {
	goto L130;
    }
    zcopy_(&km1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j + k * a_dim1;
	zdotu_(&z__1, &j, &a[j * a_dim1 + 1], &c__1, &work[1], &c__1);
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
	i__2 = j - 1;
	zaxpy_(&i__2, &work[j], &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1],
		 &c__1);
/* L120: */
    }
    i__1 = k + k * a_dim1;
    i__2 = k + k * a_dim1;
    zdotu_(&z__2, &km1, &work[1], &c__1, &a[k * a_dim1 + 1], &c__1);
    z__1.r = a[i__2].r + z__2.r, z__1.i = a[i__2].i + z__2.i;
    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
L130:
    kstep = 1;
    goto L180;
L140:

/*              2 by 2 */

    i__1 = k + (k + 1) * a_dim1;
    t.r = a[i__1].r, t.i = a[i__1].i;
    z_div(&z__1, &a[k + k * a_dim1], &t);
    ak.r = z__1.r, ak.i = z__1.i;
    z_div(&z__1, &a[k + 1 + (k + 1) * a_dim1], &t);
    akp1.r = z__1.r, akp1.i = z__1.i;
    z_div(&z__1, &a[k + (k + 1) * a_dim1], &t);
    akkp1.r = z__1.r, akkp1.i = z__1.i;
    z__3.r = ak.r * akp1.r - ak.i * akp1.i, z__3.i = ak.r * akp1.i + ak.i * 
	    akp1.r;
    z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
    z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * z__2.i + t.i * 
	    z__2.r;
    d__.r = z__1.r, d__.i = z__1.i;
    i__1 = k + k * a_dim1;
    z_div(&z__1, &akp1, &d__);
    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
    i__1 = k + 1 + (k + 1) * a_dim1;
    z_div(&z__1, &ak, &d__);
    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
    i__1 = k + (k + 1) * a_dim1;
    z__2.r = -akkp1.r, z__2.i = -akkp1.i;
    z_div(&z__1, &z__2, &d__);
    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
    if (km1 < 1) {
	goto L170;
    }
    zcopy_(&km1, &a[(k + 1) * a_dim1 + 1], &c__1, &work[1], &c__1);
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j + (k + 1) * a_dim1;
	zdotu_(&z__1, &j, &a[j * a_dim1 + 1], &c__1, &work[1], &c__1);
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
	i__2 = j - 1;
	zaxpy_(&i__2, &work[j], &a[j * a_dim1 + 1], &c__1, &a[(k + 1) * 
		a_dim1 + 1], &c__1);
/* L150: */
    }
    i__1 = k + 1 + (k + 1) * a_dim1;
    i__2 = k + 1 + (k + 1) * a_dim1;
    zdotu_(&z__2, &km1, &work[1], &c__1, &a[(k + 1) * a_dim1 + 1], &c__1);
    z__1.r = a[i__2].r + z__2.r, z__1.i = a[i__2].i + z__2.i;
    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
    i__1 = k + (k + 1) * a_dim1;
    i__2 = k + (k + 1) * a_dim1;
    zdotu_(&z__2, &km1, &a[k * a_dim1 + 1], &c__1, &a[(k + 1) * a_dim1 + 1], &
	    c__1);
    z__1.r = a[i__2].r + z__2.r, z__1.i = a[i__2].i + z__2.i;
    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
    zcopy_(&km1, &a[k * a_dim1 + 1], &c__1, &work[1], &c__1);
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j + k * a_dim1;
	zdotu_(&z__1, &j, &a[j * a_dim1 + 1], &c__1, &work[1], &c__1);
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
	i__2 = j - 1;
	zaxpy_(&i__2, &work[j], &a[j * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1],
		 &c__1);
/* L160: */
    }
    i__1 = k + k * a_dim1;
    i__2 = k + k * a_dim1;
    zdotu_(&z__2, &km1, &work[1], &c__1, &a[k * a_dim1 + 1], &c__1);
    z__1.r = a[i__2].r + z__2.r, z__1.i = a[i__2].i + z__2.i;
    a[i__1].r = z__1.r, a[i__1].i = z__1.i;
L170:
    kstep = 2;
L180:

/*           swap */

    ks = (i__1 = kpvt[k], abs(i__1));
    if (ks == k) {
	goto L210;
    }
    zswap_(&ks, &a[ks * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
    i__1 = k;
    for (jb = ks; jb <= i__1; ++jb) {
	j = k + ks - jb;
	i__2 = j + k * a_dim1;
	temp.r = a[i__2].r, temp.i = a[i__2].i;
	i__2 = j + k * a_dim1;
	i__3 = ks + j * a_dim1;
	a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
	i__2 = ks + j * a_dim1;
	a[i__2].r = temp.r, a[i__2].i = temp.i;
/* L190: */
    }
    if (kstep == 1) {
	goto L200;
    }
    i__1 = ks + (k + 1) * a_dim1;
    temp.r = a[i__1].r, temp.i = a[i__1].i;
    i__1 = ks + (k + 1) * a_dim1;
    i__2 = k + (k + 1) * a_dim1;
    a[i__1].r = a[i__2].r, a[i__1].i = a[i__2].i;
    i__1 = k + (k + 1) * a_dim1;
    a[i__1].r = temp.r, a[i__1].i = temp.i;
L200:
L210:
    k += kstep;
    goto L110;
L220:
L230:
    return 0;
} /* zsidi_ */

