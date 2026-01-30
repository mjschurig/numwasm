/* zspdi.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int zspdi_(doublecomplex *ap, integer *n, integer *kpvt, 
	doublecomplex *det, doublecomplex *work, integer *job)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    doublecomplex d__;
    integer j, k;
    doublecomplex t, ak;
    integer jb, ij, ik, jk, kk, ks, km1;
    doublereal ten;
    integer iks, ksj;
    doublecomplex akp1;
    integer ikp1, jkp1, kkp1;
    doublecomplex temp, akkp1;
    integer kskp1;
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


/*     zspdi computes the determinant and inverse */
/*     of a complex*16 symmetric matrix using the factors from zspfa, */
/*     where the matrix is stored in packed form. */

/*     on entry */

/*        ap      complex*16 (n*(n+1)/2) */
/*                the output from zspfa. */

/*        n       integer */
/*                the order of the matrix a. */

/*        kpvt    integer(n) */
/*                the pivot vector from zspfa. */

/*        work    complex*16(n) */
/*                work vector.  contents ignored. */

/*        job     integer */
/*                job has the decimal expansion  ab  where */
/*                   if  b .ne. 0, the inverse is computed, */
/*                   if  a .ne. 0, the determinant is computed, */

/*                for example, job = 11  gives both. */

/*     on return */

/*        variables not requested by job are not used. */

/*        ap     contains the upper triangle of the inverse of */
/*               the original matrix, stored in packed form. */
/*               the columns of the upper triangle are stored */
/*               sequentially in a one-dimensional array. */

/*        det    complex*16(2) */
/*               determinant of original matrix. */
/*               determinant = det(1) * 10.0**det(2) */
/*               with 1.0 .le. dabs(det(1)) .lt. 10.0 */
/*               or det(1) = 0.0. */

/*     error condition */

/*        a division by zero will occur if the inverse is requested */
/*        and  zspco  has set rcond .eq. 0.0 */
/*        or  zspfa  has set  info .ne. 0 . */

/*     linpack. this version dated 08/14/78 . */
/*     james bunch, univ. calif. san diego, argonne nat. lab. */

/*     subroutines and functions */

/*     blas zaxpy,zcopy,zdotu,zswap */
/*     fortran dabs,dcmplx,iabs,mod */

/*     internal variables. */



    /* Parameter adjustments */
    --work;
    --det;
    --kpvt;
    --ap;

    /* Function Body */
    noinv = *job % 10 == 0;
    nodet = *job % 100 / 10 == 0;

    if (nodet) {
	goto L110;
    }
    det[1].r = 1., det[1].i = 0.;
    det[2].r = 0., det[2].i = 0.;
    ten = 10.;
    t.r = 0., t.i = 0.;
    ik = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	kk = ik + k;
	i__2 = kk;
	d__.r = ap[i__2].r, d__.i = ap[i__2].i;

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
	ikp1 = ik + k;
	kkp1 = ikp1 + k;
	i__2 = kkp1;
	t.r = ap[i__2].r, t.i = ap[i__2].i;
	z_div(&z__3, &d__, &t);
	i__2 = kkp1 + 1;
	z__2.r = z__3.r * ap[i__2].r - z__3.i * ap[i__2].i, z__2.i = z__3.r * 
		ap[i__2].i + z__3.i * ap[i__2].r;
	z__1.r = z__2.r - t.r, z__1.i = z__2.i - t.i;
	d__.r = z__1.r, d__.i = z__1.i;
	goto L20;
L10:
	d__.r = t.r, d__.i = t.i;
	t.r = 0., t.i = 0.;
L20:
L30:

	if (nodet) {
	    goto L90;
	}
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
L90:
	ik += k;
/* L100: */
    }
L110:

/*     compute inverse(a) */

    if (noinv) {
	goto L240;
    }
    k = 1;
    ik = 0;
L120:
    if (k > *n) {
	goto L230;
    }
    km1 = k - 1;
    kk = ik + k;
    ikp1 = ik + k;
    if (kpvt[k] < 0) {
	goto L150;
    }

/*              1 by 1 */

    i__1 = kk;
    z_div(&z__1, &c_b3, &ap[kk]);
    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
    if (km1 < 1) {
	goto L140;
    }
    zcopy_(&km1, &ap[ik + 1], &c__1, &work[1], &c__1);
    ij = 0;
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	jk = ik + j;
	i__2 = jk;
	zdotu_(&z__1, &j, &ap[ij + 1], &c__1, &work[1], &c__1);
	ap[i__2].r = z__1.r, ap[i__2].i = z__1.i;
	i__2 = j - 1;
	zaxpy_(&i__2, &work[j], &ap[ij + 1], &c__1, &ap[ik + 1], &c__1);
	ij += j;
/* L130: */
    }
    i__1 = kk;
    i__2 = kk;
    zdotu_(&z__2, &km1, &work[1], &c__1, &ap[ik + 1], &c__1);
    z__1.r = ap[i__2].r + z__2.r, z__1.i = ap[i__2].i + z__2.i;
    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
L140:
    kstep = 1;
    goto L190;
L150:

/*              2 by 2 */

    kkp1 = ikp1 + k;
    i__1 = kkp1;
    t.r = ap[i__1].r, t.i = ap[i__1].i;
    z_div(&z__1, &ap[kk], &t);
    ak.r = z__1.r, ak.i = z__1.i;
    z_div(&z__1, &ap[kkp1 + 1], &t);
    akp1.r = z__1.r, akp1.i = z__1.i;
    z_div(&z__1, &ap[kkp1], &t);
    akkp1.r = z__1.r, akkp1.i = z__1.i;
    z__3.r = ak.r * akp1.r - ak.i * akp1.i, z__3.i = ak.r * akp1.i + ak.i * 
	    akp1.r;
    z__2.r = z__3.r - 1., z__2.i = z__3.i - 0.;
    z__1.r = t.r * z__2.r - t.i * z__2.i, z__1.i = t.r * z__2.i + t.i * 
	    z__2.r;
    d__.r = z__1.r, d__.i = z__1.i;
    i__1 = kk;
    z_div(&z__1, &akp1, &d__);
    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
    i__1 = kkp1 + 1;
    z_div(&z__1, &ak, &d__);
    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
    i__1 = kkp1;
    z__2.r = -akkp1.r, z__2.i = -akkp1.i;
    z_div(&z__1, &z__2, &d__);
    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
    if (km1 < 1) {
	goto L180;
    }
    zcopy_(&km1, &ap[ikp1 + 1], &c__1, &work[1], &c__1);
    ij = 0;
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	jkp1 = ikp1 + j;
	i__2 = jkp1;
	zdotu_(&z__1, &j, &ap[ij + 1], &c__1, &work[1], &c__1);
	ap[i__2].r = z__1.r, ap[i__2].i = z__1.i;
	i__2 = j - 1;
	zaxpy_(&i__2, &work[j], &ap[ij + 1], &c__1, &ap[ikp1 + 1], &c__1);
	ij += j;
/* L160: */
    }
    i__1 = kkp1 + 1;
    i__2 = kkp1 + 1;
    zdotu_(&z__2, &km1, &work[1], &c__1, &ap[ikp1 + 1], &c__1);
    z__1.r = ap[i__2].r + z__2.r, z__1.i = ap[i__2].i + z__2.i;
    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
    i__1 = kkp1;
    i__2 = kkp1;
    zdotu_(&z__2, &km1, &ap[ik + 1], &c__1, &ap[ikp1 + 1], &c__1);
    z__1.r = ap[i__2].r + z__2.r, z__1.i = ap[i__2].i + z__2.i;
    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
    zcopy_(&km1, &ap[ik + 1], &c__1, &work[1], &c__1);
    ij = 0;
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	jk = ik + j;
	i__2 = jk;
	zdotu_(&z__1, &j, &ap[ij + 1], &c__1, &work[1], &c__1);
	ap[i__2].r = z__1.r, ap[i__2].i = z__1.i;
	i__2 = j - 1;
	zaxpy_(&i__2, &work[j], &ap[ij + 1], &c__1, &ap[ik + 1], &c__1);
	ij += j;
/* L170: */
    }
    i__1 = kk;
    i__2 = kk;
    zdotu_(&z__2, &km1, &work[1], &c__1, &ap[ik + 1], &c__1);
    z__1.r = ap[i__2].r + z__2.r, z__1.i = ap[i__2].i + z__2.i;
    ap[i__1].r = z__1.r, ap[i__1].i = z__1.i;
L180:
    kstep = 2;
L190:

/*           swap */

    ks = (i__1 = kpvt[k], abs(i__1));
    if (ks == k) {
	goto L220;
    }
    iks = ks * (ks - 1) / 2;
    zswap_(&ks, &ap[iks + 1], &c__1, &ap[ik + 1], &c__1);
    ksj = ik + ks;
    i__1 = k;
    for (jb = ks; jb <= i__1; ++jb) {
	j = k + ks - jb;
	jk = ik + j;
	i__2 = jk;
	temp.r = ap[i__2].r, temp.i = ap[i__2].i;
	i__2 = jk;
	i__3 = ksj;
	ap[i__2].r = ap[i__3].r, ap[i__2].i = ap[i__3].i;
	i__2 = ksj;
	ap[i__2].r = temp.r, ap[i__2].i = temp.i;
	ksj -= j - 1;
/* L200: */
    }
    if (kstep == 1) {
	goto L210;
    }
    kskp1 = ikp1 + ks;
    i__1 = kskp1;
    temp.r = ap[i__1].r, temp.i = ap[i__1].i;
    i__1 = kskp1;
    i__2 = kkp1;
    ap[i__1].r = ap[i__2].r, ap[i__1].i = ap[i__2].i;
    i__1 = kkp1;
    ap[i__1].r = temp.r, ap[i__1].i = temp.i;
L210:
L220:
    ik += k;
    if (kstep == 2) {
	ik = ik + k + 1;
    }
    k += kstep;
    goto L120;
L230:
L240:
    return 0;
} /* zspdi_ */

