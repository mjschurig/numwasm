/* chpdi.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int chpdi_(complex *ap, integer *n, integer *kpvt, real *det,
	 integer *inert, complex *work, integer *job)
{
    /* System generated locals */
    integer i__1, i__2;
    real r__1;
    complex q__1, q__2, q__3;

    /* Builtin functions */
    double c_abs(complex *);
    void r_cnjg(complex *, complex *);

    /* Local variables */
    real d__;
    integer j, k;
    real t, ak;
    integer jb, ij, ik, jk, kk, ks, km1;
    real ten;
    integer iks, ksj;
    real akp1;
    integer ikp1, jkp1, kkp1;
    complex temp, akkp1;
    integer kskp1;
    extern /* Complex */ VOID cdotc_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);
    logical nodet;
    extern /* Subroutine */ int ccopy_(integer *, complex *, integer *, 
	    complex *, integer *), cswap_(integer *, complex *, integer *, 
	    complex *, integer *), caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *);
    integer kstep;
    logical noert, noinv;


/*     chpdi computes the determinant, inertia and inverse */
/*     of a complex hermitian matrix using the factors from chpfa, */
/*     where the matrix is stored in packed form. */

/*     on entry */

/*        ap      complex (n*(n+1)/2) */
/*                the output from chpfa. */

/*        n       integer */
/*                the order of the matrix a. */

/*        kpvt    integer(n) */
/*                the pivot vector from chpfa. */

/*        work    complex(n) */
/*                work vector.  contents ignored. */

/*        job     integer */
/*                job has the decimal expansion  abc  where */
/*                   if  c .ne. 0, the inverse is computed, */
/*                   if  b .ne. 0, the determinant is computed, */
/*                   if  a .ne. 0, the inertia is computed. */

/*                for example, job = 111  gives all three. */

/*     on return */

/*        variables not requested by job are not used. */

/*        ap     contains the upper triangle of the inverse of */
/*               the original matrix, stored in packed form. */
/*               the columns of the upper triangle are stored */
/*               sequentially in a one-dimensional array. */

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

/*        a division by zero will occur if the inverse is requested */
/*        and  chpco  has set rcond .eq. 0.0 */
/*        or  chpfa  has set  info .ne. 0 . */

/*     linpack. this version dated 08/14/78 . */
/*     james bunch, univ. calif. san diego, argonne nat. lab. */

/*     subroutines and functions */

/*     blas caxpy,ccopy,cdotc,cswap */
/*     fortran abs,cabs,cmplx,conjg,iabs,mod,real */

/*     internal variables. */


    /* Parameter adjustments */
    --work;
    --inert;
    --det;
    --kpvt;
    --ap;

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
    ik = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	kk = ik + k;
	i__2 = kk;
	d__ = ap[i__2].r;

/*           check if 1 by 1 */

	if (kpvt[k] > 0) {
	    goto L50;
	}

/*              2 by 2 block */
/*              use det (d  s)  =  (d/t * c - t) * t  ,  t = cabs(s) */
/*                      (s  c) */
/*              to avoid underflow/overflow troubles. */
/*              take two passes through scaling.  use  t  for flag. */

	if (t != 0.f) {
	    goto L30;
	}
	ikp1 = ik + k;
	kkp1 = ikp1 + k;
	t = c_abs(&ap[kkp1]);
	i__2 = kkp1 + 1;
	d__ = d__ / t * ap[i__2].r - t;
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
	ik += k;
/* L130: */
    }
L140:

/*     compute inverse(a) */

    if (noinv) {
	goto L270;
    }
    k = 1;
    ik = 0;
L150:
    if (k > *n) {
	goto L260;
    }
    km1 = k - 1;
    kk = ik + k;
    ikp1 = ik + k;
    kkp1 = ikp1 + k;
    if (kpvt[k] < 0) {
	goto L180;
    }

/*              1 by 1 */

    i__1 = kk;
    i__2 = kk;
    r__1 = 1.f / ap[i__2].r;
    q__1.r = r__1, q__1.i = 0.f;
    ap[i__1].r = q__1.r, ap[i__1].i = q__1.i;
    if (km1 < 1) {
	goto L170;
    }
    ccopy_(&km1, &ap[ik + 1], &c__1, &work[1], &c__1);
    ij = 0;
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	jk = ik + j;
	i__2 = jk;
	cdotc_(&q__1, &j, &ap[ij + 1], &c__1, &work[1], &c__1);
	ap[i__2].r = q__1.r, ap[i__2].i = q__1.i;
	i__2 = j - 1;
	caxpy_(&i__2, &work[j], &ap[ij + 1], &c__1, &ap[ik + 1], &c__1);
	ij += j;
/* L160: */
    }
    i__1 = kk;
    i__2 = kk;
    cdotc_(&q__3, &km1, &work[1], &c__1, &ap[ik + 1], &c__1);
    r__1 = q__3.r;
    q__2.r = r__1, q__2.i = 0.f;
    q__1.r = ap[i__2].r + q__2.r, q__1.i = ap[i__2].i + q__2.i;
    ap[i__1].r = q__1.r, ap[i__1].i = q__1.i;
L170:
    kstep = 1;
    goto L220;
L180:

/*              2 by 2 */

    t = c_abs(&ap[kkp1]);
    i__1 = kk;
    ak = ap[i__1].r / t;
    i__1 = kkp1 + 1;
    akp1 = ap[i__1].r / t;
    i__1 = kkp1;
    q__1.r = ap[i__1].r / t, q__1.i = ap[i__1].i / t;
    akkp1.r = q__1.r, akkp1.i = q__1.i;
    d__ = t * (ak * akp1 - 1.f);
    i__1 = kk;
    r__1 = akp1 / d__;
    q__1.r = r__1, q__1.i = 0.f;
    ap[i__1].r = q__1.r, ap[i__1].i = q__1.i;
    i__1 = kkp1 + 1;
    r__1 = ak / d__;
    q__1.r = r__1, q__1.i = 0.f;
    ap[i__1].r = q__1.r, ap[i__1].i = q__1.i;
    i__1 = kkp1;
    q__2.r = -akkp1.r, q__2.i = -akkp1.i;
    q__1.r = q__2.r / d__, q__1.i = q__2.i / d__;
    ap[i__1].r = q__1.r, ap[i__1].i = q__1.i;
    if (km1 < 1) {
	goto L210;
    }
    ccopy_(&km1, &ap[ikp1 + 1], &c__1, &work[1], &c__1);
    ij = 0;
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	jkp1 = ikp1 + j;
	i__2 = jkp1;
	cdotc_(&q__1, &j, &ap[ij + 1], &c__1, &work[1], &c__1);
	ap[i__2].r = q__1.r, ap[i__2].i = q__1.i;
	i__2 = j - 1;
	caxpy_(&i__2, &work[j], &ap[ij + 1], &c__1, &ap[ikp1 + 1], &c__1);
	ij += j;
/* L190: */
    }
    i__1 = kkp1 + 1;
    i__2 = kkp1 + 1;
    cdotc_(&q__3, &km1, &work[1], &c__1, &ap[ikp1 + 1], &c__1);
    r__1 = q__3.r;
    q__2.r = r__1, q__2.i = 0.f;
    q__1.r = ap[i__2].r + q__2.r, q__1.i = ap[i__2].i + q__2.i;
    ap[i__1].r = q__1.r, ap[i__1].i = q__1.i;
    i__1 = kkp1;
    i__2 = kkp1;
    cdotc_(&q__2, &km1, &ap[ik + 1], &c__1, &ap[ikp1 + 1], &c__1);
    q__1.r = ap[i__2].r + q__2.r, q__1.i = ap[i__2].i + q__2.i;
    ap[i__1].r = q__1.r, ap[i__1].i = q__1.i;
    ccopy_(&km1, &ap[ik + 1], &c__1, &work[1], &c__1);
    ij = 0;
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	jk = ik + j;
	i__2 = jk;
	cdotc_(&q__1, &j, &ap[ij + 1], &c__1, &work[1], &c__1);
	ap[i__2].r = q__1.r, ap[i__2].i = q__1.i;
	i__2 = j - 1;
	caxpy_(&i__2, &work[j], &ap[ij + 1], &c__1, &ap[ik + 1], &c__1);
	ij += j;
/* L200: */
    }
    i__1 = kk;
    i__2 = kk;
    cdotc_(&q__3, &km1, &work[1], &c__1, &ap[ik + 1], &c__1);
    r__1 = q__3.r;
    q__2.r = r__1, q__2.i = 0.f;
    q__1.r = ap[i__2].r + q__2.r, q__1.i = ap[i__2].i + q__2.i;
    ap[i__1].r = q__1.r, ap[i__1].i = q__1.i;
L210:
    kstep = 2;
L220:

/*           swap */

    ks = (i__1 = kpvt[k], abs(i__1));
    if (ks == k) {
	goto L250;
    }
    iks = ks * (ks - 1) / 2;
    cswap_(&ks, &ap[iks + 1], &c__1, &ap[ik + 1], &c__1);
    ksj = ik + ks;
    i__1 = k;
    for (jb = ks; jb <= i__1; ++jb) {
	j = k + ks - jb;
	jk = ik + j;
	r_cnjg(&q__1, &ap[jk]);
	temp.r = q__1.r, temp.i = q__1.i;
	i__2 = jk;
	r_cnjg(&q__1, &ap[ksj]);
	ap[i__2].r = q__1.r, ap[i__2].i = q__1.i;
	i__2 = ksj;
	ap[i__2].r = temp.r, ap[i__2].i = temp.i;
	ksj -= j - 1;
/* L230: */
    }
    if (kstep == 1) {
	goto L240;
    }
    kskp1 = ikp1 + ks;
    i__1 = kskp1;
    temp.r = ap[i__1].r, temp.i = ap[i__1].i;
    i__1 = kskp1;
    i__2 = kkp1;
    ap[i__1].r = ap[i__2].r, ap[i__1].i = ap[i__2].i;
    i__1 = kkp1;
    ap[i__1].r = temp.r, ap[i__1].i = temp.i;
L240:
L250:
    ik += k;
    if (kstep == 2) {
	ik = ik + k + 1;
    }
    k += kstep;
    goto L150;
L260:
L270:
    return 0;
} /* chpdi_ */

