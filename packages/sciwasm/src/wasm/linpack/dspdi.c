/* dspdi.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dspdi_(doublereal *ap, integer *n, integer *kpvt, 
	doublereal *det, integer *inert, doublereal *work, integer *job)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    doublereal d__;
    integer j, k;
    doublereal t, ak;
    integer jb, ij, ik, jk, kk, ks, km1;
    doublereal ten;
    integer iks, ksj;
    doublereal akp1;
    integer ikp1, jkp1, kkp1;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    doublereal temp, akkp1;
    integer kskp1;
    logical nodet;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dswap_(integer *, doublereal *, integer 
	    *, doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    integer kstep;
    logical noert, noinv;


/*     dspdi computes the determinant, inertia and inverse */
/*     of a double precision symmetric matrix using the factors from */
/*     dspfa, where the matrix is stored in packed form. */

/*     on entry */

/*        ap      double precision (n*(n+1)/2) */
/*                the output from dspfa. */

/*        n       integer */
/*                the order of the matrix a. */

/*        kpvt    integer(n) */
/*                the pivot vector from dspfa. */

/*        work    double precision(n) */
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

/*        a division by zero will occur if the inverse is requested */
/*        and  dspco  has set rcond .eq. 0.0 */
/*        or  dspfa  has set  info .ne. 0 . */

/*     linpack. this version dated 08/14/78 . */
/*     james bunch, univ. calif. san diego, argonne nat. lab. */

/*     subroutines and functions */

/*     blas daxpy,dcopy,ddot,dswap */
/*     fortran dabs,iabs,mod */

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
    det[1] = 1.;
    det[2] = 0.;
    ten = 10.;
L20:
    t = 0.;
    ik = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	kk = ik + k;
	d__ = ap[kk];

/*           check if 1 by 1 */

	if (kpvt[k] > 0) {
	    goto L50;
	}

/*              2 by 2 block */
/*              use det (d  s)  =  (d/t * c - t) * t  ,  t = dabs(s) */
/*                      (s  c) */
/*              to avoid underflow/overflow troubles. */
/*              take two passes through scaling.  use  t  for flag. */

	if (t != 0.) {
	    goto L30;
	}
	ikp1 = ik + k;
	kkp1 = ikp1 + k;
	t = (d__1 = ap[kkp1], abs(d__1));
	d__ = d__ / t * ap[kkp1 + 1] - t;
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

    ap[kk] = 1. / ap[kk];
    if (km1 < 1) {
	goto L170;
    }
    dcopy_(&km1, &ap[ik + 1], &c__1, &work[1], &c__1);
    ij = 0;
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	jk = ik + j;
	ap[jk] = ddot_(&j, &ap[ij + 1], &c__1, &work[1], &c__1);
	i__2 = j - 1;
	daxpy_(&i__2, &work[j], &ap[ij + 1], &c__1, &ap[ik + 1], &c__1);
	ij += j;
/* L160: */
    }
    ap[kk] += ddot_(&km1, &work[1], &c__1, &ap[ik + 1], &c__1);
L170:
    kstep = 1;
    goto L220;
L180:

/*              2 by 2 */

    t = (d__1 = ap[kkp1], abs(d__1));
    ak = ap[kk] / t;
    akp1 = ap[kkp1 + 1] / t;
    akkp1 = ap[kkp1] / t;
    d__ = t * (ak * akp1 - 1.);
    ap[kk] = akp1 / d__;
    ap[kkp1 + 1] = ak / d__;
    ap[kkp1] = -akkp1 / d__;
    if (km1 < 1) {
	goto L210;
    }
    dcopy_(&km1, &ap[ikp1 + 1], &c__1, &work[1], &c__1);
    ij = 0;
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	jkp1 = ikp1 + j;
	ap[jkp1] = ddot_(&j, &ap[ij + 1], &c__1, &work[1], &c__1);
	i__2 = j - 1;
	daxpy_(&i__2, &work[j], &ap[ij + 1], &c__1, &ap[ikp1 + 1], &c__1);
	ij += j;
/* L190: */
    }
    ap[kkp1 + 1] += ddot_(&km1, &work[1], &c__1, &ap[ikp1 + 1], &c__1);
    ap[kkp1] += ddot_(&km1, &ap[ik + 1], &c__1, &ap[ikp1 + 1], &c__1);
    dcopy_(&km1, &ap[ik + 1], &c__1, &work[1], &c__1);
    ij = 0;
    i__1 = km1;
    for (j = 1; j <= i__1; ++j) {
	jk = ik + j;
	ap[jk] = ddot_(&j, &ap[ij + 1], &c__1, &work[1], &c__1);
	i__2 = j - 1;
	daxpy_(&i__2, &work[j], &ap[ij + 1], &c__1, &ap[ik + 1], &c__1);
	ij += j;
/* L200: */
    }
    ap[kk] += ddot_(&km1, &work[1], &c__1, &ap[ik + 1], &c__1);
L210:
    kstep = 2;
L220:

/*           swap */

    ks = (i__1 = kpvt[k], abs(i__1));
    if (ks == k) {
	goto L250;
    }
    iks = ks * (ks - 1) / 2;
    dswap_(&ks, &ap[iks + 1], &c__1, &ap[ik + 1], &c__1);
    ksj = ik + ks;
    i__1 = k;
    for (jb = ks; jb <= i__1; ++jb) {
	j = k + ks - jb;
	jk = ik + j;
	temp = ap[jk];
	ap[jk] = ap[ksj];
	ap[ksj] = temp;
	ksj -= j - 1;
/* L230: */
    }
    if (kstep == 1) {
	goto L240;
    }
    kskp1 = ikp1 + ks;
    temp = ap[kskp1];
    ap[kskp1] = ap[kkp1];
    ap[kkp1] = temp;
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
} /* dspdi_ */

