/* dspsl.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dspsl_(doublereal *ap, integer *n, integer *kpvt, 
	doublereal *b)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer k;
    doublereal ak, bk;
    integer ik, kk, kp;
    doublereal akm1, bkm1;
    integer ikm1, km1k, ikp1;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    doublereal temp;
    integer km1km1;
    doublereal denom;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);


/*     dsisl solves the double precision symmetric system */
/*     a * x = b */
/*     using the factors computed by dspfa. */

/*     on entry */

/*        ap      double precision(n*(n+1)/2) */
/*                the output from dspfa. */

/*        n       integer */
/*                the order of the matrix  a . */

/*        kpvt    integer(n) */
/*                the pivot vector from dspfa. */

/*        b       double precision(n) */
/*                the right hand side vector. */

/*     on return */

/*        b       the solution vector  x . */

/*     error condition */

/*        a division by zero may occur if  dspco  has set rcond .eq. 0.0 */
/*        or  dspfa  has set info .ne. 0  . */

/*     to compute  inverse(a) * c  where  c  is a matrix */
/*     with  p  columns */
/*           call dspfa(ap,n,kpvt,info) */
/*           if (info .ne. 0) go to ... */
/*           do 10 j = 1, p */
/*              call dspsl(ap,n,kpvt,c(1,j)) */
/*        10 continue */

/*     linpack. this version dated 08/14/78 . */
/*     james bunch, univ. calif. san diego, argonne nat. lab. */

/*     subroutines and functions */

/*     blas daxpy,ddot */
/*     fortran iabs */

/*     internal variables. */


/*     loop backward applying the transformations and */
/*     d inverse to b. */

    /* Parameter adjustments */
    --b;
    --kpvt;
    --ap;

    /* Function Body */
    k = *n;
    ik = *n * (*n - 1) / 2;
L10:
    if (k == 0) {
	goto L80;
    }
    kk = ik + k;
    if (kpvt[k] < 0) {
	goto L40;
    }

/*           1 x 1 pivot block. */

    if (k == 1) {
	goto L30;
    }
    kp = kpvt[k];
    if (kp == k) {
	goto L20;
    }

/*                 interchange. */

    temp = b[k];
    b[k] = b[kp];
    b[kp] = temp;
L20:

/*              apply the transformation. */

    i__1 = k - 1;
    daxpy_(&i__1, &b[k], &ap[ik + 1], &c__1, &b[1], &c__1);
L30:

/*           apply d inverse. */

    b[k] /= ap[kk];
    --k;
    ik -= k;
    goto L70;
L40:

/*           2 x 2 pivot block. */

    ikm1 = ik - (k - 1);
    if (k == 2) {
	goto L60;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k - 1) {
	goto L50;
    }

/*                 interchange. */

    temp = b[k - 1];
    b[k - 1] = b[kp];
    b[kp] = temp;
L50:

/*              apply the transformation. */

    i__1 = k - 2;
    daxpy_(&i__1, &b[k], &ap[ik + 1], &c__1, &b[1], &c__1);
    i__1 = k - 2;
    daxpy_(&i__1, &b[k - 1], &ap[ikm1 + 1], &c__1, &b[1], &c__1);
L60:

/*           apply d inverse. */

    km1k = ik + k - 1;
    kk = ik + k;
    ak = ap[kk] / ap[km1k];
    km1km1 = ikm1 + k - 1;
    akm1 = ap[km1km1] / ap[km1k];
    bk = b[k] / ap[km1k];
    bkm1 = b[k - 1] / ap[km1k];
    denom = ak * akm1 - 1.;
    b[k] = (akm1 * bk - bkm1) / denom;
    b[k - 1] = (ak * bkm1 - bk) / denom;
    k += -2;
    ik = ik - (k + 1) - k;
L70:
    goto L10;
L80:

/*     loop forward applying the transformations. */

    k = 1;
    ik = 0;
L90:
    if (k > *n) {
	goto L160;
    }
    if (kpvt[k] < 0) {
	goto L120;
    }

/*           1 x 1 pivot block. */

    if (k == 1) {
	goto L110;
    }

/*              apply the transformation. */

    i__1 = k - 1;
    b[k] += ddot_(&i__1, &ap[ik + 1], &c__1, &b[1], &c__1);
    kp = kpvt[k];
    if (kp == k) {
	goto L100;
    }

/*                 interchange. */

    temp = b[k];
    b[k] = b[kp];
    b[kp] = temp;
L100:
L110:
    ik += k;
    ++k;
    goto L150;
L120:

/*           2 x 2 pivot block. */

    if (k == 1) {
	goto L140;
    }

/*              apply the transformation. */

    i__1 = k - 1;
    b[k] += ddot_(&i__1, &ap[ik + 1], &c__1, &b[1], &c__1);
    ikp1 = ik + k;
    i__1 = k - 1;
    b[k + 1] += ddot_(&i__1, &ap[ikp1 + 1], &c__1, &b[1], &c__1);
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k) {
	goto L130;
    }

/*                 interchange. */

    temp = b[k];
    b[k] = b[kp];
    b[kp] = temp;
L130:
L140:
    ik = ik + k + k + 1;
    k += 2;
L150:
    goto L90;
L160:
    return 0;
} /* dspsl_ */

