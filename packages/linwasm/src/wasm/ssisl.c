/* ssisl.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int ssisl_(real *a, integer *lda, integer *n, integer *kpvt, 
	real *b)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;

    /* Local variables */
    integer k;
    real ak, bk;
    integer kp;
    real akm1, bkm1, temp;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    real denom;
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *);


/*     ssisl solves the real symmetric system */
/*     a * x = b */
/*     using the factors computed by ssifa. */

/*     on entry */

/*        a       real(lda,n) */
/*                the output from ssifa. */

/*        lda     integer */
/*                the leading dimension of the array  a . */

/*        n       integer */
/*                the order of the matrix  a . */

/*        kpvt    integer(n) */
/*                the pivot vector from ssifa. */

/*        b       real(n) */
/*                the right hand side vector. */

/*     on return */

/*        b       the solution vector  x . */

/*     error condition */

/*        a division by zero may occur if  ssico  has set rcond .eq. 0.0 */
/*        or  ssifa  has set info .ne. 0  . */

/*     to compute  inverse(a) * c  where  c  is a matrix */
/*     with  p  columns */
/*           call ssifa(a,lda,n,kpvt,info) */
/*           if (info .ne. 0) go to ... */
/*           do 10 j = 1, p */
/*              call ssisl(a,lda,n,kpvt,c(1,j)) */
/*        10 continue */

/*     linpack. this version dated 08/14/78 . */
/*     james bunch, univ. calif. san diego, argonne nat. lab. */

/*     subroutines and functions */

/*     blas saxpy,sdot */
/*     fortran iabs */

/*     internal variables. */


/*     loop backward applying the transformations and */
/*     d inverse to b. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --kpvt;
    --b;

    /* Function Body */
    k = *n;
L10:
    if (k == 0) {
	goto L80;
    }
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
    saxpy_(&i__1, &b[k], &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
L30:

/*           apply d inverse. */

    b[k] /= a[k + k * a_dim1];
    --k;
    goto L70;
L40:

/*           2 x 2 pivot block. */

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
    saxpy_(&i__1, &b[k], &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
    i__1 = k - 2;
    saxpy_(&i__1, &b[k - 1], &a[(k - 1) * a_dim1 + 1], &c__1, &b[1], &c__1);
L60:

/*           apply d inverse. */

    ak = a[k + k * a_dim1] / a[k - 1 + k * a_dim1];
    akm1 = a[k - 1 + (k - 1) * a_dim1] / a[k - 1 + k * a_dim1];
    bk = b[k] / a[k - 1 + k * a_dim1];
    bkm1 = b[k - 1] / a[k - 1 + k * a_dim1];
    denom = ak * akm1 - 1.f;
    b[k] = (akm1 * bk - bkm1) / denom;
    b[k - 1] = (ak * bkm1 - bk) / denom;
    k += -2;
L70:
    goto L10;
L80:

/*     loop forward applying the transformations. */

    k = 1;
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
    b[k] += sdot_(&i__1, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
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
    ++k;
    goto L150;
L120:

/*           2 x 2 pivot block. */

    if (k == 1) {
	goto L140;
    }

/*              apply the transformation. */

    i__1 = k - 1;
    b[k] += sdot_(&i__1, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
    i__1 = k - 1;
    b[k + 1] += sdot_(&i__1, &a[(k + 1) * a_dim1 + 1], &c__1, &b[1], &c__1);
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
    k += 2;
L150:
    goto L90;
L160:
    return 0;
} /* ssisl_ */

