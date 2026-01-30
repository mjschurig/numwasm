/* cspsl.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int cspsl_(complex *ap, integer *n, integer *kpvt, complex *
	b)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    complex q__1, q__2, q__3;

    /* Builtin functions */
    void c_div(complex *, complex *, complex *);

    /* Local variables */
    integer k;
    complex ak, bk;
    integer ik, kk, kp;
    complex akm1, bkm1;
    integer ikm1, km1k, ikp1;
    complex temp;
    integer km1km1;
    complex denom;
    extern /* Complex */ VOID cdotu_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);
    extern /* Subroutine */ int caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *);


/*     csisl solves the complex symmetric system */
/*     a * x = b */
/*     using the factors computed by cspfa. */

/*     on entry */

/*        ap      complex(n*(n+1)/2) */
/*                the output from cspfa. */

/*        n       integer */
/*                the order of the matrix  a . */

/*        kpvt    integer(n) */
/*                the pivot vector from cspfa. */

/*        b       complex(n) */
/*                the right hand side vector. */

/*     on return */

/*        b       the solution vector  x . */

/*     error condition */

/*        a division by zero may occur if  cspco  has set rcond .eq. 0.0 */
/*        or  cspfa  has set info .ne. 0  . */

/*     to compute  inverse(a) * c  where  c  is a matrix */
/*     with  p  columns */
/*           call cspfa(ap,n,kpvt,info) */
/*           if (info .ne. 0) go to ... */
/*           do 10 j = 1, p */
/*              call cspsl(ap,n,kpvt,c(1,j)) */
/*        10 continue */

/*     linpack. this version dated 08/14/78 . */
/*     james bunch, univ. calif. san diego, argonne nat. lab. */

/*     subroutines and functions */

/*     blas caxpy,cdotu */
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

    i__1 = k;
    temp.r = b[i__1].r, temp.i = b[i__1].i;
    i__1 = k;
    i__2 = kp;
    b[i__1].r = b[i__2].r, b[i__1].i = b[i__2].i;
    i__1 = kp;
    b[i__1].r = temp.r, b[i__1].i = temp.i;
L20:

/*              apply the transformation. */

    i__1 = k - 1;
    caxpy_(&i__1, &b[k], &ap[ik + 1], &c__1, &b[1], &c__1);
L30:

/*           apply d inverse. */

    i__1 = k;
    c_div(&q__1, &b[k], &ap[kk]);
    b[i__1].r = q__1.r, b[i__1].i = q__1.i;
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

    i__1 = k - 1;
    temp.r = b[i__1].r, temp.i = b[i__1].i;
    i__1 = k - 1;
    i__2 = kp;
    b[i__1].r = b[i__2].r, b[i__1].i = b[i__2].i;
    i__1 = kp;
    b[i__1].r = temp.r, b[i__1].i = temp.i;
L50:

/*              apply the transformation. */

    i__1 = k - 2;
    caxpy_(&i__1, &b[k], &ap[ik + 1], &c__1, &b[1], &c__1);
    i__1 = k - 2;
    caxpy_(&i__1, &b[k - 1], &ap[ikm1 + 1], &c__1, &b[1], &c__1);
L60:

/*           apply d inverse. */

    km1k = ik + k - 1;
    kk = ik + k;
    c_div(&q__1, &ap[kk], &ap[km1k]);
    ak.r = q__1.r, ak.i = q__1.i;
    km1km1 = ikm1 + k - 1;
    c_div(&q__1, &ap[km1km1], &ap[km1k]);
    akm1.r = q__1.r, akm1.i = q__1.i;
    c_div(&q__1, &b[k], &ap[km1k]);
    bk.r = q__1.r, bk.i = q__1.i;
    c_div(&q__1, &b[k - 1], &ap[km1k]);
    bkm1.r = q__1.r, bkm1.i = q__1.i;
    q__2.r = ak.r * akm1.r - ak.i * akm1.i, q__2.i = ak.r * akm1.i + ak.i * 
	    akm1.r;
    q__1.r = q__2.r - 1.f, q__1.i = q__2.i;
    denom.r = q__1.r, denom.i = q__1.i;
    i__1 = k;
    q__3.r = akm1.r * bk.r - akm1.i * bk.i, q__3.i = akm1.r * bk.i + akm1.i * 
	    bk.r;
    q__2.r = q__3.r - bkm1.r, q__2.i = q__3.i - bkm1.i;
    c_div(&q__1, &q__2, &denom);
    b[i__1].r = q__1.r, b[i__1].i = q__1.i;
    i__1 = k - 1;
    q__3.r = ak.r * bkm1.r - ak.i * bkm1.i, q__3.i = ak.r * bkm1.i + ak.i * 
	    bkm1.r;
    q__2.r = q__3.r - bk.r, q__2.i = q__3.i - bk.i;
    c_div(&q__1, &q__2, &denom);
    b[i__1].r = q__1.r, b[i__1].i = q__1.i;
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

    i__1 = k;
    i__2 = k;
    i__3 = k - 1;
    cdotu_(&q__2, &i__3, &ap[ik + 1], &c__1, &b[1], &c__1);
    q__1.r = b[i__2].r + q__2.r, q__1.i = b[i__2].i + q__2.i;
    b[i__1].r = q__1.r, b[i__1].i = q__1.i;
    kp = kpvt[k];
    if (kp == k) {
	goto L100;
    }

/*                 interchange. */

    i__1 = k;
    temp.r = b[i__1].r, temp.i = b[i__1].i;
    i__1 = k;
    i__2 = kp;
    b[i__1].r = b[i__2].r, b[i__1].i = b[i__2].i;
    i__1 = kp;
    b[i__1].r = temp.r, b[i__1].i = temp.i;
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

    i__1 = k;
    i__2 = k;
    i__3 = k - 1;
    cdotu_(&q__2, &i__3, &ap[ik + 1], &c__1, &b[1], &c__1);
    q__1.r = b[i__2].r + q__2.r, q__1.i = b[i__2].i + q__2.i;
    b[i__1].r = q__1.r, b[i__1].i = q__1.i;
    ikp1 = ik + k;
    i__1 = k + 1;
    i__2 = k + 1;
    i__3 = k - 1;
    cdotu_(&q__2, &i__3, &ap[ikp1 + 1], &c__1, &b[1], &c__1);
    q__1.r = b[i__2].r + q__2.r, q__1.i = b[i__2].i + q__2.i;
    b[i__1].r = q__1.r, b[i__1].i = q__1.i;
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k) {
	goto L130;
    }

/*                 interchange. */

    i__1 = k;
    temp.r = b[i__1].r, temp.i = b[i__1].i;
    i__1 = k;
    i__2 = kp;
    b[i__1].r = b[i__2].r, b[i__1].i = b[i__2].i;
    i__1 = kp;
    b[i__1].r = temp.r, b[i__1].i = temp.i;
L130:
L140:
    ik = ik + k + k + 1;
    k += 2;
L150:
    goto L90;
L160:
    return 0;
} /* cspsl_ */

