/* dsifa.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dsifa_(doublereal *a, integer *lda, integer *n, integer *
	kpvt, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer j, k;
    doublereal t, ak, bk;
    integer jj, km1, km2;
    doublereal akm1, bkm1;
    integer imax, jmax;
    doublereal mulk;
    logical swap;
    doublereal alpha, denom;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    integer kstep, imaxp1;
    doublereal mulkm1, absakk;
    extern integer idamax_(integer *, doublereal *, integer *);
    doublereal colmax, rowmax;


/*     dsifa factors a double precision symmetric matrix by elimination */
/*     with symmetric pivoting. */

/*     to solve  a*x = b , follow dsifa by dsisl. */
/*     to compute  inverse(a)*c , follow dsifa by dsisl. */
/*     to compute  determinant(a) , follow dsifa by dsidi. */
/*     to compute  inertia(a) , follow dsifa by dsidi. */
/*     to compute  inverse(a) , follow dsifa by dsidi. */

/*     on entry */

/*        a       double precision(lda,n) */
/*                the symmetric matrix to be factored. */
/*                only the diagonal and upper triangle are used. */

/*        lda     integer */
/*                the leading dimension of the array  a . */

/*        n       integer */
/*                the order of the matrix  a . */

/*     on return */

/*        a       a block diagonal matrix and the multipliers which */
/*                were used to obtain it. */
/*                the factorization can be written  a = u*d*trans(u) */
/*                where  u  is a product of permutation and unit */
/*                upper triangular matrices , trans(u) is the */
/*                transpose of  u , and  d  is block diagonal */
/*                with 1 by 1 and 2 by 2 blocks. */

/*        kpvt    integer(n) */
/*                an integer vector of pivot indices. */

/*        info    integer */
/*                = 0  normal value. */
/*                = k  if the k-th pivot block is singular. this is */
/*                     not an error condition for this subroutine, */
/*                     but it does indicate that dsisl or dsidi may */
/*                     divide by zero if called. */

/*     linpack. this version dated 08/14/78 . */
/*     james bunch, univ. calif. san diego, argonne nat. lab. */

/*     subroutines and functions */

/*     blas daxpy,dswap,idamax */
/*     fortran dabs,dmax1,dsqrt */

/*     internal variables */



/*     initialize */

/*     alpha is used in choosing pivot block size. */
    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --kpvt;

    /* Function Body */
    alpha = (sqrt(17.) + 1.) / 8.;

    *info = 0;

/*     main loop on k, which goes from n to 1. */

    k = *n;
L10:

/*        leave the loop if k=0 or k=1. */

/*     ...exit */
    if (k == 0) {
	goto L200;
    }
    if (k > 1) {
	goto L20;
    }
    kpvt[1] = 1;
    if (a[a_dim1 + 1] == 0.) {
	*info = 1;
    }
/*     ......exit */
    goto L200;
L20:

/*        this section of code determines the kind of */
/*        elimination to be performed.  when it is completed, */
/*        kstep will be set to the size of the pivot block, and */
/*        swap will be set to .true. if an interchange is */
/*        required. */

    km1 = k - 1;
    absakk = (d__1 = a[k + k * a_dim1], abs(d__1));

/*        determine the largest off-diagonal element in */
/*        column k. */

    i__1 = k - 1;
    imax = idamax_(&i__1, &a[k * a_dim1 + 1], &c__1);
    colmax = (d__1 = a[imax + k * a_dim1], abs(d__1));
    if (absakk < alpha * colmax) {
	goto L30;
    }
    kstep = 1;
    swap = FALSE_;
    goto L90;
L30:

/*           determine the largest off-diagonal element in */
/*           row imax. */

    rowmax = 0.;
    imaxp1 = imax + 1;
    i__1 = k;
    for (j = imaxp1; j <= i__1; ++j) {
/* Computing MAX */
	d__2 = rowmax, d__3 = (d__1 = a[imax + j * a_dim1], abs(d__1));
	rowmax = max(d__2,d__3);
/* L40: */
    }
    if (imax == 1) {
	goto L50;
    }
    i__1 = imax - 1;
    jmax = idamax_(&i__1, &a[imax * a_dim1 + 1], &c__1);
/* Computing MAX */
    d__2 = rowmax, d__3 = (d__1 = a[jmax + imax * a_dim1], abs(d__1));
    rowmax = max(d__2,d__3);
L50:
    if ((d__1 = a[imax + imax * a_dim1], abs(d__1)) < alpha * rowmax) {
	goto L60;
    }
    kstep = 1;
    swap = TRUE_;
    goto L80;
L60:
    if (absakk < alpha * colmax * (colmax / rowmax)) {
	goto L70;
    }
    kstep = 1;
    swap = FALSE_;
    goto L80;
L70:
    kstep = 2;
    swap = imax != km1;
L80:
L90:
    if (max(absakk,colmax) != 0.) {
	goto L100;
    }

/*           column k is zero.  set info and iterate the loop. */

    kpvt[k] = k;
    *info = k;
    goto L190;
L100:
    if (kstep == 2) {
	goto L140;
    }

/*           1 x 1 pivot block. */

    if (! swap) {
	goto L120;
    }

/*              perform an interchange. */

    dswap_(&imax, &a[imax * a_dim1 + 1], &c__1, &a[k * a_dim1 + 1], &c__1);
    i__1 = k;
    for (jj = imax; jj <= i__1; ++jj) {
	j = k + imax - jj;
	t = a[j + k * a_dim1];
	a[j + k * a_dim1] = a[imax + j * a_dim1];
	a[imax + j * a_dim1] = t;
/* L110: */
    }
L120:

/*           perform the elimination. */

    i__1 = km1;
    for (jj = 1; jj <= i__1; ++jj) {
	j = k - jj;
	mulk = -a[j + k * a_dim1] / a[k + k * a_dim1];
	t = mulk;
	daxpy_(&j, &t, &a[k * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &c__1);
	a[j + k * a_dim1] = mulk;
/* L130: */
    }

/*           set the pivot array. */

    kpvt[k] = k;
    if (swap) {
	kpvt[k] = imax;
    }
    goto L190;
L140:

/*           2 x 2 pivot block. */

    if (! swap) {
	goto L160;
    }

/*              perform an interchange. */

    dswap_(&imax, &a[imax * a_dim1 + 1], &c__1, &a[(k - 1) * a_dim1 + 1], &
	    c__1);
    i__1 = km1;
    for (jj = imax; jj <= i__1; ++jj) {
	j = km1 + imax - jj;
	t = a[j + (k - 1) * a_dim1];
	a[j + (k - 1) * a_dim1] = a[imax + j * a_dim1];
	a[imax + j * a_dim1] = t;
/* L150: */
    }
    t = a[k - 1 + k * a_dim1];
    a[k - 1 + k * a_dim1] = a[imax + k * a_dim1];
    a[imax + k * a_dim1] = t;
L160:

/*           perform the elimination. */

    km2 = k - 2;
    if (km2 == 0) {
	goto L180;
    }
    ak = a[k + k * a_dim1] / a[k - 1 + k * a_dim1];
    akm1 = a[k - 1 + (k - 1) * a_dim1] / a[k - 1 + k * a_dim1];
    denom = 1. - ak * akm1;
    i__1 = km2;
    for (jj = 1; jj <= i__1; ++jj) {
	j = km1 - jj;
	bk = a[j + k * a_dim1] / a[k - 1 + k * a_dim1];
	bkm1 = a[j + (k - 1) * a_dim1] / a[k - 1 + k * a_dim1];
	mulk = (akm1 * bk - bkm1) / denom;
	mulkm1 = (ak * bkm1 - bk) / denom;
	t = mulk;
	daxpy_(&j, &t, &a[k * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &c__1);
	t = mulkm1;
	daxpy_(&j, &t, &a[(k - 1) * a_dim1 + 1], &c__1, &a[j * a_dim1 + 1], &
		c__1);
	a[j + k * a_dim1] = mulk;
	a[j + (k - 1) * a_dim1] = mulkm1;
/* L170: */
    }
L180:

/*           set the pivot array. */

    kpvt[k] = 1 - k;
    if (swap) {
	kpvt[k] = -imax;
    }
    kpvt[k - 1] = kpvt[k];
L190:
    k -= kstep;
    goto L10;
L200:
    return 0;
} /* dsifa_ */

