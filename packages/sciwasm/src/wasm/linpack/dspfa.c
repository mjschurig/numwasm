/* dspfa.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dspfa_(doublereal *ap, integer *n, integer *kpvt, 
	integer *info)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer j, k;
    doublereal t, ak, bk;
    integer ij, ik, jj, im, jk, kk, km1, km2, ijj, imj, imk;
    doublereal akm1, bkm1;
    integer ikm1, jkm1, km1k, imim, jmim, imax, jmax;
    doublereal mulk;
    logical swap;
    doublereal alpha;
    integer km1km1;
    doublereal denom;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    integer kstep, imaxp1;
    doublereal mulkm1, absakk;
    extern integer idamax_(integer *, doublereal *, integer *);
    doublereal colmax, rowmax;


/*     dspfa factors a double precision symmetric matrix stored in */
/*     packed form by elimination with symmetric pivoting. */

/*     to solve  a*x = b , follow dspfa by dspsl. */
/*     to compute  inverse(a)*c , follow dspfa by dspsl. */
/*     to compute  determinant(a) , follow dspfa by dspdi. */
/*     to compute  inertia(a) , follow dspfa by dspdi. */
/*     to compute  inverse(a) , follow dspfa by dspdi. */

/*     on entry */

/*        ap      double precision (n*(n+1)/2) */
/*                the packed form of a symmetric matrix  a .  the */
/*                columns of the upper triangle are stored sequentially */
/*                in a one-dimensional array of length  n*(n+1)/2 . */
/*                see comments below for details. */

/*        n       integer */
/*                the order of the matrix  a . */

/*     output */

/*        ap      a block diagonal matrix and the multipliers which */
/*                were used to obtain it stored in packed form. */
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
/*                     but it does indicate that dspsl or dspdi may */
/*                     divide by zero if called. */

/*     packed storage */

/*          the following program segment will pack the upper */
/*          triangle of a symmetric matrix. */

/*                k = 0 */
/*                do 20 j = 1, n */
/*                   do 10 i = 1, j */
/*                      k = k + 1 */
/*                      ap(k)  = a(i,j) */
/*             10    continue */
/*             20 continue */

/*     linpack. this version dated 08/14/78 . */
/*     james bunch, univ. calif. san diego, argonne nat. lab. */

/*     subroutines and functions */

/*     blas daxpy,dswap,idamax */
/*     fortran dabs,dmax1,dsqrt */

/*     internal variables */



/*     initialize */

/*     alpha is used in choosing pivot block size. */
    /* Parameter adjustments */
    --kpvt;
    --ap;

    /* Function Body */
    alpha = (sqrt(17.) + 1.) / 8.;

    *info = 0;

/*     main loop on k, which goes from n to 1. */

    k = *n;
    ik = *n * (*n - 1) / 2;
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
    if (ap[1] == 0.) {
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
    kk = ik + k;
    absakk = (d__1 = ap[kk], abs(d__1));

/*        determine the largest off-diagonal element in */
/*        column k. */

    i__1 = k - 1;
    imax = idamax_(&i__1, &ap[ik + 1], &c__1);
    imk = ik + imax;
    colmax = (d__1 = ap[imk], abs(d__1));
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
    im = imax * (imax - 1) / 2;
    imj = im + (imax << 1);
    i__1 = k;
    for (j = imaxp1; j <= i__1; ++j) {
/* Computing MAX */
	d__2 = rowmax, d__3 = (d__1 = ap[imj], abs(d__1));
	rowmax = max(d__2,d__3);
	imj += j;
/* L40: */
    }
    if (imax == 1) {
	goto L50;
    }
    i__1 = imax - 1;
    jmax = idamax_(&i__1, &ap[im + 1], &c__1);
    jmim = jmax + im;
/* Computing MAX */
    d__2 = rowmax, d__3 = (d__1 = ap[jmim], abs(d__1));
    rowmax = max(d__2,d__3);
L50:
    imim = imax + im;
    if ((d__1 = ap[imim], abs(d__1)) < alpha * rowmax) {
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

    dswap_(&imax, &ap[im + 1], &c__1, &ap[ik + 1], &c__1);
    imj = ik + imax;
    i__1 = k;
    for (jj = imax; jj <= i__1; ++jj) {
	j = k + imax - jj;
	jk = ik + j;
	t = ap[jk];
	ap[jk] = ap[imj];
	ap[imj] = t;
	imj -= j - 1;
/* L110: */
    }
L120:

/*           perform the elimination. */

    ij = ik - (k - 1);
    i__1 = km1;
    for (jj = 1; jj <= i__1; ++jj) {
	j = k - jj;
	jk = ik + j;
	mulk = -ap[jk] / ap[kk];
	t = mulk;
	daxpy_(&j, &t, &ap[ik + 1], &c__1, &ap[ij + 1], &c__1);
	ijj = ij + j;
	ap[jk] = mulk;
	ij -= j - 1;
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

    km1k = ik + k - 1;
    ikm1 = ik - (k - 1);
    if (! swap) {
	goto L160;
    }

/*              perform an interchange. */

    dswap_(&imax, &ap[im + 1], &c__1, &ap[ikm1 + 1], &c__1);
    imj = ikm1 + imax;
    i__1 = km1;
    for (jj = imax; jj <= i__1; ++jj) {
	j = km1 + imax - jj;
	jkm1 = ikm1 + j;
	t = ap[jkm1];
	ap[jkm1] = ap[imj];
	ap[imj] = t;
	imj -= j - 1;
/* L150: */
    }
    t = ap[km1k];
    ap[km1k] = ap[imk];
    ap[imk] = t;
L160:

/*           perform the elimination. */

    km2 = k - 2;
    if (km2 == 0) {
	goto L180;
    }
    ak = ap[kk] / ap[km1k];
    km1km1 = ikm1 + k - 1;
    akm1 = ap[km1km1] / ap[km1k];
    denom = 1. - ak * akm1;
    ij = ik - (k - 1) - (k - 2);
    i__1 = km2;
    for (jj = 1; jj <= i__1; ++jj) {
	j = km1 - jj;
	jk = ik + j;
	bk = ap[jk] / ap[km1k];
	jkm1 = ikm1 + j;
	bkm1 = ap[jkm1] / ap[km1k];
	mulk = (akm1 * bk - bkm1) / denom;
	mulkm1 = (ak * bkm1 - bk) / denom;
	t = mulk;
	daxpy_(&j, &t, &ap[ik + 1], &c__1, &ap[ij + 1], &c__1);
	t = mulkm1;
	daxpy_(&j, &t, &ap[ikm1 + 1], &c__1, &ap[ij + 1], &c__1);
	ap[jk] = mulk;
	ap[jkm1] = mulkm1;
	ijj = ij + j;
	ij -= j - 1;
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
    ik -= k - 1;
    if (kstep == 2) {
	ik -= k - 2;
    }
    k -= kstep;
    goto L10;
L200:
    return 0;
} /* dspfa_ */

