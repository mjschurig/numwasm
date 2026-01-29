/* zspfa.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int zspfa_(doublecomplex *ap, integer *n, integer *kpvt, 
	integer *info)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublereal d__1, d__2, d__3, d__4;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    double sqrt(doublereal);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    integer j, k;
    doublecomplex t, ak, bk;
    integer ij, ik, jj, im, jk, kk, km1, km2, ijj, imj, imk;
    doublecomplex akm1, bkm1;
    integer ikm1, jkm1, km1k, imim, jmim, imax, jmax;
    doublecomplex mulk;
    logical swap;
    doublereal alpha;
    integer km1km1;
    doublecomplex denom;
    integer kstep;
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    integer imaxp1;
    doublecomplex mulkm1;
    doublereal absakk, colmax;
    extern integer izamax_(integer *, doublecomplex *, integer *);
    doublereal rowmax;


/*     zspfa factors a complex*16 symmetric matrix stored in */
/*     packed form by elimination with symmetric pivoting. */

/*     to solve  a*x = b , follow zspfa by zspsl. */
/*     to compute  inverse(a)*c , follow zspfa by zspsl. */
/*     to compute  determinant(a) , follow zspfa by zspdi. */
/*     to compute  inverse(a) , follow zspfa by zspdi. */

/*     on entry */

/*        ap      complex*16 (n*(n+1)/2) */
/*                the packed form of a symmetric matrix  a .  the */
/*                columns of the upper triangle are stored sequentially */
/*                in a one-dimensional array of length  n*(n+1)/2 . */
/*                see comments below for details. */

/*        n       integer */
/*                the order of the matrix  a . */

/*     on return */

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
/*                     but it does indicate that zspsl or zspdi may */
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

/*     blas zaxpy,zswap,izamax */
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
    z__1.r = ap[1].r * 0. - ap[1].i * -1., z__1.i = ap[1].i * 0. + ap[1].r * 
	    -1.;
    if ((d__1 = ap[1].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) == 0.) {
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
    i__1 = kk;
    i__2 = kk;
    z__1.r = ap[i__2].r * 0. - ap[i__2].i * -1., z__1.i = ap[i__2].i * 0. + 
	    ap[i__2].r * -1.;
    absakk = (d__1 = ap[i__1].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2));

/*        determine the largest off-diagonal element in */
/*        column k. */

    i__1 = k - 1;
    imax = izamax_(&i__1, &ap[ik + 1], &c__1);
    imk = ik + imax;
    i__1 = imk;
    i__2 = imk;
    z__1.r = ap[i__2].r * 0. - ap[i__2].i * -1., z__1.i = ap[i__2].i * 0. + 
	    ap[i__2].r * -1.;
    colmax = (d__1 = ap[i__1].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2));
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
	i__2 = imj;
	i__3 = imj;
	z__1.r = ap[i__3].r * 0. - ap[i__3].i * -1., z__1.i = ap[i__3].i * 0. 
		+ ap[i__3].r * -1.;
	d__3 = rowmax, d__4 = (d__1 = ap[i__2].r, abs(d__1)) + (d__2 = z__1.r,
		 abs(d__2));
	rowmax = max(d__3,d__4);
	imj += j;
/* L40: */
    }
    if (imax == 1) {
	goto L50;
    }
    i__1 = imax - 1;
    jmax = izamax_(&i__1, &ap[im + 1], &c__1);
    jmim = jmax + im;
/* Computing MAX */
    i__1 = jmim;
    i__2 = jmim;
    z__1.r = ap[i__2].r * 0. - ap[i__2].i * -1., z__1.i = ap[i__2].i * 0. + 
	    ap[i__2].r * -1.;
    d__3 = rowmax, d__4 = (d__1 = ap[i__1].r, abs(d__1)) + (d__2 = z__1.r, 
	    abs(d__2));
    rowmax = max(d__3,d__4);
L50:
    imim = imax + im;
    i__1 = imim;
    i__2 = imim;
    z__1.r = ap[i__2].r * 0. - ap[i__2].i * -1., z__1.i = ap[i__2].i * 0. + 
	    ap[i__2].r * -1.;
    if ((d__1 = ap[i__1].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) < alpha * 
	    rowmax) {
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

    zswap_(&imax, &ap[im + 1], &c__1, &ap[ik + 1], &c__1);
    imj = ik + imax;
    i__1 = k;
    for (jj = imax; jj <= i__1; ++jj) {
	j = k + imax - jj;
	jk = ik + j;
	i__2 = jk;
	t.r = ap[i__2].r, t.i = ap[i__2].i;
	i__2 = jk;
	i__3 = imj;
	ap[i__2].r = ap[i__3].r, ap[i__2].i = ap[i__3].i;
	i__2 = imj;
	ap[i__2].r = t.r, ap[i__2].i = t.i;
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
	i__2 = jk;
	z__2.r = -ap[i__2].r, z__2.i = -ap[i__2].i;
	z_div(&z__1, &z__2, &ap[kk]);
	mulk.r = z__1.r, mulk.i = z__1.i;
	t.r = mulk.r, t.i = mulk.i;
	zaxpy_(&j, &t, &ap[ik + 1], &c__1, &ap[ij + 1], &c__1);
	ijj = ij + j;
	i__2 = jk;
	ap[i__2].r = mulk.r, ap[i__2].i = mulk.i;
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

    zswap_(&imax, &ap[im + 1], &c__1, &ap[ikm1 + 1], &c__1);
    imj = ikm1 + imax;
    i__1 = km1;
    for (jj = imax; jj <= i__1; ++jj) {
	j = km1 + imax - jj;
	jkm1 = ikm1 + j;
	i__2 = jkm1;
	t.r = ap[i__2].r, t.i = ap[i__2].i;
	i__2 = jkm1;
	i__3 = imj;
	ap[i__2].r = ap[i__3].r, ap[i__2].i = ap[i__3].i;
	i__2 = imj;
	ap[i__2].r = t.r, ap[i__2].i = t.i;
	imj -= j - 1;
/* L150: */
    }
    i__1 = km1k;
    t.r = ap[i__1].r, t.i = ap[i__1].i;
    i__1 = km1k;
    i__2 = imk;
    ap[i__1].r = ap[i__2].r, ap[i__1].i = ap[i__2].i;
    i__1 = imk;
    ap[i__1].r = t.r, ap[i__1].i = t.i;
L160:

/*           perform the elimination. */

    km2 = k - 2;
    if (km2 == 0) {
	goto L180;
    }
    z_div(&z__1, &ap[kk], &ap[km1k]);
    ak.r = z__1.r, ak.i = z__1.i;
    km1km1 = ikm1 + k - 1;
    z_div(&z__1, &ap[km1km1], &ap[km1k]);
    akm1.r = z__1.r, akm1.i = z__1.i;
    z__2.r = ak.r * akm1.r - ak.i * akm1.i, z__2.i = ak.r * akm1.i + ak.i * 
	    akm1.r;
    z__1.r = 1. - z__2.r, z__1.i = -z__2.i;
    denom.r = z__1.r, denom.i = z__1.i;
    ij = ik - (k - 1) - (k - 2);
    i__1 = km2;
    for (jj = 1; jj <= i__1; ++jj) {
	j = km1 - jj;
	jk = ik + j;
	z_div(&z__1, &ap[jk], &ap[km1k]);
	bk.r = z__1.r, bk.i = z__1.i;
	jkm1 = ikm1 + j;
	z_div(&z__1, &ap[jkm1], &ap[km1k]);
	bkm1.r = z__1.r, bkm1.i = z__1.i;
	z__3.r = akm1.r * bk.r - akm1.i * bk.i, z__3.i = akm1.r * bk.i + 
		akm1.i * bk.r;
	z__2.r = z__3.r - bkm1.r, z__2.i = z__3.i - bkm1.i;
	z_div(&z__1, &z__2, &denom);
	mulk.r = z__1.r, mulk.i = z__1.i;
	z__3.r = ak.r * bkm1.r - ak.i * bkm1.i, z__3.i = ak.r * bkm1.i + ak.i 
		* bkm1.r;
	z__2.r = z__3.r - bk.r, z__2.i = z__3.i - bk.i;
	z_div(&z__1, &z__2, &denom);
	mulkm1.r = z__1.r, mulkm1.i = z__1.i;
	t.r = mulk.r, t.i = mulk.i;
	zaxpy_(&j, &t, &ap[ik + 1], &c__1, &ap[ij + 1], &c__1);
	t.r = mulkm1.r, t.i = mulkm1.i;
	zaxpy_(&j, &t, &ap[ikm1 + 1], &c__1, &ap[ij + 1], &c__1);
	i__2 = jk;
	ap[i__2].r = mulk.r, ap[i__2].i = mulk.i;
	i__2 = jkm1;
	ap[i__2].r = mulkm1.r, ap[i__2].i = mulkm1.i;
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
} /* zspfa_ */

