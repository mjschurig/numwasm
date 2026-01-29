/* chpfa.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int chpfa_(complex *ap, integer *n, integer *kpvt, integer *
	info)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    real r__1, r__2, r__3, r__4;
    complex q__1, q__2, q__3;

    /* Builtin functions */
    double sqrt(doublereal), r_imag(complex *);
    void r_cnjg(complex *, complex *), c_div(complex *, complex *, complex *);

    /* Local variables */
    integer j, k;
    complex t, ak, bk;
    integer ij, ik, jj, im, jk, kk, km1, km2, ijj, imj, imk;
    complex akm1, bkm1;
    integer ikm1, jkm1, km1k, imim, jmim, imax, jmax;
    complex mulk;
    logical swap;
    real alpha;
    integer km1km1;
    complex denom;
    extern /* Subroutine */ int cswap_(integer *, complex *, integer *, 
	    complex *, integer *), caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *);
    integer kstep, imaxp1;
    complex mulkm1;
    real absakk;
    extern integer icamax_(integer *, complex *, integer *);
    real colmax, rowmax;


/*     chpfa factors a complex hermitian matrix stored in */
/*     packed form by elimination with symmetric pivoting. */

/*     to solve  a*x = b , follow chpfa by chpsl. */
/*     to compute  inverse(a)*c , follow chpfa by chpsl. */
/*     to compute  determinant(a) , follow chpfa by chpdi. */
/*     to compute  inertia(a) , follow chpfa by chpdi. */
/*     to compute  inverse(a) , follow chpfa by chpdi. */

/*     on entry */

/*        ap      complex (n*(n+1)/2) */
/*                the packed form of a hermitian matrix  a .  the */
/*                columns of the upper triangle are stored sequentially */
/*                in a one-dimensional array of length  n*(n+1)/2 . */
/*                see comments below for details. */

/*        n       integer */
/*                the order of the matrix  a . */

/*     output */

/*        ap      a block diagonal matrix and the multipliers which */
/*                were used to obtain it stored in packed form. */
/*                the factorization can be written  a = u*d*ctrans(u) */
/*                where  u  is a product of permutation and unit */
/*                upper triangular matrices , ctrans(u) is the */
/*                conjugate transpose of  u , and  d  is block diagonal */
/*                with 1 by 1 and 2 by 2 blocks. */

/*        kpvt    integer(n) */
/*                an integer vector of pivot indices. */

/*        info    integer */
/*                = 0  normal value. */
/*                = k  if the k-th pivot block is singular. this is */
/*                     not an error condition for this subroutine, */
/*                     but it does indicate that chpsl or chpdi may */
/*                     divide by zero if called. */

/*     packed storage */

/*          the following program segment will pack the upper */
/*          triangle of a hermitian matrix. */

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

/*     blas caxpy,cswap,icamax */
/*     fortran abs,aimag,amax1,cmplx,conjg,real,sqrt */

/*     internal variables */



/*     initialize */

/*     alpha is used in choosing pivot block size. */
    /* Parameter adjustments */
    --kpvt;
    --ap;

    /* Function Body */
    alpha = (sqrt(17.f) + 1.f) / 8.f;

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
    if ((r__1 = ap[1].r, dabs(r__1)) + (r__2 = r_imag(&ap[1]), dabs(r__2)) == 
	    0.f) {
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
    absakk = (r__1 = ap[i__1].r, dabs(r__1)) + (r__2 = r_imag(&ap[kk]), dabs(
	    r__2));

/*        determine the largest off-diagonal element in */
/*        column k. */

    i__1 = k - 1;
    imax = icamax_(&i__1, &ap[ik + 1], &c__1);
    imk = ik + imax;
    i__1 = imk;
    colmax = (r__1 = ap[i__1].r, dabs(r__1)) + (r__2 = r_imag(&ap[imk]), dabs(
	    r__2));
    if (absakk < alpha * colmax) {
	goto L30;
    }
    kstep = 1;
    swap = FALSE_;
    goto L90;
L30:

/*           determine the largest off-diagonal element in */
/*           row imax. */

    rowmax = 0.f;
    imaxp1 = imax + 1;
    im = imax * (imax - 1) / 2;
    imj = im + (imax << 1);
    i__1 = k;
    for (j = imaxp1; j <= i__1; ++j) {
/* Computing MAX */
	i__2 = imj;
	r__3 = rowmax, r__4 = (r__1 = ap[i__2].r, dabs(r__1)) + (r__2 = 
		r_imag(&ap[imj]), dabs(r__2));
	rowmax = dmax(r__3,r__4);
	imj += j;
/* L40: */
    }
    if (imax == 1) {
	goto L50;
    }
    i__1 = imax - 1;
    jmax = icamax_(&i__1, &ap[im + 1], &c__1);
    jmim = jmax + im;
/* Computing MAX */
    i__1 = jmim;
    r__3 = rowmax, r__4 = (r__1 = ap[i__1].r, dabs(r__1)) + (r__2 = r_imag(&
	    ap[jmim]), dabs(r__2));
    rowmax = dmax(r__3,r__4);
L50:
    imim = imax + im;
    i__1 = imim;
    if ((r__1 = ap[i__1].r, dabs(r__1)) + (r__2 = r_imag(&ap[imim]), dabs(
	    r__2)) < alpha * rowmax) {
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
    if (dmax(absakk,colmax) != 0.f) {
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

    cswap_(&imax, &ap[im + 1], &c__1, &ap[ik + 1], &c__1);
    imj = ik + imax;
    i__1 = k;
    for (jj = imax; jj <= i__1; ++jj) {
	j = k + imax - jj;
	jk = ik + j;
	r_cnjg(&q__1, &ap[jk]);
	t.r = q__1.r, t.i = q__1.i;
	i__2 = jk;
	r_cnjg(&q__1, &ap[imj]);
	ap[i__2].r = q__1.r, ap[i__2].i = q__1.i;
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
	q__2.r = -ap[i__2].r, q__2.i = -ap[i__2].i;
	c_div(&q__1, &q__2, &ap[kk]);
	mulk.r = q__1.r, mulk.i = q__1.i;
	r_cnjg(&q__1, &mulk);
	t.r = q__1.r, t.i = q__1.i;
	caxpy_(&j, &t, &ap[ik + 1], &c__1, &ap[ij + 1], &c__1);
	ijj = ij + j;
	i__2 = ijj;
	i__3 = ijj;
	r__1 = ap[i__3].r;
	q__1.r = r__1, q__1.i = 0.f;
	ap[i__2].r = q__1.r, ap[i__2].i = q__1.i;
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

    cswap_(&imax, &ap[im + 1], &c__1, &ap[ikm1 + 1], &c__1);
    imj = ikm1 + imax;
    i__1 = km1;
    for (jj = imax; jj <= i__1; ++jj) {
	j = km1 + imax - jj;
	jkm1 = ikm1 + j;
	r_cnjg(&q__1, &ap[jkm1]);
	t.r = q__1.r, t.i = q__1.i;
	i__2 = jkm1;
	r_cnjg(&q__1, &ap[imj]);
	ap[i__2].r = q__1.r, ap[i__2].i = q__1.i;
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
    c_div(&q__1, &ap[kk], &ap[km1k]);
    ak.r = q__1.r, ak.i = q__1.i;
    km1km1 = ikm1 + k - 1;
    r_cnjg(&q__2, &ap[km1k]);
    c_div(&q__1, &ap[km1km1], &q__2);
    akm1.r = q__1.r, akm1.i = q__1.i;
    q__2.r = ak.r * akm1.r - ak.i * akm1.i, q__2.i = ak.r * akm1.i + ak.i * 
	    akm1.r;
    q__1.r = 1.f - q__2.r, q__1.i = -q__2.i;
    denom.r = q__1.r, denom.i = q__1.i;
    ij = ik - (k - 1) - (k - 2);
    i__1 = km2;
    for (jj = 1; jj <= i__1; ++jj) {
	j = km1 - jj;
	jk = ik + j;
	c_div(&q__1, &ap[jk], &ap[km1k]);
	bk.r = q__1.r, bk.i = q__1.i;
	jkm1 = ikm1 + j;
	r_cnjg(&q__2, &ap[km1k]);
	c_div(&q__1, &ap[jkm1], &q__2);
	bkm1.r = q__1.r, bkm1.i = q__1.i;
	q__3.r = akm1.r * bk.r - akm1.i * bk.i, q__3.i = akm1.r * bk.i + 
		akm1.i * bk.r;
	q__2.r = q__3.r - bkm1.r, q__2.i = q__3.i - bkm1.i;
	c_div(&q__1, &q__2, &denom);
	mulk.r = q__1.r, mulk.i = q__1.i;
	q__3.r = ak.r * bkm1.r - ak.i * bkm1.i, q__3.i = ak.r * bkm1.i + ak.i 
		* bkm1.r;
	q__2.r = q__3.r - bk.r, q__2.i = q__3.i - bk.i;
	c_div(&q__1, &q__2, &denom);
	mulkm1.r = q__1.r, mulkm1.i = q__1.i;
	r_cnjg(&q__1, &mulk);
	t.r = q__1.r, t.i = q__1.i;
	caxpy_(&j, &t, &ap[ik + 1], &c__1, &ap[ij + 1], &c__1);
	r_cnjg(&q__1, &mulkm1);
	t.r = q__1.r, t.i = q__1.i;
	caxpy_(&j, &t, &ap[ikm1 + 1], &c__1, &ap[ij + 1], &c__1);
	i__2 = jk;
	ap[i__2].r = mulk.r, ap[i__2].i = mulk.i;
	i__2 = jkm1;
	ap[i__2].r = mulkm1.r, ap[i__2].i = mulkm1.i;
	ijj = ij + j;
	i__2 = ijj;
	i__3 = ijj;
	r__1 = ap[i__3].r;
	q__1.r = r__1, q__1.i = 0.f;
	ap[i__2].r = q__1.r, ap[i__2].i = q__1.i;
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
} /* chpfa_ */

