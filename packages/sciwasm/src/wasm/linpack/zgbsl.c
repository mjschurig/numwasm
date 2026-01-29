/* zgbsl.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int zgbsl_(doublecomplex *abd, integer *lda, integer *n, 
	integer *ml, integer *mu, integer *ipvt, doublecomplex *b, integer *
	job)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3;
    doublecomplex z__1, z__2, z__3;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    integer k, l, m;
    doublecomplex t;
    integer kb, la, lb, lm, nm1;
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);


/*     zgbsl solves the complex*16 band system */
/*     a * x = b  or  ctrans(a) * x = b */
/*     using the factors computed by zgbco or zgbfa. */

/*     on entry */

/*        abd     complex*16(lda, n) */
/*                the output from zgbco or zgbfa. */

/*        lda     integer */
/*                the leading dimension of the array  abd . */

/*        n       integer */
/*                the order of the original matrix. */

/*        ml      integer */
/*                number of diagonals below the main diagonal. */

/*        mu      integer */
/*                number of diagonals above the main diagonal. */

/*        ipvt    integer(n) */
/*                the pivot vector from zgbco or zgbfa. */

/*        b       complex*16(n) */
/*                the right hand side vector. */

/*        job     integer */
/*                = 0         to solve  a*x = b , */
/*                = nonzero   to solve  ctrans(a)*x = b , where */
/*                            ctrans(a)  is the conjugate transpose. */

/*     on return */

/*        b       the solution vector  x . */

/*     error condition */

/*        a division by zero will occur if the input factor contains a */
/*        zero on the diagonal.  technically this indicates singularity */
/*        but it is often caused by improper arguments or improper */
/*        setting of lda .  it will not occur if the subroutines are */
/*        called correctly and if zgbco has set rcond .gt. 0.0 */
/*        or zgbfa has set info .eq. 0 . */

/*     to compute  inverse(a) * c  where  c  is a matrix */
/*     with  p  columns */
/*           call zgbco(abd,lda,n,ml,mu,ipvt,rcond,z) */
/*           if (rcond is too small) go to ... */
/*           do 10 j = 1, p */
/*              call zgbsl(abd,lda,n,ml,mu,ipvt,c(1,j),0) */
/*        10 continue */

/*     linpack. this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas zaxpy,zdotc */
/*     fortran dconjg,min0 */

/*     internal variables */


    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --ipvt;
    --b;

    /* Function Body */
    m = *mu + *ml + 1;
    nm1 = *n - 1;
    if (*job != 0) {
	goto L50;
    }

/*        job = 0 , solve  a * x = b */
/*        first solve l*y = b */

    if (*ml == 0) {
	goto L30;
    }
    if (nm1 < 1) {
	goto L30;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	lm = min(i__2,i__3);
	l = ipvt[k];
	i__2 = l;
	t.r = b[i__2].r, t.i = b[i__2].i;
	if (l == k) {
	    goto L10;
	}
	i__2 = l;
	i__3 = k;
	b[i__2].r = b[i__3].r, b[i__2].i = b[i__3].i;
	i__2 = k;
	b[i__2].r = t.r, b[i__2].i = t.i;
L10:
	zaxpy_(&lm, &t, &abd[m + 1 + k * abd_dim1], &c__1, &b[k + 1], &c__1);
/* L20: */
    }
L30:

/*        now solve  u*x = y */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	i__2 = k;
	z_div(&z__1, &b[k], &abd[m + k * abd_dim1]);
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
	lm = min(k,m) - 1;
	la = m - lm;
	lb = k - lm;
	i__2 = k;
	z__1.r = -b[i__2].r, z__1.i = -b[i__2].i;
	t.r = z__1.r, t.i = z__1.i;
	zaxpy_(&lm, &t, &abd[la + k * abd_dim1], &c__1, &b[lb], &c__1);
/* L40: */
    }
    goto L100;
L50:

/*        job = nonzero, solve  ctrans(a) * x = b */
/*        first solve  ctrans(u)*y = b */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	lm = min(k,m) - 1;
	la = m - lm;
	lb = k - lm;
	zdotc_(&z__1, &lm, &abd[la + k * abd_dim1], &c__1, &b[lb], &c__1);
	t.r = z__1.r, t.i = z__1.i;
	i__2 = k;
	i__3 = k;
	z__2.r = b[i__3].r - t.r, z__2.i = b[i__3].i - t.i;
	d_cnjg(&z__3, &abd[m + k * abd_dim1]);
	z_div(&z__1, &z__2, &z__3);
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
/* L60: */
    }

/*        now solve ctrans(l)*x = y */

    if (*ml == 0) {
	goto L90;
    }
    if (nm1 < 1) {
	goto L90;
    }
    i__1 = nm1;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n - kb;
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	lm = min(i__2,i__3);
	i__2 = k;
	i__3 = k;
	zdotc_(&z__2, &lm, &abd[m + 1 + k * abd_dim1], &c__1, &b[k + 1], &
		c__1);
	z__1.r = b[i__3].r + z__2.r, z__1.i = b[i__3].i + z__2.i;
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
	l = ipvt[k];
	if (l == k) {
	    goto L70;
	}
	i__2 = l;
	t.r = b[i__2].r, t.i = b[i__2].i;
	i__2 = l;
	i__3 = k;
	b[i__2].r = b[i__3].r, b[i__2].i = b[i__3].i;
	i__2 = k;
	b[i__2].r = t.r, b[i__2].i = t.i;
L70:
/* L80: */
	;
    }
L90:
L100:
    return 0;
} /* zgbsl_ */

