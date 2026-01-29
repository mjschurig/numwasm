/* zchdc.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int zchdc_(doublecomplex *a, integer *lda, integer *p, 
	doublecomplex *work, integer *jpvt, integer *job, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    doublereal d__1;
    doublecomplex z__1;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *);
    double sqrt(doublereal);
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    integer j, k, l, kb, jp, pl, jt, pu, km1, kp1, plp1;
    logical negk;
    integer maxl;
    doublecomplex temp;
    logical swapk;
    extern /* Subroutine */ int zswap_(integer *, doublecomplex *, integer *, 
	    doublecomplex *, integer *), zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    doublereal maxdia;


/*     zchdc computes the cholesky decomposition of a positive definite */
/*     matrix.  a pivoting option allows the user to estimate the */
/*     condition of a positive definite matrix or determine the rank */
/*     of a positive semidefinite matrix. */

/*     on entry */

/*         a      complex*16(lda,p). */
/*                a contains the matrix whose decomposition is to */
/*                be computed.  onlt the upper half of a need be stored. */
/*                the lower part of the array a is not referenced. */

/*         lda    integer. */
/*                lda is the leading dimension of the array a. */

/*         p      integer. */
/*                p is the order of the matrix. */

/*         work   complex*16. */
/*                work is a work array. */

/*         jpvt   integer(p). */
/*                jpvt contains integers that control the selection */
/*                of the pivot elements, if pivoting has been requested. */
/*                each diagonal element a(k,k) */
/*                is placed in one of three classes according to the */
/*                value of jpvt(k). */

/*                   if jpvt(k) .gt. 0, then x(k) is an initial */
/*                                      element. */

/*                   if jpvt(k) .eq. 0, then x(k) is a free element. */

/*                   if jpvt(k) .lt. 0, then x(k) is a final element. */

/*                before the decomposition is computed, initial elements */
/*                are moved by symmetric row and column interchanges to */
/*                the beginning of the array a and final */
/*                elements to the end.  both initial and final elements */
/*                are frozen in place during the computation and only */
/*                free elements are moved.  at the k-th stage of the */
/*                reduction, if a(k,k) is occupied by a free element */
/*                it is interchanged with the largest free element */
/*                a(l,l) with l .ge. k.  jpvt is not referenced if */
/*                job .eq. 0. */

/*        job     integer. */
/*                job is an integer that initiates column pivoting. */
/*                if job .eq. 0, no pivoting is done. */
/*                if job .ne. 0, pivoting is done. */

/*     on return */

/*         a      a contains in its upper half the cholesky factor */
/*                of the matrix a as it has been permuted by pivoting. */

/*         jpvt   jpvt(j) contains the index of the diagonal element */
/*                of a that was moved into the j-th position, */
/*                provided pivoting was requested. */

/*         info   contains the index of the last positive diagonal */
/*                element of the cholesky factor. */

/*     for positive definite matrices info = p is the normal return. */
/*     for pivoting with positive semidefinite matrices info will */
/*     in general be less than p.  however, info may be greater than */
/*     the rank of a, since rounding error can cause an otherwise zero */
/*     element to be positive. indefinite systems will always cause */
/*     info to be less than p. */

/*     linpack. this version dated 03/19/79 . */
/*     j.j. dongarra and g.w. stewart, argonne national laboratory and */
/*     university of maryland. */


/*     blas zaxpy,zswap */
/*     fortran dsqrt,dconjg */

/*     internal variables */


    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --work;
    --jpvt;

    /* Function Body */
    pl = 1;
    pu = 0;
    *info = *p;
    if (*job == 0) {
	goto L160;
    }

/*        pivoting has been requested. rearrange the */
/*        the elements according to jpvt. */

    i__1 = *p;
    for (k = 1; k <= i__1; ++k) {
	swapk = jpvt[k] > 0;
	negk = jpvt[k] < 0;
	jpvt[k] = k;
	if (negk) {
	    jpvt[k] = -jpvt[k];
	}
	if (! swapk) {
	    goto L60;
	}
	if (k == pl) {
	    goto L50;
	}
	i__2 = pl - 1;
	zswap_(&i__2, &a[k * a_dim1 + 1], &c__1, &a[pl * a_dim1 + 1], &c__1);
	i__2 = k + k * a_dim1;
	temp.r = a[i__2].r, temp.i = a[i__2].i;
	i__2 = k + k * a_dim1;
	i__3 = pl + pl * a_dim1;
	a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
	i__2 = pl + pl * a_dim1;
	a[i__2].r = temp.r, a[i__2].i = temp.i;
	i__2 = pl + k * a_dim1;
	d_cnjg(&z__1, &a[pl + k * a_dim1]);
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
	plp1 = pl + 1;
	if (*p < plp1) {
	    goto L40;
	}
	i__2 = *p;
	for (j = plp1; j <= i__2; ++j) {
	    if (j >= k) {
		goto L10;
	    }
	    d_cnjg(&z__1, &a[pl + j * a_dim1]);
	    temp.r = z__1.r, temp.i = z__1.i;
	    i__3 = pl + j * a_dim1;
	    d_cnjg(&z__1, &a[j + k * a_dim1]);
	    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
	    i__3 = j + k * a_dim1;
	    a[i__3].r = temp.r, a[i__3].i = temp.i;
	    goto L20;
L10:
	    if (j == k) {
		goto L20;
	    }
	    i__3 = k + j * a_dim1;
	    temp.r = a[i__3].r, temp.i = a[i__3].i;
	    i__3 = k + j * a_dim1;
	    i__4 = pl + j * a_dim1;
	    a[i__3].r = a[i__4].r, a[i__3].i = a[i__4].i;
	    i__3 = pl + j * a_dim1;
	    a[i__3].r = temp.r, a[i__3].i = temp.i;
L20:
/* L30: */
	    ;
	}
L40:
	jpvt[k] = jpvt[pl];
	jpvt[pl] = k;
L50:
	++pl;
L60:
/* L70: */
	;
    }
    pu = *p;
    if (*p < pl) {
	goto L150;
    }
    i__1 = *p;
    for (kb = pl; kb <= i__1; ++kb) {
	k = *p - kb + pl;
	if (jpvt[k] >= 0) {
	    goto L130;
	}
	jpvt[k] = -jpvt[k];
	if (pu == k) {
	    goto L120;
	}
	i__2 = k - 1;
	zswap_(&i__2, &a[k * a_dim1 + 1], &c__1, &a[pu * a_dim1 + 1], &c__1);
	i__2 = k + k * a_dim1;
	temp.r = a[i__2].r, temp.i = a[i__2].i;
	i__2 = k + k * a_dim1;
	i__3 = pu + pu * a_dim1;
	a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
	i__2 = pu + pu * a_dim1;
	a[i__2].r = temp.r, a[i__2].i = temp.i;
	i__2 = k + pu * a_dim1;
	d_cnjg(&z__1, &a[k + pu * a_dim1]);
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
	kp1 = k + 1;
	if (*p < kp1) {
	    goto L110;
	}
	i__2 = *p;
	for (j = kp1; j <= i__2; ++j) {
	    if (j >= pu) {
		goto L80;
	    }
	    d_cnjg(&z__1, &a[k + j * a_dim1]);
	    temp.r = z__1.r, temp.i = z__1.i;
	    i__3 = k + j * a_dim1;
	    d_cnjg(&z__1, &a[j + pu * a_dim1]);
	    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
	    i__3 = j + pu * a_dim1;
	    a[i__3].r = temp.r, a[i__3].i = temp.i;
	    goto L90;
L80:
	    if (j == pu) {
		goto L90;
	    }
	    i__3 = k + j * a_dim1;
	    temp.r = a[i__3].r, temp.i = a[i__3].i;
	    i__3 = k + j * a_dim1;
	    i__4 = pu + j * a_dim1;
	    a[i__3].r = a[i__4].r, a[i__3].i = a[i__4].i;
	    i__3 = pu + j * a_dim1;
	    a[i__3].r = temp.r, a[i__3].i = temp.i;
L90:
/* L100: */
	    ;
	}
L110:
	jt = jpvt[k];
	jpvt[k] = jpvt[pu];
	jpvt[pu] = jt;
L120:
	--pu;
L130:
/* L140: */
	;
    }
L150:
L160:
    i__1 = *p;
    for (k = 1; k <= i__1; ++k) {

/*        reduction loop. */

	i__2 = k + k * a_dim1;
	maxdia = a[i__2].r;
	kp1 = k + 1;
	maxl = k;

/*        determine the pivot element. */

	if (k < pl || k >= pu) {
	    goto L190;
	}
	i__2 = pu;
	for (l = kp1; l <= i__2; ++l) {
	    i__3 = l + l * a_dim1;
	    if (a[i__3].r <= maxdia) {
		goto L170;
	    }
	    i__3 = l + l * a_dim1;
	    maxdia = a[i__3].r;
	    maxl = l;
L170:
/* L180: */
	    ;
	}
L190:

/*        quit if the pivot element is not positive. */

	if (maxdia > 0.) {
	    goto L200;
	}
	*info = k - 1;
/*     ......exit */
	goto L280;
L200:
	if (k == maxl) {
	    goto L210;
	}

/*           start the pivoting and update jpvt. */

	km1 = k - 1;
	zswap_(&km1, &a[k * a_dim1 + 1], &c__1, &a[maxl * a_dim1 + 1], &c__1);
	i__2 = maxl + maxl * a_dim1;
	i__3 = k + k * a_dim1;
	a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
	i__2 = k + k * a_dim1;
	z__1.r = maxdia, z__1.i = 0.;
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
	jp = jpvt[maxl];
	jpvt[maxl] = jpvt[k];
	jpvt[k] = jp;
	i__2 = k + maxl * a_dim1;
	d_cnjg(&z__1, &a[k + maxl * a_dim1]);
	a[i__2].r = z__1.r, a[i__2].i = z__1.i;
L210:

/*        reduction step. pivoting is contained across the rows. */

	i__2 = k;
	d__1 = sqrt((doublereal) a[k + k * a_dim1].r);
	z__1.r = d__1, z__1.i = 0.;
	work[i__2].r = z__1.r, work[i__2].i = z__1.i;
	i__2 = k + k * a_dim1;
	i__3 = k;
	a[i__2].r = work[i__3].r, a[i__2].i = work[i__3].i;
	if (*p < kp1) {
	    goto L260;
	}
	i__2 = *p;
	for (j = kp1; j <= i__2; ++j) {
	    if (k == maxl) {
		goto L240;
	    }
	    if (j >= maxl) {
		goto L220;
	    }
	    d_cnjg(&z__1, &a[k + j * a_dim1]);
	    temp.r = z__1.r, temp.i = z__1.i;
	    i__3 = k + j * a_dim1;
	    d_cnjg(&z__1, &a[j + maxl * a_dim1]);
	    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
	    i__3 = j + maxl * a_dim1;
	    a[i__3].r = temp.r, a[i__3].i = temp.i;
	    goto L230;
L220:
	    if (j == maxl) {
		goto L230;
	    }
	    i__3 = k + j * a_dim1;
	    temp.r = a[i__3].r, temp.i = a[i__3].i;
	    i__3 = k + j * a_dim1;
	    i__4 = maxl + j * a_dim1;
	    a[i__3].r = a[i__4].r, a[i__3].i = a[i__4].i;
	    i__3 = maxl + j * a_dim1;
	    a[i__3].r = temp.r, a[i__3].i = temp.i;
L230:
L240:
	    i__3 = k + j * a_dim1;
	    z_div(&z__1, &a[k + j * a_dim1], &work[k]);
	    a[i__3].r = z__1.r, a[i__3].i = z__1.i;
	    i__3 = j;
	    d_cnjg(&z__1, &a[k + j * a_dim1]);
	    work[i__3].r = z__1.r, work[i__3].i = z__1.i;
	    i__3 = k + j * a_dim1;
	    z__1.r = -a[i__3].r, z__1.i = -a[i__3].i;
	    temp.r = z__1.r, temp.i = z__1.i;
	    i__3 = j - k;
	    zaxpy_(&i__3, &temp, &work[kp1], &c__1, &a[kp1 + j * a_dim1], &
		    c__1);
/* L250: */
	}
L260:
/* L270: */
	;
    }
L280:
    return 0;
} /* zchdc_ */

