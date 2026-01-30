/* dchdc.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dchdc_(doublereal *a, integer *lda, integer *p, 
	doublereal *work, integer *jpvt, integer *job, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer j, k, l, kb, jp, pl, jt, pu, km1, kp1, plp1;
    logical negk;
    integer maxl;
    doublereal temp;
    extern /* Subroutine */ int dswap_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    logical swapk;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    doublereal maxdia;


/*     dchdc computes the cholesky decomposition of a positive definite */
/*     matrix.  a pivoting option allows the user to estimate the */
/*     condition of a positive definite matrix or determine the rank */
/*     of a positive semidefinite matrix. */

/*     on entry */

/*         a      double precision(lda,p). */
/*                a contains the matrix whose decomposition is to */
/*                be computed.  onlt the upper half of a need be stored. */
/*                the lower part of the array a is not referenced. */

/*         lda    integer. */
/*                lda is the leading dimension of the array a. */

/*         p      integer. */
/*                p is the order of the matrix. */

/*         work   double precision. */
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

/*     linpack. this version dated 08/14/78 . */
/*     j.j. dongarra and g.w. stewart, argonne national laboratory and */
/*     university of maryland. */


/*     blas daxpy,dswap */
/*     fortran dsqrt */

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
	dswap_(&i__2, &a[k * a_dim1 + 1], &c__1, &a[pl * a_dim1 + 1], &c__1);
	temp = a[k + k * a_dim1];
	a[k + k * a_dim1] = a[pl + pl * a_dim1];
	a[pl + pl * a_dim1] = temp;
	plp1 = pl + 1;
	if (*p < plp1) {
	    goto L40;
	}
	i__2 = *p;
	for (j = plp1; j <= i__2; ++j) {
	    if (j >= k) {
		goto L10;
	    }
	    temp = a[pl + j * a_dim1];
	    a[pl + j * a_dim1] = a[j + k * a_dim1];
	    a[j + k * a_dim1] = temp;
	    goto L20;
L10:
	    if (j == k) {
		goto L20;
	    }
	    temp = a[k + j * a_dim1];
	    a[k + j * a_dim1] = a[pl + j * a_dim1];
	    a[pl + j * a_dim1] = temp;
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
	dswap_(&i__2, &a[k * a_dim1 + 1], &c__1, &a[pu * a_dim1 + 1], &c__1);
	temp = a[k + k * a_dim1];
	a[k + k * a_dim1] = a[pu + pu * a_dim1];
	a[pu + pu * a_dim1] = temp;
	kp1 = k + 1;
	if (*p < kp1) {
	    goto L110;
	}
	i__2 = *p;
	for (j = kp1; j <= i__2; ++j) {
	    if (j >= pu) {
		goto L80;
	    }
	    temp = a[k + j * a_dim1];
	    a[k + j * a_dim1] = a[j + pu * a_dim1];
	    a[j + pu * a_dim1] = temp;
	    goto L90;
L80:
	    if (j == pu) {
		goto L90;
	    }
	    temp = a[k + j * a_dim1];
	    a[k + j * a_dim1] = a[pu + j * a_dim1];
	    a[pu + j * a_dim1] = temp;
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

	maxdia = a[k + k * a_dim1];
	kp1 = k + 1;
	maxl = k;

/*        determine the pivot element. */

	if (k < pl || k >= pu) {
	    goto L190;
	}
	i__2 = pu;
	for (l = kp1; l <= i__2; ++l) {
	    if (a[l + l * a_dim1] <= maxdia) {
		goto L170;
	    }
	    maxdia = a[l + l * a_dim1];
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
	dswap_(&km1, &a[k * a_dim1 + 1], &c__1, &a[maxl * a_dim1 + 1], &c__1);
	a[maxl + maxl * a_dim1] = a[k + k * a_dim1];
	a[k + k * a_dim1] = maxdia;
	jp = jpvt[maxl];
	jpvt[maxl] = jpvt[k];
	jpvt[k] = jp;
L210:

/*        reduction step. pivoting is contained across the rows. */

	work[k] = sqrt(a[k + k * a_dim1]);
	a[k + k * a_dim1] = work[k];
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
	    temp = a[k + j * a_dim1];
	    a[k + j * a_dim1] = a[j + maxl * a_dim1];
	    a[j + maxl * a_dim1] = temp;
	    goto L230;
L220:
	    if (j == maxl) {
		goto L230;
	    }
	    temp = a[k + j * a_dim1];
	    a[k + j * a_dim1] = a[maxl + j * a_dim1];
	    a[maxl + j * a_dim1] = temp;
L230:
L240:
	    a[k + j * a_dim1] /= work[k];
	    work[j] = a[k + j * a_dim1];
	    temp = -a[k + j * a_dim1];
	    i__3 = j - k;
	    daxpy_(&i__3, &temp, &work[kp1], &c__1, &a[kp1 + j * a_dim1], &
		    c__1);
/* L250: */
	}
L260:
/* L270: */
	;
    }
L280:
    return 0;
} /* dchdc_ */

