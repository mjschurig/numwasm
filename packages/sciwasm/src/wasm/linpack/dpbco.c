/* dpbco.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dpbco_(doublereal *abd, integer *lda, integer *n, 
	integer *m, doublereal *rcond, doublereal *z__, integer *info)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    integer i__, j, k, l;
    doublereal s, t;
    integer j2, kb, la, lb;
    doublereal ek;
    integer lm;
    doublereal sm, wk;
    integer mu, kp1;
    doublereal wkm;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int dpbfa_(doublereal *, integer *, integer *, 
	    integer *, integer *), dscal_(integer *, doublereal *, doublereal 
	    *, integer *);
    extern doublereal dasum_(integer *, doublereal *, integer *);
    doublereal anorm;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    doublereal ynorm;


/*     dpbco factors a double precision symmetric positive definite */
/*     matrix stored in band form and estimates the condition of the */
/*     matrix. */

/*     if  rcond  is not needed, dpbfa is slightly faster. */
/*     to solve  a*x = b , follow dpbco by dpbsl. */
/*     to compute  inverse(a)*c , follow dpbco by dpbsl. */
/*     to compute  determinant(a) , follow dpbco by dpbdi. */

/*     on entry */

/*        abd     double precision(lda, n) */
/*                the matrix to be factored.  the columns of the upper */
/*                triangle are stored in the columns of abd and the */
/*                diagonals of the upper triangle are stored in the */
/*                rows of abd .  see the comments below for details. */

/*        lda     integer */
/*                the leading dimension of the array  abd . */
/*                lda must be .ge. m + 1 . */

/*        n       integer */
/*                the order of the matrix  a . */

/*        m       integer */
/*                the number of diagonals above the main diagonal. */
/*                0 .le. m .lt. n . */

/*     on return */

/*        abd     an upper triangular matrix  r , stored in band */
/*                form, so that  a = trans(r)*r . */
/*                if  info .ne. 0 , the factorization is not complete. */

/*        rcond   double precision */
/*                an estimate of the reciprocal condition of  a . */
/*                for the system  a*x = b , relative perturbations */
/*                in  a  and  b  of size  epsilon  may cause */
/*                relative perturbations in  x  of size  epsilon/rcond . */
/*                if  rcond  is so small that the logical expression */
/*                           1.0 + rcond .eq. 1.0 */
/*                is true, then  a  may be singular to working */
/*                precision.  in particular,  rcond  is zero  if */
/*                exact singularity is detected or the estimate */
/*                underflows.  if info .ne. 0 , rcond is unchanged. */

/*        z       double precision(n) */
/*                a work vector whose contents are usually unimportant. */
/*                if  a  is singular to working precision, then  z  is */
/*                an approximate null vector in the sense that */
/*                norm(a*z) = rcond*norm(a)*norm(z) . */
/*                if  info .ne. 0 , z  is unchanged. */

/*        info    integer */
/*                = 0  for normal return. */
/*                = k  signals an error condition.  the leading minor */
/*                     of order  k  is not positive definite. */

/*     band storage */

/*           if  a  is a symmetric positive definite band matrix, */
/*           the following program segment will set up the input. */

/*                   m = (band width above diagonal) */
/*                   do 20 j = 1, n */
/*                      i1 = max0(1, j-m) */
/*                      do 10 i = i1, j */
/*                         k = i-j+m+1 */
/*                         abd(k,j) = a(i,j) */
/*                10    continue */
/*                20 continue */

/*           this uses  m + 1  rows of  a , except for the  m by m */
/*           upper left triangle, which is ignored. */

/*     example..  if the original matrix is */

/*           11 12 13  0  0  0 */
/*           12 22 23 24  0  0 */
/*           13 23 33 34 35  0 */
/*            0 24 34 44 45 46 */
/*            0  0 35 45 55 56 */
/*            0  0  0 46 56 66 */

/*     then  n = 6 , m = 2  and  abd  should contain */

/*            *  * 13 24 35 46 */
/*            * 12 23 34 45 56 */
/*           11 22 33 44 55 66 */

/*     linpack.  this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     linpack dpbfa */
/*     blas daxpy,ddot,dscal,dasum */
/*     fortran dabs,dmax1,max0,min0,dreal,dsign */

/*     internal variables */



/*     find norm of a */

    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --z__;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MIN */
	i__2 = j, i__3 = *m + 1;
	l = min(i__2,i__3);
/* Computing MAX */
	i__2 = *m + 2 - j;
	mu = max(i__2,1);
	z__[j] = dasum_(&l, &abd[mu + j * abd_dim1], &c__1);
	k = j - l;
	if (*m < mu) {
	    goto L20;
	}
	i__2 = *m;
	for (i__ = mu; i__ <= i__2; ++i__) {
	    ++k;
	    z__[k] += (d__1 = abd[i__ + j * abd_dim1], abs(d__1));
/* L10: */
	}
L20:
/* L30: */
	;
    }
    anorm = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	d__1 = anorm, d__2 = z__[j];
	anorm = max(d__1,d__2);
/* L40: */
    }

/*     factor */

    dpbfa_(&abd[abd_offset], lda, n, m, info);
    if (*info != 0) {
	goto L180;
    }

/*        rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) . */
/*        estimate = norm(z)/norm(y) where  a*z = y  and  a*y = e . */
/*        the components of  e  are chosen to cause maximum local */
/*        growth in the elements of w  where  trans(r)*w = e . */
/*        the vectors are frequently rescaled to avoid overflow. */

/*        solve trans(r)*w = e */

    ek = 1.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z__[j] = 0.;
/* L50: */
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (z__[k] != 0.) {
	    d__1 = -z__[k];
	    ek = d_sign(&ek, &d__1);
	}
	if ((d__1 = ek - z__[k], abs(d__1)) <= abd[*m + 1 + k * abd_dim1]) {
	    goto L60;
	}
	s = abd[*m + 1 + k * abd_dim1] / (d__1 = ek - z__[k], abs(d__1));
	dscal_(n, &s, &z__[1], &c__1);
	ek = s * ek;
L60:
	wk = ek - z__[k];
	wkm = -ek - z__[k];
	s = abs(wk);
	sm = abs(wkm);
	wk /= abd[*m + 1 + k * abd_dim1];
	wkm /= abd[*m + 1 + k * abd_dim1];
	kp1 = k + 1;
/* Computing MIN */
	i__2 = k + *m;
	j2 = min(i__2,*n);
	i__ = *m + 1;
	if (kp1 > j2) {
	    goto L100;
	}
	i__2 = j2;
	for (j = kp1; j <= i__2; ++j) {
	    --i__;
	    sm += (d__1 = z__[j] + wkm * abd[i__ + j * abd_dim1], abs(d__1));
	    z__[j] += wk * abd[i__ + j * abd_dim1];
	    s += (d__1 = z__[j], abs(d__1));
/* L70: */
	}
	if (s >= sm) {
	    goto L90;
	}
	t = wkm - wk;
	wk = wkm;
	i__ = *m + 1;
	i__2 = j2;
	for (j = kp1; j <= i__2; ++j) {
	    --i__;
	    z__[j] += t * abd[i__ + j * abd_dim1];
/* L80: */
	}
L90:
L100:
	z__[k] = wk;
/* L110: */
    }
    s = 1. / dasum_(n, &z__[1], &c__1);
    dscal_(n, &s, &z__[1], &c__1);

/*        solve  r*y = w */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if ((d__1 = z__[k], abs(d__1)) <= abd[*m + 1 + k * abd_dim1]) {
	    goto L120;
	}
	s = abd[*m + 1 + k * abd_dim1] / (d__1 = z__[k], abs(d__1));
	dscal_(n, &s, &z__[1], &c__1);
L120:
	z__[k] /= abd[*m + 1 + k * abd_dim1];
/* Computing MIN */
	i__2 = k - 1;
	lm = min(i__2,*m);
	la = *m + 1 - lm;
	lb = k - lm;
	t = -z__[k];
	daxpy_(&lm, &t, &abd[la + k * abd_dim1], &c__1, &z__[lb], &c__1);
/* L130: */
    }
    s = 1. / dasum_(n, &z__[1], &c__1);
    dscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.;

/*        solve trans(r)*v = y */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* Computing MIN */
	i__2 = k - 1;
	lm = min(i__2,*m);
	la = *m + 1 - lm;
	lb = k - lm;
	z__[k] -= ddot_(&lm, &abd[la + k * abd_dim1], &c__1, &z__[lb], &c__1);
	if ((d__1 = z__[k], abs(d__1)) <= abd[*m + 1 + k * abd_dim1]) {
	    goto L140;
	}
	s = abd[*m + 1 + k * abd_dim1] / (d__1 = z__[k], abs(d__1));
	dscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L140:
	z__[k] /= abd[*m + 1 + k * abd_dim1];
/* L150: */
    }
    s = 1. / dasum_(n, &z__[1], &c__1);
    dscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*        solve  r*z = w */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if ((d__1 = z__[k], abs(d__1)) <= abd[*m + 1 + k * abd_dim1]) {
	    goto L160;
	}
	s = abd[*m + 1 + k * abd_dim1] / (d__1 = z__[k], abs(d__1));
	dscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L160:
	z__[k] /= abd[*m + 1 + k * abd_dim1];
/* Computing MIN */
	i__2 = k - 1;
	lm = min(i__2,*m);
	la = *m + 1 - lm;
	lb = k - lm;
	t = -z__[k];
	daxpy_(&lm, &t, &abd[la + k * abd_dim1], &c__1, &z__[lb], &c__1);
/* L170: */
    }
/*        make znorm = 1.0 */
    s = 1. / dasum_(n, &z__[1], &c__1);
    dscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.) {
	*rcond = 0.;
    }
L180:
    return 0;
} /* dpbco_ */

