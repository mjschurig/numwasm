/* spbco.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int spbco_(real *abd, integer *lda, integer *n, integer *m, 
	real *rcond, real *z__, integer *info)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3;
    real r__1, r__2;

    /* Builtin functions */
    double r_sign(real *, real *);

    /* Local variables */
    integer i__, j, k, l;
    real s, t;
    integer j2, kb, la, lb;
    real ek;
    integer lm;
    real sm, wk;
    integer mu, kp1;
    real wkm;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    extern /* Subroutine */ int spbfa_(real *, integer *, integer *, integer *
	    , integer *), sscal_(integer *, real *, real *, integer *);
    real anorm;
    extern doublereal sasum_(integer *, real *, integer *);
    real ynorm;
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *);


/*     spbco factors a real symmetric positive definite matrix */
/*     stored in band form and estimates the condition of the matrix. */

/*     if  rcond  is not needed, spbfa is slightly faster. */
/*     to solve  a*x = b , follow spbco by spbsl. */
/*     to compute  inverse(a)*c , follow spbco by spbsl. */
/*     to compute  determinant(a) , follow spbco by spbdi. */

/*     on entry */

/*        abd     real(lda, n) */
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

/*        rcond   real */
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

/*        z       real(n) */
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

/*     linpack spbfa */
/*     blas saxpy,sdot,sscal,sasum */
/*     fortran abs,amax1,max0,min0,real,sign */

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
	z__[j] = sasum_(&l, &abd[mu + j * abd_dim1], &c__1);
	k = j - l;
	if (*m < mu) {
	    goto L20;
	}
	i__2 = *m;
	for (i__ = mu; i__ <= i__2; ++i__) {
	    ++k;
	    z__[k] += (r__1 = abd[i__ + j * abd_dim1], dabs(r__1));
/* L10: */
	}
L20:
/* L30: */
	;
    }
    anorm = 0.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	r__1 = anorm, r__2 = z__[j];
	anorm = dmax(r__1,r__2);
/* L40: */
    }

/*     factor */

    spbfa_(&abd[abd_offset], lda, n, m, info);
    if (*info != 0) {
	goto L180;
    }

/*        rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) . */
/*        estimate = norm(z)/norm(y) where  a*z = y  and  a*y = e . */
/*        the components of  e  are chosen to cause maximum local */
/*        growth in the elements of w  where  trans(r)*w = e . */
/*        the vectors are frequently rescaled to avoid overflow. */

/*        solve trans(r)*w = e */

    ek = 1.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z__[j] = 0.f;
/* L50: */
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (z__[k] != 0.f) {
	    r__1 = -z__[k];
	    ek = r_sign(&ek, &r__1);
	}
	if ((r__1 = ek - z__[k], dabs(r__1)) <= abd[*m + 1 + k * abd_dim1]) {
	    goto L60;
	}
	s = abd[*m + 1 + k * abd_dim1] / (r__1 = ek - z__[k], dabs(r__1));
	sscal_(n, &s, &z__[1], &c__1);
	ek = s * ek;
L60:
	wk = ek - z__[k];
	wkm = -ek - z__[k];
	s = dabs(wk);
	sm = dabs(wkm);
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
	    sm += (r__1 = z__[j] + wkm * abd[i__ + j * abd_dim1], dabs(r__1));
	    z__[j] += wk * abd[i__ + j * abd_dim1];
	    s += (r__1 = z__[j], dabs(r__1));
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
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

/*        solve  r*y = w */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if ((r__1 = z__[k], dabs(r__1)) <= abd[*m + 1 + k * abd_dim1]) {
	    goto L120;
	}
	s = abd[*m + 1 + k * abd_dim1] / (r__1 = z__[k], dabs(r__1));
	sscal_(n, &s, &z__[1], &c__1);
L120:
	z__[k] /= abd[*m + 1 + k * abd_dim1];
/* Computing MIN */
	i__2 = k - 1;
	lm = min(i__2,*m);
	la = *m + 1 - lm;
	lb = k - lm;
	t = -z__[k];
	saxpy_(&lm, &t, &abd[la + k * abd_dim1], &c__1, &z__[lb], &c__1);
/* L130: */
    }
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.f;

/*        solve trans(r)*v = y */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* Computing MIN */
	i__2 = k - 1;
	lm = min(i__2,*m);
	la = *m + 1 - lm;
	lb = k - lm;
	z__[k] -= sdot_(&lm, &abd[la + k * abd_dim1], &c__1, &z__[lb], &c__1);
	if ((r__1 = z__[k], dabs(r__1)) <= abd[*m + 1 + k * abd_dim1]) {
	    goto L140;
	}
	s = abd[*m + 1 + k * abd_dim1] / (r__1 = z__[k], dabs(r__1));
	sscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L140:
	z__[k] /= abd[*m + 1 + k * abd_dim1];
/* L150: */
    }
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*        solve  r*z = w */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if ((r__1 = z__[k], dabs(r__1)) <= abd[*m + 1 + k * abd_dim1]) {
	    goto L160;
	}
	s = abd[*m + 1 + k * abd_dim1] / (r__1 = z__[k], dabs(r__1));
	sscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L160:
	z__[k] /= abd[*m + 1 + k * abd_dim1];
/* Computing MIN */
	i__2 = k - 1;
	lm = min(i__2,*m);
	la = *m + 1 - lm;
	lb = k - lm;
	t = -z__[k];
	saxpy_(&lm, &t, &abd[la + k * abd_dim1], &c__1, &z__[lb], &c__1);
/* L170: */
    }
/*        make znorm = 1.0 */
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.f) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.f) {
	*rcond = 0.f;
    }
L180:
    return 0;
} /* spbco_ */

