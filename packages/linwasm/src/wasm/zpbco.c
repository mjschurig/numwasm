/* zpbco.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int zpbco_(doublecomplex *abd, integer *lda, integer *n, 
	integer *m, doublereal *rcond, doublecomplex *z__, integer *info)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    integer i__, j, k, l;
    doublereal s;
    doublecomplex t;
    integer j2, kb, la, lb;
    doublecomplex ek;
    integer lm;
    doublereal sm;
    doublecomplex wk;
    integer mu, kp1;
    doublecomplex wkm;
    extern /* Subroutine */ int zpbfa_(doublecomplex *, integer *, integer *, 
	    integer *, integer *);
    doublereal anorm;
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    doublereal ynorm;
    extern /* Subroutine */ int zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), zdscal_(
	    integer *, doublereal *, doublecomplex *, integer *);
    extern doublereal dzasum_(integer *, doublecomplex *, integer *);


/*     zpbco factors a complex*16 hermitian positive definite matrix */
/*     stored in band form and estimates the condition of the matrix. */

/*     if  rcond  is not needed, zpbfa is slightly faster. */
/*     to solve  a*x = b , follow zpbco by zpbsl. */
/*     to compute  inverse(a)*c , follow zpbco by zpbsl. */
/*     to compute  determinant(a) , follow zpbco by zpbdi. */

/*     on entry */

/*        abd     complex*16(lda, n) */
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
/*                form, so that  a = ctrans(r)*r . */
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

/*        z       complex*16(n) */
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

/*           if  a  is a hermitian positive definite band matrix, */
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

/*     linpack zpbfa */
/*     blas zaxpy,zdotc,zdscal,dzasum */
/*     fortran dabs,dmax1,dcmplx,dconjg,max0,min0 */

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
	i__2 = j;
	d__1 = dzasum_(&l, &abd[mu + j * abd_dim1], &c__1);
	z__1.r = d__1, z__1.i = 0.;
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
	k = j - l;
	if (*m < mu) {
	    goto L20;
	}
	i__2 = *m;
	for (i__ = mu; i__ <= i__2; ++i__) {
	    ++k;
	    i__3 = k;
	    i__4 = k;
	    i__5 = i__ + j * abd_dim1;
	    i__6 = i__ + j * abd_dim1;
	    z__2.r = abd[i__6].r * 0. - abd[i__6].i * -1., z__2.i = abd[i__6]
		    .i * 0. + abd[i__6].r * -1.;
	    d__3 = z__[i__4].r + ((d__1 = abd[i__5].r, abs(d__1)) + (d__2 = 
		    z__2.r, abs(d__2)));
	    z__1.r = d__3, z__1.i = 0.;
	    z__[i__3].r = z__1.r, z__[i__3].i = z__1.i;
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
	i__2 = j;
	d__1 = anorm, d__2 = z__[i__2].r;
	anorm = max(d__1,d__2);
/* L40: */
    }

/*     factor */

    zpbfa_(&abd[abd_offset], lda, n, m, info);
    if (*info != 0) {
	goto L180;
    }

/*        rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) . */
/*        estimate = norm(z)/norm(y) where  a*z = y  and  a*y = e . */
/*        the components of  e  are chosen to cause maximum local */
/*        growth in the elements of w  where  ctrans(r)*w = e . */
/*        the vectors are frequently rescaled to avoid overflow. */

/*        solve ctrans(r)*w = e */

    ek.r = 1., ek.i = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	z__[i__2].r = 0., z__[i__2].i = 0.;
/* L50: */
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k;
	i__3 = k;
	z__1.r = z__[i__3].r * 0. - z__[i__3].i * -1., z__1.i = z__[i__3].i * 
		0. + z__[i__3].r * -1.;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) != 
		0.) {
	    i__4 = k;
	    z__3.r = -z__[i__4].r, z__3.i = -z__[i__4].i;
	    z__2.r = z__3.r, z__2.i = z__3.i;
	    z__5.r = ek.r * 0. - ek.i * -1., z__5.i = ek.i * 0. + ek.r * -1.;
	    d__7 = (d__3 = ek.r, abs(d__3)) + (d__4 = z__5.r, abs(d__4));
	    z__7.r = z__2.r * 0. - z__2.i * -1., z__7.i = z__2.i * 0. + 
		    z__2.r * -1.;
	    d__8 = (d__5 = z__2.r, abs(d__5)) + (d__6 = z__7.r, abs(d__6));
	    z__6.r = z__2.r / d__8, z__6.i = z__2.i / d__8;
	    z__4.r = d__7 * z__6.r, z__4.i = d__7 * z__6.i;
	    ek.r = z__4.r, ek.i = z__4.i;
	}
	i__2 = k;
	z__2.r = ek.r - z__[i__2].r, z__2.i = ek.i - z__[i__2].i;
	z__1.r = z__2.r, z__1.i = z__2.i;
	z__3.r = z__1.r * 0. - z__1.i * -1., z__3.i = z__1.i * 0. + z__1.r * 
		-1.;
	i__3 = *m + 1 + k * abd_dim1;
	if ((d__1 = z__1.r, abs(d__1)) + (d__2 = z__3.r, abs(d__2)) <= abd[
		i__3].r) {
	    goto L60;
	}
	i__2 = k;
	z__2.r = ek.r - z__[i__2].r, z__2.i = ek.i - z__[i__2].i;
	z__1.r = z__2.r, z__1.i = z__2.i;
	i__3 = *m + 1 + k * abd_dim1;
	z__3.r = z__1.r * 0. - z__1.i * -1., z__3.i = z__1.i * 0. + z__1.r * 
		-1.;
	s = abd[i__3].r / ((d__1 = z__1.r, abs(d__1)) + (d__2 = z__3.r, abs(
		d__2)));
	zdscal_(n, &s, &z__[1], &c__1);
	z__2.r = s, z__2.i = 0.;
	z__1.r = z__2.r * ek.r - z__2.i * ek.i, z__1.i = z__2.r * ek.i + 
		z__2.i * ek.r;
	ek.r = z__1.r, ek.i = z__1.i;
L60:
	i__2 = k;
	z__1.r = ek.r - z__[i__2].r, z__1.i = ek.i - z__[i__2].i;
	wk.r = z__1.r, wk.i = z__1.i;
	z__2.r = -ek.r, z__2.i = -ek.i;
	i__2 = k;
	z__1.r = z__2.r - z__[i__2].r, z__1.i = z__2.i - z__[i__2].i;
	wkm.r = z__1.r, wkm.i = z__1.i;
	z__1.r = wk.r * 0. - wk.i * -1., z__1.i = wk.i * 0. + wk.r * -1.;
	s = (d__1 = wk.r, abs(d__1)) + (d__2 = z__1.r, abs(d__2));
	z__1.r = wkm.r * 0. - wkm.i * -1., z__1.i = wkm.i * 0. + wkm.r * -1.;
	sm = (d__1 = wkm.r, abs(d__1)) + (d__2 = z__1.r, abs(d__2));
	z_div(&z__1, &wk, &abd[*m + 1 + k * abd_dim1]);
	wk.r = z__1.r, wk.i = z__1.i;
	z_div(&z__1, &wkm, &abd[*m + 1 + k * abd_dim1]);
	wkm.r = z__1.r, wkm.i = z__1.i;
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
	    i__3 = j;
	    d_cnjg(&z__4, &abd[i__ + j * abd_dim1]);
	    z__3.r = wkm.r * z__4.r - wkm.i * z__4.i, z__3.i = wkm.r * z__4.i 
		    + wkm.i * z__4.r;
	    z__2.r = z__[i__3].r + z__3.r, z__2.i = z__[i__3].i + z__3.i;
	    z__1.r = z__2.r, z__1.i = z__2.i;
	    z__5.r = z__1.r * 0. - z__1.i * -1., z__5.i = z__1.i * 0. + 
		    z__1.r * -1.;
	    sm += (d__1 = z__1.r, abs(d__1)) + (d__2 = z__5.r, abs(d__2));
	    i__3 = j;
	    i__4 = j;
	    d_cnjg(&z__3, &abd[i__ + j * abd_dim1]);
	    z__2.r = wk.r * z__3.r - wk.i * z__3.i, z__2.i = wk.r * z__3.i + 
		    wk.i * z__3.r;
	    z__1.r = z__[i__4].r + z__2.r, z__1.i = z__[i__4].i + z__2.i;
	    z__[i__3].r = z__1.r, z__[i__3].i = z__1.i;
	    i__3 = j;
	    i__4 = j;
	    z__1.r = z__[i__4].r * 0. - z__[i__4].i * -1., z__1.i = z__[i__4]
		    .i * 0. + z__[i__4].r * -1.;
	    s += (d__1 = z__[i__3].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2));
/* L70: */
	}
	if (s >= sm) {
	    goto L90;
	}
	z__1.r = wkm.r - wk.r, z__1.i = wkm.i - wk.i;
	t.r = z__1.r, t.i = z__1.i;
	wk.r = wkm.r, wk.i = wkm.i;
	i__ = *m + 1;
	i__2 = j2;
	for (j = kp1; j <= i__2; ++j) {
	    --i__;
	    i__3 = j;
	    i__4 = j;
	    d_cnjg(&z__3, &abd[i__ + j * abd_dim1]);
	    z__2.r = t.r * z__3.r - t.i * z__3.i, z__2.i = t.r * z__3.i + t.i 
		    * z__3.r;
	    z__1.r = z__[i__4].r + z__2.r, z__1.i = z__[i__4].i + z__2.i;
	    z__[i__3].r = z__1.r, z__[i__3].i = z__1.i;
/* L80: */
	}
L90:
L100:
	i__2 = k;
	z__[i__2].r = wk.r, z__[i__2].i = wk.i;
/* L110: */
    }
    s = 1. / dzasum_(n, &z__[1], &c__1);
    zdscal_(n, &s, &z__[1], &c__1);

/*        solve  r*y = w */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	i__2 = k;
	i__3 = k;
	z__1.r = z__[i__3].r * 0. - z__[i__3].i * -1., z__1.i = z__[i__3].i * 
		0. + z__[i__3].r * -1.;
	i__4 = *m + 1 + k * abd_dim1;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) <= 
		abd[i__4].r) {
	    goto L120;
	}
	i__2 = *m + 1 + k * abd_dim1;
	i__3 = k;
	i__4 = k;
	z__1.r = z__[i__4].r * 0. - z__[i__4].i * -1., z__1.i = z__[i__4].i * 
		0. + z__[i__4].r * -1.;
	s = abd[i__2].r / ((d__1 = z__[i__3].r, abs(d__1)) + (d__2 = z__1.r, 
		abs(d__2)));
	zdscal_(n, &s, &z__[1], &c__1);
L120:
	i__2 = k;
	z_div(&z__1, &z__[k], &abd[*m + 1 + k * abd_dim1]);
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
/* Computing MIN */
	i__2 = k - 1;
	lm = min(i__2,*m);
	la = *m + 1 - lm;
	lb = k - lm;
	i__2 = k;
	z__1.r = -z__[i__2].r, z__1.i = -z__[i__2].i;
	t.r = z__1.r, t.i = z__1.i;
	zaxpy_(&lm, &t, &abd[la + k * abd_dim1], &c__1, &z__[lb], &c__1);
/* L130: */
    }
    s = 1. / dzasum_(n, &z__[1], &c__1);
    zdscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.;

/*        solve ctrans(r)*v = y */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* Computing MIN */
	i__2 = k - 1;
	lm = min(i__2,*m);
	la = *m + 1 - lm;
	lb = k - lm;
	i__2 = k;
	i__3 = k;
	zdotc_(&z__2, &lm, &abd[la + k * abd_dim1], &c__1, &z__[lb], &c__1);
	z__1.r = z__[i__3].r - z__2.r, z__1.i = z__[i__3].i - z__2.i;
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
	i__2 = k;
	i__3 = k;
	z__1.r = z__[i__3].r * 0. - z__[i__3].i * -1., z__1.i = z__[i__3].i * 
		0. + z__[i__3].r * -1.;
	i__4 = *m + 1 + k * abd_dim1;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) <= 
		abd[i__4].r) {
	    goto L140;
	}
	i__2 = *m + 1 + k * abd_dim1;
	i__3 = k;
	i__4 = k;
	z__1.r = z__[i__4].r * 0. - z__[i__4].i * -1., z__1.i = z__[i__4].i * 
		0. + z__[i__4].r * -1.;
	s = abd[i__2].r / ((d__1 = z__[i__3].r, abs(d__1)) + (d__2 = z__1.r, 
		abs(d__2)));
	zdscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L140:
	i__2 = k;
	z_div(&z__1, &z__[k], &abd[*m + 1 + k * abd_dim1]);
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
/* L150: */
    }
    s = 1. / dzasum_(n, &z__[1], &c__1);
    zdscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*        solve  r*z = w */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	i__2 = k;
	i__3 = k;
	z__1.r = z__[i__3].r * 0. - z__[i__3].i * -1., z__1.i = z__[i__3].i * 
		0. + z__[i__3].r * -1.;
	i__4 = *m + 1 + k * abd_dim1;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) <= 
		abd[i__4].r) {
	    goto L160;
	}
	i__2 = *m + 1 + k * abd_dim1;
	i__3 = k;
	i__4 = k;
	z__1.r = z__[i__4].r * 0. - z__[i__4].i * -1., z__1.i = z__[i__4].i * 
		0. + z__[i__4].r * -1.;
	s = abd[i__2].r / ((d__1 = z__[i__3].r, abs(d__1)) + (d__2 = z__1.r, 
		abs(d__2)));
	zdscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L160:
	i__2 = k;
	z_div(&z__1, &z__[k], &abd[*m + 1 + k * abd_dim1]);
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
/* Computing MIN */
	i__2 = k - 1;
	lm = min(i__2,*m);
	la = *m + 1 - lm;
	lb = k - lm;
	i__2 = k;
	z__1.r = -z__[i__2].r, z__1.i = -z__[i__2].i;
	t.r = z__1.r, t.i = z__1.i;
	zaxpy_(&lm, &t, &abd[la + k * abd_dim1], &c__1, &z__[lb], &c__1);
/* L170: */
    }
/*        make znorm = 1.0 */
    s = 1. / dzasum_(n, &z__[1], &c__1);
    zdscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.) {
	*rcond = 0.;
    }
L180:
    return 0;
} /* zpbco_ */

