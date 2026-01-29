/* cpbco.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int cpbco_(complex *abd, integer *lda, integer *n, integer *
	m, real *rcond, complex *z__, integer *info)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3, i__4, i__5;
    real r__1, r__2, r__3, r__4, r__5, r__6, r__7, r__8;
    complex q__1, q__2, q__3, q__4;

    /* Builtin functions */
    double r_imag(complex *);
    void c_div(complex *, complex *, complex *), r_cnjg(complex *, complex *);

    /* Local variables */
    integer i__, j, k, l;
    real s;
    complex t;
    integer j2, kb, la, lb;
    complex ek;
    integer lm;
    real sm;
    complex wk;
    integer mu, kp1;
    complex wkm;
    extern /* Subroutine */ int cpbfa_(complex *, integer *, integer *, 
	    integer *, integer *);
    extern /* Complex */ VOID cdotc_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);
    real anorm;
    extern /* Subroutine */ int caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *);
    real ynorm;
    extern /* Subroutine */ int csscal_(integer *, real *, complex *, integer 
	    *);
    extern doublereal scasum_(integer *, complex *, integer *);


/*     cpbco factors a complex hermitian positive definite matrix */
/*     stored in band form and estimates the condition of the matrix. */

/*     if  rcond  is not needed, cpbfa is slightly faster. */
/*     to solve  a*x = b , follow cpbco by cpbsl. */
/*     to compute  inverse(a)*c , follow cpbco by cpbsl. */
/*     to compute  determinant(a) , follow cpbco by cpbdi. */

/*     on entry */

/*        abd     complex(lda, n) */
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

/*        z       complex(n) */
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

/*     linpack cpbfa */
/*     blas caxpy,cdotc,csscal,scasum */
/*     fortran abs,aimag,amax1,cmplx,conjg,max0,min0,real */

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
	r__1 = scasum_(&l, &abd[mu + j * abd_dim1], &c__1);
	q__1.r = r__1, q__1.i = 0.f;
	z__[i__2].r = q__1.r, z__[i__2].i = q__1.i;
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
	    r__3 = z__[i__4].r + ((r__1 = abd[i__5].r, dabs(r__1)) + (r__2 = 
		    r_imag(&abd[i__ + j * abd_dim1]), dabs(r__2)));
	    q__1.r = r__3, q__1.i = 0.f;
	    z__[i__3].r = q__1.r, z__[i__3].i = q__1.i;
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
	i__2 = j;
	r__1 = anorm, r__2 = z__[i__2].r;
	anorm = dmax(r__1,r__2);
/* L40: */
    }

/*     factor */

    cpbfa_(&abd[abd_offset], lda, n, m, info);
    if (*info != 0) {
	goto L180;
    }

/*        rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) . */
/*        estimate = norm(z)/norm(y) where  a*z = y  and  a*y = e . */
/*        the components of  e  are chosen to cause maximum local */
/*        growth in the elements of w  where  ctrans(r)*w = e . */
/*        the vectors are frequently rescaled to avoid overflow. */

/*        solve ctrans(r)*w = e */

    ek.r = 1.f, ek.i = 0.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	z__[i__2].r = 0.f, z__[i__2].i = 0.f;
/* L50: */
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k;
	if ((r__1 = z__[i__2].r, dabs(r__1)) + (r__2 = r_imag(&z__[k]), dabs(
		r__2)) != 0.f) {
	    i__3 = k;
	    q__2.r = -z__[i__3].r, q__2.i = -z__[i__3].i;
	    q__1.r = q__2.r, q__1.i = q__2.i;
	    r__7 = (r__3 = ek.r, dabs(r__3)) + (r__4 = r_imag(&ek), dabs(r__4)
		    );
	    r__8 = (r__5 = q__1.r, dabs(r__5)) + (r__6 = r_imag(&q__1), dabs(
		    r__6));
	    q__4.r = q__1.r / r__8, q__4.i = q__1.i / r__8;
	    q__3.r = r__7 * q__4.r, q__3.i = r__7 * q__4.i;
	    ek.r = q__3.r, ek.i = q__3.i;
	}
	i__2 = k;
	q__2.r = ek.r - z__[i__2].r, q__2.i = ek.i - z__[i__2].i;
	q__1.r = q__2.r, q__1.i = q__2.i;
	i__3 = *m + 1 + k * abd_dim1;
	if ((r__1 = q__1.r, dabs(r__1)) + (r__2 = r_imag(&q__1), dabs(r__2)) 
		<= abd[i__3].r) {
	    goto L60;
	}
	i__2 = k;
	q__2.r = ek.r - z__[i__2].r, q__2.i = ek.i - z__[i__2].i;
	q__1.r = q__2.r, q__1.i = q__2.i;
	i__3 = *m + 1 + k * abd_dim1;
	s = abd[i__3].r / ((r__1 = q__1.r, dabs(r__1)) + (r__2 = r_imag(&q__1)
		, dabs(r__2)));
	csscal_(n, &s, &z__[1], &c__1);
	q__2.r = s, q__2.i = 0.f;
	q__1.r = q__2.r * ek.r - q__2.i * ek.i, q__1.i = q__2.r * ek.i + 
		q__2.i * ek.r;
	ek.r = q__1.r, ek.i = q__1.i;
L60:
	i__2 = k;
	q__1.r = ek.r - z__[i__2].r, q__1.i = ek.i - z__[i__2].i;
	wk.r = q__1.r, wk.i = q__1.i;
	q__2.r = -ek.r, q__2.i = -ek.i;
	i__2 = k;
	q__1.r = q__2.r - z__[i__2].r, q__1.i = q__2.i - z__[i__2].i;
	wkm.r = q__1.r, wkm.i = q__1.i;
	s = (r__1 = wk.r, dabs(r__1)) + (r__2 = r_imag(&wk), dabs(r__2));
	sm = (r__1 = wkm.r, dabs(r__1)) + (r__2 = r_imag(&wkm), dabs(r__2));
	c_div(&q__1, &wk, &abd[*m + 1 + k * abd_dim1]);
	wk.r = q__1.r, wk.i = q__1.i;
	c_div(&q__1, &wkm, &abd[*m + 1 + k * abd_dim1]);
	wkm.r = q__1.r, wkm.i = q__1.i;
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
	    r_cnjg(&q__4, &abd[i__ + j * abd_dim1]);
	    q__3.r = wkm.r * q__4.r - wkm.i * q__4.i, q__3.i = wkm.r * q__4.i 
		    + wkm.i * q__4.r;
	    q__2.r = z__[i__3].r + q__3.r, q__2.i = z__[i__3].i + q__3.i;
	    q__1.r = q__2.r, q__1.i = q__2.i;
	    sm += (r__1 = q__1.r, dabs(r__1)) + (r__2 = r_imag(&q__1), dabs(
		    r__2));
	    i__3 = j;
	    i__4 = j;
	    r_cnjg(&q__3, &abd[i__ + j * abd_dim1]);
	    q__2.r = wk.r * q__3.r - wk.i * q__3.i, q__2.i = wk.r * q__3.i + 
		    wk.i * q__3.r;
	    q__1.r = z__[i__4].r + q__2.r, q__1.i = z__[i__4].i + q__2.i;
	    z__[i__3].r = q__1.r, z__[i__3].i = q__1.i;
	    i__3 = j;
	    s += (r__1 = z__[i__3].r, dabs(r__1)) + (r__2 = r_imag(&z__[j]), 
		    dabs(r__2));
/* L70: */
	}
	if (s >= sm) {
	    goto L90;
	}
	q__1.r = wkm.r - wk.r, q__1.i = wkm.i - wk.i;
	t.r = q__1.r, t.i = q__1.i;
	wk.r = wkm.r, wk.i = wkm.i;
	i__ = *m + 1;
	i__2 = j2;
	for (j = kp1; j <= i__2; ++j) {
	    --i__;
	    i__3 = j;
	    i__4 = j;
	    r_cnjg(&q__3, &abd[i__ + j * abd_dim1]);
	    q__2.r = t.r * q__3.r - t.i * q__3.i, q__2.i = t.r * q__3.i + t.i 
		    * q__3.r;
	    q__1.r = z__[i__4].r + q__2.r, q__1.i = z__[i__4].i + q__2.i;
	    z__[i__3].r = q__1.r, z__[i__3].i = q__1.i;
/* L80: */
	}
L90:
L100:
	i__2 = k;
	z__[i__2].r = wk.r, z__[i__2].i = wk.i;
/* L110: */
    }
    s = 1.f / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

/*        solve  r*y = w */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	i__2 = k;
	i__3 = *m + 1 + k * abd_dim1;
	if ((r__1 = z__[i__2].r, dabs(r__1)) + (r__2 = r_imag(&z__[k]), dabs(
		r__2)) <= abd[i__3].r) {
	    goto L120;
	}
	i__2 = *m + 1 + k * abd_dim1;
	i__3 = k;
	s = abd[i__2].r / ((r__1 = z__[i__3].r, dabs(r__1)) + (r__2 = r_imag(&
		z__[k]), dabs(r__2)));
	csscal_(n, &s, &z__[1], &c__1);
L120:
	i__2 = k;
	c_div(&q__1, &z__[k], &abd[*m + 1 + k * abd_dim1]);
	z__[i__2].r = q__1.r, z__[i__2].i = q__1.i;
/* Computing MIN */
	i__2 = k - 1;
	lm = min(i__2,*m);
	la = *m + 1 - lm;
	lb = k - lm;
	i__2 = k;
	q__1.r = -z__[i__2].r, q__1.i = -z__[i__2].i;
	t.r = q__1.r, t.i = q__1.i;
	caxpy_(&lm, &t, &abd[la + k * abd_dim1], &c__1, &z__[lb], &c__1);
/* L130: */
    }
    s = 1.f / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.f;

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
	cdotc_(&q__2, &lm, &abd[la + k * abd_dim1], &c__1, &z__[lb], &c__1);
	q__1.r = z__[i__3].r - q__2.r, q__1.i = z__[i__3].i - q__2.i;
	z__[i__2].r = q__1.r, z__[i__2].i = q__1.i;
	i__2 = k;
	i__3 = *m + 1 + k * abd_dim1;
	if ((r__1 = z__[i__2].r, dabs(r__1)) + (r__2 = r_imag(&z__[k]), dabs(
		r__2)) <= abd[i__3].r) {
	    goto L140;
	}
	i__2 = *m + 1 + k * abd_dim1;
	i__3 = k;
	s = abd[i__2].r / ((r__1 = z__[i__3].r, dabs(r__1)) + (r__2 = r_imag(&
		z__[k]), dabs(r__2)));
	csscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L140:
	i__2 = k;
	c_div(&q__1, &z__[k], &abd[*m + 1 + k * abd_dim1]);
	z__[i__2].r = q__1.r, z__[i__2].i = q__1.i;
/* L150: */
    }
    s = 1.f / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*        solve  r*z = w */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	i__2 = k;
	i__3 = *m + 1 + k * abd_dim1;
	if ((r__1 = z__[i__2].r, dabs(r__1)) + (r__2 = r_imag(&z__[k]), dabs(
		r__2)) <= abd[i__3].r) {
	    goto L160;
	}
	i__2 = *m + 1 + k * abd_dim1;
	i__3 = k;
	s = abd[i__2].r / ((r__1 = z__[i__3].r, dabs(r__1)) + (r__2 = r_imag(&
		z__[k]), dabs(r__2)));
	csscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L160:
	i__2 = k;
	c_div(&q__1, &z__[k], &abd[*m + 1 + k * abd_dim1]);
	z__[i__2].r = q__1.r, z__[i__2].i = q__1.i;
/* Computing MIN */
	i__2 = k - 1;
	lm = min(i__2,*m);
	la = *m + 1 - lm;
	lb = k - lm;
	i__2 = k;
	q__1.r = -z__[i__2].r, q__1.i = -z__[i__2].i;
	t.r = q__1.r, t.i = q__1.i;
	caxpy_(&lm, &t, &abd[la + k * abd_dim1], &c__1, &z__[lb], &c__1);
/* L170: */
    }
/*        make znorm = 1.0 */
    s = 1.f / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.f) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.f) {
	*rcond = 0.f;
    }
L180:
    return 0;
} /* cpbco_ */

