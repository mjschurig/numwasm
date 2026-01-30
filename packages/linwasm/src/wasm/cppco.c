/* cppco.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int cppco_(complex *ap, integer *n, real *rcond, complex *
	z__, integer *info)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    real r__1, r__2, r__3, r__4, r__5, r__6, r__7, r__8;
    complex q__1, q__2, q__3, q__4;

    /* Builtin functions */
    double r_imag(complex *);
    void c_div(complex *, complex *, complex *), r_cnjg(complex *, complex *);

    /* Local variables */
    integer i__, j, k;
    real s;
    complex t;
    integer j1, kb;
    complex ek;
    integer ij, kj, kk;
    real sm;
    complex wk;
    integer jm1, kp1;
    complex wkm;
    extern /* Subroutine */ int cppfa_(complex *, integer *, integer *);
    extern /* Complex */ VOID cdotc_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);
    real anorm;
    extern /* Subroutine */ int caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *);
    real ynorm;
    extern /* Subroutine */ int csscal_(integer *, real *, complex *, integer 
	    *);
    extern doublereal scasum_(integer *, complex *, integer *);


/*     cppco factors a complex hermitian positive definite matrix */
/*     stored in packed form */
/*     and estimates the condition of the matrix. */

/*     if  rcond  is not needed, cppfa is slightly faster. */
/*     to solve  a*x = b , follow cppco by cppsl. */
/*     to compute  inverse(a)*c , follow cppco by cppsl. */
/*     to compute  determinant(a) , follow cppco by cppdi. */
/*     to compute  inverse(a) , follow cppco by cppdi. */

/*     on entry */

/*        ap      complex (n*(n+1)/2) */
/*                the packed form of a hermitian matrix  a .  the */
/*                columns of the upper triangle are stored sequentially */
/*                in a one-dimensional array of length  n*(n+1)/2 . */
/*                see comments below for details. */

/*        n       integer */
/*                the order of the matrix  a . */

/*     on return */

/*        ap      an upper triangular matrix  r , stored in packed */
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

/*     packed storage */

/*          the following program segment will pack the upper */
/*          triangle of a hermitian matrix. */

/*                k = 0 */
/*                do 20 j = 1, n */
/*                   do 10 i = 1, j */
/*                      k = k + 1 */
/*                      ap(k) = a(i,j) */
/*             10    continue */
/*             20 continue */

/*     linpack.  this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     linpack cppfa */
/*     blas caxpy,cdotc,csscal,scasum */
/*     fortran abs,aimag,amax1,cmplx,conjg,real */

/*     internal variables */



/*     find norm of a */

    /* Parameter adjustments */
    --z__;
    --ap;

    /* Function Body */
    j1 = 1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	r__1 = scasum_(&j, &ap[j1], &c__1);
	q__1.r = r__1, q__1.i = 0.f;
	z__[i__2].r = q__1.r, z__[i__2].i = q__1.i;
	ij = j1;
	j1 += j;
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    i__4 = i__;
	    i__5 = ij;
	    r__3 = z__[i__4].r + ((r__1 = ap[i__5].r, dabs(r__1)) + (r__2 = 
		    r_imag(&ap[ij]), dabs(r__2)));
	    q__1.r = r__3, q__1.i = 0.f;
	    z__[i__3].r = q__1.r, z__[i__3].i = q__1.i;
	    ++ij;
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

    cppfa_(&ap[1], n, info);
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
    kk = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	kk += k;
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
	i__3 = kk;
	if ((r__1 = q__1.r, dabs(r__1)) + (r__2 = r_imag(&q__1), dabs(r__2)) 
		<= ap[i__3].r) {
	    goto L60;
	}
	i__2 = k;
	q__2.r = ek.r - z__[i__2].r, q__2.i = ek.i - z__[i__2].i;
	q__1.r = q__2.r, q__1.i = q__2.i;
	i__3 = kk;
	s = ap[i__3].r / ((r__1 = q__1.r, dabs(r__1)) + (r__2 = r_imag(&q__1),
		 dabs(r__2)));
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
	c_div(&q__1, &wk, &ap[kk]);
	wk.r = q__1.r, wk.i = q__1.i;
	c_div(&q__1, &wkm, &ap[kk]);
	wkm.r = q__1.r, wkm.i = q__1.i;
	kp1 = k + 1;
	kj = kk + k;
	if (kp1 > *n) {
	    goto L100;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    i__3 = j;
	    r_cnjg(&q__4, &ap[kj]);
	    q__3.r = wkm.r * q__4.r - wkm.i * q__4.i, q__3.i = wkm.r * q__4.i 
		    + wkm.i * q__4.r;
	    q__2.r = z__[i__3].r + q__3.r, q__2.i = z__[i__3].i + q__3.i;
	    q__1.r = q__2.r, q__1.i = q__2.i;
	    sm += (r__1 = q__1.r, dabs(r__1)) + (r__2 = r_imag(&q__1), dabs(
		    r__2));
	    i__3 = j;
	    i__4 = j;
	    r_cnjg(&q__3, &ap[kj]);
	    q__2.r = wk.r * q__3.r - wk.i * q__3.i, q__2.i = wk.r * q__3.i + 
		    wk.i * q__3.r;
	    q__1.r = z__[i__4].r + q__2.r, q__1.i = z__[i__4].i + q__2.i;
	    z__[i__3].r = q__1.r, z__[i__3].i = q__1.i;
	    i__3 = j;
	    s += (r__1 = z__[i__3].r, dabs(r__1)) + (r__2 = r_imag(&z__[j]), 
		    dabs(r__2));
	    kj += j;
/* L70: */
	}
	if (s >= sm) {
	    goto L90;
	}
	q__1.r = wkm.r - wk.r, q__1.i = wkm.i - wk.i;
	t.r = q__1.r, t.i = q__1.i;
	wk.r = wkm.r, wk.i = wkm.i;
	kj = kk + k;
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    i__3 = j;
	    i__4 = j;
	    r_cnjg(&q__3, &ap[kj]);
	    q__2.r = t.r * q__3.r - t.i * q__3.i, q__2.i = t.r * q__3.i + t.i 
		    * q__3.r;
	    q__1.r = z__[i__4].r + q__2.r, q__1.i = z__[i__4].i + q__2.i;
	    z__[i__3].r = q__1.r, z__[i__3].i = q__1.i;
	    kj += j;
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

/*        solve r*y = w */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	i__2 = k;
	i__3 = kk;
	if ((r__1 = z__[i__2].r, dabs(r__1)) + (r__2 = r_imag(&z__[k]), dabs(
		r__2)) <= ap[i__3].r) {
	    goto L120;
	}
	i__2 = kk;
	i__3 = k;
	s = ap[i__2].r / ((r__1 = z__[i__3].r, dabs(r__1)) + (r__2 = r_imag(&
		z__[k]), dabs(r__2)));
	csscal_(n, &s, &z__[1], &c__1);
L120:
	i__2 = k;
	c_div(&q__1, &z__[k], &ap[kk]);
	z__[i__2].r = q__1.r, z__[i__2].i = q__1.i;
	kk -= k;
	i__2 = k;
	q__1.r = -z__[i__2].r, q__1.i = -z__[i__2].i;
	t.r = q__1.r, t.i = q__1.i;
	i__2 = k - 1;
	caxpy_(&i__2, &t, &ap[kk + 1], &c__1, &z__[1], &c__1);
/* L130: */
    }
    s = 1.f / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.f;

/*        solve ctrans(r)*v = y */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k;
	i__3 = k;
	i__4 = k - 1;
	cdotc_(&q__2, &i__4, &ap[kk + 1], &c__1, &z__[1], &c__1);
	q__1.r = z__[i__3].r - q__2.r, q__1.i = z__[i__3].i - q__2.i;
	z__[i__2].r = q__1.r, z__[i__2].i = q__1.i;
	kk += k;
	i__2 = k;
	i__3 = kk;
	if ((r__1 = z__[i__2].r, dabs(r__1)) + (r__2 = r_imag(&z__[k]), dabs(
		r__2)) <= ap[i__3].r) {
	    goto L140;
	}
	i__2 = kk;
	i__3 = k;
	s = ap[i__2].r / ((r__1 = z__[i__3].r, dabs(r__1)) + (r__2 = r_imag(&
		z__[k]), dabs(r__2)));
	csscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L140:
	i__2 = k;
	c_div(&q__1, &z__[k], &ap[kk]);
	z__[i__2].r = q__1.r, z__[i__2].i = q__1.i;
/* L150: */
    }
    s = 1.f / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*        solve r*z = v */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	i__2 = k;
	i__3 = kk;
	if ((r__1 = z__[i__2].r, dabs(r__1)) + (r__2 = r_imag(&z__[k]), dabs(
		r__2)) <= ap[i__3].r) {
	    goto L160;
	}
	i__2 = kk;
	i__3 = k;
	s = ap[i__2].r / ((r__1 = z__[i__3].r, dabs(r__1)) + (r__2 = r_imag(&
		z__[k]), dabs(r__2)));
	csscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L160:
	i__2 = k;
	c_div(&q__1, &z__[k], &ap[kk]);
	z__[i__2].r = q__1.r, z__[i__2].i = q__1.i;
	kk -= k;
	i__2 = k;
	q__1.r = -z__[i__2].r, q__1.i = -z__[i__2].i;
	t.r = q__1.r, t.i = q__1.i;
	i__2 = k - 1;
	caxpy_(&i__2, &t, &ap[kk + 1], &c__1, &z__[1], &c__1);
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
} /* cppco_ */

