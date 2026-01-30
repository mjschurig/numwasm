/* zpoco.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int zpoco_(doublecomplex *a, integer *lda, integer *n, 
	doublereal *rcond, doublecomplex *z__, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4, i__5, i__6;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    integer i__, j, k;
    doublereal s;
    doublecomplex t;
    integer kb;
    doublecomplex ek;
    doublereal sm;
    doublecomplex wk;
    integer jm1, kp1;
    doublecomplex wkm;
    doublereal anorm;
    extern /* Subroutine */ int zpofa_(doublecomplex *, integer *, integer *, 
	    integer *);
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    doublereal ynorm;
    extern /* Subroutine */ int zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), zdscal_(
	    integer *, doublereal *, doublecomplex *, integer *);
    extern doublereal dzasum_(integer *, doublecomplex *, integer *);


/*     zpoco factors a complex*16 hermitian positive definite matrix */
/*     and estimates the condition of the matrix. */

/*     if  rcond  is not needed, zpofa is slightly faster. */
/*     to solve  a*x = b , follow zpoco by zposl. */
/*     to compute  inverse(a)*c , follow zpoco by zposl. */
/*     to compute  determinant(a) , follow zpoco by zpodi. */
/*     to compute  inverse(a) , follow zpoco by zpodi. */

/*     on entry */

/*        a       complex*16(lda, n) */
/*                the hermitian matrix to be factored.  only the */
/*                diagonal and upper triangle are used. */

/*        lda     integer */
/*                the leading dimension of the array  a . */

/*        n       integer */
/*                the order of the matrix  a . */

/*     on return */

/*        a       an upper triangular matrix  r  so that  a = */
/*                ctrans(r)*r where  ctrans(r)  is the conjugate */
/*                transpose.  the strict lower triangle is unaltered. */
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
/*                if  a  is close to a singular matrix, then  z  is */
/*                an approximate null vector in the sense that */
/*                norm(a*z) = rcond*norm(a)*norm(z) . */
/*                if  info .ne. 0 , z  is unchanged. */

/*        info    integer */
/*                = 0  for normal return. */
/*                = k  signals an error condition.  the leading minor */
/*                     of order  k  is not positive definite. */

/*     linpack.  this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     linpack zpofa */
/*     blas zaxpy,zdotc,zdscal,dzasum */
/*     fortran dabs,dmax1,dcmplx,dconjg */

/*     internal variables */



/*     find norm of a using only upper half */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --z__;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	d__1 = dzasum_(&j, &a[j * a_dim1 + 1], &c__1);
	z__1.r = d__1, z__1.i = 0.;
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    i__4 = i__;
	    i__5 = i__ + j * a_dim1;
	    i__6 = i__ + j * a_dim1;
	    z__2.r = a[i__6].r * 0. - a[i__6].i * -1., z__2.i = a[i__6].i * 
		    0. + a[i__6].r * -1.;
	    d__3 = z__[i__4].r + ((d__1 = a[i__5].r, abs(d__1)) + (d__2 = 
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

    zpofa_(&a[a_offset], lda, n, info);
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
	i__3 = k + k * a_dim1;
	if ((d__1 = z__1.r, abs(d__1)) + (d__2 = z__3.r, abs(d__2)) <= a[i__3]
		.r) {
	    goto L60;
	}
	i__2 = k;
	z__2.r = ek.r - z__[i__2].r, z__2.i = ek.i - z__[i__2].i;
	z__1.r = z__2.r, z__1.i = z__2.i;
	i__3 = k + k * a_dim1;
	z__3.r = z__1.r * 0. - z__1.i * -1., z__3.i = z__1.i * 0. + z__1.r * 
		-1.;
	s = a[i__3].r / ((d__1 = z__1.r, abs(d__1)) + (d__2 = z__3.r, abs(
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
	z_div(&z__1, &wk, &a[k + k * a_dim1]);
	wk.r = z__1.r, wk.i = z__1.i;
	z_div(&z__1, &wkm, &a[k + k * a_dim1]);
	wkm.r = z__1.r, wkm.i = z__1.i;
	kp1 = k + 1;
	if (kp1 > *n) {
	    goto L100;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    i__3 = j;
	    d_cnjg(&z__4, &a[k + j * a_dim1]);
	    z__3.r = wkm.r * z__4.r - wkm.i * z__4.i, z__3.i = wkm.r * z__4.i 
		    + wkm.i * z__4.r;
	    z__2.r = z__[i__3].r + z__3.r, z__2.i = z__[i__3].i + z__3.i;
	    z__1.r = z__2.r, z__1.i = z__2.i;
	    z__5.r = z__1.r * 0. - z__1.i * -1., z__5.i = z__1.i * 0. + 
		    z__1.r * -1.;
	    sm += (d__1 = z__1.r, abs(d__1)) + (d__2 = z__5.r, abs(d__2));
	    i__3 = j;
	    i__4 = j;
	    d_cnjg(&z__3, &a[k + j * a_dim1]);
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
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    i__3 = j;
	    i__4 = j;
	    d_cnjg(&z__3, &a[k + j * a_dim1]);
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

/*        solve r*y = w */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	i__2 = k;
	i__3 = k;
	z__1.r = z__[i__3].r * 0. - z__[i__3].i * -1., z__1.i = z__[i__3].i * 
		0. + z__[i__3].r * -1.;
	i__4 = k + k * a_dim1;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) <= a[
		i__4].r) {
	    goto L120;
	}
	i__2 = k + k * a_dim1;
	i__3 = k;
	i__4 = k;
	z__1.r = z__[i__4].r * 0. - z__[i__4].i * -1., z__1.i = z__[i__4].i * 
		0. + z__[i__4].r * -1.;
	s = a[i__2].r / ((d__1 = z__[i__3].r, abs(d__1)) + (d__2 = z__1.r, 
		abs(d__2)));
	zdscal_(n, &s, &z__[1], &c__1);
L120:
	i__2 = k;
	z_div(&z__1, &z__[k], &a[k + k * a_dim1]);
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
	i__2 = k;
	z__1.r = -z__[i__2].r, z__1.i = -z__[i__2].i;
	t.r = z__1.r, t.i = z__1.i;
	i__2 = k - 1;
	zaxpy_(&i__2, &t, &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
/* L130: */
    }
    s = 1. / dzasum_(n, &z__[1], &c__1);
    zdscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.;

/*        solve ctrans(r)*v = y */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k;
	i__3 = k;
	i__4 = k - 1;
	zdotc_(&z__2, &i__4, &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
	z__1.r = z__[i__3].r - z__2.r, z__1.i = z__[i__3].i - z__2.i;
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
	i__2 = k;
	i__3 = k;
	z__1.r = z__[i__3].r * 0. - z__[i__3].i * -1., z__1.i = z__[i__3].i * 
		0. + z__[i__3].r * -1.;
	i__4 = k + k * a_dim1;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) <= a[
		i__4].r) {
	    goto L140;
	}
	i__2 = k + k * a_dim1;
	i__3 = k;
	i__4 = k;
	z__1.r = z__[i__4].r * 0. - z__[i__4].i * -1., z__1.i = z__[i__4].i * 
		0. + z__[i__4].r * -1.;
	s = a[i__2].r / ((d__1 = z__[i__3].r, abs(d__1)) + (d__2 = z__1.r, 
		abs(d__2)));
	zdscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L140:
	i__2 = k;
	z_div(&z__1, &z__[k], &a[k + k * a_dim1]);
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
/* L150: */
    }
    s = 1. / dzasum_(n, &z__[1], &c__1);
    zdscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*        solve r*z = v */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	i__2 = k;
	i__3 = k;
	z__1.r = z__[i__3].r * 0. - z__[i__3].i * -1., z__1.i = z__[i__3].i * 
		0. + z__[i__3].r * -1.;
	i__4 = k + k * a_dim1;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) <= a[
		i__4].r) {
	    goto L160;
	}
	i__2 = k + k * a_dim1;
	i__3 = k;
	i__4 = k;
	z__1.r = z__[i__4].r * 0. - z__[i__4].i * -1., z__1.i = z__[i__4].i * 
		0. + z__[i__4].r * -1.;
	s = a[i__2].r / ((d__1 = z__[i__3].r, abs(d__1)) + (d__2 = z__1.r, 
		abs(d__2)));
	zdscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L160:
	i__2 = k;
	z_div(&z__1, &z__[k], &a[k + k * a_dim1]);
	z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
	i__2 = k;
	z__1.r = -z__[i__2].r, z__1.i = -z__[i__2].i;
	t.r = z__1.r, t.i = z__1.i;
	i__2 = k - 1;
	zaxpy_(&i__2, &t, &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
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
} /* zpoco_ */

