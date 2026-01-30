/* sppco.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int sppco_(real *ap, integer *n, real *rcond, real *z__, 
	integer *info)
{
    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2;

    /* Builtin functions */
    double r_sign(real *, real *);

    /* Local variables */
    integer i__, j, k;
    real s, t;
    integer j1, kb;
    real ek;
    integer ij, kj, kk;
    real sm, wk;
    integer jm1, kp1;
    real wkm;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *), 
	    sppfa_(real *, integer *, integer *);
    real anorm;
    extern doublereal sasum_(integer *, real *, integer *);
    real ynorm;
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *);


/*     sppco factors a real symmetric positive definite matrix */
/*     stored in packed form */
/*     and estimates the condition of the matrix. */

/*     if  rcond  is not needed, sppfa is slightly faster. */
/*     to solve  a*x = b , follow sppco by sppsl. */
/*     to compute  inverse(a)*c , follow sppco by sppsl. */
/*     to compute  determinant(a) , follow sppco by sppdi. */
/*     to compute  inverse(a) , follow sppco by sppdi. */

/*     on entry */

/*        ap      real (n*(n+1)/2) */
/*                the packed form of a symmetric matrix  a .  the */
/*                columns of the upper triangle are stored sequentially */
/*                in a one-dimensional array of length  n*(n+1)/2 . */
/*                see comments below for details. */

/*        n       integer */
/*                the order of the matrix  a . */

/*     on return */

/*        ap      an upper triangular matrix  r , stored in packed */
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

/*     packed storage */

/*          the following program segment will pack the upper */
/*          triangle of a symmetric matrix. */

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

/*     linpack sppfa */
/*     blas saxpy,sdot,sscal,sasum */
/*     fortran abs,amax1,real,sign */

/*     internal variables */



/*     find norm of a */

    /* Parameter adjustments */
    --z__;
    --ap;

    /* Function Body */
    j1 = 1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z__[j] = sasum_(&j, &ap[j1], &c__1);
	ij = j1;
	j1 += j;
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    z__[i__] += (r__1 = ap[ij], dabs(r__1));
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
	r__1 = anorm, r__2 = z__[j];
	anorm = dmax(r__1,r__2);
/* L40: */
    }

/*     factor */

    sppfa_(&ap[1], n, info);
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
    kk = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	kk += k;
	if (z__[k] != 0.f) {
	    r__1 = -z__[k];
	    ek = r_sign(&ek, &r__1);
	}
	if ((r__1 = ek - z__[k], dabs(r__1)) <= ap[kk]) {
	    goto L60;
	}
	s = ap[kk] / (r__1 = ek - z__[k], dabs(r__1));
	sscal_(n, &s, &z__[1], &c__1);
	ek = s * ek;
L60:
	wk = ek - z__[k];
	wkm = -ek - z__[k];
	s = dabs(wk);
	sm = dabs(wkm);
	wk /= ap[kk];
	wkm /= ap[kk];
	kp1 = k + 1;
	kj = kk + k;
	if (kp1 > *n) {
	    goto L100;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    sm += (r__1 = z__[j] + wkm * ap[kj], dabs(r__1));
	    z__[j] += wk * ap[kj];
	    s += (r__1 = z__[j], dabs(r__1));
	    kj += j;
/* L70: */
	}
	if (s >= sm) {
	    goto L90;
	}
	t = wkm - wk;
	wk = wkm;
	kj = kk + k;
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    z__[j] += t * ap[kj];
	    kj += j;
/* L80: */
	}
L90:
L100:
	z__[k] = wk;
/* L110: */
    }
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

/*        solve r*y = w */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if ((r__1 = z__[k], dabs(r__1)) <= ap[kk]) {
	    goto L120;
	}
	s = ap[kk] / (r__1 = z__[k], dabs(r__1));
	sscal_(n, &s, &z__[1], &c__1);
L120:
	z__[k] /= ap[kk];
	kk -= k;
	t = -z__[k];
	i__2 = k - 1;
	saxpy_(&i__2, &t, &ap[kk + 1], &c__1, &z__[1], &c__1);
/* L130: */
    }
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.f;

/*        solve trans(r)*v = y */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k - 1;
	z__[k] -= sdot_(&i__2, &ap[kk + 1], &c__1, &z__[1], &c__1);
	kk += k;
	if ((r__1 = z__[k], dabs(r__1)) <= ap[kk]) {
	    goto L140;
	}
	s = ap[kk] / (r__1 = z__[k], dabs(r__1));
	sscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L140:
	z__[k] /= ap[kk];
/* L150: */
    }
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*        solve r*z = v */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if ((r__1 = z__[k], dabs(r__1)) <= ap[kk]) {
	    goto L160;
	}
	s = ap[kk] / (r__1 = z__[k], dabs(r__1));
	sscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L160:
	z__[k] /= ap[kk];
	kk -= k;
	t = -z__[k];
	i__2 = k - 1;
	saxpy_(&i__2, &t, &ap[kk + 1], &c__1, &z__[1], &c__1);
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
} /* sppco_ */

