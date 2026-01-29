/* dppco.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dppco_(doublereal *ap, integer *n, doublereal *rcond, 
	doublereal *z__, integer *info)
{
    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *);

    /* Local variables */
    integer i__, j, k;
    doublereal s, t;
    integer j1, kb;
    doublereal ek;
    integer ij, kj, kk;
    doublereal sm, wk;
    integer jm1, kp1;
    doublereal wkm;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dppfa_(doublereal *, integer *, integer *);
    extern doublereal dasum_(integer *, doublereal *, integer *);
    doublereal anorm;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    doublereal ynorm;


/*     dppco factors a double precision symmetric positive definite */
/*     matrix stored in packed form */
/*     and estimates the condition of the matrix. */

/*     if  rcond  is not needed, dppfa is slightly faster. */
/*     to solve  a*x = b , follow dppco by dppsl. */
/*     to compute  inverse(a)*c , follow dppco by dppsl. */
/*     to compute  determinant(a) , follow dppco by dppdi. */
/*     to compute  inverse(a) , follow dppco by dppdi. */

/*     on entry */

/*        ap      double precision (n*(n+1)/2) */
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

/*     linpack dppfa */
/*     blas daxpy,ddot,dscal,dasum */
/*     fortran dabs,dmax1,dreal,dsign */

/*     internal variables */



/*     find norm of a */

    /* Parameter adjustments */
    --z__;
    --ap;

    /* Function Body */
    j1 = 1;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z__[j] = dasum_(&j, &ap[j1], &c__1);
	ij = j1;
	j1 += j;
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    z__[i__] += (d__1 = ap[ij], abs(d__1));
	    ++ij;
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

    dppfa_(&ap[1], n, info);
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
    kk = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	kk += k;
	if (z__[k] != 0.) {
	    d__1 = -z__[k];
	    ek = d_sign(&ek, &d__1);
	}
	if ((d__1 = ek - z__[k], abs(d__1)) <= ap[kk]) {
	    goto L60;
	}
	s = ap[kk] / (d__1 = ek - z__[k], abs(d__1));
	dscal_(n, &s, &z__[1], &c__1);
	ek = s * ek;
L60:
	wk = ek - z__[k];
	wkm = -ek - z__[k];
	s = abs(wk);
	sm = abs(wkm);
	wk /= ap[kk];
	wkm /= ap[kk];
	kp1 = k + 1;
	kj = kk + k;
	if (kp1 > *n) {
	    goto L100;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    sm += (d__1 = z__[j] + wkm * ap[kj], abs(d__1));
	    z__[j] += wk * ap[kj];
	    s += (d__1 = z__[j], abs(d__1));
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
    s = 1. / dasum_(n, &z__[1], &c__1);
    dscal_(n, &s, &z__[1], &c__1);

/*        solve r*y = w */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if ((d__1 = z__[k], abs(d__1)) <= ap[kk]) {
	    goto L120;
	}
	s = ap[kk] / (d__1 = z__[k], abs(d__1));
	dscal_(n, &s, &z__[1], &c__1);
L120:
	z__[k] /= ap[kk];
	kk -= k;
	t = -z__[k];
	i__2 = k - 1;
	daxpy_(&i__2, &t, &ap[kk + 1], &c__1, &z__[1], &c__1);
/* L130: */
    }
    s = 1. / dasum_(n, &z__[1], &c__1);
    dscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.;

/*        solve trans(r)*v = y */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k - 1;
	z__[k] -= ddot_(&i__2, &ap[kk + 1], &c__1, &z__[1], &c__1);
	kk += k;
	if ((d__1 = z__[k], abs(d__1)) <= ap[kk]) {
	    goto L140;
	}
	s = ap[kk] / (d__1 = z__[k], abs(d__1));
	dscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L140:
	z__[k] /= ap[kk];
/* L150: */
    }
    s = 1. / dasum_(n, &z__[1], &c__1);
    dscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*        solve r*z = v */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if ((d__1 = z__[k], abs(d__1)) <= ap[kk]) {
	    goto L160;
	}
	s = ap[kk] / (d__1 = z__[k], abs(d__1));
	dscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L160:
	z__[k] /= ap[kk];
	kk -= k;
	t = -z__[k];
	i__2 = k - 1;
	daxpy_(&i__2, &t, &ap[kk + 1], &c__1, &z__[1], &c__1);
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
} /* dppco_ */

