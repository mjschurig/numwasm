/* ssico.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int ssico_(real *a, integer *lda, integer *n, integer *kpvt, 
	real *rcond, real *z__)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    real r__1, r__2;

    /* Builtin functions */
    double r_sign(real *, real *);

    /* Local variables */
    integer i__, j, k;
    real s, t, ak, bk, ek;
    integer kp, ks, jm1, kps;
    real akm1, bkm1;
    integer info;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    real denom;
    extern /* Subroutine */ int ssifa_(real *, integer *, integer *, integer *
	    , integer *), sscal_(integer *, real *, real *, integer *);
    real anorm;
    extern doublereal sasum_(integer *, real *, integer *);
    real ynorm;
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *);


/*     ssico factors a real symmetric matrix by elimination with */
/*     symmetric pivoting and estimates the condition of the matrix. */

/*     if  rcond  is not needed, ssifa is slightly faster. */
/*     to solve  a*x = b , follow ssico by ssisl. */
/*     to compute  inverse(a)*c , follow ssico by ssisl. */
/*     to compute  inverse(a) , follow ssico by ssidi. */
/*     to compute  determinant(a) , follow ssico by ssidi. */
/*     to compute  inertia(a), follow ssico by ssidi. */

/*     on entry */

/*        a       real(lda, n) */
/*                the symmetric matrix to be factored. */
/*                only the diagonal and upper triangle are used. */

/*        lda     integer */
/*                the leading dimension of the array  a . */

/*        n       integer */
/*                the order of the matrix  a . */

/*     output */

/*        a       a block diagonal matrix and the multipliers which */
/*                were used to obtain it. */
/*                the factorization can be written  a = u*d*trans(u) */
/*                where  u  is a product of permutation and unit */
/*                upper triangular matrices , trans(u) is the */
/*                transpose of  u , and  d  is block diagonal */
/*                with 1 by 1 and 2 by 2 blocks. */

/*        kpvt    integer(n) */
/*                an integer vector of pivot indices. */

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
/*                underflows. */

/*        z       real(n) */
/*                a work vector whose contents are usually unimportant. */
/*                if  a  is close to a singular matrix, then  z  is */
/*                an approximate null vector in the sense that */
/*                norm(a*z) = rcond*norm(a)*norm(z) . */

/*     linpack. this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     linpack ssifa */
/*     blas saxpy,sdot,sscal,sasum */
/*     fortran abs,amax1,iabs,sign */

/*     internal variables */



/*     find norm of a using only upper half */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --kpvt;
    --z__;

    /* Function Body */
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z__[j] = sasum_(&j, &a[j * a_dim1 + 1], &c__1);
	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    z__[i__] += (r__1 = a[i__ + j * a_dim1], dabs(r__1));
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

    ssifa_(&a[a_offset], lda, n, &kpvt[1], &info);

/*     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) . */
/*     estimate = norm(z)/norm(y) where  a*z = y  and  a*y = e . */
/*     the components of  e  are chosen to cause maximum local */
/*     growth in the elements of w  where  u*d*w = e . */
/*     the vectors are frequently rescaled to avoid overflow. */

/*     solve u*d*w = e */

    ek = 1.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z__[j] = 0.f;
/* L50: */
    }
    k = *n;
L60:
    if (k == 0) {
	goto L120;
    }
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    kps = k + 1 - ks;
    if (kp == kps) {
	goto L70;
    }
    t = z__[kps];
    z__[kps] = z__[kp];
    z__[kp] = t;
L70:
    if (z__[k] != 0.f) {
	ek = r_sign(&ek, &z__[k]);
    }
    z__[k] += ek;
    i__1 = k - ks;
    saxpy_(&i__1, &z__[k], &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
    if (ks == 1) {
	goto L80;
    }
    if (z__[k - 1] != 0.f) {
	ek = r_sign(&ek, &z__[k - 1]);
    }
    z__[k - 1] += ek;
    i__1 = k - ks;
    saxpy_(&i__1, &z__[k - 1], &a[(k - 1) * a_dim1 + 1], &c__1, &z__[1], &
	    c__1);
L80:
    if (ks == 2) {
	goto L100;
    }
    if ((r__1 = z__[k], dabs(r__1)) <= (r__2 = a[k + k * a_dim1], dabs(r__2)))
	     {
	goto L90;
    }
    s = (r__1 = a[k + k * a_dim1], dabs(r__1)) / (r__2 = z__[k], dabs(r__2));
    sscal_(n, &s, &z__[1], &c__1);
    ek = s * ek;
L90:
    if (a[k + k * a_dim1] != 0.f) {
	z__[k] /= a[k + k * a_dim1];
    }
    if (a[k + k * a_dim1] == 0.f) {
	z__[k] = 1.f;
    }
    goto L110;
L100:
    ak = a[k + k * a_dim1] / a[k - 1 + k * a_dim1];
    akm1 = a[k - 1 + (k - 1) * a_dim1] / a[k - 1 + k * a_dim1];
    bk = z__[k] / a[k - 1 + k * a_dim1];
    bkm1 = z__[k - 1] / a[k - 1 + k * a_dim1];
    denom = ak * akm1 - 1.f;
    z__[k] = (akm1 * bk - bkm1) / denom;
    z__[k - 1] = (ak * bkm1 - bk) / denom;
L110:
    k -= ks;
    goto L60;
L120:
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

/*     solve trans(u)*y = w */

    k = 1;
L130:
    if (k > *n) {
	goto L160;
    }
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    if (k == 1) {
	goto L150;
    }
    i__1 = k - 1;
    z__[k] += sdot_(&i__1, &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
    if (ks == 2) {
	i__1 = k - 1;
	z__[k + 1] += sdot_(&i__1, &a[(k + 1) * a_dim1 + 1], &c__1, &z__[1], &
		c__1);
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k) {
	goto L140;
    }
    t = z__[k];
    z__[k] = z__[kp];
    z__[kp] = t;
L140:
L150:
    k += ks;
    goto L130;
L160:
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.f;

/*     solve u*d*v = y */

    k = *n;
L170:
    if (k == 0) {
	goto L230;
    }
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    if (k == ks) {
	goto L190;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    kps = k + 1 - ks;
    if (kp == kps) {
	goto L180;
    }
    t = z__[kps];
    z__[kps] = z__[kp];
    z__[kp] = t;
L180:
    i__1 = k - ks;
    saxpy_(&i__1, &z__[k], &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
    if (ks == 2) {
	i__1 = k - ks;
	saxpy_(&i__1, &z__[k - 1], &a[(k - 1) * a_dim1 + 1], &c__1, &z__[1], &
		c__1);
    }
L190:
    if (ks == 2) {
	goto L210;
    }
    if ((r__1 = z__[k], dabs(r__1)) <= (r__2 = a[k + k * a_dim1], dabs(r__2)))
	     {
	goto L200;
    }
    s = (r__1 = a[k + k * a_dim1], dabs(r__1)) / (r__2 = z__[k], dabs(r__2));
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;
L200:
    if (a[k + k * a_dim1] != 0.f) {
	z__[k] /= a[k + k * a_dim1];
    }
    if (a[k + k * a_dim1] == 0.f) {
	z__[k] = 1.f;
    }
    goto L220;
L210:
    ak = a[k + k * a_dim1] / a[k - 1 + k * a_dim1];
    akm1 = a[k - 1 + (k - 1) * a_dim1] / a[k - 1 + k * a_dim1];
    bk = z__[k] / a[k - 1 + k * a_dim1];
    bkm1 = z__[k - 1] / a[k - 1 + k * a_dim1];
    denom = ak * akm1 - 1.f;
    z__[k] = (akm1 * bk - bkm1) / denom;
    z__[k - 1] = (ak * bkm1 - bk) / denom;
L220:
    k -= ks;
    goto L170;
L230:
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*     solve trans(u)*z = v */

    k = 1;
L240:
    if (k > *n) {
	goto L270;
    }
    ks = 1;
    if (kpvt[k] < 0) {
	ks = 2;
    }
    if (k == 1) {
	goto L260;
    }
    i__1 = k - 1;
    z__[k] += sdot_(&i__1, &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
    if (ks == 2) {
	i__1 = k - 1;
	z__[k + 1] += sdot_(&i__1, &a[(k + 1) * a_dim1 + 1], &c__1, &z__[1], &
		c__1);
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k) {
	goto L250;
    }
    t = z__[k];
    z__[k] = z__[kp];
    z__[kp] = t;
L250:
L260:
    k += ks;
    goto L240;
L270:
/*     make znorm = 1.0 */
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.f) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.f) {
	*rcond = 0.f;
    }
    return 0;
} /* ssico_ */

