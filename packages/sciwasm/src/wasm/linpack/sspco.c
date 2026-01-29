/* sspco.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int sspco_(real *ap, integer *n, integer *kpvt, real *rcond, 
	real *z__)
{
    /* System generated locals */
    integer i__1, i__2;
    real r__1, r__2;

    /* Builtin functions */
    double r_sign(real *, real *);

    /* Local variables */
    integer i__, j, k;
    real s, t;
    integer j1;
    real ak, bk, ek;
    integer ij, ik, kk, kp, ks, jm1, kps;
    real akm1, bkm1;
    integer ikm1, km1k, ikp1, info;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    integer km1km1;
    real denom;
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *);
    real anorm;
    extern /* Subroutine */ int sspfa_(real *, integer *, integer *, integer *
	    );
    extern doublereal sasum_(integer *, real *, integer *);
    real ynorm;
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *);


/*     sspco factors a real symmetric matrix stored in packed */
/*     form by elimination with symmetric pivoting and estimates */
/*     the condition of the matrix. */

/*     if  rcond  is not needed, sspfa is slightly faster. */
/*     to solve  a*x = b , follow sspco by sspsl. */
/*     to compute  inverse(a)*c , follow sspco by sspsl. */
/*     to compute  inverse(a) , follow sspco by sspdi. */
/*     to compute  determinant(a) , follow sspco by sspdi. */
/*     to compute  inertia(a), follow sspco by sspdi. */

/*     on entry */

/*        ap      real (n*(n+1)/2) */
/*                the packed form of a symmetric matrix  a .  the */
/*                columns of the upper triangle are stored sequentially */
/*                in a one-dimensional array of length  n*(n+1)/2 . */
/*                see comments below for details. */

/*        n       integer */
/*                the order of the matrix  a . */

/*     output */

/*        ap      a block diagonal matrix and the multipliers which */
/*                were used to obtain it stored in packed form. */
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

/*     linpack. this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     linpack sspfa */
/*     blas saxpy,sdot,sscal,sasum */
/*     fortran abs,amax1,iabs,sign */

/*     internal variables */



/*     find norm of a using only upper half */

    /* Parameter adjustments */
    --z__;
    --kpvt;
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

    sspfa_(&ap[1], n, &kpvt[1], &info);

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
    ik = *n * (*n - 1) / 2;
L60:
    if (k == 0) {
	goto L120;
    }
    kk = ik + k;
    ikm1 = ik - (k - 1);
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
    saxpy_(&i__1, &z__[k], &ap[ik + 1], &c__1, &z__[1], &c__1);
    if (ks == 1) {
	goto L80;
    }
    if (z__[k - 1] != 0.f) {
	ek = r_sign(&ek, &z__[k - 1]);
    }
    z__[k - 1] += ek;
    i__1 = k - ks;
    saxpy_(&i__1, &z__[k - 1], &ap[ikm1 + 1], &c__1, &z__[1], &c__1);
L80:
    if (ks == 2) {
	goto L100;
    }
    if ((r__1 = z__[k], dabs(r__1)) <= (r__2 = ap[kk], dabs(r__2))) {
	goto L90;
    }
    s = (r__1 = ap[kk], dabs(r__1)) / (r__2 = z__[k], dabs(r__2));
    sscal_(n, &s, &z__[1], &c__1);
    ek = s * ek;
L90:
    if (ap[kk] != 0.f) {
	z__[k] /= ap[kk];
    }
    if (ap[kk] == 0.f) {
	z__[k] = 1.f;
    }
    goto L110;
L100:
    km1k = ik + k - 1;
    km1km1 = ikm1 + k - 1;
    ak = ap[kk] / ap[km1k];
    akm1 = ap[km1km1] / ap[km1k];
    bk = z__[k] / ap[km1k];
    bkm1 = z__[k - 1] / ap[km1k];
    denom = ak * akm1 - 1.f;
    z__[k] = (akm1 * bk - bkm1) / denom;
    z__[k - 1] = (ak * bkm1 - bk) / denom;
L110:
    k -= ks;
    ik -= k;
    if (ks == 2) {
	ik -= k + 1;
    }
    goto L60;
L120:
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

/*     solve trans(u)*y = w */

    k = 1;
    ik = 0;
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
    z__[k] += sdot_(&i__1, &ap[ik + 1], &c__1, &z__[1], &c__1);
    ikp1 = ik + k;
    if (ks == 2) {
	i__1 = k - 1;
	z__[k + 1] += sdot_(&i__1, &ap[ikp1 + 1], &c__1, &z__[1], &c__1);
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
    ik += k;
    if (ks == 2) {
	ik += k + 1;
    }
    k += ks;
    goto L130;
L160:
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.f;

/*     solve u*d*v = y */

    k = *n;
    ik = *n * (*n - 1) / 2;
L170:
    if (k == 0) {
	goto L230;
    }
    kk = ik + k;
    ikm1 = ik - (k - 1);
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
    saxpy_(&i__1, &z__[k], &ap[ik + 1], &c__1, &z__[1], &c__1);
    if (ks == 2) {
	i__1 = k - ks;
	saxpy_(&i__1, &z__[k - 1], &ap[ikm1 + 1], &c__1, &z__[1], &c__1);
    }
L190:
    if (ks == 2) {
	goto L210;
    }
    if ((r__1 = z__[k], dabs(r__1)) <= (r__2 = ap[kk], dabs(r__2))) {
	goto L200;
    }
    s = (r__1 = ap[kk], dabs(r__1)) / (r__2 = z__[k], dabs(r__2));
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;
L200:
    if (ap[kk] != 0.f) {
	z__[k] /= ap[kk];
    }
    if (ap[kk] == 0.f) {
	z__[k] = 1.f;
    }
    goto L220;
L210:
    km1k = ik + k - 1;
    km1km1 = ikm1 + k - 1;
    ak = ap[kk] / ap[km1k];
    akm1 = ap[km1km1] / ap[km1k];
    bk = z__[k] / ap[km1k];
    bkm1 = z__[k - 1] / ap[km1k];
    denom = ak * akm1 - 1.f;
    z__[k] = (akm1 * bk - bkm1) / denom;
    z__[k - 1] = (ak * bkm1 - bk) / denom;
L220:
    k -= ks;
    ik -= k;
    if (ks == 2) {
	ik -= k + 1;
    }
    goto L170;
L230:
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*     solve trans(u)*z = v */

    k = 1;
    ik = 0;
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
    z__[k] += sdot_(&i__1, &ap[ik + 1], &c__1, &z__[1], &c__1);
    ikp1 = ik + k;
    if (ks == 2) {
	i__1 = k - 1;
	z__[k + 1] += sdot_(&i__1, &ap[ikp1 + 1], &c__1, &z__[1], &c__1);
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
    ik += k;
    if (ks == 2) {
	ik += k + 1;
    }
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
} /* sspco_ */

