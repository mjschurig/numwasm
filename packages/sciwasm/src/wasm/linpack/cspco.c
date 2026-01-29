/* cspco.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int cspco_(complex *ap, integer *n, integer *kpvt, real *
	rcond, complex *z__)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5;
    real r__1, r__2, r__3, r__4, r__5, r__6, r__7, r__8;
    complex q__1, q__2, q__3;

    /* Builtin functions */
    double r_imag(complex *);
    void c_div(complex *, complex *, complex *);

    /* Local variables */
    integer i__, j, k;
    real s;
    complex t;
    integer j1;
    complex ak, bk, ek;
    integer ij, ik, kk, kp, ks, jm1, kps;
    complex akm1, bkm1;
    integer ikm1, km1k, ikp1, info;
    extern /* Subroutine */ int cspfa_(complex *, integer *, integer *, 
	    integer *);
    integer km1km1;
    complex denom;
    real anorm;
    extern /* Complex */ VOID cdotu_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);
    extern /* Subroutine */ int caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *);
    real ynorm;
    extern /* Subroutine */ int csscal_(integer *, real *, complex *, integer 
	    *);
    extern doublereal scasum_(integer *, complex *, integer *);


/*     cspco factors a complex symmetric matrix stored in packed */
/*     form by elimination with symmetric pivoting and estimates */
/*     the condition of the matrix. */

/*     if  rcond  is not needed, cspfa is slightly faster. */
/*     to solve  a*x = b , follow cspco by cspsl. */
/*     to compute  inverse(a)*c , follow cspco by cspsl. */
/*     to compute  inverse(a) , follow cspco by cspdi. */
/*     to compute  determinant(a) , follow cspco by cspdi. */

/*     on entry */

/*        ap      complex (n*(n+1)/2) */
/*                the packed form of a symmetric matrix  a .  the */
/*                columns of the upper triangle are stored sequentially */
/*                in a one-dimensional array of length  n*(n+1)/2 . */
/*                see comments below for details. */

/*        n       integer */
/*                the order of the matrix  a . */

/*     on return */

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

/*        z       complex(n) */
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

/*     linpack cspfa */
/*     blas caxpy,cdotu,csscal,scasum */
/*     fortran abs,aimag,amax1,cmplx,iabs,real */

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

    cspfa_(&ap[1], n, &kpvt[1], &info);

/*     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) . */
/*     estimate = norm(z)/norm(y) where  a*z = y  and  a*y = e . */
/*     the components of  e  are chosen to cause maximum local */
/*     growth in the elements of w  where  u*d*w = e . */
/*     the vectors are frequently rescaled to avoid overflow. */

/*     solve u*d*w = e */

    ek.r = 1.f, ek.i = 0.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	z__[i__2].r = 0.f, z__[i__2].i = 0.f;
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
    i__1 = kps;
    t.r = z__[i__1].r, t.i = z__[i__1].i;
    i__1 = kps;
    i__2 = kp;
    z__[i__1].r = z__[i__2].r, z__[i__1].i = z__[i__2].i;
    i__1 = kp;
    z__[i__1].r = t.r, z__[i__1].i = t.i;
L70:
    i__1 = k;
    if ((r__1 = z__[i__1].r, dabs(r__1)) + (r__2 = r_imag(&z__[k]), dabs(r__2)
	    ) != 0.f) {
	r__7 = (r__3 = ek.r, dabs(r__3)) + (r__4 = r_imag(&ek), dabs(r__4));
	i__2 = k;
	i__3 = k;
	r__8 = (r__5 = z__[i__3].r, dabs(r__5)) + (r__6 = r_imag(&z__[k]), 
		dabs(r__6));
	q__2.r = z__[i__2].r / r__8, q__2.i = z__[i__2].i / r__8;
	q__1.r = r__7 * q__2.r, q__1.i = r__7 * q__2.i;
	ek.r = q__1.r, ek.i = q__1.i;
    }
    i__1 = k;
    i__2 = k;
    q__1.r = z__[i__2].r + ek.r, q__1.i = z__[i__2].i + ek.i;
    z__[i__1].r = q__1.r, z__[i__1].i = q__1.i;
    i__1 = k - ks;
    caxpy_(&i__1, &z__[k], &ap[ik + 1], &c__1, &z__[1], &c__1);
    if (ks == 1) {
	goto L80;
    }
    i__1 = k - 1;
    if ((r__1 = z__[i__1].r, dabs(r__1)) + (r__2 = r_imag(&z__[k - 1]), dabs(
	    r__2)) != 0.f) {
	r__7 = (r__3 = ek.r, dabs(r__3)) + (r__4 = r_imag(&ek), dabs(r__4));
	i__2 = k - 1;
	i__3 = k - 1;
	r__8 = (r__5 = z__[i__3].r, dabs(r__5)) + (r__6 = r_imag(&z__[k - 1]),
		 dabs(r__6));
	q__2.r = z__[i__2].r / r__8, q__2.i = z__[i__2].i / r__8;
	q__1.r = r__7 * q__2.r, q__1.i = r__7 * q__2.i;
	ek.r = q__1.r, ek.i = q__1.i;
    }
    i__1 = k - 1;
    i__2 = k - 1;
    q__1.r = z__[i__2].r + ek.r, q__1.i = z__[i__2].i + ek.i;
    z__[i__1].r = q__1.r, z__[i__1].i = q__1.i;
    i__1 = k - ks;
    caxpy_(&i__1, &z__[k - 1], &ap[ikm1 + 1], &c__1, &z__[1], &c__1);
L80:
    if (ks == 2) {
	goto L100;
    }
    i__1 = k;
    i__2 = kk;
    if ((r__1 = z__[i__1].r, dabs(r__1)) + (r__2 = r_imag(&z__[k]), dabs(r__2)
	    ) <= (r__3 = ap[i__2].r, dabs(r__3)) + (r__4 = r_imag(&ap[kk]), 
	    dabs(r__4))) {
	goto L90;
    }
    i__1 = kk;
    i__2 = k;
    s = ((r__1 = ap[i__1].r, dabs(r__1)) + (r__2 = r_imag(&ap[kk]), dabs(r__2)
	    )) / ((r__3 = z__[i__2].r, dabs(r__3)) + (r__4 = r_imag(&z__[k]), 
	    dabs(r__4)));
    csscal_(n, &s, &z__[1], &c__1);
    q__2.r = s, q__2.i = 0.f;
    q__1.r = q__2.r * ek.r - q__2.i * ek.i, q__1.i = q__2.r * ek.i + q__2.i * 
	    ek.r;
    ek.r = q__1.r, ek.i = q__1.i;
L90:
    i__1 = kk;
    if ((r__1 = ap[i__1].r, dabs(r__1)) + (r__2 = r_imag(&ap[kk]), dabs(r__2))
	     != 0.f) {
	i__2 = k;
	c_div(&q__1, &z__[k], &ap[kk]);
	z__[i__2].r = q__1.r, z__[i__2].i = q__1.i;
    }
    i__1 = kk;
    if ((r__1 = ap[i__1].r, dabs(r__1)) + (r__2 = r_imag(&ap[kk]), dabs(r__2))
	     == 0.f) {
	i__2 = k;
	z__[i__2].r = 1.f, z__[i__2].i = 0.f;
    }
    goto L110;
L100:
    km1k = ik + k - 1;
    km1km1 = ikm1 + k - 1;
    c_div(&q__1, &ap[kk], &ap[km1k]);
    ak.r = q__1.r, ak.i = q__1.i;
    c_div(&q__1, &ap[km1km1], &ap[km1k]);
    akm1.r = q__1.r, akm1.i = q__1.i;
    c_div(&q__1, &z__[k], &ap[km1k]);
    bk.r = q__1.r, bk.i = q__1.i;
    c_div(&q__1, &z__[k - 1], &ap[km1k]);
    bkm1.r = q__1.r, bkm1.i = q__1.i;
    q__2.r = ak.r * akm1.r - ak.i * akm1.i, q__2.i = ak.r * akm1.i + ak.i * 
	    akm1.r;
    q__1.r = q__2.r - 1.f, q__1.i = q__2.i;
    denom.r = q__1.r, denom.i = q__1.i;
    i__1 = k;
    q__3.r = akm1.r * bk.r - akm1.i * bk.i, q__3.i = akm1.r * bk.i + akm1.i * 
	    bk.r;
    q__2.r = q__3.r - bkm1.r, q__2.i = q__3.i - bkm1.i;
    c_div(&q__1, &q__2, &denom);
    z__[i__1].r = q__1.r, z__[i__1].i = q__1.i;
    i__1 = k - 1;
    q__3.r = ak.r * bkm1.r - ak.i * bkm1.i, q__3.i = ak.r * bkm1.i + ak.i * 
	    bkm1.r;
    q__2.r = q__3.r - bk.r, q__2.i = q__3.i - bk.i;
    c_div(&q__1, &q__2, &denom);
    z__[i__1].r = q__1.r, z__[i__1].i = q__1.i;
L110:
    k -= ks;
    ik -= k;
    if (ks == 2) {
	ik -= k + 1;
    }
    goto L60;
L120:
    s = 1.f / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

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
    i__1 = k;
    i__2 = k;
    i__3 = k - 1;
    cdotu_(&q__2, &i__3, &ap[ik + 1], &c__1, &z__[1], &c__1);
    q__1.r = z__[i__2].r + q__2.r, q__1.i = z__[i__2].i + q__2.i;
    z__[i__1].r = q__1.r, z__[i__1].i = q__1.i;
    ikp1 = ik + k;
    if (ks == 2) {
	i__1 = k + 1;
	i__2 = k + 1;
	i__3 = k - 1;
	cdotu_(&q__2, &i__3, &ap[ikp1 + 1], &c__1, &z__[1], &c__1);
	q__1.r = z__[i__2].r + q__2.r, q__1.i = z__[i__2].i + q__2.i;
	z__[i__1].r = q__1.r, z__[i__1].i = q__1.i;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k) {
	goto L140;
    }
    i__1 = k;
    t.r = z__[i__1].r, t.i = z__[i__1].i;
    i__1 = k;
    i__2 = kp;
    z__[i__1].r = z__[i__2].r, z__[i__1].i = z__[i__2].i;
    i__1 = kp;
    z__[i__1].r = t.r, z__[i__1].i = t.i;
L140:
L150:
    ik += k;
    if (ks == 2) {
	ik += k + 1;
    }
    k += ks;
    goto L130;
L160:
    s = 1.f / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);

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
    i__1 = kps;
    t.r = z__[i__1].r, t.i = z__[i__1].i;
    i__1 = kps;
    i__2 = kp;
    z__[i__1].r = z__[i__2].r, z__[i__1].i = z__[i__2].i;
    i__1 = kp;
    z__[i__1].r = t.r, z__[i__1].i = t.i;
L180:
    i__1 = k - ks;
    caxpy_(&i__1, &z__[k], &ap[ik + 1], &c__1, &z__[1], &c__1);
    if (ks == 2) {
	i__1 = k - ks;
	caxpy_(&i__1, &z__[k - 1], &ap[ikm1 + 1], &c__1, &z__[1], &c__1);
    }
L190:
    if (ks == 2) {
	goto L210;
    }
    i__1 = k;
    i__2 = kk;
    if ((r__1 = z__[i__1].r, dabs(r__1)) + (r__2 = r_imag(&z__[k]), dabs(r__2)
	    ) <= (r__3 = ap[i__2].r, dabs(r__3)) + (r__4 = r_imag(&ap[kk]), 
	    dabs(r__4))) {
	goto L200;
    }
    i__1 = kk;
    i__2 = k;
    s = ((r__1 = ap[i__1].r, dabs(r__1)) + (r__2 = r_imag(&ap[kk]), dabs(r__2)
	    )) / ((r__3 = z__[i__2].r, dabs(r__3)) + (r__4 = r_imag(&z__[k]), 
	    dabs(r__4)));
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;
L200:
    i__1 = kk;
    if ((r__1 = ap[i__1].r, dabs(r__1)) + (r__2 = r_imag(&ap[kk]), dabs(r__2))
	     != 0.f) {
	i__2 = k;
	c_div(&q__1, &z__[k], &ap[kk]);
	z__[i__2].r = q__1.r, z__[i__2].i = q__1.i;
    }
    i__1 = kk;
    if ((r__1 = ap[i__1].r, dabs(r__1)) + (r__2 = r_imag(&ap[kk]), dabs(r__2))
	     == 0.f) {
	i__2 = k;
	z__[i__2].r = 1.f, z__[i__2].i = 0.f;
    }
    goto L220;
L210:
    km1k = ik + k - 1;
    km1km1 = ikm1 + k - 1;
    c_div(&q__1, &ap[kk], &ap[km1k]);
    ak.r = q__1.r, ak.i = q__1.i;
    c_div(&q__1, &ap[km1km1], &ap[km1k]);
    akm1.r = q__1.r, akm1.i = q__1.i;
    c_div(&q__1, &z__[k], &ap[km1k]);
    bk.r = q__1.r, bk.i = q__1.i;
    c_div(&q__1, &z__[k - 1], &ap[km1k]);
    bkm1.r = q__1.r, bkm1.i = q__1.i;
    q__2.r = ak.r * akm1.r - ak.i * akm1.i, q__2.i = ak.r * akm1.i + ak.i * 
	    akm1.r;
    q__1.r = q__2.r - 1.f, q__1.i = q__2.i;
    denom.r = q__1.r, denom.i = q__1.i;
    i__1 = k;
    q__3.r = akm1.r * bk.r - akm1.i * bk.i, q__3.i = akm1.r * bk.i + akm1.i * 
	    bk.r;
    q__2.r = q__3.r - bkm1.r, q__2.i = q__3.i - bkm1.i;
    c_div(&q__1, &q__2, &denom);
    z__[i__1].r = q__1.r, z__[i__1].i = q__1.i;
    i__1 = k - 1;
    q__3.r = ak.r * bkm1.r - ak.i * bkm1.i, q__3.i = ak.r * bkm1.i + ak.i * 
	    bkm1.r;
    q__2.r = q__3.r - bk.r, q__2.i = q__3.i - bk.i;
    c_div(&q__1, &q__2, &denom);
    z__[i__1].r = q__1.r, z__[i__1].i = q__1.i;
L220:
    k -= ks;
    ik -= k;
    if (ks == 2) {
	ik -= k + 1;
    }
    goto L170;
L230:
    s = 1.f / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
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
    i__1 = k;
    i__2 = k;
    i__3 = k - 1;
    cdotu_(&q__2, &i__3, &ap[ik + 1], &c__1, &z__[1], &c__1);
    q__1.r = z__[i__2].r + q__2.r, q__1.i = z__[i__2].i + q__2.i;
    z__[i__1].r = q__1.r, z__[i__1].i = q__1.i;
    ikp1 = ik + k;
    if (ks == 2) {
	i__1 = k + 1;
	i__2 = k + 1;
	i__3 = k - 1;
	cdotu_(&q__2, &i__3, &ap[ikp1 + 1], &c__1, &z__[1], &c__1);
	q__1.r = z__[i__2].r + q__2.r, q__1.i = z__[i__2].i + q__2.i;
	z__[i__1].r = q__1.r, z__[i__1].i = q__1.i;
    }
    kp = (i__1 = kpvt[k], abs(i__1));
    if (kp == k) {
	goto L250;
    }
    i__1 = k;
    t.r = z__[i__1].r, t.i = z__[i__1].i;
    i__1 = k;
    i__2 = kp;
    z__[i__1].r = z__[i__2].r, z__[i__1].i = z__[i__2].i;
    i__1 = kp;
    z__[i__1].r = t.r, z__[i__1].i = t.i;
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
    s = 1.f / scasum_(n, &z__[1], &c__1);
    csscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.f) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.f) {
	*rcond = 0.f;
    }
    return 0;
} /* cspco_ */

