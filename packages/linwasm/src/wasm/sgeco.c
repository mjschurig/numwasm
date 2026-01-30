/* sgeco.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int sgeco_(real *a, integer *lda, integer *n, integer *ipvt, 
	real *rcond, real *z__)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;
    real r__1, r__2;

    /* Builtin functions */
    double r_sign(real *, real *);

    /* Local variables */
    integer j, k, l;
    real s, t;
    integer kb;
    real ek, sm, wk;
    integer kp1;
    real wkm;
    integer info;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    extern /* Subroutine */ int sgefa_(real *, integer *, integer *, integer *
	    , integer *), sscal_(integer *, real *, real *, integer *);
    real anorm;
    extern doublereal sasum_(integer *, real *, integer *);
    real ynorm;
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *);


/*     sgeco factors a real matrix by gaussian elimination */
/*     and estimates the condition of the matrix. */

/*     if  rcond  is not needed, sgefa is slightly faster. */
/*     to solve  a*x = b , follow sgeco by sgesl. */
/*     to compute  inverse(a)*c , follow sgeco by sgesl. */
/*     to compute  determinant(a) , follow sgeco by sgedi. */
/*     to compute  inverse(a) , follow sgeco by sgedi. */

/*     on entry */

/*        a       real(lda, n) */
/*                the matrix to be factored. */

/*        lda     integer */
/*                the leading dimension of the array  a . */

/*        n       integer */
/*                the order of the matrix  a . */

/*     on return */

/*        a       an upper triangular matrix and the multipliers */
/*                which were used to obtain it. */
/*                the factorization can be written  a = l*u  where */
/*                l  is a product of permutation and unit lower */
/*                triangular matrices and  u  is upper triangular. */

/*        ipvt    integer(n) */
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

/*     linpack sgefa */
/*     blas saxpy,sdot,sscal,sasum */
/*     fortran abs,amax1,sign */

/*     internal variables */



/*     compute 1-norm of a */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipvt;
    --z__;

    /* Function Body */
    anorm = 0.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	r__1 = anorm, r__2 = sasum_(n, &a[j * a_dim1 + 1], &c__1);
	anorm = dmax(r__1,r__2);
/* L10: */
    }

/*     factor */

    sgefa_(&a[a_offset], lda, n, &ipvt[1], &info);

/*     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) . */
/*     estimate = norm(z)/norm(y) where  a*z = y  and  trans(a)*y = e . */
/*     trans(a)  is the transpose of a .  the components of  e  are */
/*     chosen to cause maximum local growth in the elements of w  where */
/*     trans(u)*w = e .  the vectors are frequently rescaled to avoid */
/*     overflow. */

/*     solve trans(u)*w = e */

    ek = 1.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z__[j] = 0.f;
/* L20: */
    }
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (z__[k] != 0.f) {
	    r__1 = -z__[k];
	    ek = r_sign(&ek, &r__1);
	}
	if ((r__1 = ek - z__[k], dabs(r__1)) <= (r__2 = a[k + k * a_dim1], 
		dabs(r__2))) {
	    goto L30;
	}
	s = (r__1 = a[k + k * a_dim1], dabs(r__1)) / (r__2 = ek - z__[k], 
		dabs(r__2));
	sscal_(n, &s, &z__[1], &c__1);
	ek = s * ek;
L30:
	wk = ek - z__[k];
	wkm = -ek - z__[k];
	s = dabs(wk);
	sm = dabs(wkm);
	if (a[k + k * a_dim1] == 0.f) {
	    goto L40;
	}
	wk /= a[k + k * a_dim1];
	wkm /= a[k + k * a_dim1];
	goto L50;
L40:
	wk = 1.f;
	wkm = 1.f;
L50:
	kp1 = k + 1;
	if (kp1 > *n) {
	    goto L90;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    sm += (r__1 = z__[j] + wkm * a[k + j * a_dim1], dabs(r__1));
	    z__[j] += wk * a[k + j * a_dim1];
	    s += (r__1 = z__[j], dabs(r__1));
/* L60: */
	}
	if (s >= sm) {
	    goto L80;
	}
	t = wkm - wk;
	wk = wkm;
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    z__[j] += t * a[k + j * a_dim1];
/* L70: */
	}
L80:
L90:
	z__[k] = wk;
/* L100: */
    }
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

/*     solve trans(l)*y = w */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if (k < *n) {
	    i__2 = *n - k;
	    z__[k] += sdot_(&i__2, &a[k + 1 + k * a_dim1], &c__1, &z__[k + 1],
		     &c__1);
	}
	if ((r__1 = z__[k], dabs(r__1)) <= 1.f) {
	    goto L110;
	}
	s = 1.f / (r__1 = z__[k], dabs(r__1));
	sscal_(n, &s, &z__[1], &c__1);
L110:
	l = ipvt[k];
	t = z__[l];
	z__[l] = z__[k];
	z__[k] = t;
/* L120: */
    }
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.f;

/*     solve l*v = y */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	l = ipvt[k];
	t = z__[l];
	z__[l] = z__[k];
	z__[k] = t;
	if (k < *n) {
	    i__2 = *n - k;
	    saxpy_(&i__2, &t, &a[k + 1 + k * a_dim1], &c__1, &z__[k + 1], &
		    c__1);
	}
	if ((r__1 = z__[k], dabs(r__1)) <= 1.f) {
	    goto L130;
	}
	s = 1.f / (r__1 = z__[k], dabs(r__1));
	sscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L130:
/* L140: */
	;
    }
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*     solve  u*z = v */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if ((r__1 = z__[k], dabs(r__1)) <= (r__2 = a[k + k * a_dim1], dabs(
		r__2))) {
	    goto L150;
	}
	s = (r__1 = a[k + k * a_dim1], dabs(r__1)) / (r__2 = z__[k], dabs(
		r__2));
	sscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L150:
	if (a[k + k * a_dim1] != 0.f) {
	    z__[k] /= a[k + k * a_dim1];
	}
	if (a[k + k * a_dim1] == 0.f) {
	    z__[k] = 1.f;
	}
	t = -z__[k];
	i__2 = k - 1;
	saxpy_(&i__2, &t, &a[k * a_dim1 + 1], &c__1, &z__[1], &c__1);
/* L160: */
    }
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
} /* sgeco_ */

