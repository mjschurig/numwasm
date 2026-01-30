/* strco.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int strco_(real *t, integer *ldt, integer *n, real *rcond, 
	real *z__, integer *job)
{
    /* System generated locals */
    integer t_dim1, t_offset, i__1, i__2;
    real r__1, r__2;

    /* Builtin functions */
    double r_sign(real *, real *);

    /* Local variables */
    integer j, k, l;
    real s, w;
    integer i1, j1, j2;
    real ek;
    integer kk;
    real sm, wk, wkm;
    extern /* Subroutine */ int sscal_(integer *, real *, real *, integer *);
    extern doublereal sasum_(integer *, real *, integer *);
    logical lower;
    real tnorm, ynorm;
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *);


/*     strco estimates the condition of a real triangular matrix. */

/*     on entry */

/*        t       real(ldt,n) */
/*                t contains the triangular matrix. the zero */
/*                elements of the matrix are not referenced, and */
/*                the corresponding elements of the array can be */
/*                used to store other information. */

/*        ldt     integer */
/*                ldt is the leading dimension of the array t. */

/*        n       integer */
/*                n is the order of the system. */

/*        job     integer */
/*                = 0         t  is lower triangular. */
/*                = nonzero   t  is upper triangular. */

/*     on return */

/*        rcond   real */
/*                an estimate of the reciprocal condition of  t . */
/*                for the system  t*x = b , relative perturbations */
/*                in  t  and  b  of size  epsilon  may cause */
/*                relative perturbations in  x  of size  epsilon/rcond . */
/*                if  rcond  is so small that the logical expression */
/*                           1.0 + rcond .eq. 1.0 */
/*                is true, then  t  may be singular to working */
/*                precision.  in particular,  rcond  is zero  if */
/*                exact singularity is detected or the estimate */
/*                underflows. */

/*        z       real(n) */
/*                a work vector whose contents are usually unimportant. */
/*                if  t  is close to a singular matrix, then  z  is */
/*                an approximate null vector in the sense that */
/*                norm(a*z) = rcond*norm(a)*norm(z) . */

/*     linpack. this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas saxpy,sscal,sasum */
/*     fortran abs,amax1,sign */

/*     internal variables */


    /* Parameter adjustments */
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --z__;

    /* Function Body */
    lower = *job == 0;

/*     compute 1-norm of t */

    tnorm = 0.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	l = j;
	if (lower) {
	    l = *n + 1 - j;
	}
	i1 = 1;
	if (lower) {
	    i1 = j;
	}
/* Computing MAX */
	r__1 = tnorm, r__2 = sasum_(&l, &t[i1 + j * t_dim1], &c__1);
	tnorm = dmax(r__1,r__2);
/* L10: */
    }

/*     rcond = 1/(norm(t)*(estimate of norm(inverse(t)))) . */
/*     estimate = norm(z)/norm(y) where  t*z = y  and  trans(t)*y = e . */
/*     trans(t)  is the transpose of t . */
/*     the components of  e  are chosen to cause maximum local */
/*     growth in the elements of y . */
/*     the vectors are frequently rescaled to avoid overflow. */

/*     solve trans(t)*y = e */

    ek = 1.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z__[j] = 0.f;
/* L20: */
    }
    i__1 = *n;
    for (kk = 1; kk <= i__1; ++kk) {
	k = kk;
	if (lower) {
	    k = *n + 1 - kk;
	}
	if (z__[k] != 0.f) {
	    r__1 = -z__[k];
	    ek = r_sign(&ek, &r__1);
	}
	if ((r__1 = ek - z__[k], dabs(r__1)) <= (r__2 = t[k + k * t_dim1], 
		dabs(r__2))) {
	    goto L30;
	}
	s = (r__1 = t[k + k * t_dim1], dabs(r__1)) / (r__2 = ek - z__[k], 
		dabs(r__2));
	sscal_(n, &s, &z__[1], &c__1);
	ek = s * ek;
L30:
	wk = ek - z__[k];
	wkm = -ek - z__[k];
	s = dabs(wk);
	sm = dabs(wkm);
	if (t[k + k * t_dim1] == 0.f) {
	    goto L40;
	}
	wk /= t[k + k * t_dim1];
	wkm /= t[k + k * t_dim1];
	goto L50;
L40:
	wk = 1.f;
	wkm = 1.f;
L50:
	if (kk == *n) {
	    goto L90;
	}
	j1 = k + 1;
	if (lower) {
	    j1 = 1;
	}
	j2 = *n;
	if (lower) {
	    j2 = k - 1;
	}
	i__2 = j2;
	for (j = j1; j <= i__2; ++j) {
	    sm += (r__1 = z__[j] + wkm * t[k + j * t_dim1], dabs(r__1));
	    z__[j] += wk * t[k + j * t_dim1];
	    s += (r__1 = z__[j], dabs(r__1));
/* L60: */
	}
	if (s >= sm) {
	    goto L80;
	}
	w = wkm - wk;
	wk = wkm;
	i__2 = j2;
	for (j = j1; j <= i__2; ++j) {
	    z__[j] += w * t[k + j * t_dim1];
/* L70: */
	}
L80:
L90:
	z__[k] = wk;
/* L100: */
    }
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.f;

/*     solve t*z = y */

    i__1 = *n;
    for (kk = 1; kk <= i__1; ++kk) {
	k = *n + 1 - kk;
	if (lower) {
	    k = kk;
	}
	if ((r__1 = z__[k], dabs(r__1)) <= (r__2 = t[k + k * t_dim1], dabs(
		r__2))) {
	    goto L110;
	}
	s = (r__1 = t[k + k * t_dim1], dabs(r__1)) / (r__2 = z__[k], dabs(
		r__2));
	sscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L110:
	if (t[k + k * t_dim1] != 0.f) {
	    z__[k] /= t[k + k * t_dim1];
	}
	if (t[k + k * t_dim1] == 0.f) {
	    z__[k] = 1.f;
	}
	i1 = 1;
	if (lower) {
	    i1 = k + 1;
	}
	if (kk >= *n) {
	    goto L120;
	}
	w = -z__[k];
	i__2 = *n - kk;
	saxpy_(&i__2, &w, &t[i1 + k * t_dim1], &c__1, &z__[i1], &c__1);
L120:
/* L130: */
	;
    }
/*     make znorm = 1.0 */
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (tnorm != 0.f) {
	*rcond = ynorm / tnorm;
    }
    if (tnorm == 0.f) {
	*rcond = 0.f;
    }
    return 0;
} /* strco_ */

