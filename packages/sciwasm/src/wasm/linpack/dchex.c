/* dchex.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dchex_(doublereal *r__, integer *ldr, integer *p, 
	integer *k, integer *l, doublereal *z__, integer *ldz, integer *nz, 
	doublereal *c__, doublereal *s, integer *job)
{
    /* System generated locals */
    integer r_dim1, r_offset, z_dim1, z_offset, i__1, i__2, i__3;

    /* Local variables */
    integer i__, j;
    doublereal t;
    integer ii, jj, il, iu, km1, lm1, kp1, lmk;
    extern /* Subroutine */ int drotg_(doublereal *, doublereal *, doublereal 
	    *, doublereal *);


/*     dchex updates the cholesky factorization */

/*                   a = trans(r)*r */

/*     of a positive definite matrix a of order p under diagonal */
/*     permutations of the form */

/*                   trans(e)*a*e */

/*     where e is a permutation matrix.  specifically, given */
/*     an upper triangular matrix r and a permutation matrix */
/*     e (which is specified by k, l, and job), dchex determines */
/*     a orthogonal matrix u such that */

/*                           u*r*e = rr, */

/*     where rr is upper triangular.  at the users option, the */
/*     transformation u will be multiplied into the array z. */
/*     if a = trans(x)*x, so that r is the triangular part of the */
/*     qr factorization of x, then rr is the triangular part of the */
/*     qr factorization of x*e, i.e. x with its columns permuted. */
/*     for a less terse description of what dchex does and how */
/*     it may be applied, see the linpack guide. */

/*     the matrix q is determined as the product u(l-k)*...*u(1) */
/*     of plane rotations of the form */

/*                           (    c(i)       s(i) ) */
/*                           (                    ) , */
/*                           (    -s(i)      c(i) ) */

/*     where c(i) is double precision, the rows these rotations operate */
/*     on are described below. */

/*     there are two types of permutations, which are determined */
/*     by the value of job. */

/*     1. right circular shift (job = 1). */

/*         the columns are rearranged in the following order. */

/*                1,...,k-1,l,k,k+1,...,l-1,l+1,...,p. */

/*         u is the product of l-k rotations u(i), where u(i) */
/*         acts in the (l-i,l-i+1)-plane. */

/*     2. left circular shift (job = 2). */
/*         the columns are rearranged in the following order */

/*                1,...,k-1,k+1,k+2,...,l,k,l+1,...,p. */

/*         u is the product of l-k rotations u(i), where u(i) */
/*         acts in the (k+i-1,k+i)-plane. */

/*     on entry */

/*         r      double precision(ldr,p), where ldr.ge.p. */
/*                r contains the upper triangular factor */
/*                that is to be updated.  elements of r */
/*                below the diagonal are not referenced. */

/*         ldr    integer. */
/*                ldr is the leading dimension of the array r. */

/*         p      integer. */
/*                p is the order of the matrix r. */

/*         k      integer. */
/*                k is the first column to be permuted. */

/*         l      integer. */
/*                l is the last column to be permuted. */
/*                l must be strictly greater than k. */

/*         z      double precision(ldz,nz), where ldz.ge.p. */
/*                z is an array of nz p-vectors into which the */
/*                transformation u is multiplied.  z is */
/*                not referenced if nz = 0. */

/*         ldz    integer. */
/*                ldz is the leading dimension of the array z. */

/*         nz     integer. */
/*                nz is the number of columns of the matrix z. */

/*         job    integer. */
/*                job determines the type of permutation. */
/*                       job = 1  right circular shift. */
/*                       job = 2  left circular shift. */

/*     on return */

/*         r      contains the updated factor. */

/*         z      contains the updated matrix z. */

/*         c      double precision(p). */
/*                c contains the cosines of the transforming rotations. */

/*         s      double precision(p). */
/*                s contains the sines of the transforming rotations. */

/*     linpack. this version dated 08/14/78 . */
/*     g.w. stewart, university of maryland, argonne national lab. */

/*     dchex uses the following functions and subroutines. */

/*     blas drotg */
/*     fortran min0 */


/*     initialize */

    /* Parameter adjustments */
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --c__;
    --s;

    /* Function Body */
    km1 = *k - 1;
    kp1 = *k + 1;
    lmk = *l - *k;
    lm1 = *l - 1;

/*     perform the appropriate task. */

    switch (*job) {
	case 1:  goto L10;
	case 2:  goto L130;
    }

/*     right circular shift. */

L10:

/*        reorder the columns. */

    i__1 = *l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = *l - i__ + 1;
	s[i__] = r__[ii + *l * r_dim1];
/* L20: */
    }
    i__1 = lm1;
    for (jj = *k; jj <= i__1; ++jj) {
	j = lm1 - jj + *k;
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    r__[i__ + (j + 1) * r_dim1] = r__[i__ + j * r_dim1];
/* L30: */
	}
	r__[j + 1 + (j + 1) * r_dim1] = 0.;
/* L40: */
    }
    if (*k == 1) {
	goto L60;
    }
    i__1 = km1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = *l - i__ + 1;
	r__[i__ + *k * r_dim1] = s[ii];
/* L50: */
    }
L60:

/*        calculate the rotations. */

    t = s[1];
    i__1 = lmk;
    for (i__ = 1; i__ <= i__1; ++i__) {
	drotg_(&s[i__ + 1], &t, &c__[i__], &s[i__]);
	t = s[i__ + 1];
/* L70: */
    }
    r__[*k + *k * r_dim1] = t;
    i__1 = *p;
    for (j = kp1; j <= i__1; ++j) {
/* Computing MAX */
	i__2 = 1, i__3 = *l - j + 1;
	il = max(i__2,i__3);
	i__2 = lmk;
	for (ii = il; ii <= i__2; ++ii) {
	    i__ = *l - ii;
	    t = c__[ii] * r__[i__ + j * r_dim1] + s[ii] * r__[i__ + 1 + j * 
		    r_dim1];
	    r__[i__ + 1 + j * r_dim1] = c__[ii] * r__[i__ + 1 + j * r_dim1] - 
		    s[ii] * r__[i__ + j * r_dim1];
	    r__[i__ + j * r_dim1] = t;
/* L80: */
	}
/* L90: */
    }

/*        if required, apply the transformations to z. */

    if (*nz < 1) {
	goto L120;
    }
    i__1 = *nz;
    for (j = 1; j <= i__1; ++j) {
	i__2 = lmk;
	for (ii = 1; ii <= i__2; ++ii) {
	    i__ = *l - ii;
	    t = c__[ii] * z__[i__ + j * z_dim1] + s[ii] * z__[i__ + 1 + j * 
		    z_dim1];
	    z__[i__ + 1 + j * z_dim1] = c__[ii] * z__[i__ + 1 + j * z_dim1] - 
		    s[ii] * z__[i__ + j * z_dim1];
	    z__[i__ + j * z_dim1] = t;
/* L100: */
	}
/* L110: */
    }
L120:
    goto L260;

/*     left circular shift */

L130:

/*        reorder the columns */

    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = lmk + i__;
	s[ii] = r__[i__ + *k * r_dim1];
/* L140: */
    }
    i__1 = lm1;
    for (j = *k; j <= i__1; ++j) {
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    r__[i__ + j * r_dim1] = r__[i__ + (j + 1) * r_dim1];
/* L150: */
	}
	jj = j - km1;
	s[jj] = r__[j + 1 + (j + 1) * r_dim1];
/* L160: */
    }
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = lmk + i__;
	r__[i__ + *l * r_dim1] = s[ii];
/* L170: */
    }
    i__1 = *l;
    for (i__ = kp1; i__ <= i__1; ++i__) {
	r__[i__ + *l * r_dim1] = 0.;
/* L180: */
    }

/*        reduction loop. */

    i__1 = *p;
    for (j = *k; j <= i__1; ++j) {
	if (j == *k) {
	    goto L200;
	}

/*              apply the rotations. */

/* Computing MIN */
	i__2 = j - 1, i__3 = *l - 1;
	iu = min(i__2,i__3);
	i__2 = iu;
	for (i__ = *k; i__ <= i__2; ++i__) {
	    ii = i__ - *k + 1;
	    t = c__[ii] * r__[i__ + j * r_dim1] + s[ii] * r__[i__ + 1 + j * 
		    r_dim1];
	    r__[i__ + 1 + j * r_dim1] = c__[ii] * r__[i__ + 1 + j * r_dim1] - 
		    s[ii] * r__[i__ + j * r_dim1];
	    r__[i__ + j * r_dim1] = t;
/* L190: */
	}
L200:
	if (j >= *l) {
	    goto L210;
	}
	jj = j - *k + 1;
	t = s[jj];
	drotg_(&r__[j + j * r_dim1], &t, &c__[jj], &s[jj]);
L210:
/* L220: */
	;
    }

/*        apply the rotations to z. */

    if (*nz < 1) {
	goto L250;
    }
    i__1 = *nz;
    for (j = 1; j <= i__1; ++j) {
	i__2 = lm1;
	for (i__ = *k; i__ <= i__2; ++i__) {
	    ii = i__ - km1;
	    t = c__[ii] * z__[i__ + j * z_dim1] + s[ii] * z__[i__ + 1 + j * 
		    z_dim1];
	    z__[i__ + 1 + j * z_dim1] = c__[ii] * z__[i__ + 1 + j * z_dim1] - 
		    s[ii] * z__[i__ + j * z_dim1];
	    z__[i__ + j * z_dim1] = t;
/* L230: */
	}
/* L240: */
    }
L250:
L260:
    return 0;
} /* dchex_ */

