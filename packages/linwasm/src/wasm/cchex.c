/* cchex.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int cchex_(complex *r__, integer *ldr, integer *p, integer *
	k, integer *l, complex *z__, integer *ldz, integer *nz, real *c__, 
	complex *s, integer *job)
{
    /* System generated locals */
    integer r_dim1, r_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5, 
	    i__6;
    complex q__1, q__2, q__3, q__4;

    /* Builtin functions */
    void r_cnjg(complex *, complex *);

    /* Local variables */
    integer i__, j;
    complex t;
    integer ii, jj, il, iu, km1, lm1, kp1, lmk;
    extern /* Subroutine */ int crotg_(complex *, complex *, real *, complex *
	    );


/*     cchex updates the cholesky factorization */

/*                   a = ctrans(r)*r */

/*     of a positive definite matrix a of order p under diagonal */
/*     permutations of the form */

/*                   trans(e)*a*e */

/*     where e is a permutation matrix.  specifically, given */
/*     an upper triangular matrix r and a permutation matrix */
/*     e (which is specified by k, l, and job), cchex determines */
/*     a unitary matrix u such that */

/*                           u*r*e = rr, */

/*     where rr is upper triangular.  at the users option, the */
/*     transformation u will be multiplied into the array z. */
/*     if a = ctrans(x)*x, so that r is the triangular part of the */
/*     qr factorization of x, then rr is the triangular part of the */
/*     qr factorization of x*e, i.e. x with its columns permuted. */
/*     for a less terse description of what cchex does and how */
/*     it may be applied, see the linpack guide. */

/*     the matrix q is determined as the product u(l-k)*...*u(1) */
/*     of plane rotations of the form */

/*                           (    c(i)       s(i) ) */
/*                           (                    ) , */
/*                           ( -conjg(s(i))  c(i) ) */

/*     where c(i) is real, the rows these rotations operate on */
/*     are described below. */

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

/*         r      complex(ldr,p), where ldr.ge.p. */
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

/*         z      complex(ldz,nz), where ldz.ge.p. */
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

/*         c      real(p). */
/*                c contains the cosines of the transforming rotations. */

/*         s      complex(p). */
/*                s contains the sines of the transforming rotations. */

/*     linpack. this version dated 08/14/78 . */
/*     g.w. stewart, university of maryland, argonne national lab. */

/*     cchex uses the following functions and subroutines. */

/*     extended blas crotg */
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
	i__2 = i__;
	i__3 = ii + *l * r_dim1;
	s[i__2].r = r__[i__3].r, s[i__2].i = r__[i__3].i;
/* L20: */
    }
    i__1 = lm1;
    for (jj = *k; jj <= i__1; ++jj) {
	j = lm1 - jj + *k;
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + (j + 1) * r_dim1;
	    i__4 = i__ + j * r_dim1;
	    r__[i__3].r = r__[i__4].r, r__[i__3].i = r__[i__4].i;
/* L30: */
	}
	i__2 = j + 1 + (j + 1) * r_dim1;
	r__[i__2].r = 0.f, r__[i__2].i = 0.f;
/* L40: */
    }
    if (*k == 1) {
	goto L60;
    }
    i__1 = km1;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = *l - i__ + 1;
	i__2 = i__ + *k * r_dim1;
	i__3 = ii;
	r__[i__2].r = s[i__3].r, r__[i__2].i = s[i__3].i;
/* L50: */
    }
L60:

/*        calculate the rotations. */

    t.r = s[1].r, t.i = s[1].i;
    i__1 = lmk;
    for (i__ = 1; i__ <= i__1; ++i__) {
	crotg_(&s[i__ + 1], &t, &c__[i__], &s[i__]);
	i__2 = i__ + 1;
	t.r = s[i__2].r, t.i = s[i__2].i;
/* L70: */
    }
    i__1 = *k + *k * r_dim1;
    r__[i__1].r = t.r, r__[i__1].i = t.i;
    i__1 = *p;
    for (j = kp1; j <= i__1; ++j) {
/* Computing MAX */
	i__2 = 1, i__3 = *l - j + 1;
	il = max(i__2,i__3);
	i__2 = lmk;
	for (ii = il; ii <= i__2; ++ii) {
	    i__ = *l - ii;
	    i__3 = ii;
	    i__4 = i__ + j * r_dim1;
	    q__2.r = c__[i__3] * r__[i__4].r, q__2.i = c__[i__3] * r__[i__4]
		    .i;
	    i__5 = ii;
	    i__6 = i__ + 1 + j * r_dim1;
	    q__3.r = s[i__5].r * r__[i__6].r - s[i__5].i * r__[i__6].i, 
		    q__3.i = s[i__5].r * r__[i__6].i + s[i__5].i * r__[i__6]
		    .r;
	    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
	    t.r = q__1.r, t.i = q__1.i;
	    i__3 = i__ + 1 + j * r_dim1;
	    i__4 = ii;
	    i__5 = i__ + 1 + j * r_dim1;
	    q__2.r = c__[i__4] * r__[i__5].r, q__2.i = c__[i__4] * r__[i__5]
		    .i;
	    r_cnjg(&q__4, &s[ii]);
	    i__6 = i__ + j * r_dim1;
	    q__3.r = q__4.r * r__[i__6].r - q__4.i * r__[i__6].i, q__3.i = 
		    q__4.r * r__[i__6].i + q__4.i * r__[i__6].r;
	    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
	    r__[i__3].r = q__1.r, r__[i__3].i = q__1.i;
	    i__3 = i__ + j * r_dim1;
	    r__[i__3].r = t.r, r__[i__3].i = t.i;
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
	    i__3 = ii;
	    i__4 = i__ + j * z_dim1;
	    q__2.r = c__[i__3] * z__[i__4].r, q__2.i = c__[i__3] * z__[i__4]
		    .i;
	    i__5 = ii;
	    i__6 = i__ + 1 + j * z_dim1;
	    q__3.r = s[i__5].r * z__[i__6].r - s[i__5].i * z__[i__6].i, 
		    q__3.i = s[i__5].r * z__[i__6].i + s[i__5].i * z__[i__6]
		    .r;
	    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
	    t.r = q__1.r, t.i = q__1.i;
	    i__3 = i__ + 1 + j * z_dim1;
	    i__4 = ii;
	    i__5 = i__ + 1 + j * z_dim1;
	    q__2.r = c__[i__4] * z__[i__5].r, q__2.i = c__[i__4] * z__[i__5]
		    .i;
	    r_cnjg(&q__4, &s[ii]);
	    i__6 = i__ + j * z_dim1;
	    q__3.r = q__4.r * z__[i__6].r - q__4.i * z__[i__6].i, q__3.i = 
		    q__4.r * z__[i__6].i + q__4.i * z__[i__6].r;
	    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
	    z__[i__3].r = q__1.r, z__[i__3].i = q__1.i;
	    i__3 = i__ + j * z_dim1;
	    z__[i__3].r = t.r, z__[i__3].i = t.i;
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
	i__2 = ii;
	i__3 = i__ + *k * r_dim1;
	s[i__2].r = r__[i__3].r, s[i__2].i = r__[i__3].i;
/* L140: */
    }
    i__1 = lm1;
    for (j = *k; j <= i__1; ++j) {
	i__2 = j;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + j * r_dim1;
	    i__4 = i__ + (j + 1) * r_dim1;
	    r__[i__3].r = r__[i__4].r, r__[i__3].i = r__[i__4].i;
/* L150: */
	}
	jj = j - km1;
	i__2 = jj;
	i__3 = j + 1 + (j + 1) * r_dim1;
	s[i__2].r = r__[i__3].r, s[i__2].i = r__[i__3].i;
/* L160: */
    }
    i__1 = *k;
    for (i__ = 1; i__ <= i__1; ++i__) {
	ii = lmk + i__;
	i__2 = i__ + *l * r_dim1;
	i__3 = ii;
	r__[i__2].r = s[i__3].r, r__[i__2].i = s[i__3].i;
/* L170: */
    }
    i__1 = *l;
    for (i__ = kp1; i__ <= i__1; ++i__) {
	i__2 = i__ + *l * r_dim1;
	r__[i__2].r = 0.f, r__[i__2].i = 0.f;
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
	    i__3 = ii;
	    i__4 = i__ + j * r_dim1;
	    q__2.r = c__[i__3] * r__[i__4].r, q__2.i = c__[i__3] * r__[i__4]
		    .i;
	    i__5 = ii;
	    i__6 = i__ + 1 + j * r_dim1;
	    q__3.r = s[i__5].r * r__[i__6].r - s[i__5].i * r__[i__6].i, 
		    q__3.i = s[i__5].r * r__[i__6].i + s[i__5].i * r__[i__6]
		    .r;
	    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
	    t.r = q__1.r, t.i = q__1.i;
	    i__3 = i__ + 1 + j * r_dim1;
	    i__4 = ii;
	    i__5 = i__ + 1 + j * r_dim1;
	    q__2.r = c__[i__4] * r__[i__5].r, q__2.i = c__[i__4] * r__[i__5]
		    .i;
	    r_cnjg(&q__4, &s[ii]);
	    i__6 = i__ + j * r_dim1;
	    q__3.r = q__4.r * r__[i__6].r - q__4.i * r__[i__6].i, q__3.i = 
		    q__4.r * r__[i__6].i + q__4.i * r__[i__6].r;
	    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
	    r__[i__3].r = q__1.r, r__[i__3].i = q__1.i;
	    i__3 = i__ + j * r_dim1;
	    r__[i__3].r = t.r, r__[i__3].i = t.i;
/* L190: */
	}
L200:
	if (j >= *l) {
	    goto L210;
	}
	jj = j - *k + 1;
	i__2 = jj;
	t.r = s[i__2].r, t.i = s[i__2].i;
	crotg_(&r__[j + j * r_dim1], &t, &c__[jj], &s[jj]);
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
	    i__3 = ii;
	    i__4 = i__ + j * z_dim1;
	    q__2.r = c__[i__3] * z__[i__4].r, q__2.i = c__[i__3] * z__[i__4]
		    .i;
	    i__5 = ii;
	    i__6 = i__ + 1 + j * z_dim1;
	    q__3.r = s[i__5].r * z__[i__6].r - s[i__5].i * z__[i__6].i, 
		    q__3.i = s[i__5].r * z__[i__6].i + s[i__5].i * z__[i__6]
		    .r;
	    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
	    t.r = q__1.r, t.i = q__1.i;
	    i__3 = i__ + 1 + j * z_dim1;
	    i__4 = ii;
	    i__5 = i__ + 1 + j * z_dim1;
	    q__2.r = c__[i__4] * z__[i__5].r, q__2.i = c__[i__4] * z__[i__5]
		    .i;
	    r_cnjg(&q__4, &s[ii]);
	    i__6 = i__ + j * z_dim1;
	    q__3.r = q__4.r * z__[i__6].r - q__4.i * z__[i__6].i, q__3.i = 
		    q__4.r * z__[i__6].i + q__4.i * z__[i__6].r;
	    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
	    z__[i__3].r = q__1.r, z__[i__3].i = q__1.i;
	    i__3 = i__ + j * z_dim1;
	    z__[i__3].r = t.r, z__[i__3].i = t.i;
/* L230: */
	}
/* L240: */
    }
L250:
L260:
    return 0;
} /* cchex_ */

