/* cchud.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int cchud_(complex *r__, integer *ldr, integer *p, complex *
	x, complex *z__, integer *ldz, integer *nz, complex *y, real *rho, 
	real *c__, complex *s)
{
    /* System generated locals */
    integer r_dim1, r_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5;
    real r__1, r__2;
    complex q__1, q__2, q__3, q__4;

    /* Builtin functions */
    void r_cnjg(complex *, complex *);
    double c_abs(complex *), sqrt(doublereal);

    /* Local variables */
    integer i__, j;
    complex t, xj;
    integer jm1;
    complex zeta;
    real scale, azeta;
    extern /* Subroutine */ int crotg_(complex *, complex *, real *, complex *
	    );


/*     cchud updates an augmented cholesky decomposition of the */
/*     triangular part of an augmented qr decomposition.  specifically, */
/*     given an upper triangular matrix r of order p, a row vector */
/*     x, a column vector z, and a scalar y, cchud determines a */
/*     untiary matrix u and a scalar zeta such that */


/*                              (r  z)     (rr   zz ) */
/*                         u  * (    )  =  (        ) , */
/*                              (x  y)     ( 0  zeta) */

/*     where rr is upper triangular.  if r and z have been */
/*     obtained from the factorization of a least squares */
/*     problem, then rr and zz are the factors corresponding to */
/*     the problem with the observation (x,y) appended.  in this */
/*     case, if rho is the norm of the residual vector, then the */
/*     norm of the residual vector of the updated problem is */
/*     sqrt(rho**2 + zeta**2).  cchud will simultaneously update */
/*     several triplets (z,y,rho). */
/*     for a less terse description of what cchud does and how */
/*     it may be applied see the linpack guide. */

/*     the matrix u is determined as the product u(p)*...*u(1), */
/*     where u(i) is a rotation in the (i,p+1) plane of the */
/*     form */

/*                       (     c(i)      s(i) ) */
/*                       (                    ) . */
/*                       ( -conjg(s(i))  c(i) ) */

/*     the rotations are chosen so that c(i) is real. */

/*     on entry */

/*         r      complex(ldr,p), where ldr .ge. p. */
/*                r contains the upper triangular matrix */
/*                that is to be updated.  the part of r */
/*                below the diagonal is not referenced. */

/*         ldr    integer. */
/*                ldr is the leading dimension of the array r. */

/*         p      integer. */
/*                p is the order of the matrix r. */

/*         x      complex(p). */
/*                x contains the row to be added to r.  x is */
/*                not altered by cchud. */

/*         z      complex(ldz,nz), where ldz .ge. p. */
/*                z is an array containing nz p-vectors to */
/*                be updated with r. */

/*         ldz    integer. */
/*                ldz is the leading dimension of the array z. */

/*         nz     integer. */
/*                nz is the number of vectors to be updated */
/*                nz may be zero, in which case z, y, and rho */
/*                are not referenced. */

/*         y      complex(nz). */
/*                y contains the scalars for updating the vectors */
/*                z.  y is not altered by cchud. */

/*         rho    real(nz). */
/*                rho contains the norms of the residual */
/*                vectors that are to be updated.  if rho(j) */
/*                is negative, it is left unaltered. */

/*     on return */

/*         rc */
/*         rho    contain the updated quantities. */
/*         z */

/*         c      real(p). */
/*                c contains the cosines of the transforming */
/*                rotations. */

/*         s      complex(p). */
/*                s contains the sines of the transforming */
/*                rotations. */

/*     linpack. this version dated 08/14/78 . */
/*     g.w. stewart, university of maryland, argonne national lab. */

/*     cchud uses the following functions and subroutines. */

/*     extended blas crotg */
/*     fortran conjg,sqrt */


/*     update r. */

    /* Parameter adjustments */
    r_dim1 = *ldr;
    r_offset = 1 + r_dim1;
    r__ -= r_offset;
    --x;
    z_dim1 = *ldz;
    z_offset = 1 + z_dim1;
    z__ -= z_offset;
    --y;
    --rho;
    --c__;
    --s;

    /* Function Body */
    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	xj.r = x[i__2].r, xj.i = x[i__2].i;

/*        apply the previous rotations. */

	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    i__4 = i__ + j * r_dim1;
	    q__2.r = c__[i__3] * r__[i__4].r, q__2.i = c__[i__3] * r__[i__4]
		    .i;
	    i__5 = i__;
	    q__3.r = s[i__5].r * xj.r - s[i__5].i * xj.i, q__3.i = s[i__5].r *
		     xj.i + s[i__5].i * xj.r;
	    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
	    t.r = q__1.r, t.i = q__1.i;
	    i__3 = i__;
	    q__2.r = c__[i__3] * xj.r, q__2.i = c__[i__3] * xj.i;
	    r_cnjg(&q__4, &s[i__]);
	    i__4 = i__ + j * r_dim1;
	    q__3.r = q__4.r * r__[i__4].r - q__4.i * r__[i__4].i, q__3.i = 
		    q__4.r * r__[i__4].i + q__4.i * r__[i__4].r;
	    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
	    xj.r = q__1.r, xj.i = q__1.i;
	    i__3 = i__ + j * r_dim1;
	    r__[i__3].r = t.r, r__[i__3].i = t.i;
/* L10: */
	}
L20:

/*        compute the next rotation. */

	crotg_(&r__[j + j * r_dim1], &xj, &c__[j], &s[j]);
/* L30: */
    }

/*     if required, update z and rho. */

    if (*nz < 1) {
	goto L70;
    }
    i__1 = *nz;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	zeta.r = y[i__2].r, zeta.i = y[i__2].i;
	i__2 = *p;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    i__4 = i__ + j * z_dim1;
	    q__2.r = c__[i__3] * z__[i__4].r, q__2.i = c__[i__3] * z__[i__4]
		    .i;
	    i__5 = i__;
	    q__3.r = s[i__5].r * zeta.r - s[i__5].i * zeta.i, q__3.i = s[i__5]
		    .r * zeta.i + s[i__5].i * zeta.r;
	    q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
	    t.r = q__1.r, t.i = q__1.i;
	    i__3 = i__;
	    q__2.r = c__[i__3] * zeta.r, q__2.i = c__[i__3] * zeta.i;
	    r_cnjg(&q__4, &s[i__]);
	    i__4 = i__ + j * z_dim1;
	    q__3.r = q__4.r * z__[i__4].r - q__4.i * z__[i__4].i, q__3.i = 
		    q__4.r * z__[i__4].i + q__4.i * z__[i__4].r;
	    q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
	    zeta.r = q__1.r, zeta.i = q__1.i;
	    i__3 = i__ + j * z_dim1;
	    z__[i__3].r = t.r, z__[i__3].i = t.i;
/* L40: */
	}
	azeta = c_abs(&zeta);
	if (azeta == 0.f || rho[j] < 0.f) {
	    goto L50;
	}
	scale = azeta + rho[j];
/* Computing 2nd power */
	r__1 = azeta / scale;
/* Computing 2nd power */
	r__2 = rho[j] / scale;
	rho[j] = scale * sqrt(r__1 * r__1 + r__2 * r__2);
L50:
/* L60: */
	;
    }
L70:
    return 0;
} /* cchud_ */

