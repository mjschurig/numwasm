/* dchud.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dchud_(doublereal *r__, integer *ldr, integer *p, 
	doublereal *x, doublereal *z__, integer *ldz, integer *nz, doublereal 
	*y, doublereal *rho, doublereal *c__, doublereal *s)
{
    /* System generated locals */
    integer r_dim1, r_offset, z_dim1, z_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer i__, j;
    doublereal t, xj;
    integer jm1;
    doublereal zeta, scale, azeta;
    extern /* Subroutine */ int drotg_(doublereal *, doublereal *, doublereal 
	    *, doublereal *);


/*     dchud updates an augmented cholesky decomposition of the */
/*     triangular part of an augmented qr decomposition.  specifically, */
/*     given an upper triangular matrix r of order p, a row vector */
/*     x, a column vector z, and a scalar y, dchud determines a */
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
/*     dsqrt(rho**2 + zeta**2).  dchud will simultaneously update */
/*     several triplets (z,y,rho). */
/*     for a less terse description of what dchud does and how */
/*     it may be applied, see the linpack guide. */

/*     the matrix u is determined as the product u(p)*...*u(1), */
/*     where u(i) is a rotation in the (i,p+1) plane of the */
/*     form */

/*                       (     c(i)      s(i) ) */
/*                       (                    ) . */
/*                       (    -s(i)      c(i) ) */

/*     the rotations are chosen so that c(i) is double precision. */

/*     on entry */

/*         r      double precision(ldr,p), where ldr .ge. p. */
/*                r contains the upper triangular matrix */
/*                that is to be updated.  the part of r */
/*                below the diagonal is not referenced. */

/*         ldr    integer. */
/*                ldr is the leading dimension of the array r. */

/*         p      integer. */
/*                p is the order of the matrix r. */

/*         x      double precision(p). */
/*                x contains the row to be added to r.  x is */
/*                not altered by dchud. */

/*         z      double precision(ldz,nz), where ldz .ge. p. */
/*                z is an array containing nz p-vectors to */
/*                be updated with r. */

/*         ldz    integer. */
/*                ldz is the leading dimension of the array z. */

/*         nz     integer. */
/*                nz is the number of vectors to be updated */
/*                nz may be zero, in which case z, y, and rho */
/*                are not referenced. */

/*         y      double precision(nz). */
/*                y contains the scalars for updating the vectors */
/*                z.  y is not altered by dchud. */

/*         rho    double precision(nz). */
/*                rho contains the norms of the residual */
/*                vectors that are to be updated.  if rho(j) */
/*                is negative, it is left unaltered. */

/*     on return */

/*         rc */
/*         rho    contain the updated quantities. */
/*         z */

/*         c      double precision(p). */
/*                c contains the cosines of the transforming */
/*                rotations. */

/*         s      double precision(p). */
/*                s contains the sines of the transforming */
/*                rotations. */

/*     linpack. this version dated 08/14/78 . */
/*     g.w. stewart, university of maryland, argonne national lab. */

/*     dchud uses the following functions and subroutines. */

/*     extended blas drotg */
/*     fortran dsqrt */


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
	xj = x[j];

/*        apply the previous rotations. */

	jm1 = j - 1;
	if (jm1 < 1) {
	    goto L20;
	}
	i__2 = jm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    t = c__[i__] * r__[i__ + j * r_dim1] + s[i__] * xj;
	    xj = c__[i__] * xj - s[i__] * r__[i__ + j * r_dim1];
	    r__[i__ + j * r_dim1] = t;
/* L10: */
	}
L20:

/*        compute the next rotation. */

	drotg_(&r__[j + j * r_dim1], &xj, &c__[j], &s[j]);
/* L30: */
    }

/*     if required, update z and rho. */

    if (*nz < 1) {
	goto L70;
    }
    i__1 = *nz;
    for (j = 1; j <= i__1; ++j) {
	zeta = y[j];
	i__2 = *p;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    t = c__[i__] * z__[i__ + j * z_dim1] + s[i__] * zeta;
	    zeta = c__[i__] * zeta - s[i__] * z__[i__ + j * z_dim1];
	    z__[i__ + j * z_dim1] = t;
/* L40: */
	}
	azeta = abs(zeta);
	if (azeta == 0. || rho[j] < 0.) {
	    goto L50;
	}
	scale = azeta + rho[j];
/* Computing 2nd power */
	d__1 = azeta / scale;
/* Computing 2nd power */
	d__2 = rho[j] / scale;
	rho[j] = scale * sqrt(d__1 * d__1 + d__2 * d__2);
L50:
/* L60: */
	;
    }
L70:
    return 0;
} /* dchud_ */

