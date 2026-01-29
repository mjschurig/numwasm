/* dchdd.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dchdd_(doublereal *r__, integer *ldr, integer *p, 
	doublereal *x, doublereal *z__, integer *ldz, integer *nz, doublereal 
	*y, doublereal *rho, doublereal *c__, doublereal *s, integer *info)
{
    /* System generated locals */
    integer r_dim1, r_offset, z_dim1, z_offset, i__1, i__2;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    doublereal a, b;
    integer i__, j;
    doublereal t;
    integer ii;
    doublereal xx;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    doublereal zeta, norm;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    doublereal alpha;
    real scale;
    doublereal azeta;


/*     dchdd downdates an augmented cholesky decomposition or the */
/*     triangular factor of an augmented qr decomposition. */
/*     specifically, given an upper triangular matrix r of order p,  a */
/*     row vector x, a column vector z, and a scalar y, dchdd */
/*     determineds a orthogonal matrix u and a scalar zeta such that */

/*                        (r   z )     (rr  zz) */
/*                    u * (      )  =  (      ) , */
/*                        (0 zeta)     ( x   y) */

/*     where rr is upper triangular.  if r and z have been obtained */
/*     from the factorization of a least squares problem, then */
/*     rr and zz are the factors corresponding to the problem */
/*     with the observation (x,y) removed.  in this case, if rho */
/*     is the norm of the residual vector, then the norm of */
/*     the residual vector of the downdated problem is */
/*     dsqrt(rho**2 - zeta**2). dchdd will simultaneously downdate */
/*     several triplets (z,y,rho) along with r. */
/*     for a less terse description of what dchdd does and how */
/*     it may be applied, see the linpack guide. */

/*     the matrix u is determined as the product u(1)*...*u(p) */
/*     where u(i) is a rotation in the (p+1,i)-plane of the */
/*     form */

/*                       ( c(i)     -s(i)     ) */
/*                       (                    ) . */
/*                       ( s(i)       c(i)    ) */

/*     the rotations are chosen so that c(i) is double precision. */

/*     the user is warned that a given downdating problem may */
/*     be impossible to accomplish or may produce */
/*     inaccurate results.  for example, this can happen */
/*     if x is near a vector whose removal will reduce the */
/*     rank of r.  beware. */

/*     on entry */

/*         r      double precision(ldr,p), where ldr .ge. p. */
/*                r contains the upper triangular matrix */
/*                that is to be downdated.  the part of  r */
/*                below the diagonal is not referenced. */

/*         ldr    integer. */
/*                ldr is the leading dimension fo the array r. */

/*         p      integer. */
/*                p is the order of the matrix r. */

/*         x      double precision(p). */
/*                x contains the row vector that is to */
/*                be removed from r.  x is not altered by dchdd. */

/*         z      double precision(ldz,nz), where ldz .ge. p. */
/*                z is an array of nz p-vectors which */
/*                are to be downdated along with r. */

/*         ldz    integer. */
/*                ldz is the leading dimension of the array z. */

/*         nz     integer. */
/*                nz is the number of vectors to be downdated */
/*                nz may be zero, in which case z, y, and rho */
/*                are not referenced. */

/*         y      double precision(nz). */
/*                y contains the scalars for the downdating */
/*                of the vectors z.  y is not altered by dchdd. */

/*         rho    double precision(nz). */
/*                rho contains the norms of the residual */
/*                vectors that are to be downdated. */

/*     on return */

/*         r */
/*         z      contain the downdated quantities. */
/*         rho */

/*         c      double precision(p). */
/*                c contains the cosines of the transforming */
/*                rotations. */

/*         s      double precision(p). */
/*                s contains the sines of the transforming */
/*                rotations. */

/*         info   integer. */
/*                info is set as follows. */

/*                   info = 0  if the entire downdating */
/*                             was successful. */

/*                   info =-1  if r could not be downdated. */
/*                             in this case, all quantities */
/*                             are left unaltered. */

/*                   info = 1  if some rho could not be */
/*                             downdated.  the offending rhos are */
/*                             set to -1. */

/*     linpack. this version dated 08/14/78 . */
/*     g.w. stewart, university of maryland, argonne national lab. */

/*     dchdd uses the following functions and subprograms. */

/*     fortran dabs */
/*     blas ddot, dnrm2 */


/*     solve the system trans(r)*a = x, placing the result */
/*     in the array s. */

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
    *info = 0;
    s[1] = x[1] / r__[r_dim1 + 1];
    if (*p < 2) {
	goto L20;
    }
    i__1 = *p;
    for (j = 2; j <= i__1; ++j) {
	i__2 = j - 1;
	s[j] = x[j] - ddot_(&i__2, &r__[j * r_dim1 + 1], &c__1, &s[1], &c__1);
	s[j] /= r__[j + j * r_dim1];
/* L10: */
    }
L20:
    norm = dnrm2_(p, &s[1], &c__1);
    if (norm < 1.) {
	goto L30;
    }
    *info = -1;
    goto L120;
L30:
/* Computing 2nd power */
    d__1 = norm;
    alpha = sqrt(1. - d__1 * d__1);

/*        determine the transformations. */

    i__1 = *p;
    for (ii = 1; ii <= i__1; ++ii) {
	i__ = *p - ii + 1;
	scale = alpha + (d__1 = s[i__], abs(d__1));
	a = alpha / scale;
	b = s[i__] / scale;
/* Computing 2nd power */
	d__1 = a;
/* Computing 2nd power */
	d__2 = b;
	norm = sqrt(d__1 * d__1 + d__2 * d__2 + 0.);
	c__[i__] = a / norm;
	s[i__] = b / norm;
	alpha = scale * norm;
/* L40: */
    }

/*        apply the transformations to r. */

    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
	xx = 0.;
	i__2 = j;
	for (ii = 1; ii <= i__2; ++ii) {
	    i__ = j - ii + 1;
	    t = c__[i__] * xx + s[i__] * r__[i__ + j * r_dim1];
	    r__[i__ + j * r_dim1] = c__[i__] * r__[i__ + j * r_dim1] - s[i__] 
		    * xx;
	    xx = t;
/* L50: */
	}
/* L60: */
    }

/*        if required, downdate z and rho. */

    if (*nz < 1) {
	goto L110;
    }
    i__1 = *nz;
    for (j = 1; j <= i__1; ++j) {
	zeta = y[j];
	i__2 = *p;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    z__[i__ + j * z_dim1] = (z__[i__ + j * z_dim1] - s[i__] * zeta) / 
		    c__[i__];
	    zeta = c__[i__] * zeta - s[i__] * z__[i__ + j * z_dim1];
/* L70: */
	}
	azeta = abs(zeta);
	if (azeta <= rho[j]) {
	    goto L80;
	}
	*info = 1;
	rho[j] = -1.;
	goto L90;
L80:
/* Computing 2nd power */
	d__1 = azeta / rho[j];
	rho[j] *= sqrt(1. - d__1 * d__1);
L90:
/* L100: */
	;
    }
L110:
L120:
    return 0;
} /* dchdd_ */

