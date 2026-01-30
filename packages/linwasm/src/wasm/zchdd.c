/* zchdd.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int zchdd_(doublecomplex *r__, integer *ldr, integer *p, 
	doublecomplex *x, doublecomplex *z__, integer *ldz, integer *nz, 
	doublecomplex *y, doublereal *rho, doublereal *c__, doublecomplex *s, 
	integer *info)
{
    /* System generated locals */
    integer r_dim1, r_offset, z_dim1, z_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3;
    doublecomplex z__1, z__2, z__3, z__4;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *), z_div(doublecomplex *, 
	    doublecomplex *, doublecomplex *);
    double sqrt(doublereal), z_abs(doublecomplex *);

    /* Local variables */
    doublereal a;
    doublecomplex b;
    integer i__, j;
    doublecomplex t;
    integer ii;
    doublecomplex xx, zeta;
    doublereal norm, alpha;
    real scale;
    doublereal azeta;
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern doublereal dznrm2_(integer *, doublecomplex *, integer *);


/*     zchdd downdates an augmented cholesky decomposition or the */
/*     triangular factor of an augmented qr decomposition. */
/*     specifically, given an upper triangular matrix r of order p,  a */
/*     row vector x, a column vector z, and a scalar y, zchdd */
/*     determineds a unitary matrix u and a scalar zeta such that */

/*                        (r   z )     (rr  zz) */
/*                    u * (      )  =  (      ) , */
/*                        (0 zeta)     ( x   y) */

/*     where rr is upper triangular.  if r and z have been obtained */
/*     from the factorization of a least squares problem, then */
/*     rr and zz are the factors corresponding to the problem */
/*     with the observation (x,y) removed.  in this case, if rho */
/*     is the norm of the residual vector, then the norm of */
/*     the residual vector of the downdated problem is */
/*     dsqrt(rho**2 - zeta**2). zchdd will simultaneously downdate */
/*     several triplets (z,y,rho) along with r. */
/*     for a less terse description of what zchdd does and how */
/*     it may be applied, see the linpack guide. */

/*     the matrix u is determined as the product u(1)*...*u(p) */
/*     where u(i) is a rotation in the (p+1,i)-plane of the */
/*     form */

/*                       ( c(i)  -dconjg(s(i)) ) */
/*                       (                    ) . */
/*                       ( s(i)       c(i)    ) */

/*     the rotations are chosen so that c(i) is double precision. */

/*     the user is warned that a given downdating problem may */
/*     be impossible to accomplish or may produce */
/*     inaccurate results.  for example, this can happen */
/*     if x is near a vector whose removal will reduce the */
/*     rank of r.  beware. */

/*     on entry */

/*         r      complex*16(ldr,p), where ldr .ge. p. */
/*                r contains the upper triangular matrix */
/*                that is to be downdated.  the part of  r */
/*                below the diagonal is not referenced. */

/*         ldr    integer. */
/*                ldr is the leading dimension fo the array r. */

/*         p      integer. */
/*                p is the order of the matrix r. */

/*         x      complex*16(p). */
/*                x contains the row vector that is to */
/*                be removed from r.  x is not altered by zchdd. */

/*         z      complex*16(ldz,nz), where ldz .ge. p. */
/*                z is an array of nz p-vectors which */
/*                are to be downdated along with r. */

/*         ldz    integer. */
/*                ldz is the leading dimension of the array z. */

/*         nz     integer. */
/*                nz is the number of vectors to be downdated */
/*                nz may be zero, in which case z, y, and rho */
/*                are not referenced. */

/*         y      complex*16(nz). */
/*                y contains the scalars for the downdating */
/*                of the vectors z.  y is not altered by zchdd. */

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

/*         s      complex*16(p). */
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

/*     zchdd uses the following functions and subprograms. */

/*     fortran cdabs,dconjg */
/*     blas zdotc, dznrm2 */


/*     solve the system ctrans(r)*a = x, placing the result */
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
    d_cnjg(&z__2, &x[1]);
    d_cnjg(&z__3, &r__[r_dim1 + 1]);
    z_div(&z__1, &z__2, &z__3);
    s[1].r = z__1.r, s[1].i = z__1.i;
    if (*p < 2) {
	goto L20;
    }
    i__1 = *p;
    for (j = 2; j <= i__1; ++j) {
	i__2 = j;
	d_cnjg(&z__2, &x[j]);
	i__3 = j - 1;
	zdotc_(&z__3, &i__3, &r__[j * r_dim1 + 1], &c__1, &s[1], &c__1);
	z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	s[i__2].r = z__1.r, s[i__2].i = z__1.i;
	i__2 = j;
	d_cnjg(&z__2, &r__[j + j * r_dim1]);
	z_div(&z__1, &s[j], &z__2);
	s[i__2].r = z__1.r, s[i__2].i = z__1.i;
/* L10: */
    }
L20:
    norm = dznrm2_(p, &s[1], &c__1);
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
	scale = alpha + z_abs(&s[i__]);
	a = alpha / scale;
	i__2 = i__;
	z__1.r = s[i__2].r / scale, z__1.i = s[i__2].i / scale;
	b.r = z__1.r, b.i = z__1.i;
/* Computing 2nd power */
	d__1 = a;
/* Computing 2nd power */
	d__2 = b.r;
	z__1.r = b.r * 0. - b.i * -1., z__1.i = b.i * 0. + b.r * -1.;
/* Computing 2nd power */
	d__3 = z__1.r;
	norm = sqrt(d__1 * d__1 + d__2 * d__2 + d__3 * d__3);
	c__[i__] = a / norm;
	i__2 = i__;
	d_cnjg(&z__2, &b);
	z__1.r = z__2.r / norm, z__1.i = z__2.i / norm;
	s[i__2].r = z__1.r, s[i__2].i = z__1.i;
	alpha = scale * norm;
/* L40: */
    }

/*        apply the transformations to r. */

    i__1 = *p;
    for (j = 1; j <= i__1; ++j) {
	xx.r = 0., xx.i = 0.;
	i__2 = j;
	for (ii = 1; ii <= i__2; ++ii) {
	    i__ = j - ii + 1;
	    i__3 = i__;
	    z__2.r = c__[i__3] * xx.r, z__2.i = c__[i__3] * xx.i;
	    i__4 = i__;
	    i__5 = i__ + j * r_dim1;
	    z__3.r = s[i__4].r * r__[i__5].r - s[i__4].i * r__[i__5].i, 
		    z__3.i = s[i__4].r * r__[i__5].i + s[i__4].i * r__[i__5]
		    .r;
	    z__1.r = z__2.r + z__3.r, z__1.i = z__2.i + z__3.i;
	    t.r = z__1.r, t.i = z__1.i;
	    i__3 = i__ + j * r_dim1;
	    i__4 = i__;
	    i__5 = i__ + j * r_dim1;
	    z__2.r = c__[i__4] * r__[i__5].r, z__2.i = c__[i__4] * r__[i__5]
		    .i;
	    d_cnjg(&z__4, &s[i__]);
	    z__3.r = z__4.r * xx.r - z__4.i * xx.i, z__3.i = z__4.r * xx.i + 
		    z__4.i * xx.r;
	    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	    r__[i__3].r = z__1.r, r__[i__3].i = z__1.i;
	    xx.r = t.r, xx.i = t.i;
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
	i__2 = j;
	zeta.r = y[i__2].r, zeta.i = y[i__2].i;
	i__2 = *p;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + j * z_dim1;
	    i__4 = i__ + j * z_dim1;
	    d_cnjg(&z__4, &s[i__]);
	    z__3.r = z__4.r * zeta.r - z__4.i * zeta.i, z__3.i = z__4.r * 
		    zeta.i + z__4.i * zeta.r;
	    z__2.r = z__[i__4].r - z__3.r, z__2.i = z__[i__4].i - z__3.i;
	    i__5 = i__;
	    z__1.r = z__2.r / c__[i__5], z__1.i = z__2.i / c__[i__5];
	    z__[i__3].r = z__1.r, z__[i__3].i = z__1.i;
	    i__3 = i__;
	    z__2.r = c__[i__3] * zeta.r, z__2.i = c__[i__3] * zeta.i;
	    i__4 = i__;
	    i__5 = i__ + j * z_dim1;
	    z__3.r = s[i__4].r * z__[i__5].r - s[i__4].i * z__[i__5].i, 
		    z__3.i = s[i__4].r * z__[i__5].i + s[i__4].i * z__[i__5]
		    .r;
	    z__1.r = z__2.r - z__3.r, z__1.i = z__2.i - z__3.i;
	    zeta.r = z__1.r, zeta.i = z__1.i;
/* L70: */
	}
	azeta = z_abs(&zeta);
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
} /* zchdd_ */

