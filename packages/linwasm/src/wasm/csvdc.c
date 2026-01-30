/* csvdc.f -- translated by f2c (version 20240504).
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
static complex c_b8 = {1.f,0.f};
static complex c_b53 = {-1.f,0.f};

/* Subroutine */ int csvdc_(complex *x, integer *ldx, integer *n, integer *p, 
	complex *s, complex *e, complex *u, integer *ldu, complex *v, integer 
	*ldv, complex *work, integer *job, integer *info)
{
    /* System generated locals */
    integer x_dim1, x_offset, u_dim1, u_offset, v_dim1, v_offset, i__1, i__2, 
	    i__3, i__4;
    real r__1, r__2, r__3, r__4;
    complex q__1, q__2, q__3;

    /* Builtin functions */
    double r_imag(complex *), c_abs(complex *);
    void c_div(complex *, complex *, complex *), r_cnjg(complex *, complex *);
    double sqrt(doublereal);

    /* Local variables */
    real b, c__, f, g;
    integer i__, j, k, l, m;
    complex r__, t;
    real t1, el;
    integer kk;
    real cs;
    integer ll, mm, ls;
    real sl;
    integer lu;
    real sm, sn;
    integer lm1, mm1, lp1, mp1, nct, ncu, lls, nrt;
    real emm1, smm1;
    integer kase, jobu, iter;
    real test;
    integer nctp1, nrtp1;
    extern /* Subroutine */ int cscal_(integer *, complex *, complex *, 
	    integer *);
    real scale;
    extern /* Complex */ VOID cdotc_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);
    real shift;
    extern /* Subroutine */ int cswap_(integer *, complex *, integer *, 
	    complex *, integer *);
    integer maxit;
    extern /* Subroutine */ int caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *), csrot_(integer *, complex *, 
	    integer *, complex *, integer *, real *, real *);
    logical wantu, wantv;
    extern /* Subroutine */ int srotg_(real *, real *, real *, real *);
    real ztest;
    extern doublereal scnrm2_(integer *, complex *, integer *);



/*     csvdc is a subroutine to reduce a complex nxp matrix x by */
/*     unitary transformations u and v to diagonal form.  the */
/*     diagonal elements s(i) are the singular values of x.  the */
/*     columns of u are the corresponding left singular vectors, */
/*     and the columns of v the right singular vectors. */

/*     on entry */

/*         x         complex(ldx,p), where ldx.ge.n. */
/*                   x contains the matrix whose singular value */
/*                   decomposition is to be computed.  x is */
/*                   destroyed by csvdc. */

/*         ldx       integer. */
/*                   ldx is the leading dimension of the array x. */

/*         n         integer. */
/*                   n is the number of rows of the matrix x. */

/*         p         integer. */
/*                   p is the number of columns of the matrix x. */

/*         ldu       integer. */
/*                   ldu is the leading dimension of the array u */
/*                   (see below). */

/*         ldv       integer. */
/*                   ldv is the leading dimension of the array v */
/*                   (see below). */

/*         work      complex(n). */
/*                   work is a scratch array. */

/*         job       integer. */
/*                   job controls the computation of the singular */
/*                   vectors.  it has the decimal expansion ab */
/*                   with the following meaning */

/*                        a.eq.0    do not compute the left singular */
/*                                  vectors. */
/*                        a.eq.1    return the n left singular vectors */
/*                                  in u. */
/*                        a.ge.2    returns the first min(n,p) */
/*                                  left singular vectors in u. */
/*                        b.eq.0    do not compute the right singular */
/*                                  vectors. */
/*                        b.eq.1    return the right singular vectors */
/*                                  in v. */

/*     on return */

/*         s         complex(mm), where mm=min(n+1,p). */
/*                   the first min(n,p) entries of s contain the */
/*                   singular values of x arranged in descending */
/*                   order of magnitude. */

/*         e         complex(p). */
/*                   e ordinarily contains zeros.  however see the */
/*                   discussion of info for exceptions. */

/*         u         complex(ldu,k), where ldu.ge.n.  if joba.eq.1 then */
/*                                   k.eq.n, if joba.ge.2 then */
/*                                   k.eq.min(n,p). */
/*                   u contains the matrix of left singular vectors. */
/*                   u is not referenced if joba.eq.0.  if n.le.p */
/*                   or if joba.gt.2, then u may be identified with x */
/*                   in the subroutine call. */

/*         v         complex(ldv,p), where ldv.ge.p. */
/*                   v contains the matrix of right singular vectors. */
/*                   v is not referenced if jobb.eq.0.  if p.le.n, */
/*                   then v may be identified whth x in the */
/*                   subroutine call. */

/*         info      integer. */
/*                   the singular values (and their corresponding */
/*                   singular vectors) s(info+1),s(info+2),...,s(m) */
/*                   are correct (here m=min(n,p)).  thus if */
/*                   info.eq.0, all the singular values and their */
/*                   vectors are correct.  in any event, the matrix */
/*                   b = ctrans(u)*x*v is the bidiagonal matrix */
/*                   with the elements of s on its diagonal and the */
/*                   elements of e on its super-diagonal (ctrans(u) */
/*                   is the conjugate-transpose of u).  thus the */
/*                   singular values of x and b are the same. */

/*     linpack. this version dated 03/19/79 . */
/*              correction to shift calculation made 2/85. */
/*     g.w. stewart, university of maryland, argonne national lab. */

/*     csvdc uses the following functions and subprograms. */

/*     external csrot */
/*     blas caxpy,cdotc,cscal,cswap,scnrm2,srotg */
/*     fortran abs,aimag,amax1,cabs,cmplx */
/*     fortran conjg,max0,min0,mod,real,sqrt */

/*     internal variables */



/*     set the maximum number of iterations. */

    /* Parameter adjustments */
    x_dim1 = *ldx;
    x_offset = 1 + x_dim1;
    x -= x_offset;
    --s;
    --e;
    u_dim1 = *ldu;
    u_offset = 1 + u_dim1;
    u -= u_offset;
    v_dim1 = *ldv;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    --work;

    /* Function Body */
    maxit = 30;

/*     determine what is to be computed. */

    wantu = FALSE_;
    wantv = FALSE_;
    jobu = *job % 100 / 10;
    ncu = *n;
    if (jobu > 1) {
	ncu = min(*n,*p);
    }
    if (jobu != 0) {
	wantu = TRUE_;
    }
    if (*job % 10 != 0) {
	wantv = TRUE_;
    }

/*     reduce x to bidiagonal form, storing the diagonal elements */
/*     in s and the super-diagonal elements in e. */

    *info = 0;
/* Computing MIN */
    i__1 = *n - 1;
    nct = min(i__1,*p);
/* Computing MAX */
/* Computing MIN */
    i__3 = *p - 2;
    i__1 = 0, i__2 = min(i__3,*n);
    nrt = max(i__1,i__2);
    lu = max(nct,nrt);
    if (lu < 1) {
	goto L170;
    }
    i__1 = lu;
    for (l = 1; l <= i__1; ++l) {
	lp1 = l + 1;
	if (l > nct) {
	    goto L20;
	}

/*           compute the transformation for the l-th column and */
/*           place the l-th diagonal in s(l). */

	i__2 = l;
	i__3 = *n - l + 1;
	r__1 = scnrm2_(&i__3, &x[l + l * x_dim1], &c__1);
	q__1.r = r__1, q__1.i = 0.f;
	s[i__2].r = q__1.r, s[i__2].i = q__1.i;
	i__2 = l;
	if ((r__1 = s[i__2].r, dabs(r__1)) + (r__2 = r_imag(&s[l]), dabs(r__2)
		) == 0.f) {
	    goto L10;
	}
	i__2 = l + l * x_dim1;
	if ((r__1 = x[i__2].r, dabs(r__1)) + (r__2 = r_imag(&x[l + l * x_dim1]
		), dabs(r__2)) != 0.f) {
	    i__3 = l;
	    r__3 = c_abs(&s[l]);
	    i__4 = l + l * x_dim1;
	    r__4 = c_abs(&x[l + l * x_dim1]);
	    q__2.r = x[i__4].r / r__4, q__2.i = x[i__4].i / r__4;
	    q__1.r = r__3 * q__2.r, q__1.i = r__3 * q__2.i;
	    s[i__3].r = q__1.r, s[i__3].i = q__1.i;
	}
	i__2 = *n - l + 1;
	c_div(&q__1, &c_b8, &s[l]);
	cscal_(&i__2, &q__1, &x[l + l * x_dim1], &c__1);
	i__2 = l + l * x_dim1;
	i__3 = l + l * x_dim1;
	q__1.r = x[i__3].r + 1.f, q__1.i = x[i__3].i + 0.f;
	x[i__2].r = q__1.r, x[i__2].i = q__1.i;
L10:
	i__2 = l;
	i__3 = l;
	q__1.r = -s[i__3].r, q__1.i = -s[i__3].i;
	s[i__2].r = q__1.r, s[i__2].i = q__1.i;
L20:
	if (*p < lp1) {
	    goto L50;
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    if (l > nct) {
		goto L30;
	    }
	    i__3 = l;
	    if ((r__1 = s[i__3].r, dabs(r__1)) + (r__2 = r_imag(&s[l]), dabs(
		    r__2)) == 0.f) {
		goto L30;
	    }

/*              apply the transformation. */

	    i__3 = *n - l + 1;
	    cdotc_(&q__3, &i__3, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1]
		    , &c__1);
	    q__2.r = -q__3.r, q__2.i = -q__3.i;
	    c_div(&q__1, &q__2, &x[l + l * x_dim1]);
	    t.r = q__1.r, t.i = q__1.i;
	    i__3 = *n - l + 1;
	    caxpy_(&i__3, &t, &x[l + l * x_dim1], &c__1, &x[l + j * x_dim1], &
		    c__1);
L30:

/*           place the l-th row of x into  e for the */
/*           subsequent calculation of the row transformation. */

	    i__3 = j;
	    r_cnjg(&q__1, &x[l + j * x_dim1]);
	    e[i__3].r = q__1.r, e[i__3].i = q__1.i;
/* L40: */
	}
L50:
	if (! wantu || l > nct) {
	    goto L70;
	}

/*           place the transformation in u for subsequent back */
/*           multiplication. */

	i__2 = *n;
	for (i__ = l; i__ <= i__2; ++i__) {
	    i__3 = i__ + l * u_dim1;
	    i__4 = i__ + l * x_dim1;
	    u[i__3].r = x[i__4].r, u[i__3].i = x[i__4].i;
/* L60: */
	}
L70:
	if (l > nrt) {
	    goto L150;
	}

/*           compute the l-th row transformation and place the */
/*           l-th super-diagonal in e(l). */

	i__2 = l;
	i__3 = *p - l;
	r__1 = scnrm2_(&i__3, &e[lp1], &c__1);
	q__1.r = r__1, q__1.i = 0.f;
	e[i__2].r = q__1.r, e[i__2].i = q__1.i;
	i__2 = l;
	if ((r__1 = e[i__2].r, dabs(r__1)) + (r__2 = r_imag(&e[l]), dabs(r__2)
		) == 0.f) {
	    goto L80;
	}
	i__2 = lp1;
	if ((r__1 = e[i__2].r, dabs(r__1)) + (r__2 = r_imag(&e[lp1]), dabs(
		r__2)) != 0.f) {
	    i__3 = l;
	    r__3 = c_abs(&e[l]);
	    i__4 = lp1;
	    r__4 = c_abs(&e[lp1]);
	    q__2.r = e[i__4].r / r__4, q__2.i = e[i__4].i / r__4;
	    q__1.r = r__3 * q__2.r, q__1.i = r__3 * q__2.i;
	    e[i__3].r = q__1.r, e[i__3].i = q__1.i;
	}
	i__2 = *p - l;
	c_div(&q__1, &c_b8, &e[l]);
	cscal_(&i__2, &q__1, &e[lp1], &c__1);
	i__2 = lp1;
	i__3 = lp1;
	q__1.r = e[i__3].r + 1.f, q__1.i = e[i__3].i + 0.f;
	e[i__2].r = q__1.r, e[i__2].i = q__1.i;
L80:
	i__2 = l;
	r_cnjg(&q__2, &e[l]);
	q__1.r = -q__2.r, q__1.i = -q__2.i;
	e[i__2].r = q__1.r, e[i__2].i = q__1.i;
	i__2 = l;
	if (lp1 > *n || (r__1 = e[i__2].r, dabs(r__1)) + (r__2 = r_imag(&e[l])
		, dabs(r__2)) == 0.f) {
	    goto L120;
	}

/*              apply the transformation. */

	i__2 = *n;
	for (i__ = lp1; i__ <= i__2; ++i__) {
	    i__3 = i__;
	    work[i__3].r = 0.f, work[i__3].i = 0.f;
/* L90: */
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *n - l;
	    caxpy_(&i__3, &e[j], &x[lp1 + j * x_dim1], &c__1, &work[lp1], &
		    c__1);
/* L100: */
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *n - l;
	    i__4 = j;
	    q__3.r = -e[i__4].r, q__3.i = -e[i__4].i;
	    c_div(&q__2, &q__3, &e[lp1]);
	    r_cnjg(&q__1, &q__2);
	    caxpy_(&i__3, &q__1, &work[lp1], &c__1, &x[lp1 + j * x_dim1], &
		    c__1);
/* L110: */
	}
L120:
	if (! wantv) {
	    goto L140;
	}

/*              place the transformation in v for subsequent */
/*              back multiplication. */

	i__2 = *p;
	for (i__ = lp1; i__ <= i__2; ++i__) {
	    i__3 = i__ + l * v_dim1;
	    i__4 = i__;
	    v[i__3].r = e[i__4].r, v[i__3].i = e[i__4].i;
/* L130: */
	}
L140:
L150:
/* L160: */
	;
    }
L170:

/*     set up the final bidiagonal matrix or order m. */

/* Computing MIN */
    i__1 = *p, i__2 = *n + 1;
    m = min(i__1,i__2);
    nctp1 = nct + 1;
    nrtp1 = nrt + 1;
    if (nct < *p) {
	i__1 = nctp1;
	i__2 = nctp1 + nctp1 * x_dim1;
	s[i__1].r = x[i__2].r, s[i__1].i = x[i__2].i;
    }
    if (*n < m) {
	i__1 = m;
	s[i__1].r = 0.f, s[i__1].i = 0.f;
    }
    if (nrtp1 < m) {
	i__1 = nrtp1;
	i__2 = nrtp1 + m * x_dim1;
	e[i__1].r = x[i__2].r, e[i__1].i = x[i__2].i;
    }
    i__1 = m;
    e[i__1].r = 0.f, e[i__1].i = 0.f;

/*     if required, generate u. */

    if (! wantu) {
	goto L300;
    }
    if (ncu < nctp1) {
	goto L200;
    }
    i__1 = ncu;
    for (j = nctp1; j <= i__1; ++j) {
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + j * u_dim1;
	    u[i__3].r = 0.f, u[i__3].i = 0.f;
/* L180: */
	}
	i__2 = j + j * u_dim1;
	u[i__2].r = 1.f, u[i__2].i = 0.f;
/* L190: */
    }
L200:
    if (nct < 1) {
	goto L290;
    }
    i__1 = nct;
    for (ll = 1; ll <= i__1; ++ll) {
	l = nct - ll + 1;
	i__2 = l;
	if ((r__1 = s[i__2].r, dabs(r__1)) + (r__2 = r_imag(&s[l]), dabs(r__2)
		) == 0.f) {
	    goto L250;
	}
	lp1 = l + 1;
	if (ncu < lp1) {
	    goto L220;
	}
	i__2 = ncu;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *n - l + 1;
	    cdotc_(&q__3, &i__3, &u[l + l * u_dim1], &c__1, &u[l + j * u_dim1]
		    , &c__1);
	    q__2.r = -q__3.r, q__2.i = -q__3.i;
	    c_div(&q__1, &q__2, &u[l + l * u_dim1]);
	    t.r = q__1.r, t.i = q__1.i;
	    i__3 = *n - l + 1;
	    caxpy_(&i__3, &t, &u[l + l * u_dim1], &c__1, &u[l + j * u_dim1], &
		    c__1);
/* L210: */
	}
L220:
	i__2 = *n - l + 1;
	cscal_(&i__2, &c_b53, &u[l + l * u_dim1], &c__1);
	i__2 = l + l * u_dim1;
	i__3 = l + l * u_dim1;
	q__1.r = u[i__3].r + 1.f, q__1.i = u[i__3].i + 0.f;
	u[i__2].r = q__1.r, u[i__2].i = q__1.i;
	lm1 = l - 1;
	if (lm1 < 1) {
	    goto L240;
	}
	i__2 = lm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + l * u_dim1;
	    u[i__3].r = 0.f, u[i__3].i = 0.f;
/* L230: */
	}
L240:
	goto L270;
L250:
	i__2 = *n;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + l * u_dim1;
	    u[i__3].r = 0.f, u[i__3].i = 0.f;
/* L260: */
	}
	i__2 = l + l * u_dim1;
	u[i__2].r = 1.f, u[i__2].i = 0.f;
L270:
/* L280: */
	;
    }
L290:
L300:

/*     if it is required, generate v. */

    if (! wantv) {
	goto L350;
    }
    i__1 = *p;
    for (ll = 1; ll <= i__1; ++ll) {
	l = *p - ll + 1;
	lp1 = l + 1;
	if (l > nrt) {
	    goto L320;
	}
	i__2 = l;
	if ((r__1 = e[i__2].r, dabs(r__1)) + (r__2 = r_imag(&e[l]), dabs(r__2)
		) == 0.f) {
	    goto L320;
	}
	i__2 = *p;
	for (j = lp1; j <= i__2; ++j) {
	    i__3 = *p - l;
	    cdotc_(&q__3, &i__3, &v[lp1 + l * v_dim1], &c__1, &v[lp1 + j * 
		    v_dim1], &c__1);
	    q__2.r = -q__3.r, q__2.i = -q__3.i;
	    c_div(&q__1, &q__2, &v[lp1 + l * v_dim1]);
	    t.r = q__1.r, t.i = q__1.i;
	    i__3 = *p - l;
	    caxpy_(&i__3, &t, &v[lp1 + l * v_dim1], &c__1, &v[lp1 + j * 
		    v_dim1], &c__1);
/* L310: */
	}
L320:
	i__2 = *p;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + l * v_dim1;
	    v[i__3].r = 0.f, v[i__3].i = 0.f;
/* L330: */
	}
	i__2 = l + l * v_dim1;
	v[i__2].r = 1.f, v[i__2].i = 0.f;
/* L340: */
    }
L350:

/*     transform s and e so that they are real. */

    i__1 = m;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__;
	if ((r__1 = s[i__2].r, dabs(r__1)) + (r__2 = r_imag(&s[i__]), dabs(
		r__2)) == 0.f) {
	    goto L360;
	}
	r__1 = c_abs(&s[i__]);
	q__1.r = r__1, q__1.i = 0.f;
	t.r = q__1.r, t.i = q__1.i;
	c_div(&q__1, &s[i__], &t);
	r__.r = q__1.r, r__.i = q__1.i;
	i__2 = i__;
	s[i__2].r = t.r, s[i__2].i = t.i;
	if (i__ < m) {
	    i__2 = i__;
	    c_div(&q__1, &e[i__], &r__);
	    e[i__2].r = q__1.r, e[i__2].i = q__1.i;
	}
	if (wantu) {
	    cscal_(n, &r__, &u[i__ * u_dim1 + 1], &c__1);
	}
L360:
/*     ...exit */
	if (i__ == m) {
	    goto L390;
	}
	i__2 = i__;
	if ((r__1 = e[i__2].r, dabs(r__1)) + (r__2 = r_imag(&e[i__]), dabs(
		r__2)) == 0.f) {
	    goto L370;
	}
	r__1 = c_abs(&e[i__]);
	q__1.r = r__1, q__1.i = 0.f;
	t.r = q__1.r, t.i = q__1.i;
	c_div(&q__1, &t, &e[i__]);
	r__.r = q__1.r, r__.i = q__1.i;
	i__2 = i__;
	e[i__2].r = t.r, e[i__2].i = t.i;
	i__2 = i__ + 1;
	i__3 = i__ + 1;
	q__1.r = s[i__3].r * r__.r - s[i__3].i * r__.i, q__1.i = s[i__3].r * 
		r__.i + s[i__3].i * r__.r;
	s[i__2].r = q__1.r, s[i__2].i = q__1.i;
	if (wantv) {
	    cscal_(p, &r__, &v[(i__ + 1) * v_dim1 + 1], &c__1);
	}
L370:
/* L380: */
	;
    }
L390:

/*     main iteration loop for the singular values. */

    mm = m;
    iter = 0;
L400:

/*        quit if all the singular values have been found. */

/*     ...exit */
    if (m == 0) {
	goto L660;
    }

/*        if too many iterations have been performed, set */
/*        flag and return. */

    if (iter < maxit) {
	goto L410;
    }
    *info = m;
/*     ......exit */
    goto L660;
L410:

/*        this section of the program inspects for */
/*        negligible elements in the s and e arrays.  on */
/*        completion the variables kase and l are set as follows. */

/*           kase = 1     if s(m) and e(l-1) are negligible and l.lt.m */
/*           kase = 2     if s(l) is negligible and l.lt.m */
/*           kase = 3     if e(l-1) is negligible, l.lt.m, and */
/*                        s(l), ..., s(m) are not negligible (qr step). */
/*           kase = 4     if e(m-1) is negligible (convergence). */

    i__1 = m;
    for (ll = 1; ll <= i__1; ++ll) {
	l = m - ll;
/*        ...exit */
	if (l == 0) {
	    goto L440;
	}
	test = c_abs(&s[l]) + c_abs(&s[l + 1]);
	ztest = test + c_abs(&e[l]);
	if (ztest != test) {
	    goto L420;
	}
	i__2 = l;
	e[i__2].r = 0.f, e[i__2].i = 0.f;
/*        ......exit */
	goto L440;
L420:
/* L430: */
	;
    }
L440:
    if (l != m - 1) {
	goto L450;
    }
    kase = 4;
    goto L520;
L450:
    lp1 = l + 1;
    mp1 = m + 1;
    i__1 = mp1;
    for (lls = lp1; lls <= i__1; ++lls) {
	ls = m - lls + lp1;
/*           ...exit */
	if (ls == l) {
	    goto L480;
	}
	test = 0.f;
	if (ls != m) {
	    test += c_abs(&e[ls]);
	}
	if (ls != l + 1) {
	    test += c_abs(&e[ls - 1]);
	}
	ztest = test + c_abs(&s[ls]);
	if (ztest != test) {
	    goto L460;
	}
	i__2 = ls;
	s[i__2].r = 0.f, s[i__2].i = 0.f;
/*           ......exit */
	goto L480;
L460:
/* L470: */
	;
    }
L480:
    if (ls != l) {
	goto L490;
    }
    kase = 3;
    goto L510;
L490:
    if (ls != m) {
	goto L500;
    }
    kase = 1;
    goto L510;
L500:
    kase = 2;
    l = ls;
L510:
L520:
    ++l;

/*        perform the task indicated by kase. */

    switch (kase) {
	case 1:  goto L530;
	case 2:  goto L560;
	case 3:  goto L580;
	case 4:  goto L610;
    }

/*        deflate negligible s(m). */

L530:
    mm1 = m - 1;
    i__1 = m - 1;
    f = e[i__1].r;
    i__1 = m - 1;
    e[i__1].r = 0.f, e[i__1].i = 0.f;
    i__1 = mm1;
    for (kk = l; kk <= i__1; ++kk) {
	k = mm1 - kk + l;
	i__2 = k;
	t1 = s[i__2].r;
	srotg_(&t1, &f, &cs, &sn);
	i__2 = k;
	q__1.r = t1, q__1.i = 0.f;
	s[i__2].r = q__1.r, s[i__2].i = q__1.i;
	if (k == l) {
	    goto L540;
	}
	i__2 = k - 1;
	f = -sn * e[i__2].r;
	i__2 = k - 1;
	i__3 = k - 1;
	q__1.r = cs * e[i__3].r, q__1.i = cs * e[i__3].i;
	e[i__2].r = q__1.r, e[i__2].i = q__1.i;
L540:
	if (wantv) {
	    csrot_(p, &v[k * v_dim1 + 1], &c__1, &v[m * v_dim1 + 1], &c__1, &
		    cs, &sn);
	}
/* L550: */
    }
    goto L650;

/*        split at negligible s(l). */

L560:
    i__1 = l - 1;
    f = e[i__1].r;
    i__1 = l - 1;
    e[i__1].r = 0.f, e[i__1].i = 0.f;
    i__1 = m;
    for (k = l; k <= i__1; ++k) {
	i__2 = k;
	t1 = s[i__2].r;
	srotg_(&t1, &f, &cs, &sn);
	i__2 = k;
	q__1.r = t1, q__1.i = 0.f;
	s[i__2].r = q__1.r, s[i__2].i = q__1.i;
	i__2 = k;
	f = -sn * e[i__2].r;
	i__2 = k;
	i__3 = k;
	q__1.r = cs * e[i__3].r, q__1.i = cs * e[i__3].i;
	e[i__2].r = q__1.r, e[i__2].i = q__1.i;
	if (wantu) {
	    csrot_(n, &u[k * u_dim1 + 1], &c__1, &u[(l - 1) * u_dim1 + 1], &
		    c__1, &cs, &sn);
	}
/* L570: */
    }
    goto L650;

/*        perform one qr step. */

L580:

/*           calculate the shift. */

/* Computing MAX */
    r__1 = c_abs(&s[m]), r__2 = c_abs(&s[m - 1]), r__1 = max(r__1,r__2), r__2 
	    = c_abs(&e[m - 1]), r__1 = max(r__1,r__2), r__2 = c_abs(&s[l]), 
	    r__1 = max(r__1,r__2), r__2 = c_abs(&e[l]);
    scale = dmax(r__1,r__2);
    i__1 = m;
    sm = s[i__1].r / scale;
    i__1 = m - 1;
    smm1 = s[i__1].r / scale;
    i__1 = m - 1;
    emm1 = e[i__1].r / scale;
    i__1 = l;
    sl = s[i__1].r / scale;
    i__1 = l;
    el = e[i__1].r / scale;
/* Computing 2nd power */
    r__1 = emm1;
    b = ((smm1 + sm) * (smm1 - sm) + r__1 * r__1) / 2.f;
/* Computing 2nd power */
    r__1 = sm * emm1;
    c__ = r__1 * r__1;
    shift = 0.f;
    if (b == 0.f && c__ == 0.f) {
	goto L590;
    }
/* Computing 2nd power */
    r__1 = b;
    shift = sqrt(r__1 * r__1 + c__);
    if (b < 0.f) {
	shift = -shift;
    }
    shift = c__ / (b + shift);
L590:
    f = (sl + sm) * (sl - sm) + shift;
    g = sl * el;

/*           chase zeros. */

    mm1 = m - 1;
    i__1 = mm1;
    for (k = l; k <= i__1; ++k) {
	srotg_(&f, &g, &cs, &sn);
	if (k != l) {
	    i__2 = k - 1;
	    q__1.r = f, q__1.i = 0.f;
	    e[i__2].r = q__1.r, e[i__2].i = q__1.i;
	}
	i__2 = k;
	i__3 = k;
	f = cs * s[i__2].r + sn * e[i__3].r;
	i__2 = k;
	i__3 = k;
	q__2.r = cs * e[i__3].r, q__2.i = cs * e[i__3].i;
	i__4 = k;
	q__3.r = sn * s[i__4].r, q__3.i = sn * s[i__4].i;
	q__1.r = q__2.r - q__3.r, q__1.i = q__2.i - q__3.i;
	e[i__2].r = q__1.r, e[i__2].i = q__1.i;
	i__2 = k + 1;
	g = sn * s[i__2].r;
	i__2 = k + 1;
	i__3 = k + 1;
	q__1.r = cs * s[i__3].r, q__1.i = cs * s[i__3].i;
	s[i__2].r = q__1.r, s[i__2].i = q__1.i;
	if (wantv) {
	    csrot_(p, &v[k * v_dim1 + 1], &c__1, &v[(k + 1) * v_dim1 + 1], &
		    c__1, &cs, &sn);
	}
	srotg_(&f, &g, &cs, &sn);
	i__2 = k;
	q__1.r = f, q__1.i = 0.f;
	s[i__2].r = q__1.r, s[i__2].i = q__1.i;
	i__2 = k;
	i__3 = k + 1;
	f = cs * e[i__2].r + sn * s[i__3].r;
	i__2 = k + 1;
	r__1 = -sn;
	i__3 = k;
	q__2.r = r__1 * e[i__3].r, q__2.i = r__1 * e[i__3].i;
	i__4 = k + 1;
	q__3.r = cs * s[i__4].r, q__3.i = cs * s[i__4].i;
	q__1.r = q__2.r + q__3.r, q__1.i = q__2.i + q__3.i;
	s[i__2].r = q__1.r, s[i__2].i = q__1.i;
	i__2 = k + 1;
	g = sn * e[i__2].r;
	i__2 = k + 1;
	i__3 = k + 1;
	q__1.r = cs * e[i__3].r, q__1.i = cs * e[i__3].i;
	e[i__2].r = q__1.r, e[i__2].i = q__1.i;
	if (wantu && k < *n) {
	    csrot_(n, &u[k * u_dim1 + 1], &c__1, &u[(k + 1) * u_dim1 + 1], &
		    c__1, &cs, &sn);
	}
/* L600: */
    }
    i__1 = m - 1;
    q__1.r = f, q__1.i = 0.f;
    e[i__1].r = q__1.r, e[i__1].i = q__1.i;
    ++iter;
    goto L650;

/*        convergence. */

L610:

/*           make the singular value  positive */

    i__1 = l;
    if (s[i__1].r >= 0.f) {
	goto L620;
    }
    i__1 = l;
    i__2 = l;
    q__1.r = -s[i__2].r, q__1.i = -s[i__2].i;
    s[i__1].r = q__1.r, s[i__1].i = q__1.i;
    if (wantv) {
	cscal_(p, &c_b53, &v[l * v_dim1 + 1], &c__1);
    }
L620:

/*           order the singular value. */

L630:
    if (l == mm) {
	goto L640;
    }
/*           ...exit */
    i__1 = l;
    i__2 = l + 1;
    if (s[i__1].r >= s[i__2].r) {
	goto L640;
    }
    i__1 = l;
    t.r = s[i__1].r, t.i = s[i__1].i;
    i__1 = l;
    i__2 = l + 1;
    s[i__1].r = s[i__2].r, s[i__1].i = s[i__2].i;
    i__1 = l + 1;
    s[i__1].r = t.r, s[i__1].i = t.i;
    if (wantv && l < *p) {
	cswap_(p, &v[l * v_dim1 + 1], &c__1, &v[(l + 1) * v_dim1 + 1], &c__1);
    }
    if (wantu && l < *n) {
	cswap_(n, &u[l * u_dim1 + 1], &c__1, &u[(l + 1) * u_dim1 + 1], &c__1);
    }
    ++l;
    goto L630;
L640:
    iter = 0;
    --m;
L650:
    goto L400;
L660:
    return 0;
} /* csvdc_ */

