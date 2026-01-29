/* ztrco.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int ztrco_(doublecomplex *t, integer *ldt, integer *n, 
	doublereal *rcond, doublecomplex *z__, integer *job)
{
    /* System generated locals */
    integer t_dim1, t_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *), z_div(doublecomplex *, 
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    integer j, k, l;
    doublereal s;
    doublecomplex w;
    integer i1, j1, j2;
    doublecomplex ek;
    integer kk;
    doublereal sm;
    doublecomplex wk, wkm;
    logical lower;
    doublereal tnorm, ynorm;
    extern /* Subroutine */ int zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), zdscal_(
	    integer *, doublereal *, doublecomplex *, integer *);
    extern doublereal dzasum_(integer *, doublecomplex *, integer *);


/*     ztrco estimates the condition of a complex*16 triangular matrix. */

/*     on entry */

/*        t       complex*16(ldt,n) */
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

/*        rcond   double precision */
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

/*        z       complex*16(n) */
/*                a work vector whose contents are usually unimportant. */
/*                if  t  is close to a singular matrix, then  z  is */
/*                an approximate null vector in the sense that */
/*                norm(a*z) = rcond*norm(a)*norm(z) . */

/*     linpack. this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas zaxpy,zdscal,dzasum */
/*     fortran dabs,dmax1,dcmplx,dconjg */

/*     internal variables */


    /* Parameter adjustments */
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --z__;

    /* Function Body */
    lower = *job == 0;

/*     compute 1-norm of t */

    tnorm = 0.;
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
	d__1 = tnorm, d__2 = dzasum_(&l, &t[i1 + j * t_dim1], &c__1);
	tnorm = max(d__1,d__2);
/* L10: */
    }

/*     rcond = 1/(norm(t)*(estimate of norm(inverse(t)))) . */
/*     estimate = norm(z)/norm(y) where  t*z = y  and  ctrans(t)*y = e . */
/*     ctrans(t)  is the conjugate transpose of t . */
/*     the components of  e  are chosen to cause maximum local */
/*     growth in the elements of y . */
/*     the vectors are frequently rescaled to avoid overflow. */

/*     solve ctrans(t)*y = e */

    ek.r = 1., ek.i = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	z__[i__2].r = 0., z__[i__2].i = 0.;
/* L20: */
    }
    i__1 = *n;
    for (kk = 1; kk <= i__1; ++kk) {
	k = kk;
	if (lower) {
	    k = *n + 1 - kk;
	}
	i__2 = k;
	i__3 = k;
	z__1.r = z__[i__3].r * 0. - z__[i__3].i * -1., z__1.i = z__[i__3].i * 
		0. + z__[i__3].r * -1.;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) != 
		0.) {
	    i__4 = k;
	    z__3.r = -z__[i__4].r, z__3.i = -z__[i__4].i;
	    z__2.r = z__3.r, z__2.i = z__3.i;
	    z__5.r = ek.r * 0. - ek.i * -1., z__5.i = ek.i * 0. + ek.r * -1.;
	    d__7 = (d__3 = ek.r, abs(d__3)) + (d__4 = z__5.r, abs(d__4));
	    z__7.r = z__2.r * 0. - z__2.i * -1., z__7.i = z__2.i * 0. + 
		    z__2.r * -1.;
	    d__8 = (d__5 = z__2.r, abs(d__5)) + (d__6 = z__7.r, abs(d__6));
	    z__6.r = z__2.r / d__8, z__6.i = z__2.i / d__8;
	    z__4.r = d__7 * z__6.r, z__4.i = d__7 * z__6.i;
	    ek.r = z__4.r, ek.i = z__4.i;
	}
	i__2 = k;
	z__2.r = ek.r - z__[i__2].r, z__2.i = ek.i - z__[i__2].i;
	z__1.r = z__2.r, z__1.i = z__2.i;
	z__3.r = z__1.r * 0. - z__1.i * -1., z__3.i = z__1.i * 0. + z__1.r * 
		-1.;
	i__3 = k + k * t_dim1;
	i__4 = k + k * t_dim1;
	z__4.r = t[i__4].r * 0. - t[i__4].i * -1., z__4.i = t[i__4].i * 0. + 
		t[i__4].r * -1.;
	if ((d__1 = z__1.r, abs(d__1)) + (d__2 = z__3.r, abs(d__2)) <= (d__3 =
		 t[i__3].r, abs(d__3)) + (d__4 = z__4.r, abs(d__4))) {
	    goto L30;
	}
	i__2 = k;
	z__2.r = ek.r - z__[i__2].r, z__2.i = ek.i - z__[i__2].i;
	z__1.r = z__2.r, z__1.i = z__2.i;
	i__3 = k + k * t_dim1;
	i__4 = k + k * t_dim1;
	z__3.r = t[i__4].r * 0. - t[i__4].i * -1., z__3.i = t[i__4].i * 0. + 
		t[i__4].r * -1.;
	z__4.r = z__1.r * 0. - z__1.i * -1., z__4.i = z__1.i * 0. + z__1.r * 
		-1.;
	s = ((d__1 = t[i__3].r, abs(d__1)) + (d__2 = z__3.r, abs(d__2))) / ((
		d__3 = z__1.r, abs(d__3)) + (d__4 = z__4.r, abs(d__4)));
	zdscal_(n, &s, &z__[1], &c__1);
	z__2.r = s, z__2.i = 0.;
	z__1.r = z__2.r * ek.r - z__2.i * ek.i, z__1.i = z__2.r * ek.i + 
		z__2.i * ek.r;
	ek.r = z__1.r, ek.i = z__1.i;
L30:
	i__2 = k;
	z__1.r = ek.r - z__[i__2].r, z__1.i = ek.i - z__[i__2].i;
	wk.r = z__1.r, wk.i = z__1.i;
	z__2.r = -ek.r, z__2.i = -ek.i;
	i__2 = k;
	z__1.r = z__2.r - z__[i__2].r, z__1.i = z__2.i - z__[i__2].i;
	wkm.r = z__1.r, wkm.i = z__1.i;
	z__1.r = wk.r * 0. - wk.i * -1., z__1.i = wk.i * 0. + wk.r * -1.;
	s = (d__1 = wk.r, abs(d__1)) + (d__2 = z__1.r, abs(d__2));
	z__1.r = wkm.r * 0. - wkm.i * -1., z__1.i = wkm.i * 0. + wkm.r * -1.;
	sm = (d__1 = wkm.r, abs(d__1)) + (d__2 = z__1.r, abs(d__2));
	i__2 = k + k * t_dim1;
	i__3 = k + k * t_dim1;
	z__1.r = t[i__3].r * 0. - t[i__3].i * -1., z__1.i = t[i__3].i * 0. + 
		t[i__3].r * -1.;
	if ((d__1 = t[i__2].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) == 0.) 
		{
	    goto L40;
	}
	d_cnjg(&z__2, &t[k + k * t_dim1]);
	z_div(&z__1, &wk, &z__2);
	wk.r = z__1.r, wk.i = z__1.i;
	d_cnjg(&z__2, &t[k + k * t_dim1]);
	z_div(&z__1, &wkm, &z__2);
	wkm.r = z__1.r, wkm.i = z__1.i;
	goto L50;
L40:
	wk.r = 1., wk.i = 0.;
	wkm.r = 1., wkm.i = 0.;
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
	    i__3 = j;
	    d_cnjg(&z__4, &t[k + j * t_dim1]);
	    z__3.r = wkm.r * z__4.r - wkm.i * z__4.i, z__3.i = wkm.r * z__4.i 
		    + wkm.i * z__4.r;
	    z__2.r = z__[i__3].r + z__3.r, z__2.i = z__[i__3].i + z__3.i;
	    z__1.r = z__2.r, z__1.i = z__2.i;
	    z__5.r = z__1.r * 0. - z__1.i * -1., z__5.i = z__1.i * 0. + 
		    z__1.r * -1.;
	    sm += (d__1 = z__1.r, abs(d__1)) + (d__2 = z__5.r, abs(d__2));
	    i__3 = j;
	    i__4 = j;
	    d_cnjg(&z__3, &t[k + j * t_dim1]);
	    z__2.r = wk.r * z__3.r - wk.i * z__3.i, z__2.i = wk.r * z__3.i + 
		    wk.i * z__3.r;
	    z__1.r = z__[i__4].r + z__2.r, z__1.i = z__[i__4].i + z__2.i;
	    z__[i__3].r = z__1.r, z__[i__3].i = z__1.i;
	    i__3 = j;
	    i__4 = j;
	    z__1.r = z__[i__4].r * 0. - z__[i__4].i * -1., z__1.i = z__[i__4]
		    .i * 0. + z__[i__4].r * -1.;
	    s += (d__1 = z__[i__3].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2));
/* L60: */
	}
	if (s >= sm) {
	    goto L80;
	}
	z__1.r = wkm.r - wk.r, z__1.i = wkm.i - wk.i;
	w.r = z__1.r, w.i = z__1.i;
	wk.r = wkm.r, wk.i = wkm.i;
	i__2 = j2;
	for (j = j1; j <= i__2; ++j) {
	    i__3 = j;
	    i__4 = j;
	    d_cnjg(&z__3, &t[k + j * t_dim1]);
	    z__2.r = w.r * z__3.r - w.i * z__3.i, z__2.i = w.r * z__3.i + w.i 
		    * z__3.r;
	    z__1.r = z__[i__4].r + z__2.r, z__1.i = z__[i__4].i + z__2.i;
	    z__[i__3].r = z__1.r, z__[i__3].i = z__1.i;
/* L70: */
	}
L80:
L90:
	i__2 = k;
	z__[i__2].r = wk.r, z__[i__2].i = wk.i;
/* L100: */
    }
    s = 1. / dzasum_(n, &z__[1], &c__1);
    zdscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.;

/*     solve t*z = y */

    i__1 = *n;
    for (kk = 1; kk <= i__1; ++kk) {
	k = *n + 1 - kk;
	if (lower) {
	    k = kk;
	}
	i__2 = k;
	i__3 = k;
	z__1.r = z__[i__3].r * 0. - z__[i__3].i * -1., z__1.i = z__[i__3].i * 
		0. + z__[i__3].r * -1.;
	i__4 = k + k * t_dim1;
	i__5 = k + k * t_dim1;
	z__2.r = t[i__5].r * 0. - t[i__5].i * -1., z__2.i = t[i__5].i * 0. + 
		t[i__5].r * -1.;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) <= (
		d__3 = t[i__4].r, abs(d__3)) + (d__4 = z__2.r, abs(d__4))) {
	    goto L110;
	}
	i__2 = k + k * t_dim1;
	i__3 = k + k * t_dim1;
	z__1.r = t[i__3].r * 0. - t[i__3].i * -1., z__1.i = t[i__3].i * 0. + 
		t[i__3].r * -1.;
	i__4 = k;
	i__5 = k;
	z__2.r = z__[i__5].r * 0. - z__[i__5].i * -1., z__2.i = z__[i__5].i * 
		0. + z__[i__5].r * -1.;
	s = ((d__1 = t[i__2].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2))) / ((
		d__3 = z__[i__4].r, abs(d__3)) + (d__4 = z__2.r, abs(d__4)));
	zdscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L110:
	i__2 = k + k * t_dim1;
	i__3 = k + k * t_dim1;
	z__1.r = t[i__3].r * 0. - t[i__3].i * -1., z__1.i = t[i__3].i * 0. + 
		t[i__3].r * -1.;
	if ((d__1 = t[i__2].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) != 0.) 
		{
	    i__4 = k;
	    z_div(&z__2, &z__[k], &t[k + k * t_dim1]);
	    z__[i__4].r = z__2.r, z__[i__4].i = z__2.i;
	}
	i__2 = k + k * t_dim1;
	i__3 = k + k * t_dim1;
	z__1.r = t[i__3].r * 0. - t[i__3].i * -1., z__1.i = t[i__3].i * 0. + 
		t[i__3].r * -1.;
	if ((d__1 = t[i__2].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) == 0.) 
		{
	    i__4 = k;
	    z__[i__4].r = 1., z__[i__4].i = 0.;
	}
	i1 = 1;
	if (lower) {
	    i1 = k + 1;
	}
	if (kk >= *n) {
	    goto L120;
	}
	i__2 = k;
	z__1.r = -z__[i__2].r, z__1.i = -z__[i__2].i;
	w.r = z__1.r, w.i = z__1.i;
	i__2 = *n - kk;
	zaxpy_(&i__2, &w, &t[i1 + k * t_dim1], &c__1, &z__[i1], &c__1);
L120:
/* L130: */
	;
    }
/*     make znorm = 1.0 */
    s = 1. / dzasum_(n, &z__[1], &c__1);
    zdscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (tnorm != 0.) {
	*rcond = ynorm / tnorm;
    }
    if (tnorm == 0.) {
	*rcond = 0.;
    }
    return 0;
} /* ztrco_ */

