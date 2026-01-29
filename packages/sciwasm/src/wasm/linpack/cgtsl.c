/* cgtsl.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int cgtsl_(integer *n, complex *c__, complex *d__, complex *
	e, complex *b, integer *info)
{
    /* System generated locals */
    integer i__1, i__2, i__3, i__4, i__5, i__6, i__7;
    real r__1, r__2, r__3, r__4;
    complex q__1, q__2, q__3, q__4, q__5;

    /* Builtin functions */
    double r_imag(complex *);
    void c_div(complex *, complex *, complex *);

    /* Local variables */
    integer k;
    complex t;
    integer kb, kp1, nm1, nm2;


/*     cgtsl given a general tridiagonal matrix and a right hand */
/*     side will find the solution. */

/*     on entry */

/*        n       integer */
/*                is the order of the tridiagonal matrix. */

/*        c       complex(n) */
/*                is the subdiagonal of the tridiagonal matrix. */
/*                c(2) through c(n) should contain the subdiagonal. */
/*                on output c is destroyed. */

/*        d       complex(n) */
/*                is the diagonal of the tridiagonal matrix. */
/*                on output d is destroyed. */

/*        e       complex(n) */
/*                is the superdiagonal of the tridiagonal matrix. */
/*                e(1) through e(n-1) should contain the superdiagonal. */
/*                on output e is destroyed. */

/*        b       complex(n) */
/*                is the right hand side vector. */

/*     on return */

/*        b       is the solution vector. */

/*        info    integer */
/*                = 0 normal value. */
/*                = k if the k-th element of the diagonal becomes */
/*                    exactly zero.  the subroutine returns when */
/*                    this is detected. */

/*     linpack. this version dated 08/14/78 . */
/*     jack dongarra, argonne national laboratory. */

/*     no externals */
/*     fortran abs,aimag,real */

/*     internal variables */

/*     begin block permitting ...exits to 100 */

    /* Parameter adjustments */
    --b;
    --e;
    --d__;
    --c__;

    /* Function Body */
    *info = 0;
    c__[1].r = d__[1].r, c__[1].i = d__[1].i;
    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L40;
    }
    d__[1].r = e[1].r, d__[1].i = e[1].i;
    e[1].r = 0.f, e[1].i = 0.f;
    i__1 = *n;
    e[i__1].r = 0.f, e[i__1].i = 0.f;

    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;

/*              find the largest of the two rows */

	i__2 = kp1;
	i__3 = k;
	if ((r__1 = c__[i__2].r, dabs(r__1)) + (r__2 = r_imag(&c__[kp1]), 
		dabs(r__2)) < (r__3 = c__[i__3].r, dabs(r__3)) + (r__4 = 
		r_imag(&c__[k]), dabs(r__4))) {
	    goto L10;
	}

/*                 interchange row */

	i__2 = kp1;
	t.r = c__[i__2].r, t.i = c__[i__2].i;
	i__2 = kp1;
	i__3 = k;
	c__[i__2].r = c__[i__3].r, c__[i__2].i = c__[i__3].i;
	i__2 = k;
	c__[i__2].r = t.r, c__[i__2].i = t.i;
	i__2 = kp1;
	t.r = d__[i__2].r, t.i = d__[i__2].i;
	i__2 = kp1;
	i__3 = k;
	d__[i__2].r = d__[i__3].r, d__[i__2].i = d__[i__3].i;
	i__2 = k;
	d__[i__2].r = t.r, d__[i__2].i = t.i;
	i__2 = kp1;
	t.r = e[i__2].r, t.i = e[i__2].i;
	i__2 = kp1;
	i__3 = k;
	e[i__2].r = e[i__3].r, e[i__2].i = e[i__3].i;
	i__2 = k;
	e[i__2].r = t.r, e[i__2].i = t.i;
	i__2 = kp1;
	t.r = b[i__2].r, t.i = b[i__2].i;
	i__2 = kp1;
	i__3 = k;
	b[i__2].r = b[i__3].r, b[i__2].i = b[i__3].i;
	i__2 = k;
	b[i__2].r = t.r, b[i__2].i = t.i;
L10:

/*              zero elements */

	i__2 = k;
	if ((r__1 = c__[i__2].r, dabs(r__1)) + (r__2 = r_imag(&c__[k]), dabs(
		r__2)) != 0.f) {
	    goto L20;
	}
	*info = k;
/*     ............exit */
	goto L100;
L20:
	i__2 = kp1;
	q__2.r = -c__[i__2].r, q__2.i = -c__[i__2].i;
	c_div(&q__1, &q__2, &c__[k]);
	t.r = q__1.r, t.i = q__1.i;
	i__2 = kp1;
	i__3 = kp1;
	i__4 = k;
	q__2.r = t.r * d__[i__4].r - t.i * d__[i__4].i, q__2.i = t.r * d__[
		i__4].i + t.i * d__[i__4].r;
	q__1.r = d__[i__3].r + q__2.r, q__1.i = d__[i__3].i + q__2.i;
	c__[i__2].r = q__1.r, c__[i__2].i = q__1.i;
	i__2 = kp1;
	i__3 = kp1;
	i__4 = k;
	q__2.r = t.r * e[i__4].r - t.i * e[i__4].i, q__2.i = t.r * e[i__4].i 
		+ t.i * e[i__4].r;
	q__1.r = e[i__3].r + q__2.r, q__1.i = e[i__3].i + q__2.i;
	d__[i__2].r = q__1.r, d__[i__2].i = q__1.i;
	i__2 = kp1;
	e[i__2].r = 0.f, e[i__2].i = 0.f;
	i__2 = kp1;
	i__3 = kp1;
	i__4 = k;
	q__2.r = t.r * b[i__4].r - t.i * b[i__4].i, q__2.i = t.r * b[i__4].i 
		+ t.i * b[i__4].r;
	q__1.r = b[i__3].r + q__2.r, q__1.i = b[i__3].i + q__2.i;
	b[i__2].r = q__1.r, b[i__2].i = q__1.i;
/* L30: */
    }
L40:
    i__1 = *n;
    if ((r__1 = c__[i__1].r, dabs(r__1)) + (r__2 = r_imag(&c__[*n]), dabs(
	    r__2)) != 0.f) {
	goto L50;
    }
    *info = *n;
    goto L90;
L50:

/*           back solve */

    nm2 = *n - 2;
    i__1 = *n;
    c_div(&q__1, &b[*n], &c__[*n]);
    b[i__1].r = q__1.r, b[i__1].i = q__1.i;
    if (*n == 1) {
	goto L80;
    }
    i__1 = nm1;
    i__2 = nm1;
    i__3 = nm1;
    i__4 = *n;
    q__3.r = d__[i__3].r * b[i__4].r - d__[i__3].i * b[i__4].i, q__3.i = d__[
	    i__3].r * b[i__4].i + d__[i__3].i * b[i__4].r;
    q__2.r = b[i__2].r - q__3.r, q__2.i = b[i__2].i - q__3.i;
    c_div(&q__1, &q__2, &c__[nm1]);
    b[i__1].r = q__1.r, b[i__1].i = q__1.i;
    if (nm2 < 1) {
	goto L70;
    }
    i__1 = nm2;
    for (kb = 1; kb <= i__1; ++kb) {
	k = nm2 - kb + 1;
	i__2 = k;
	i__3 = k;
	i__4 = k;
	i__5 = k + 1;
	q__4.r = d__[i__4].r * b[i__5].r - d__[i__4].i * b[i__5].i, q__4.i = 
		d__[i__4].r * b[i__5].i + d__[i__4].i * b[i__5].r;
	q__3.r = b[i__3].r - q__4.r, q__3.i = b[i__3].i - q__4.i;
	i__6 = k;
	i__7 = k + 2;
	q__5.r = e[i__6].r * b[i__7].r - e[i__6].i * b[i__7].i, q__5.i = e[
		i__6].r * b[i__7].i + e[i__6].i * b[i__7].r;
	q__2.r = q__3.r - q__5.r, q__2.i = q__3.i - q__5.i;
	c_div(&q__1, &q__2, &c__[k]);
	b[i__2].r = q__1.r, b[i__2].i = q__1.i;
/* L60: */
    }
L70:
L80:
L90:
L100:

    return 0;
} /* cgtsl_ */

