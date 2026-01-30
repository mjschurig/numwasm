/* ztrsl.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int ztrsl_(doublecomplex *t, integer *ldt, integer *n, 
	doublecomplex *b, integer *job, integer *info)
{
    /* System generated locals */
    integer t_dim1, t_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *), d_cnjg(
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    integer j, jj, case__;
    doublecomplex temp;
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);



/*     ztrsl solves systems of the form */

/*                   t * x = b */
/*     or */
/*                   ctrans(t) * x = b */

/*     where t is a triangular matrix of order n. here ctrans(t) */
/*     denotes the conjugate transpose of the matrix t. */

/*     on entry */

/*         t         complex*16(ldt,n) */
/*                   t contains the matrix of the system. the zero */
/*                   elements of the matrix are not referenced, and */
/*                   the corresponding elements of the array can be */
/*                   used to store other information. */

/*         ldt       integer */
/*                   ldt is the leading dimension of the array t. */

/*         n         integer */
/*                   n is the order of the system. */

/*         b         complex*16(n). */
/*                   b contains the right hand side of the system. */

/*         job       integer */
/*                   job specifies what kind of system is to be solved. */
/*                   if job is */

/*                        00   solve t*x=b, t lower triangular, */
/*                        01   solve t*x=b, t upper triangular, */
/*                        10   solve ctrans(t)*x=b, t lower triangular, */
/*                        11   solve ctrans(t)*x=b, t upper triangular. */

/*     on return */

/*         b         b contains the solution, if info .eq. 0. */
/*                   otherwise b is unaltered. */

/*         info      integer */
/*                   info contains zero if the system is nonsingular. */
/*                   otherwise info contains the index of */
/*                   the first zero diagonal element of t. */

/*     linpack. this version dated 08/14/78 . */
/*     g. w. stewart, university of maryland, argonne national lab. */

/*     subroutines and functions */

/*     blas zaxpy,zdotc */
/*     fortran dabs,dconjg,mod */

/*     internal variables */


/*     begin block permitting ...exits to 150 */

/*        check for zero diagonal elements. */

    /* Parameter adjustments */
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --b;

    /* Function Body */
    i__1 = *n;
    for (*info = 1; *info <= i__1; ++(*info)) {
/*     ......exit */
	i__2 = *info + *info * t_dim1;
	i__3 = *info + *info * t_dim1;
	z__1.r = t[i__3].r * 0. - t[i__3].i * -1., z__1.i = t[i__3].i * 0. + 
		t[i__3].r * -1.;
	if ((d__1 = t[i__2].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) == 0.) 
		{
	    goto L150;
	}
/* L10: */
    }
    *info = 0;

/*        determine the task and go to it. */

    case__ = 1;
    if (*job % 10 != 0) {
	case__ = 2;
    }
    if (*job % 100 / 10 != 0) {
	case__ += 2;
    }
    switch (case__) {
	case 1:  goto L20;
	case 2:  goto L50;
	case 3:  goto L80;
	case 4:  goto L110;
    }

/*        solve t*x=b for t lower triangular */

L20:
    z_div(&z__1, &b[1], &t[t_dim1 + 1]);
    b[1].r = z__1.r, b[1].i = z__1.i;
    if (*n < 2) {
	goto L40;
    }
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
	i__2 = j - 1;
	z__1.r = -b[i__2].r, z__1.i = -b[i__2].i;
	temp.r = z__1.r, temp.i = z__1.i;
	i__2 = *n - j + 1;
	zaxpy_(&i__2, &temp, &t[j + (j - 1) * t_dim1], &c__1, &b[j], &c__1);
	i__2 = j;
	z_div(&z__1, &b[j], &t[j + j * t_dim1]);
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
/* L30: */
    }
L40:
    goto L140;

/*        solve t*x=b for t upper triangular. */

L50:
    i__1 = *n;
    z_div(&z__1, &b[*n], &t[*n + *n * t_dim1]);
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    if (*n < 2) {
	goto L70;
    }
    i__1 = *n;
    for (jj = 2; jj <= i__1; ++jj) {
	j = *n - jj + 1;
	i__2 = j + 1;
	z__1.r = -b[i__2].r, z__1.i = -b[i__2].i;
	temp.r = z__1.r, temp.i = z__1.i;
	zaxpy_(&j, &temp, &t[(j + 1) * t_dim1 + 1], &c__1, &b[1], &c__1);
	i__2 = j;
	z_div(&z__1, &b[j], &t[j + j * t_dim1]);
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
/* L60: */
    }
L70:
    goto L140;

/*        solve ctrans(t)*x=b for t lower triangular. */

L80:
    i__1 = *n;
    d_cnjg(&z__2, &t[*n + *n * t_dim1]);
    z_div(&z__1, &b[*n], &z__2);
    b[i__1].r = z__1.r, b[i__1].i = z__1.i;
    if (*n < 2) {
	goto L100;
    }
    i__1 = *n;
    for (jj = 2; jj <= i__1; ++jj) {
	j = *n - jj + 1;
	i__2 = j;
	i__3 = j;
	i__4 = jj - 1;
	zdotc_(&z__2, &i__4, &t[j + 1 + j * t_dim1], &c__1, &b[j + 1], &c__1);
	z__1.r = b[i__3].r - z__2.r, z__1.i = b[i__3].i - z__2.i;
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
	i__2 = j;
	d_cnjg(&z__2, &t[j + j * t_dim1]);
	z_div(&z__1, &b[j], &z__2);
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
/* L90: */
    }
L100:
    goto L140;

/*        solve ctrans(t)*x=b for t upper triangular. */

L110:
    d_cnjg(&z__2, &t[t_dim1 + 1]);
    z_div(&z__1, &b[1], &z__2);
    b[1].r = z__1.r, b[1].i = z__1.i;
    if (*n < 2) {
	goto L130;
    }
    i__1 = *n;
    for (j = 2; j <= i__1; ++j) {
	i__2 = j;
	i__3 = j;
	i__4 = j - 1;
	zdotc_(&z__2, &i__4, &t[j * t_dim1 + 1], &c__1, &b[1], &c__1);
	z__1.r = b[i__3].r - z__2.r, z__1.i = b[i__3].i - z__2.i;
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
	i__2 = j;
	d_cnjg(&z__2, &t[j + j * t_dim1]);
	z_div(&z__1, &b[j], &z__2);
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
/* L120: */
    }
L130:
L140:
L150:
    return 0;
} /* ztrsl_ */

