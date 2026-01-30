/* zppsl.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int zppsl_(doublecomplex *ap, integer *n, doublecomplex *b)
{
    /* System generated locals */
    integer i__1, i__2, i__3;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    integer k;
    doublecomplex t;
    integer kb, kk;
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern /* Subroutine */ int zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);


/*     zppsl solves the complex*16 hermitian positive definite system */
/*     a * x = b */
/*     using the factors computed by zppco or zppfa. */

/*     on entry */

/*        ap      complex*16 (n*(n+1)/2) */
/*                the output from zppco or zppfa. */

/*        n       integer */
/*                the order of the matrix  a . */

/*        b       complex*16(n) */
/*                the right hand side vector. */

/*     on return */

/*        b       the solution vector  x . */

/*     error condition */

/*        a division by zero will occur if the input factor contains */
/*        a zero on the diagonal.  technically this indicates */
/*        singularity but it is usually caused by improper subroutine */
/*        arguments.  it will not occur if the subroutines are called */
/*        correctly and  info .eq. 0 . */

/*     to compute  inverse(a) * c  where  c  is a matrix */
/*     with  p  columns */
/*           call zppco(ap,n,rcond,z,info) */
/*           if (rcond is too small .or. info .ne. 0) go to ... */
/*           do 10 j = 1, p */
/*              call zppsl(ap,n,c(1,j)) */
/*        10 continue */

/*     linpack.  this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas zaxpy,zdotc */

/*     internal variables */


    /* Parameter adjustments */
    --b;
    --ap;

    /* Function Body */
    kk = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k - 1;
	zdotc_(&z__1, &i__2, &ap[kk + 1], &c__1, &b[1], &c__1);
	t.r = z__1.r, t.i = z__1.i;
	kk += k;
	i__2 = k;
	i__3 = k;
	z__2.r = b[i__3].r - t.r, z__2.i = b[i__3].i - t.i;
	z_div(&z__1, &z__2, &ap[kk]);
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
/* L10: */
    }
    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	i__2 = k;
	z_div(&z__1, &b[k], &ap[kk]);
	b[i__2].r = z__1.r, b[i__2].i = z__1.i;
	kk -= k;
	i__2 = k;
	z__1.r = -b[i__2].r, z__1.i = -b[i__2].i;
	t.r = z__1.r, t.i = z__1.i;
	i__2 = k - 1;
	zaxpy_(&i__2, &t, &ap[kk + 1], &c__1, &b[1], &c__1);
/* L20: */
    }
    return 0;
} /* zppsl_ */

