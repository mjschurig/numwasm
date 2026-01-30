/* cpbsl.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int cpbsl_(complex *abd, integer *lda, integer *n, integer *
	m, complex *b)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3;
    complex q__1, q__2;

    /* Builtin functions */
    void c_div(complex *, complex *, complex *);

    /* Local variables */
    integer k;
    complex t;
    integer kb, la, lb, lm;
    extern /* Complex */ VOID cdotc_(complex *, integer *, complex *, integer 
	    *, complex *, integer *);
    extern /* Subroutine */ int caxpy_(integer *, complex *, complex *, 
	    integer *, complex *, integer *);


/*     cpbsl solves the complex hermitian positive definite band */
/*     system  a*x = b */
/*     using the factors computed by cpbco or cpbfa. */

/*     on entry */

/*        abd     complex(lda, n) */
/*                the output from cpbco or cpbfa. */

/*        lda     integer */
/*                the leading dimension of the array  abd . */

/*        n       integer */
/*                the order of the matrix  a . */

/*        m       integer */
/*                the number of diagonals above the main diagonal. */

/*        b       complex(n) */
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
/*           call cpbco(abd,lda,n,rcond,z,info) */
/*           if (rcond is too small .or. info .ne. 0) go to ... */
/*           do 10 j = 1, p */
/*              call cpbsl(abd,lda,n,c(1,j)) */
/*        10 continue */

/*     linpack.  this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas caxpy,cdotc */
/*     fortran min0 */

/*     internal variables */


/*     solve ctrans(r)*y = b */

    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --b;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* Computing MIN */
	i__2 = k - 1;
	lm = min(i__2,*m);
	la = *m + 1 - lm;
	lb = k - lm;
	cdotc_(&q__1, &lm, &abd[la + k * abd_dim1], &c__1, &b[lb], &c__1);
	t.r = q__1.r, t.i = q__1.i;
	i__2 = k;
	i__3 = k;
	q__2.r = b[i__3].r - t.r, q__2.i = b[i__3].i - t.i;
	c_div(&q__1, &q__2, &abd[*m + 1 + k * abd_dim1]);
	b[i__2].r = q__1.r, b[i__2].i = q__1.i;
/* L10: */
    }

/*     solve r*x = y */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
/* Computing MIN */
	i__2 = k - 1;
	lm = min(i__2,*m);
	la = *m + 1 - lm;
	lb = k - lm;
	i__2 = k;
	c_div(&q__1, &b[k], &abd[*m + 1 + k * abd_dim1]);
	b[i__2].r = q__1.r, b[i__2].i = q__1.i;
	i__2 = k;
	q__1.r = -b[i__2].r, q__1.i = -b[i__2].i;
	t.r = q__1.r, t.i = q__1.i;
	caxpy_(&lm, &t, &abd[la + k * abd_dim1], &c__1, &b[lb], &c__1);
/* L20: */
    }
    return 0;
} /* cpbsl_ */

