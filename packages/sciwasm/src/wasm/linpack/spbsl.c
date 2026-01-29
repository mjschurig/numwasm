/* spbsl.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int spbsl_(real *abd, integer *lda, integer *n, integer *m, 
	real *b)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2;

    /* Local variables */
    integer k;
    real t;
    integer kb, la, lb, lm;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *);


/*     spbsl solves the real symmetric positive definite band */
/*     system  a*x = b */
/*     using the factors computed by spbco or spbfa. */

/*     on entry */

/*        abd     real(lda, n) */
/*                the output from spbco or spbfa. */

/*        lda     integer */
/*                the leading dimension of the array  abd . */

/*        n       integer */
/*                the order of the matrix  a . */

/*        m       integer */
/*                the number of diagonals above the main diagonal. */

/*        b       real(n) */
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
/*           call spbco(abd,lda,n,rcond,z,info) */
/*           if (rcond is too small .or. info .ne. 0) go to ... */
/*           do 10 j = 1, p */
/*              call spbsl(abd,lda,n,c(1,j)) */
/*        10 continue */

/*     linpack.  this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas saxpy,sdot */
/*     fortran min0 */

/*     internal variables */


/*     solve trans(r)*y = b */

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
	t = sdot_(&lm, &abd[la + k * abd_dim1], &c__1, &b[lb], &c__1);
	b[k] = (b[k] - t) / abd[*m + 1 + k * abd_dim1];
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
	b[k] /= abd[*m + 1 + k * abd_dim1];
	t = -b[k];
	saxpy_(&lm, &t, &abd[la + k * abd_dim1], &c__1, &b[lb], &c__1);
/* L20: */
    }
    return 0;
} /* spbsl_ */

