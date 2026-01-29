/* sppsl.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int sppsl_(real *ap, integer *n, real *b)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer k;
    real t;
    integer kb, kk;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *);


/*     sppsl solves the real symmetric positive definite system */
/*     a * x = b */
/*     using the factors computed by sppco or sppfa. */

/*     on entry */

/*        ap      real (n*(n+1)/2) */
/*                the output from sppco or sppfa. */

/*        n       integer */
/*                the order of the matrix  a . */

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
/*           call sppco(ap,n,rcond,z,info) */
/*           if (rcond is too small .or. info .ne. 0) go to ... */
/*           do 10 j = 1, p */
/*              call sppsl(ap,n,c(1,j)) */
/*        10 continue */

/*     linpack.  this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas saxpy,sdot */

/*     internal variables */


    /* Parameter adjustments */
    --b;
    --ap;

    /* Function Body */
    kk = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	i__2 = k - 1;
	t = sdot_(&i__2, &ap[kk + 1], &c__1, &b[1], &c__1);
	kk += k;
	b[k] = (b[k] - t) / ap[kk];
/* L10: */
    }
    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	b[k] /= ap[kk];
	kk -= k;
	t = -b[k];
	i__2 = k - 1;
	saxpy_(&i__2, &t, &ap[kk + 1], &c__1, &b[1], &c__1);
/* L20: */
    }
    return 0;
} /* sppsl_ */

