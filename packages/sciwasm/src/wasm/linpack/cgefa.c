/* cgefa.f -- translated by f2c (version 20240504).
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
static complex c_b7 = {-1.f,-0.f};

/* Subroutine */ int cgefa_(complex *a, integer *lda, integer *n, integer *
	ipvt, integer *info)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2;
    complex q__1;

    /* Builtin functions */
    double r_imag(complex *);
    void c_div(complex *, complex *, complex *);

    /* Local variables */
    integer j, k, l;
    complex t;
    integer kp1, nm1;
    extern /* Subroutine */ int cscal_(integer *, complex *, complex *, 
	    integer *), caxpy_(integer *, complex *, complex *, integer *, 
	    complex *, integer *);
    extern integer icamax_(integer *, complex *, integer *);


/*     cgefa factors a complex matrix by gaussian elimination. */

/*     cgefa is usually called by cgeco, but it can be called */
/*     directly with a saving in time if  rcond  is not needed. */
/*     (time for cgeco) = (1 + 9/n)*(time for cgefa) . */

/*     on entry */

/*        a       complex(lda, n) */
/*                the matrix to be factored. */

/*        lda     integer */
/*                the leading dimension of the array  a . */

/*        n       integer */
/*                the order of the matrix  a . */

/*     on return */

/*        a       an upper triangular matrix and the multipliers */
/*                which were used to obtain it. */
/*                the factorization can be written  a = l*u  where */
/*                l  is a product of permutation and unit lower */
/*                triangular matrices and  u  is upper triangular. */

/*        ipvt    integer(n) */
/*                an integer vector of pivot indices. */

/*        info    integer */
/*                = 0  normal value. */
/*                = k  if  u(k,k) .eq. 0.0 .  this is not an error */
/*                     condition for this subroutine, but it does */
/*                     indicate that cgesl or cgedi will divide by zero */
/*                     if called.  use  rcond  in cgeco for a reliable */
/*                     indication of singularity. */

/*     linpack. this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas caxpy,cscal,icamax */
/*     fortran abs,aimag,real */

/*     internal variables */



/*     gaussian elimination with partial pivoting */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --ipvt;

    /* Function Body */
    *info = 0;
    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L70;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;

/*        find l = pivot index */

	i__2 = *n - k + 1;
	l = icamax_(&i__2, &a[k + k * a_dim1], &c__1) + k - 1;
	ipvt[k] = l;

/*        zero pivot implies this column already triangularized */

	i__2 = l + k * a_dim1;
	if ((r__1 = a[i__2].r, dabs(r__1)) + (r__2 = r_imag(&a[l + k * a_dim1]
		), dabs(r__2)) == 0.f) {
	    goto L40;
	}

/*           interchange if necessary */

	if (l == k) {
	    goto L10;
	}
	i__2 = l + k * a_dim1;
	t.r = a[i__2].r, t.i = a[i__2].i;
	i__2 = l + k * a_dim1;
	i__3 = k + k * a_dim1;
	a[i__2].r = a[i__3].r, a[i__2].i = a[i__3].i;
	i__2 = k + k * a_dim1;
	a[i__2].r = t.r, a[i__2].i = t.i;
L10:

/*           compute multipliers */

	c_div(&q__1, &c_b7, &a[k + k * a_dim1]);
	t.r = q__1.r, t.i = q__1.i;
	i__2 = *n - k;
	cscal_(&i__2, &t, &a[k + 1 + k * a_dim1], &c__1);

/*           row elimination with column indexing */

	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    i__3 = l + j * a_dim1;
	    t.r = a[i__3].r, t.i = a[i__3].i;
	    if (l == k) {
		goto L20;
	    }
	    i__3 = l + j * a_dim1;
	    i__4 = k + j * a_dim1;
	    a[i__3].r = a[i__4].r, a[i__3].i = a[i__4].i;
	    i__3 = k + j * a_dim1;
	    a[i__3].r = t.r, a[i__3].i = t.i;
L20:
	    i__3 = *n - k;
	    caxpy_(&i__3, &t, &a[k + 1 + k * a_dim1], &c__1, &a[k + 1 + j * 
		    a_dim1], &c__1);
/* L30: */
	}
	goto L50;
L40:
	*info = k;
L50:
/* L60: */
	;
    }
L70:
    ipvt[*n] = *n;
    i__1 = *n + *n * a_dim1;
    if ((r__1 = a[i__1].r, dabs(r__1)) + (r__2 = r_imag(&a[*n + *n * a_dim1]),
	     dabs(r__2)) == 0.f) {
	*info = *n;
    }
    return 0;
} /* cgefa_ */

