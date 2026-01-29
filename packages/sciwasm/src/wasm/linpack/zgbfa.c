/* zgbfa.f -- translated by f2c (version 20240504).
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
static doublecomplex c_b15 = {-1.,-0.};

/* Subroutine */ int zgbfa_(doublecomplex *abd, integer *lda, integer *n, 
	integer *ml, integer *mu, integer *ipvt, integer *info)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3, i__4;
    doublereal d__1, d__2;
    doublecomplex z__1;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    integer i__, j, k, l, m;
    doublecomplex t;
    integer i0, j0, j1, lm, mm, ju, jz, kp1, nm1;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    extern integer izamax_(integer *, doublecomplex *, integer *);


/*     zgbfa factors a complex*16 band matrix by elimination. */

/*     zgbfa is usually called by zgbco, but it can be called */
/*     directly with a saving in time if  rcond  is not needed. */

/*     on entry */

/*        abd     complex*16(lda, n) */
/*                contains the matrix in band storage.  the columns */
/*                of the matrix are stored in the columns of  abd  and */
/*                the diagonals of the matrix are stored in rows */
/*                ml+1 through 2*ml+mu+1 of  abd . */
/*                see the comments below for details. */

/*        lda     integer */
/*                the leading dimension of the array  abd . */
/*                lda must be .ge. 2*ml + mu + 1 . */

/*        n       integer */
/*                the order of the original matrix. */

/*        ml      integer */
/*                number of diagonals below the main diagonal. */
/*                0 .le. ml .lt. n . */

/*        mu      integer */
/*                number of diagonals above the main diagonal. */
/*                0 .le. mu .lt. n . */
/*                more efficient if  ml .le. mu . */
/*     on return */

/*        abd     an upper triangular matrix in band storage and */
/*                the multipliers which were used to obtain it. */
/*                the factorization can be written  a = l*u  where */
/*                l  is a product of permutation and unit lower */
/*                triangular matrices and  u  is upper triangular. */

/*        ipvt    integer(n) */
/*                an integer vector of pivot indices. */

/*        info    integer */
/*                = 0  normal value. */
/*                = k  if  u(k,k) .eq. 0.0 .  this is not an error */
/*                     condition for this subroutine, but it does */
/*                     indicate that zgbsl will divide by zero if */
/*                     called.  use  rcond  in zgbco for a reliable */
/*                     indication of singularity. */

/*     band storage */

/*           if  a  is a band matrix, the following program segment */
/*           will set up the input. */

/*                   ml = (band width below the diagonal) */
/*                   mu = (band width above the diagonal) */
/*                   m = ml + mu + 1 */
/*                   do 20 j = 1, n */
/*                      i1 = max0(1, j-mu) */
/*                      i2 = min0(n, j+ml) */
/*                      do 10 i = i1, i2 */
/*                         k = i - j + m */
/*                         abd(k,j) = a(i,j) */
/*                10    continue */
/*                20 continue */

/*           this uses rows  ml+1  through  2*ml+mu+1  of  abd . */
/*           in addition, the first  ml  rows in  abd  are used for */
/*           elements generated during the triangularization. */
/*           the total number of rows needed in  abd  is  2*ml+mu+1 . */
/*           the  ml+mu by ml+mu  upper left triangle and the */
/*           ml by ml  lower right triangle are not referenced. */

/*     linpack. this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas zaxpy,zscal,izamax */
/*     fortran dabs,max0,min0 */

/*     internal variables */



    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --ipvt;

    /* Function Body */
    m = *ml + *mu + 1;
    *info = 0;

/*     zero initial fill-in columns */

    j0 = *mu + 2;
    j1 = min(*n,m) - 1;
    if (j1 < j0) {
	goto L30;
    }
    i__1 = j1;
    for (jz = j0; jz <= i__1; ++jz) {
	i0 = m + 1 - jz;
	i__2 = *ml;
	for (i__ = i0; i__ <= i__2; ++i__) {
	    i__3 = i__ + jz * abd_dim1;
	    abd[i__3].r = 0., abd[i__3].i = 0.;
/* L10: */
	}
/* L20: */
    }
L30:
    jz = j1;
    ju = 0;

/*     gaussian elimination with partial pivoting */

    nm1 = *n - 1;
    if (nm1 < 1) {
	goto L130;
    }
    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;

/*        zero next fill-in column */

	++jz;
	if (jz > *n) {
	    goto L50;
	}
	if (*ml < 1) {
	    goto L50;
	}
	i__2 = *ml;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    i__3 = i__ + jz * abd_dim1;
	    abd[i__3].r = 0., abd[i__3].i = 0.;
/* L40: */
	}
L50:

/*        find l = pivot index */

/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	lm = min(i__2,i__3);
	i__2 = lm + 1;
	l = izamax_(&i__2, &abd[m + k * abd_dim1], &c__1) + m - 1;
	ipvt[k] = l + k - m;

/*        zero pivot implies this column already triangularized */

	i__2 = l + k * abd_dim1;
	i__3 = l + k * abd_dim1;
	z__1.r = abd[i__3].r * 0. - abd[i__3].i * -1., z__1.i = abd[i__3].i * 
		0. + abd[i__3].r * -1.;
	if ((d__1 = abd[i__2].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) == 
		0.) {
	    goto L100;
	}

/*           interchange if necessary */

	if (l == m) {
	    goto L60;
	}
	i__2 = l + k * abd_dim1;
	t.r = abd[i__2].r, t.i = abd[i__2].i;
	i__2 = l + k * abd_dim1;
	i__3 = m + k * abd_dim1;
	abd[i__2].r = abd[i__3].r, abd[i__2].i = abd[i__3].i;
	i__2 = m + k * abd_dim1;
	abd[i__2].r = t.r, abd[i__2].i = t.i;
L60:

/*           compute multipliers */

	z_div(&z__1, &c_b15, &abd[m + k * abd_dim1]);
	t.r = z__1.r, t.i = z__1.i;
	zscal_(&lm, &t, &abd[m + 1 + k * abd_dim1], &c__1);

/*           row elimination with column indexing */

/* Computing MIN */
/* Computing MAX */
	i__3 = ju, i__4 = *mu + ipvt[k];
	i__2 = max(i__3,i__4);
	ju = min(i__2,*n);
	mm = m;
	if (ju < kp1) {
	    goto L90;
	}
	i__2 = ju;
	for (j = kp1; j <= i__2; ++j) {
	    --l;
	    --mm;
	    i__3 = l + j * abd_dim1;
	    t.r = abd[i__3].r, t.i = abd[i__3].i;
	    if (l == mm) {
		goto L70;
	    }
	    i__3 = l + j * abd_dim1;
	    i__4 = mm + j * abd_dim1;
	    abd[i__3].r = abd[i__4].r, abd[i__3].i = abd[i__4].i;
	    i__3 = mm + j * abd_dim1;
	    abd[i__3].r = t.r, abd[i__3].i = t.i;
L70:
	    zaxpy_(&lm, &t, &abd[m + 1 + k * abd_dim1], &c__1, &abd[mm + 1 + 
		    j * abd_dim1], &c__1);
/* L80: */
	}
L90:
	goto L110;
L100:
	*info = k;
L110:
/* L120: */
	;
    }
L130:
    ipvt[*n] = *n;
    i__1 = m + *n * abd_dim1;
    i__2 = m + *n * abd_dim1;
    z__1.r = abd[i__2].r * 0. - abd[i__2].i * -1., z__1.i = abd[i__2].i * 0. 
	    + abd[i__2].r * -1.;
    if ((d__1 = abd[i__1].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) == 0.) {
	*info = *n;
    }
    return 0;
} /* zgbfa_ */

