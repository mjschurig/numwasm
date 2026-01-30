/* zgbco.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int zgbco_(doublecomplex *abd, integer *lda, integer *n, 
	integer *ml, integer *mu, integer *ipvt, doublereal *rcond, 
	doublecomplex *z__)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3, i__4, i__5;
    doublereal d__1, d__2, d__3, d__4, d__5, d__6, d__7, d__8;
    doublecomplex z__1, z__2, z__3, z__4, z__5, z__6, z__7;

    /* Builtin functions */
    void d_cnjg(doublecomplex *, doublecomplex *), z_div(doublecomplex *, 
	    doublecomplex *, doublecomplex *);

    /* Local variables */
    integer j, k, l, m;
    doublereal s;
    doublecomplex t;
    integer kb, la;
    doublecomplex ek;
    integer lm, mm, is, ju;
    doublereal sm;
    doublecomplex wk;
    integer lz, kp1;
    doublecomplex wkm;
    integer info;
    extern /* Subroutine */ int zgbfa_(doublecomplex *, integer *, integer *, 
	    integer *, integer *, integer *, integer *);
    doublereal anorm;
    extern /* Double Complex */ VOID zdotc_(doublecomplex *, integer *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);
    doublereal ynorm;
    extern /* Subroutine */ int zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *), zdscal_(
	    integer *, doublereal *, doublecomplex *, integer *);
    extern doublereal dzasum_(integer *, doublecomplex *, integer *);


/*     zgbco factors a complex*16 band matrix by gaussian */
/*     elimination and estimates the condition of the matrix. */

/*     if  rcond  is not needed, zgbfa is slightly faster. */
/*     to solve  a*x = b , follow zgbco by zgbsl. */
/*     to compute  inverse(a)*c , follow zgbco by zgbsl. */
/*     to compute  determinant(a) , follow zgbco by zgbdi. */

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

/*        rcond   double precision */
/*                an estimate of the reciprocal condition of  a . */
/*                for the system  a*x = b , relative perturbations */
/*                in  a  and  b  of size  epsilon  may cause */
/*                relative perturbations in  x  of size  epsilon/rcond . */
/*                if  rcond  is so small that the logical expression */
/*                           1.0 + rcond .eq. 1.0 */
/*                is true, then  a  may be singular to working */
/*                precision.  in particular,  rcond  is zero  if */
/*                exact singularity is detected or the estimate */
/*                underflows. */

/*        z       complex*16(n) */
/*                a work vector whose contents are usually unimportant. */
/*                if  a  is close to a singular matrix, then  z  is */
/*                an approximate null vector in the sense that */
/*                norm(a*z) = rcond*norm(a)*norm(z) . */

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

/*     example..  if the original matrix is */

/*           11 12 13  0  0  0 */
/*           21 22 23 24  0  0 */
/*            0 32 33 34 35  0 */
/*            0  0 43 44 45 46 */
/*            0  0  0 54 55 56 */
/*            0  0  0  0 65 66 */

/*      then  n = 6, ml = 1, mu = 2, lda .ge. 5  and abd should contain */

/*            *  *  *  +  +  +  , * = not used */
/*            *  * 13 24 35 46  , + = used for pivoting */
/*            * 12 23 34 45 56 */
/*           11 22 33 44 55 66 */
/*           21 32 43 54 65  * */

/*     linpack. this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     linpack zgbfa */
/*     blas zaxpy,zdotc,zdscal,dzasum */
/*     fortran dabs,dmax1,dcmplx,dconjg,max0,min0 */

/*     internal variables */



/*     compute 1-norm of a */

    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --ipvt;
    --z__;

    /* Function Body */
    anorm = 0.;
    l = *ml + 1;
    is = l + *mu;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	d__1 = anorm, d__2 = dzasum_(&l, &abd[is + j * abd_dim1], &c__1);
	anorm = max(d__1,d__2);
	if (is > *ml + 1) {
	    --is;
	}
	if (j <= *mu) {
	    ++l;
	}
	if (j >= *n - *ml) {
	    --l;
	}
/* L10: */
    }

/*     factor */

    zgbfa_(&abd[abd_offset], lda, n, ml, mu, &ipvt[1], &info);

/*     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) . */
/*     estimate = norm(z)/norm(y) where  a*z = y  and  ctrans(a)*y = e . */
/*     ctrans(a)  is the conjugate transpose of a . */
/*     the components of  e  are chosen to cause maximum local */
/*     growth in the elements of w  where  ctrans(u)*w = e . */
/*     the vectors are frequently rescaled to avoid overflow. */

/*     solve ctrans(u)*w = e */

    ek.r = 1., ek.i = 0.;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	i__2 = j;
	z__[i__2].r = 0., z__[i__2].i = 0.;
/* L20: */
    }
    m = *ml + *mu + 1;
    ju = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
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
	i__3 = m + k * abd_dim1;
	i__4 = m + k * abd_dim1;
	z__4.r = abd[i__4].r * 0. - abd[i__4].i * -1., z__4.i = abd[i__4].i * 
		0. + abd[i__4].r * -1.;
	if ((d__1 = z__1.r, abs(d__1)) + (d__2 = z__3.r, abs(d__2)) <= (d__3 =
		 abd[i__3].r, abs(d__3)) + (d__4 = z__4.r, abs(d__4))) {
	    goto L30;
	}
	i__2 = k;
	z__2.r = ek.r - z__[i__2].r, z__2.i = ek.i - z__[i__2].i;
	z__1.r = z__2.r, z__1.i = z__2.i;
	i__3 = m + k * abd_dim1;
	i__4 = m + k * abd_dim1;
	z__3.r = abd[i__4].r * 0. - abd[i__4].i * -1., z__3.i = abd[i__4].i * 
		0. + abd[i__4].r * -1.;
	z__4.r = z__1.r * 0. - z__1.i * -1., z__4.i = z__1.i * 0. + z__1.r * 
		-1.;
	s = ((d__1 = abd[i__3].r, abs(d__1)) + (d__2 = z__3.r, abs(d__2))) / (
		(d__3 = z__1.r, abs(d__3)) + (d__4 = z__4.r, abs(d__4)));
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
	i__2 = m + k * abd_dim1;
	i__3 = m + k * abd_dim1;
	z__1.r = abd[i__3].r * 0. - abd[i__3].i * -1., z__1.i = abd[i__3].i * 
		0. + abd[i__3].r * -1.;
	if ((d__1 = abd[i__2].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) == 
		0.) {
	    goto L40;
	}
	d_cnjg(&z__2, &abd[m + k * abd_dim1]);
	z_div(&z__1, &wk, &z__2);
	wk.r = z__1.r, wk.i = z__1.i;
	d_cnjg(&z__2, &abd[m + k * abd_dim1]);
	z_div(&z__1, &wkm, &z__2);
	wkm.r = z__1.r, wkm.i = z__1.i;
	goto L50;
L40:
	wk.r = 1., wk.i = 0.;
	wkm.r = 1., wkm.i = 0.;
L50:
	kp1 = k + 1;
/* Computing MIN */
/* Computing MAX */
	i__3 = ju, i__4 = *mu + ipvt[k];
	i__2 = max(i__3,i__4);
	ju = min(i__2,*n);
	mm = m;
	if (kp1 > ju) {
	    goto L90;
	}
	i__2 = ju;
	for (j = kp1; j <= i__2; ++j) {
	    --mm;
	    i__3 = j;
	    d_cnjg(&z__4, &abd[mm + j * abd_dim1]);
	    z__3.r = wkm.r * z__4.r - wkm.i * z__4.i, z__3.i = wkm.r * z__4.i 
		    + wkm.i * z__4.r;
	    z__2.r = z__[i__3].r + z__3.r, z__2.i = z__[i__3].i + z__3.i;
	    z__1.r = z__2.r, z__1.i = z__2.i;
	    z__5.r = z__1.r * 0. - z__1.i * -1., z__5.i = z__1.i * 0. + 
		    z__1.r * -1.;
	    sm += (d__1 = z__1.r, abs(d__1)) + (d__2 = z__5.r, abs(d__2));
	    i__3 = j;
	    i__4 = j;
	    d_cnjg(&z__3, &abd[mm + j * abd_dim1]);
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
	t.r = z__1.r, t.i = z__1.i;
	wk.r = wkm.r, wk.i = wkm.i;
	mm = m;
	i__2 = ju;
	for (j = kp1; j <= i__2; ++j) {
	    --mm;
	    i__3 = j;
	    i__4 = j;
	    d_cnjg(&z__3, &abd[mm + j * abd_dim1]);
	    z__2.r = t.r * z__3.r - t.i * z__3.i, z__2.i = t.r * z__3.i + t.i 
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

/*     solve ctrans(l)*y = w */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	lm = min(i__2,i__3);
	if (k < *n) {
	    i__2 = k;
	    i__3 = k;
	    zdotc_(&z__2, &lm, &abd[m + 1 + k * abd_dim1], &c__1, &z__[k + 1],
		     &c__1);
	    z__1.r = z__[i__3].r + z__2.r, z__1.i = z__[i__3].i + z__2.i;
	    z__[i__2].r = z__1.r, z__[i__2].i = z__1.i;
	}
	i__2 = k;
	i__3 = k;
	z__1.r = z__[i__3].r * 0. - z__[i__3].i * -1., z__1.i = z__[i__3].i * 
		0. + z__[i__3].r * -1.;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) <= 
		1.) {
	    goto L110;
	}
	i__2 = k;
	i__3 = k;
	z__1.r = z__[i__3].r * 0. - z__[i__3].i * -1., z__1.i = z__[i__3].i * 
		0. + z__[i__3].r * -1.;
	s = 1. / ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2))
		);
	zdscal_(n, &s, &z__[1], &c__1);
L110:
	l = ipvt[k];
	i__2 = l;
	t.r = z__[i__2].r, t.i = z__[i__2].i;
	i__2 = l;
	i__3 = k;
	z__[i__2].r = z__[i__3].r, z__[i__2].i = z__[i__3].i;
	i__2 = k;
	z__[i__2].r = t.r, z__[i__2].i = t.i;
/* L120: */
    }
    s = 1. / dzasum_(n, &z__[1], &c__1);
    zdscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.;

/*     solve l*v = y */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	l = ipvt[k];
	i__2 = l;
	t.r = z__[i__2].r, t.i = z__[i__2].i;
	i__2 = l;
	i__3 = k;
	z__[i__2].r = z__[i__3].r, z__[i__2].i = z__[i__3].i;
	i__2 = k;
	z__[i__2].r = t.r, z__[i__2].i = t.i;
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	lm = min(i__2,i__3);
	if (k < *n) {
	    zaxpy_(&lm, &t, &abd[m + 1 + k * abd_dim1], &c__1, &z__[k + 1], &
		    c__1);
	}
	i__2 = k;
	i__3 = k;
	z__1.r = z__[i__3].r * 0. - z__[i__3].i * -1., z__1.i = z__[i__3].i * 
		0. + z__[i__3].r * -1.;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) <= 
		1.) {
	    goto L130;
	}
	i__2 = k;
	i__3 = k;
	z__1.r = z__[i__3].r * 0. - z__[i__3].i * -1., z__1.i = z__[i__3].i * 
		0. + z__[i__3].r * -1.;
	s = 1. / ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2))
		);
	zdscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L130:
/* L140: */
	;
    }
    s = 1. / dzasum_(n, &z__[1], &c__1);
    zdscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*     solve  u*z = w */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	i__2 = k;
	i__3 = k;
	z__1.r = z__[i__3].r * 0. - z__[i__3].i * -1., z__1.i = z__[i__3].i * 
		0. + z__[i__3].r * -1.;
	i__4 = m + k * abd_dim1;
	i__5 = m + k * abd_dim1;
	z__2.r = abd[i__5].r * 0. - abd[i__5].i * -1., z__2.i = abd[i__5].i * 
		0. + abd[i__5].r * -1.;
	if ((d__1 = z__[i__2].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) <= (
		d__3 = abd[i__4].r, abs(d__3)) + (d__4 = z__2.r, abs(d__4))) {
	    goto L150;
	}
	i__2 = m + k * abd_dim1;
	i__3 = m + k * abd_dim1;
	z__1.r = abd[i__3].r * 0. - abd[i__3].i * -1., z__1.i = abd[i__3].i * 
		0. + abd[i__3].r * -1.;
	i__4 = k;
	i__5 = k;
	z__2.r = z__[i__5].r * 0. - z__[i__5].i * -1., z__2.i = z__[i__5].i * 
		0. + z__[i__5].r * -1.;
	s = ((d__1 = abd[i__2].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2))) / (
		(d__3 = z__[i__4].r, abs(d__3)) + (d__4 = z__2.r, abs(d__4)));
	zdscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L150:
	i__2 = m + k * abd_dim1;
	i__3 = m + k * abd_dim1;
	z__1.r = abd[i__3].r * 0. - abd[i__3].i * -1., z__1.i = abd[i__3].i * 
		0. + abd[i__3].r * -1.;
	if ((d__1 = abd[i__2].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) != 
		0.) {
	    i__4 = k;
	    z_div(&z__2, &z__[k], &abd[m + k * abd_dim1]);
	    z__[i__4].r = z__2.r, z__[i__4].i = z__2.i;
	}
	i__2 = m + k * abd_dim1;
	i__3 = m + k * abd_dim1;
	z__1.r = abd[i__3].r * 0. - abd[i__3].i * -1., z__1.i = abd[i__3].i * 
		0. + abd[i__3].r * -1.;
	if ((d__1 = abd[i__2].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) == 
		0.) {
	    i__4 = k;
	    z__[i__4].r = 1., z__[i__4].i = 0.;
	}
	lm = min(k,m) - 1;
	la = m - lm;
	lz = k - lm;
	i__2 = k;
	z__1.r = -z__[i__2].r, z__1.i = -z__[i__2].i;
	t.r = z__1.r, t.i = z__1.i;
	zaxpy_(&lm, &t, &abd[la + k * abd_dim1], &c__1, &z__[lz], &c__1);
/* L160: */
    }
/*     make znorm = 1.0 */
    s = 1. / dzasum_(n, &z__[1], &c__1);
    zdscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.) {
	*rcond = 0.;
    }
    return 0;
} /* zgbco_ */

