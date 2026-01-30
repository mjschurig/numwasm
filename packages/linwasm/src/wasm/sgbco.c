/* sgbco.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int sgbco_(real *abd, integer *lda, integer *n, integer *ml, 
	integer *mu, integer *ipvt, real *rcond, real *z__)
{
    /* System generated locals */
    integer abd_dim1, abd_offset, i__1, i__2, i__3, i__4;
    real r__1, r__2;

    /* Builtin functions */
    double r_sign(real *, real *);

    /* Local variables */
    integer j, k, l, m;
    real s, t;
    integer kb, la;
    real ek;
    integer lm, mm, is, ju;
    real sm, wk;
    integer lz, kp1;
    real wkm;
    integer info;
    extern doublereal sdot_(integer *, real *, integer *, real *, integer *);
    extern /* Subroutine */ int sgbfa_(real *, integer *, integer *, integer *
	    , integer *, integer *, integer *), sscal_(integer *, real *, 
	    real *, integer *);
    real anorm;
    extern doublereal sasum_(integer *, real *, integer *);
    real ynorm;
    extern /* Subroutine */ int saxpy_(integer *, real *, real *, integer *, 
	    real *, integer *);


/*     sgbco factors a real band matrix by gaussian */
/*     elimination and estimates the condition of the matrix. */

/*     if  rcond  is not needed, sgbfa is slightly faster. */
/*     to solve  a*x = b , follow sgbco by sgbsl. */
/*     to compute  inverse(a)*c , follow sgbco by sgbsl. */
/*     to compute  determinant(a) , follow sgbco by sgbdi. */

/*     on entry */

/*        abd     real(lda, n) */
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

/*        rcond   real */
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

/*        z       real(n) */
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

/*     linpack sgbfa */
/*     blas saxpy,sdot,sscal,sasum */
/*     fortran abs,amax1,max0,min0,sign */

/*     internal variables */



/*     compute 1-norm of a */

    /* Parameter adjustments */
    abd_dim1 = *lda;
    abd_offset = 1 + abd_dim1;
    abd -= abd_offset;
    --ipvt;
    --z__;

    /* Function Body */
    anorm = 0.f;
    l = *ml + 1;
    is = l + *mu;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
/* Computing MAX */
	r__1 = anorm, r__2 = sasum_(&l, &abd[is + j * abd_dim1], &c__1);
	anorm = dmax(r__1,r__2);
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

    sgbfa_(&abd[abd_offset], lda, n, ml, mu, &ipvt[1], &info);

/*     rcond = 1/(norm(a)*(estimate of norm(inverse(a)))) . */
/*     estimate = norm(z)/norm(y) where  a*z = y  and  trans(a)*y = e . */
/*     trans(a)  is the transpose of a .  the components of  e  are */
/*     chosen to cause maximum local growth in the elements of w  where */
/*     trans(u)*w = e .  the vectors are frequently rescaled to avoid */
/*     overflow. */

/*     solve trans(u)*w = e */

    ek = 1.f;
    i__1 = *n;
    for (j = 1; j <= i__1; ++j) {
	z__[j] = 0.f;
/* L20: */
    }
    m = *ml + *mu + 1;
    ju = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	if (z__[k] != 0.f) {
	    r__1 = -z__[k];
	    ek = r_sign(&ek, &r__1);
	}
	if ((r__1 = ek - z__[k], dabs(r__1)) <= (r__2 = abd[m + k * abd_dim1],
		 dabs(r__2))) {
	    goto L30;
	}
	s = (r__1 = abd[m + k * abd_dim1], dabs(r__1)) / (r__2 = ek - z__[k], 
		dabs(r__2));
	sscal_(n, &s, &z__[1], &c__1);
	ek = s * ek;
L30:
	wk = ek - z__[k];
	wkm = -ek - z__[k];
	s = dabs(wk);
	sm = dabs(wkm);
	if (abd[m + k * abd_dim1] == 0.f) {
	    goto L40;
	}
	wk /= abd[m + k * abd_dim1];
	wkm /= abd[m + k * abd_dim1];
	goto L50;
L40:
	wk = 1.f;
	wkm = 1.f;
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
	    sm += (r__1 = z__[j] + wkm * abd[mm + j * abd_dim1], dabs(r__1));
	    z__[j] += wk * abd[mm + j * abd_dim1];
	    s += (r__1 = z__[j], dabs(r__1));
/* L60: */
	}
	if (s >= sm) {
	    goto L80;
	}
	t = wkm - wk;
	wk = wkm;
	mm = m;
	i__2 = ju;
	for (j = kp1; j <= i__2; ++j) {
	    --mm;
	    z__[j] += t * abd[mm + j * abd_dim1];
/* L70: */
	}
L80:
L90:
	z__[k] = wk;
/* L100: */
    }
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

/*     solve trans(l)*y = w */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	lm = min(i__2,i__3);
	if (k < *n) {
	    z__[k] += sdot_(&lm, &abd[m + 1 + k * abd_dim1], &c__1, &z__[k + 
		    1], &c__1);
	}
	if ((r__1 = z__[k], dabs(r__1)) <= 1.f) {
	    goto L110;
	}
	s = 1.f / (r__1 = z__[k], dabs(r__1));
	sscal_(n, &s, &z__[1], &c__1);
L110:
	l = ipvt[k];
	t = z__[l];
	z__[l] = z__[k];
	z__[k] = t;
/* L120: */
    }
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);

    ynorm = 1.f;

/*     solve l*v = y */

    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	l = ipvt[k];
	t = z__[l];
	z__[l] = z__[k];
	z__[k] = t;
/* Computing MIN */
	i__2 = *ml, i__3 = *n - k;
	lm = min(i__2,i__3);
	if (k < *n) {
	    saxpy_(&lm, &t, &abd[m + 1 + k * abd_dim1], &c__1, &z__[k + 1], &
		    c__1);
	}
	if ((r__1 = z__[k], dabs(r__1)) <= 1.f) {
	    goto L130;
	}
	s = 1.f / (r__1 = z__[k], dabs(r__1));
	sscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L130:
/* L140: */
	;
    }
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

/*     solve  u*z = w */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	if ((r__1 = z__[k], dabs(r__1)) <= (r__2 = abd[m + k * abd_dim1], 
		dabs(r__2))) {
	    goto L150;
	}
	s = (r__1 = abd[m + k * abd_dim1], dabs(r__1)) / (r__2 = z__[k], dabs(
		r__2));
	sscal_(n, &s, &z__[1], &c__1);
	ynorm = s * ynorm;
L150:
	if (abd[m + k * abd_dim1] != 0.f) {
	    z__[k] /= abd[m + k * abd_dim1];
	}
	if (abd[m + k * abd_dim1] == 0.f) {
	    z__[k] = 1.f;
	}
	lm = min(k,m) - 1;
	la = m - lm;
	lz = k - lm;
	t = -z__[k];
	saxpy_(&lm, &t, &abd[la + k * abd_dim1], &c__1, &z__[lz], &c__1);
/* L160: */
    }
/*     make znorm = 1.0 */
    s = 1.f / sasum_(n, &z__[1], &c__1);
    sscal_(n, &s, &z__[1], &c__1);
    ynorm = s * ynorm;

    if (anorm != 0.f) {
	*rcond = ynorm / anorm;
    }
    if (anorm == 0.f) {
	*rcond = 0.f;
    }
    return 0;
} /* sgbco_ */

