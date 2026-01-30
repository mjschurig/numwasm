/* dtrdi.f -- translated by f2c (version 20240504).
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

/* Subroutine */ int dtrdi_(doublereal *t, integer *ldt, integer *n, 
	doublereal *det, integer *job, integer *info)
{
    /* System generated locals */
    integer t_dim1, t_offset, i__1, i__2, i__3;

    /* Local variables */
    integer i__, j, k, kb, km1, kp1;
    doublereal ten, temp;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), daxpy_(integer *, doublereal *, doublereal *, integer 
	    *, doublereal *, integer *);


/*     dtrdi computes the determinant and inverse of a double precision */
/*     triangular matrix. */

/*     on entry */

/*        t       double precision(ldt,n) */
/*                t contains the triangular matrix. the zero */
/*                elements of the matrix are not referenced, and */
/*                the corresponding elements of the array can be */
/*                used to store other information. */

/*        ldt     integer */
/*                ldt is the leading dimension of the array t. */

/*        n       integer */
/*                n is the order of the system. */

/*        job     integer */
/*                = 010       no det, inverse of lower triangular. */
/*                = 011       no det, inverse of upper triangular. */
/*                = 100       det, no inverse. */
/*                = 110       det, inverse of lower triangular. */
/*                = 111       det, inverse of upper triangular. */

/*     on return */

/*        t       inverse of original matrix if requested. */
/*                otherwise unchanged. */

/*        det     double precision(2) */
/*                determinant of original matrix if requested. */
/*                otherwise not referenced. */
/*                determinant = det(1) * 10.0**det(2) */
/*                with  1.0 .le. dabs(det(1)) .lt. 10.0 */
/*                or  det(1) .eq. 0.0 . */

/*        info    integer */
/*                info contains zero if the system is nonsingular */
/*                and the inverse is requested. */
/*                otherwise info contains the index of */
/*                a zero diagonal element of t. */


/*     linpack. this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas daxpy,dscal */
/*     fortran dabs,mod */

/*     internal variables */


/*     begin block permitting ...exits to 180 */

/*        compute determinant */

    /* Parameter adjustments */
    t_dim1 = *ldt;
    t_offset = 1 + t_dim1;
    t -= t_offset;
    --det;

    /* Function Body */
    if (*job / 100 == 0) {
	goto L70;
    }
    det[1] = 1.;
    det[2] = 0.;
    ten = 10.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	det[1] = t[i__ + i__ * t_dim1] * det[1];
/*           ...exit */
	if (det[1] == 0.) {
	    goto L60;
	}
L10:
	if (abs(det[1]) >= 1.) {
	    goto L20;
	}
	det[1] = ten * det[1];
	det[2] += -1.;
	goto L10;
L20:
L30:
	if (abs(det[1]) < ten) {
	    goto L40;
	}
	det[1] /= ten;
	det[2] += 1.;
	goto L30;
L40:
/* L50: */
	;
    }
L60:
L70:

/*        compute inverse of upper triangular */

    if (*job / 10 % 10 == 0) {
	goto L170;
    }
    if (*job % 10 == 0) {
	goto L120;
    }
/*              begin block permitting ...exits to 110 */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	*info = k;
/*              ......exit */
	if (t[k + k * t_dim1] == 0.) {
	    goto L110;
	}
	t[k + k * t_dim1] = 1. / t[k + k * t_dim1];
	temp = -t[k + k * t_dim1];
	i__2 = k - 1;
	dscal_(&i__2, &temp, &t[k * t_dim1 + 1], &c__1);
	kp1 = k + 1;
	if (*n < kp1) {
	    goto L90;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    temp = t[k + j * t_dim1];
	    t[k + j * t_dim1] = 0.;
	    daxpy_(&k, &temp, &t[k * t_dim1 + 1], &c__1, &t[j * t_dim1 + 1], &
		    c__1);
/* L80: */
	}
L90:
/* L100: */
	;
    }
    *info = 0;
L110:
    goto L160;
L120:

/*              compute inverse of lower triangular */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	*info = k;
/*     ............exit */
	if (t[k + k * t_dim1] == 0.) {
	    goto L180;
	}
	t[k + k * t_dim1] = 1. / t[k + k * t_dim1];
	temp = -t[k + k * t_dim1];
	if (k != *n) {
	    i__2 = *n - k;
	    dscal_(&i__2, &temp, &t[k + 1 + k * t_dim1], &c__1);
	}
	km1 = k - 1;
	if (km1 < 1) {
	    goto L140;
	}
	i__2 = km1;
	for (j = 1; j <= i__2; ++j) {
	    temp = t[k + j * t_dim1];
	    t[k + j * t_dim1] = 0.;
	    i__3 = *n - k + 1;
	    daxpy_(&i__3, &temp, &t[k + k * t_dim1], &c__1, &t[k + j * t_dim1]
		    , &c__1);
/* L130: */
	}
L140:
/* L150: */
	;
    }
    *info = 0;
L160:
L170:
L180:
    return 0;
} /* dtrdi_ */

