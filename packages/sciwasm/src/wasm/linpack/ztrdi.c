/* ztrdi.f -- translated by f2c (version 20240504).
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

static doublecomplex c_b3 = {1.,0.};
static integer c__1 = 1;

/* Subroutine */ int ztrdi_(doublecomplex *t, integer *ldt, integer *n, 
	doublecomplex *det, integer *job, integer *info)
{
    /* System generated locals */
    integer t_dim1, t_offset, i__1, i__2, i__3;
    doublereal d__1, d__2;
    doublecomplex z__1, z__2;

    /* Builtin functions */
    void z_div(doublecomplex *, doublecomplex *, doublecomplex *);

    /* Local variables */
    integer i__, j, k, kb, km1, kp1;
    doublereal ten;
    doublecomplex temp;
    extern /* Subroutine */ int zscal_(integer *, doublecomplex *, 
	    doublecomplex *, integer *), zaxpy_(integer *, doublecomplex *, 
	    doublecomplex *, integer *, doublecomplex *, integer *);


/*     ztrdi computes the determinant and inverse of a complex*16 */
/*     triangular matrix. */

/*     on entry */

/*        t       complex*16(ldt,n) */
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

/*        det     complex*16(2) */
/*                determinant of original matrix if requested. */
/*                otherwise not referenced. */
/*                determinant = det(1) * 10.0**det(2) */
/*                with  1.0 .le. cabs1(det(1)) .lt. 10.0 */
/*                or  det(1) .eq. 0.0 . */

/*        info    integer */
/*                info contains zero if the system is nonsingular */
/*                and the inverse is requested. */
/*                otherwise info contains the index of */
/*                a zero diagonal element of t. */


/*     linpack. this version dated 08/14/78 . */
/*     cleve moler, university of new mexico, argonne national lab. */

/*     subroutines and functions */

/*     blas zaxpy,zscal */
/*     fortran dabs,dcmplx,mod */

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
    det[1].r = 1., det[1].i = 0.;
    det[2].r = 0., det[2].i = 0.;
    ten = 10.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	i__2 = i__ + i__ * t_dim1;
	z__1.r = t[i__2].r * det[1].r - t[i__2].i * det[1].i, z__1.i = t[i__2]
		.r * det[1].i + t[i__2].i * det[1].r;
	det[1].r = z__1.r, det[1].i = z__1.i;
/*           ...exit */
	z__1.r = det[1].r * 0. - det[1].i * -1., z__1.i = det[1].i * 0. + det[
		1].r * -1.;
	if ((d__1 = det[1].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) == 0.) {
	    goto L60;
	}
L10:
	z__1.r = det[1].r * 0. - det[1].i * -1., z__1.i = det[1].i * 0. + det[
		1].r * -1.;
	if ((d__1 = det[1].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) >= 1.) {
	    goto L20;
	}
	z__2.r = ten, z__2.i = 0.;
	z__1.r = z__2.r * det[1].r - z__2.i * det[1].i, z__1.i = z__2.r * det[
		1].i + z__2.i * det[1].r;
	det[1].r = z__1.r, det[1].i = z__1.i;
	z__1.r = det[2].r - 1., z__1.i = det[2].i - 0.;
	det[2].r = z__1.r, det[2].i = z__1.i;
	goto L10;
L20:
L30:
	z__1.r = det[1].r * 0. - det[1].i * -1., z__1.i = det[1].i * 0. + det[
		1].r * -1.;
	if ((d__1 = det[1].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) < ten) {
	    goto L40;
	}
	z__2.r = ten, z__2.i = 0.;
	z_div(&z__1, &det[1], &z__2);
	det[1].r = z__1.r, det[1].i = z__1.i;
	z__1.r = det[2].r + 1., z__1.i = det[2].i + 0.;
	det[2].r = z__1.r, det[2].i = z__1.i;
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
	i__2 = k + k * t_dim1;
	i__3 = k + k * t_dim1;
	z__1.r = t[i__3].r * 0. - t[i__3].i * -1., z__1.i = t[i__3].i * 0. + 
		t[i__3].r * -1.;
	if ((d__1 = t[i__2].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) == 0.) 
		{
	    goto L110;
	}
	i__2 = k + k * t_dim1;
	z_div(&z__1, &c_b3, &t[k + k * t_dim1]);
	t[i__2].r = z__1.r, t[i__2].i = z__1.i;
	i__2 = k + k * t_dim1;
	z__1.r = -t[i__2].r, z__1.i = -t[i__2].i;
	temp.r = z__1.r, temp.i = z__1.i;
	i__2 = k - 1;
	zscal_(&i__2, &temp, &t[k * t_dim1 + 1], &c__1);
	kp1 = k + 1;
	if (*n < kp1) {
	    goto L90;
	}
	i__2 = *n;
	for (j = kp1; j <= i__2; ++j) {
	    i__3 = k + j * t_dim1;
	    temp.r = t[i__3].r, temp.i = t[i__3].i;
	    i__3 = k + j * t_dim1;
	    t[i__3].r = 0., t[i__3].i = 0.;
	    zaxpy_(&k, &temp, &t[k * t_dim1 + 1], &c__1, &t[j * t_dim1 + 1], &
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
	i__2 = k + k * t_dim1;
	i__3 = k + k * t_dim1;
	z__1.r = t[i__3].r * 0. - t[i__3].i * -1., z__1.i = t[i__3].i * 0. + 
		t[i__3].r * -1.;
	if ((d__1 = t[i__2].r, abs(d__1)) + (d__2 = z__1.r, abs(d__2)) == 0.) 
		{
	    goto L180;
	}
	i__2 = k + k * t_dim1;
	z_div(&z__1, &c_b3, &t[k + k * t_dim1]);
	t[i__2].r = z__1.r, t[i__2].i = z__1.i;
	i__2 = k + k * t_dim1;
	z__1.r = -t[i__2].r, z__1.i = -t[i__2].i;
	temp.r = z__1.r, temp.i = z__1.i;
	if (k != *n) {
	    i__2 = *n - k;
	    zscal_(&i__2, &temp, &t[k + 1 + k * t_dim1], &c__1);
	}
	km1 = k - 1;
	if (km1 < 1) {
	    goto L140;
	}
	i__2 = km1;
	for (j = 1; j <= i__2; ++j) {
	    i__3 = k + j * t_dim1;
	    temp.r = t[i__3].r, temp.i = t[i__3].i;
	    i__3 = k + j * t_dim1;
	    t[i__3].r = 0., t[i__3].i = 0.;
	    i__3 = *n - k + 1;
	    zaxpy_(&i__3, &temp, &t[k + k * t_dim1], &c__1, &t[k + j * t_dim1]
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
} /* ztrdi_ */

