/* ../reference/netlib/vodpk.f -- translated by f2c (version 20240504).
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

#ifdef VODE_COMMON_DEFINED
/* Common blocks and shared functions are defined in vode_common.c */
#include "vode_common.h"
/* Note: VODPK uses dvod02_3.npe (defined in vode_common.h) instead of dvod02_1.nje */
#else
/* Common Block Declarations (when compiling standalone) */

union {
    struct {
	doublereal acnrm, ccmxj, conp, crate, drc, el[13], eta, etamax, h__,
		hmin, hmxi, hnew, hscal, prl1, rc, rl1, tau[13], tq[5], tn,
		uround;
	integer icf, init, ipup, jcur, jstart, jsv, kflag, kuth, l, lmax, lyh,
		 lewt, lacor, lsavf, lwm, liwm, locjs, maxord, meth, miter,
		msbj, mxhnil, mxstep, n, newh, newq, nhnil, nq, nqnyh, nqwait,
		 nslj, nslp, nyh;
    } _1;
    struct {
	doublereal rvod1[48];
	integer ivod1[33];
    } _2;
} dvod01_;

#define dvod01_1 (dvod01_._1)
#define dvod01_2 (dvod01_._2)

union {
    struct {
	doublereal hu;
	integer ncfn, netf, nfe, npe, nlu, nni, nqu, nst;
    } _1;
    struct {
	doublereal rvod2[1];
	integer ivod2[8];
    } _2;
    struct {
	doublereal hu;
	integer ncfn, netf, nfe, nje, nlu, nni, nqu, nst;
    } _3;
} dvod02_;

#define dvod02_1 (dvod02_._1)
#define dvod02_2 (dvod02_._2)
#define dvod02_3 (dvod02_._3)
#endif

/* VODPK-specific common block */
union {
    struct {
	doublereal delt, sqrtn, rsqrtn;
	integer jpre, jacflg, lociwp, locwp, lvsav, kmp, maxl, mnewt, nli, 
		nps, ncfl;
    } _1;
    struct {
	doublereal rvpk1[3];
	integer ivpk1[11];
    } _2;
} dvpk01_;

#define dvpk01_1 (dvpk01_._1)
#define dvpk01_2 (dvpk01_._2)

/* Table of constant values */

static integer c__1 = 1;
static integer c__0 = 0;
static integer c__56 = 56;
static integer c__111 = 111;
static integer c__2 = 2;
static integer c__112 = 112;
static integer c__113 = 113;
static integer c__50 = 50;
static integer c__101 = 101;
static integer c__60 = 60;
static integer c__102 = 102;
static integer c__201 = 201;
static integer c__202 = 202;
static integer c__203 = 203;
static integer c__204 = 204;
static integer c__205 = 205;
static integer c__30 = 30;
static integer c__206 = 206;
static integer c__40 = 40;
static integer c__207 = 207;
static integer c__3 = 3;
static integer c__4 = 4;
static integer c__5 = 5;
static integer c__6 = 6;
static integer c__7 = 7;
static integer c__8 = 8;
static integer c__11 = 11;
static integer c__12 = 12;
static integer c__13 = 13;
static integer c__14 = 14;
static integer c__15 = 15;
static integer c__16 = 16;
static integer c__17 = 17;
static integer c__18 = 18;
static integer c__19 = 19;
static integer c__20 = 20;
static integer c__21 = 21;
static integer c__22 = 22;
static integer c__23 = 23;
static integer c__24 = 24;
static integer c__25 = 25;
static integer c__26 = 26;
static integer c__27 = 27;
static integer c__303 = 303;
static integer c__51 = 51;
static integer c__52 = 52;
static integer c_n1 = -1;
static logical c_false = FALSE_;
static logical c_true = TRUE_;
static doublereal c_b811 = 1.;

/* DECK DVODPK */
/* Subroutine */ int dvodpk_(S_fp f, integer *neq, doublereal *y, doublereal *
	t, doublereal *tout, integer *itol, doublereal *rtol, doublereal *
	atol, integer *itask, integer *istate, integer *iopt, doublereal *
	rwork, integer *lrw, integer *iwork, integer *liw, U_fp jac, U_fp 
	psol, integer *mf, doublereal *rpar, integer *ipar)
{
    /* Initialized data */

    static integer mord[2] = { 12,5 };
    static doublereal pt9 = .9;
    static doublereal hun = 100.;
    static integer mxstp0 = 500;
    static integer mxhnl0 = 10;
    static doublereal zero = 0.;
    static doublereal one = 1.;
    static doublereal two = 2.;
    static doublereal four = 4.;
    static doublereal pt05 = .05;
    static doublereal pt2 = .2;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double sqrt(doublereal);
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    integer i__;
    doublereal h0, rh, tp;
    integer lf0;
    doublereal big;
    integer ier, kgo;
    char msg[80];
    doublereal hmx;
    integer lwp, nli0, nni0;
    logical lcfl, lcfn, lavd;
    doublereal rcfl, rcfn;
    integer nnid;
    logical ihit;
    doublereal hmax, ewti;
    integer nstd;
    doublereal size;
    integer liwp, ncfl0, ncfn0, iflag;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    doublereal avdim, atoli;
    extern /* Subroutine */ int dvhin_(integer *, doublereal *, doublereal *, 
	    doublereal *, S_fp, doublereal *, integer *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, doublereal *,
	     doublereal *, doublereal *, integer *, integer *);
    integer leniw;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    integer lenwk, niter, lenwm;
    logical lwarn;
    integer imxer;
    doublereal tcrit;
    integer nwarn;
    doublereal tolsf;
    integer lenrw;
    doublereal rtoli, tnext;
    extern doublereal dumach_(void);
    integer leniwk;
    extern /* Subroutine */ int dewset_(integer *, integer *, doublereal *, 
	    doublereal *, doublereal *, doublereal *), dvindy_(doublereal *, 
	    integer *, doublereal *, integer *, doublereal *, integer *);
    extern /* Subroutine */ int dvnlsk_();
    integer nslast;
    extern doublereal dvnorm_(integer *, doublereal *, doublereal *);
    extern /* Subroutine */ int dvstep_(doublereal *, doublereal *, integer *,
	     doublereal *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, S_fp, U_fp, U_fp, U_fp, 
	    doublereal *, integer *), xerrwd_(char *, integer *, integer *, 
	    integer *, integer *, integer *, integer *, integer *, doublereal 
	    *, doublereal *, ftnlen);

/* ----------------------------------------------------------------------- */
/* This is the 26 April 2002 version of */
/* DVODPK: Variable-coefficient Ordinary Differential equation solver */
/*         with the Preconditioned Krylov method GMRES for the solution */
/*         of linear systems. */

/* This version is in double precision. */

/* DVODPK solves the initial value problem for stiff or nonstiff */
/* systems of first order ODEs, */
/*     dy/dt = f(t,y) ,  or, in component form, */
/*     dy(i)/dt = f(i) = f(i,t,y(1),y(2),...,y(neq)) (i = 1,...,neq). */
/* DVODPK is a package based on the VODE and LSODPK packages, and on */
/* the October 23, 1978 version of the ODEPACK user interface standard, */
/* with minor modifications. */
/* ----------------------------------------------------------------------- */
/* Authors: */
/*               Alan C. Hindmarsh and Peter N. Brown */
/*               Center for Applied Scientific Computing, L-561 */
/*               Lawrence Livermore National Laboratory */
/*               Livermore, CA 94551 */
/* and */
/*               George D. Byrne */
/*               Illinois Institute of Technology */
/*               Chicago, IL 60616 */
/* ----------------------------------------------------------------------- */
/* References: */
/* 1. P. N. Brown, G. D. Byrne, and A. C. Hindmarsh, "VODE, A Variable- */
/*    Coefficient ODE Solver," SIAM J. Sci. Stat. Comput., 10  (1989), */
/*    pp., 1038-1051.  Also LLNL report UCRL-98412, June 1988. */
/* 2. P. N. Brown and A. C. Hindmarsh, "Reduced Storage Matrix Methods */
/*    in Stiff ODE Systems," J. Appl. Math. & Comp., 31 (1989), pp.40-91. */
/*    Also LLNL report UCRL-95088, Rev. 1, June 1987. */
/* 3. G. D. Byrne, "Pragmatic Experiments with Krylov Methods in the */
/*    Stiff ODE Setting," Computational Ordinary Differential Equations, */
/*    J. Cash and I. Gladwell, eds., Oxford Univ. Press, Oxford, 1992, */
/*    pp. 323-356. */
/* ----------------------------------------------------------------------- */
/* Introduction. */

/* This is a modification of the VODE package which incorporates */
/* the preconditioned Krylov subspace iterative method SPIGMR for the */
/* linear algebraic systems that arise in the case of stiff systems. */
/* SPIGMR denotes a scaled preconditioned incomplete version of the */
/* GMRES (Generalized Minimum Residual) method. */

/* The linear systems that are solved have the form */
/*   A * x  = b ,  where  A = I - hrl1 * (df/dy) . */
/* here hrl1 is a scalar, I is the identity matrix, and df/dy is the */
/* Jacobian matrix of partial derivatives of f with respect to y */
/* (an NEQ by NEQ matrix). */

/* The particular Krylov method is chosen by setting the second digit, */
/* MITER, in the method flag MF. */
/* Currently, the values of MITER have the following meanings: */

/*          1 means SPIGMR, a scaled, preconditioned, incomplete version */
/*            of GMRES, a generalized minimum residual method. */
/*            This is the best choice in general. */

/*          9 means that only a user-supplied matrix P (approximating A) */
/*            will be used, with no Krylov iteration done internally to */
/*            DVODPK.  This option allows the user to provide the */
/*            complete linear system solution algorithm, if desired. */

/* The user can apply preconditioning to the linear system A*x = b, */
/* by means of arbitrary matrices (the preconditioners). */

/*     In the case of SPIGMR, one can apply left and right */
/* preconditioners P1 and P2, and the basic iterative method is then */
/* applied to the matrix (P1-inverse)*A*(P2-inverse) instead of to the */
/* matrix A.  The product P1*P2 should be an approximation to A */
/* such that linear systems with P1 or P2 are easier to solve than with */
/* A alone.  Preconditioning from the left only or right only means using */
/* P2 = I  or  P1 = I, respectively. */

/*     If the Jacobian  J = df/dy  splits in a natural way into a sum */
/* J = J1 + J2, then one possible choice of preconditioners is */
/*            P1 = I - hrl1 * J1  and  P2 = I - hrl1 * J2 */
/* provided each of these is easy to solve (or to approximately solve). */

/* NOTE:  To achieve an efficient solution, the preconditioned Krylov */
/* methods in DVODPK generally require a thoughtful choice of */
/* preconditioners.  If the ODE system produces linear systems that are */
/* not amenable to solution by such iterative methods, the cost can be */
/* higher than with a solver that uses sparse direct methods.  However, */
/* for many systems, careful use of DVODPK can be highly effective. */

/* See Ref. 2 for more details on the methods and applications. */
/* ----------------------------------------------------------------------- */
/* Summary of usage. */

/* Communication between the user and the DVODPK package, for normal */
/* situations, is summarized here.  This summary describes only a subset */
/* of the full set of options available.  See full description (below) */
/* for details, including optional communication, nonstandard options, */
/* and instructions for special situations.  See also the example */
/* program embedded in the comments below. */

/* A. First provide a subroutine of the form */
/*     SUBROUTINE F (NEQ, T, Y, YDOT, RPAR, IPAR) */
/*     DOUBLE PRECISION T, Y(NEQ), YDOT(NEQ), RPAR(*) */
/*     INTEGER IPAR(*) */
/* which supplies the vector function f by loading YDOT(i) with f(i). */

/* B. Next determine (or guess) whether or not the problem is stiff. */
/* Stiffness occurs when the Jacobian matrix df/dy has an eigenvalue */
/* whose real part is negative and large in magnitude, compared to the */
/* reciprocal of the t span of interest.  If the problem is nonstiff, */
/* use method flag MF = 10.  If it is stiff, MF should be 21. */

/* The following four parameters must also be set. */
/*  IWORK(1) = LWP  = length of real array WP for preconditioning. */
/*  IWORK(2) = LIWP = length of integer array IWP for preconditioning. */
/*  IWORK(3) = JPRE = preconditioner type flag: */
/*                  = 0 for no preconditioning (P1 = P2 = I) */
/*                  = 1 for left-only preconditioning (P2 = I) */
/*                  = 2 for right-only preconditioning (P1 = I) */
/*                  = 3 for two-sided preconditioning */
/*  IWORK(4) = JACFLG = flag for whether JAC is called. */
/*                    = 0 if JAC is not to be called, */
/*                    = 1 if JAC is to be called. */
/*  Use JACFLG = 1 if JAC computes any nonconstant data for use in */
/*  preconditioning, such as Jacobian elements.  See next paragraph. */
/*  The arrays WP and IWP are work arrays under the user's control, */
/*  for use in the routines that perform preconditioning operations. */

/* C. If the problem is stiff, you must supply two routines that deal */
/* with the preconditioning of the linear systems to be solved. */
/* These are as follows: */

/*     SUBROUTINE JAC (F, NEQ, T, Y, YSV, REWT, FTY, V, HRL1, WP, IWP, */
/*    1                IER, RPAR, IPAR) */
/*     DOUBLE PRECISION T, Y(NEQ), YSV(NEQ), REWT(NEQ), FTY(NEQ), V(NEQ), */
/*    1                 HRL1, WP(*), RPAR(*) */
/*     INTEGER IWP(*), IPAR(*) */

/*        This routine is optional, and is to evaluate and preprocess */
/*     any parts of the Jacobian matrix df/dy involved in the */
/*     preconditioners P1 and P2. */
/*     The Y and FTY arrays contain the current values of y and f(t,y), */
/*     respectively, and YSV also contains the current value of y. */
/*     The array V is work space of length NEQ. */
/*     JAC must multiply all computed Jacobian elements by the scalar */
/*     -hrl1, add the identity matrix I, and do any factorization */
/*     operations called for, in preparation for solving linear systems */
/*     with a coefficient matrix of P1 or P2.  The matrix P1*P2 should */
/*     be an approximation to  I - hrl1 * (df/dy), where hrl1 is a */
/*     scalar stored in HRL1. */
/*     JAC should return IER = 0 if successful, and IER .ne. 0 if not. */
/*     (If IER .ne. 0, a smaller time step will be tried.) */

/*     SUBROUTINE PSOL (NEQ, T, Y, FTY, WK, HRL1, WP, IWP, B, LR, IER, */
/*    1                 RPAR, IPAR) */
/*     DOUBLE PRECISION T, Y(NEQ), FTY(NEQ), WK(NEQ), HRL1, WP(*), */
/*    1                 B(NEQ), RPAR(*) */
/*     INTEGER IWP(*), IPAR(*) */

/*        This routine must solve a linear system with b (stored in B) */
/*     as right-hand side and one of the preconditioning matrices, P1 or */
/*     P2, as coefficient matrix, and return the solution vector in B. */
/*     LR is a flag concerning left vs. right preconditioning, input */
/*     to PSOL.  PSOL is to use P1 if LR = 1, and P2 if LR = 2. */

/*        PSOL can use data generated in the JAC routine and stored in */
/*     WP and IWP.  WK is a work array of length NEQ. */
/*     The argument HRL1 is the current value of the scalar appearing */
/*     in the linear system.  If the old value, at the time of the last */
/*     JAC call, is needed, it must have been saved by JAC in WP. */
/*     on return, PSOL should set the error flag  IER as follows: */
/*        IER = 0 if PSOL was successful, */
/*        IER .gt. 0 if a recoverable error occurred, meaning that the */
/*              time step will be retried, */
/*        IER .lt. 0 if an unrecoverable error occurred, meaning that the */
/*              solver is to stop immediately. */

/* D. Write a main program which calls subroutine DVODPK once for */
/* each point at which answers are desired.  This should also provide */
/* for possible use of logical unit 6 for output of error messages */
/* by DVODPK.  on the first call to DVODPK, supply arguments as follows: */
/* F      = name of subroutine for right-hand side vector f. */
/*          This name must be declared EXTERNAL in calling program. */
/* NEQ    = number of first order ODEs. */
/* Y      = array of initial values, of length NEQ. */
/* T      = the initial value of the independent variable. */
/* TOUT   = first point where output is desired (.ne. T). */
/* ITOL   = 1 or 2 according as ATOL (below) is a scalar or array. */
/* RTOL   = relative tolerance parameter (scalar). */
/* ATOL   = absolute tolerance parameter (scalar or array). */
/*          The estimated local error in Y(i) will be controlled so as */
/*          to be roughly less (in magnitude) than */
/*             EWT(i) = RTOL*abs(Y(i)) + ATOL     if ITOL = 1, or */
/*             EWT(i) = RTOL*abs(Y(i)) + ATOL(i)  if ITOL = 2. */
/*          Thus the local error test passes if, in each component, */
/*          either the absolute error is less than ATOL (or ATOL(i)), */
/*          or the relative error is less than RTOL. */
/*          Use RTOL = 0.0 for pure absolute error control, and */
/*          use ATOL = 0.0 (or ATOL(i) = 0.0) for pure relative error */
/*          control.  Caution: Actual (global) errors may exceed these */
/*          local tolerances, so choose them conservatively. */
/* ITASK  = 1 for normal computation of output values of Y at t = TOUT. */
/* ISTATE = integer flag (input and output).  Set ISTATE = 1. */
/* IOPT   = 0 to indicate no optional input used. */
/* RWORK  = real work array of length at least: */
/*             20 + 16*NEQ           for MF = 10, */
/*             61 + 17*NEQ + LWP     for MF = 21. */
/* LRW    = declared length of RWORK (in user's DIMENSION statement). */
/* IWORK  = integer work array of length at least: */
/*             30            for MF = 10, */
/*             30 + LIWP     for MF = 21. */
/* LIW    = declared length of IWORK (in user's DIMENSION statement). */
/* JAC,PSOL = names of subroutines for preconditioning.  These names */
/*            must be declared EXTERNAL in the user's calling program. */
/* MF     = method flag.  Standard values are: */
/*          10 for nonstiff (Adams) method. */
/*          21 for stiff (BDF) method, with SPIGMR. */

/* RPAR, IPAR  User-specified arrays used to communicate real and integer */
/*             parameters (respectively) to user-supplied subroutines. */
/*             to user-supplied subroutines.  If RPAR is a vector, then */
/*             it must be dimensioned in the user's main program.  If it */
/*             is unused or a scalar, then it need not be dimensioned. */

/* IPAR     User-specified array used to communicate integer parameter */
/*          to user-supplied subroutines.  The comments on dimensioning */
/*          RPAR apply to IPAR. */

/* Note that the user's main (calling) program must declare arrays */
/* Y, RWORK, IWORK, and possibly ATOL, RPAR, and IPAR. */

/* E. The output from the first call (or any call) is: */
/*      Y = array of computed values of y(t) vector. */
/*      T = corresponding value of independent variable (normally TOUT). */
/* ISTATE = 2  if DVODPK was successful, negative otherwise. */
/*         -1 means excess work done on this call (perhaps wrong MF). */
/*         -2 means excess accuracy requested (tolerances too small). */
/*         -3 means illegal input detected (see printed message). */
/*         -4 means repeated error test failures (check all input). */
/*         -5 means repeated convergence failures (perhaps bad JAC */
/*            or PSOL routine supplied or wrong choice of MF or */
/*            tolerances, or this solver is inappropriate). */
/*         -6 means error weight became zero during problem. (Solution */
/*            component i vanished, and ATOL or ATOL(i) = 0.) */
/*         -7 means an unrecoverable error occurred in JAC or PSOL. */

/* F. To continue the integration after a successful return, simply */
/* reset TOUT and call DVODPK again.  No other parameters need be reset. */

/* ----------------------------------------------------------------------- */
/* Example problem. */
/* An ODE system is generated from the following 2-species diurnal */
/* kinetics advection-diffusion PDE system in 2 space dimensions: */

/* dc(i)/dt = Kh*(d/dx)**2 c(i) + V*dc(i)/dx + (d/dz)(Kv(z)*dc(i)/dz) */
/*                 + Ri(c1,c2,t)      for i = 1,2,   where */
/*   R1(c1,c2,t) = -q1*c1*c3 - q2*c1*c2 + 2*q3(t)*c3 + q4(t)*c2 , */
/*   R2(c1,c2,t) =  q1*c1*c3 - q2*c1*c2 - q4(t)*c2 , */
/*   Kv(z) = Kv0*exp(z/5) , */
/* Kh, V, Kv0, q1, q2, and c3 are constants, and q3(t) and q4(t) */
/* vary diurnally.   The problem is posed on the square */
/*   0 .le. x .le. 20,    30 .le. z .le. 50   (all in km), */
/* with homogeneous Neumann boundary conditions, and for time t in */
/*   0 .le. t .le. 86400 sec (1 day). */
/* The PDE system is treated by central differences on a uniform */
/* 10 x 10 mesh, with simple polynomial initial profiles. */
/* The problem is solved with DVODPK, with the BDF/GMRES method and */
/* the block-diagonal part of the Jacobian as a left preconditioner. */
/* ----------------------------------------------------------------------- */
/*      EXTERNAL FEX, JACBD, SOLBD */
/*      DOUBLE PRECISION Q1,Q2,Q3,Q4, A3,A4, OM, C3, DZ, HDCO,VDCO,HACO */
/*      COMMON /PCOM/ Q1,Q2,Q3,Q4,A3,A4,OM,C3,DZ,HDCO,VDCO,HACO,MX,MZ,MM */
/*      DOUBLE PRECISION ATOL, AVDIM, CX, CZ, DKH, DKV0, DX, FLOOR, */
/*     1     HALFDA, PI, RPAR, RTOL, RWORK, T, TOUT, TWOHR, VEL, X, Y, Z */
/*      DIMENSION Y(2,10,10), RWORK(3861), IWORK(230) */
/*      DATA DKH/4.0D-6/, VEL/0.001D0/, DKV0/1.0D-8/, HALFDA/4.32D4/, */
/*     1  PI/3.1415926535898D0/, TWOHR/7200.0D0/, RTOL/1.0D-5/, */
/*     2  FLOOR/100.0D0/, LRW/3861/, LIW/230/, MF/21/, JPRE/1/, JACFLG/1/ */

/* Load Common block of problem parameters. */
/*      MX = 10 */
/*      MZ = 10 */
/*      MM = MX*MZ */
/*      Q1 = 1.63D-16 */
/*      Q2 = 4.66D-16 */
/*      A3 = 22.62D0 */
/*      A4 = 7.601D0 */
/*      OM = PI/HALFDA */
/*      C3 = 3.7D16 */
/*      DX = 20.0D0/(MX - 1.0D0) */
/*      DZ = 20.0D0/(MZ - 1.0D0) */
/*      HDCO = DKH/DX**2 */
/*      HACO = VEL/(2.0D0*DX) */
/*      VDCO = (1.0D0/DZ**2)*DKV0 */
/* Set other input arguments. */
/*      ATOL = RTOL*FLOOR */
/*      NEQ = 2*MX*MZ */
/*      IWORK(1) = 4*MX*MZ */
/*      IWORK(2) = NEQ */
/*      IWORK(3) = JPRE */
/*      IWORK(4) = JACFLG */
/*      T = 0.0D0 */
/*      TOUT = TWOHR */
/*      ISTATE = 1 */
/* Set initial profiles. */
/*      DO 20 JZ = 1,MZ */
/*        Z = 30.0D0 + (JZ - 1.0D0)*DZ */
/*        CZ = (0.1D0*(Z - 40.0D0))**2 */
/*        CZ = 1.0D0 - CZ + 0.5D0*CZ**2 */
/*        DO 10 JX = 1,MX */
/*          X = (JX - 1.0D0)*DX */
/*          CX = (0.1D0*(X - 10.0D0))**2 */
/*          CX = 1.0D0 - CX + 0.5D0*CX**2 */
/*          Y(1,JX,JZ) = 1.0D6*CX*CZ */
/*          Y(2,JX,JZ) = 1.0D12*CX*CZ */
/* 10       CONTINUE */
/* 20     CONTINUE */

/* Loop over output points, call DVODPK, print sample solution values. */
/*      DO 70 IOUT = 1,12 */
/*        CALL DVODPK (FEX, NEQ, Y, T, TOUT, 1, RTOL, ATOL, 1, ISTATE, 0, */
/*     1            RWORK, LRW, IWORK, LIW, JACBD, SOLBD, MF, RPAR, IPAR) */
/*        WRITE(6,50) T,IWORK(11),IWORK(14),RWORK(11) */
/* 50     FORMAT(/' t =',D10.2,5X,'no. steps =',I5, */
/*     1                      '   order =',I3,'   stepsize =',D10.2) */
/*        WRITE(6,60) Y(1,1,1), Y(1,5,5), Y(1,10,10), */
/*     1              Y(2,1,1), Y(2,5,5), Y(2,10,10) */
/* 60     FORMAT('  c1 (bot.left/middle/top rt.) =',3D12.3/ */
/*     1         '  c2 (bot.left/middle/top rt.) =',3D12.3) */
/*        IF (ISTATE .NE. 2) STOP */
/*        TOUT = TOUT + TWOHR */
/* 70     CONTINUE */

/* Print final statistics. */
/*      LENRW = IWORK(17) */
/*      LENIW = IWORK(18) */
/*      NST = IWORK(11) */
/*      NFE = IWORK(12) */
/*      NPE = IWORK(13) */
/*      NPS = IWORK(24) */
/*      NNI = IWORK(20) */
/*      NLI = IWORK(23) */
/*      AVDIM = REAL(NLI)/REAL(NNI) */
/*      NCFN = IWORK(21) */
/*      NCFL = IWORK(25) */
/*      WRITE (6,80) LENRW,LENIW,NST,NFE,NPE,NPS,NNI,NLI,AVDIM,NCFN,NCFL */
/* 80   FORMAT(//' Final statistics:'/ */
/*     1 ' RWORK size =',I5,5X,' IWORK size =',I4/ */
/*     2 ' Number of steps        =',I5,5X,'Number of f evals.     =',I5/ */
/*     3 ' Number of prec. evals. =',I5,5X,'Number of prec. solves =',I5/ */
/*     4 ' Number of nonl. iters. =',I5,5X,'Number of lin. iters.  =',I5/ */
/*     5 ' Average Krylov subspace dimension (NLI/NNI)  =',F8.4/ */
/*     6 ' Number of conv. failures:  nonlinear =',I3,'  linear =',I3) */
/*      STOP */
/*      END */

/*      SUBROUTINE FEX (NEQ, T, Y, YDOT, RPAR, IPAR) */
/*      DOUBLE PRECISION T, Y, YDOT, RPAR */
/*      DIMENSION Y(2,*), YDOT(2,*) */
/*      DOUBLE PRECISION Q1,Q2,Q3,Q4, A3,A4, OM, C3, DZ, HDCO,VDCO,HACO */
/*      COMMON /PCOM/ Q1,Q2,Q3,Q4,A3,A4,OM,C3,DZ,HDCO,VDCO,HACO,MX,MZ,MM */
/*      DOUBLE PRECISION C1, C2, C1DN, C2DN, C1UP, C2UP, C1LT, C2LT, */
/*     1    C1RT, C2RT, CZDN, CZUP, HORD1, HORD2, HORAD1, HORAD2, */
/*     2    QQ1, QQ2, QQ3, QQ4, RKIN1, RKIN2, S, VERTD1, VERTD2, ZDN, ZUP */

/* Set diurnal rate coefficients. */
/*      S = SIN(OM*T) */
/*      IF (S .GT. 0.0D0) THEN */
/*        Q3 = EXP(-A3/S) */
/*        Q4 = EXP(-A4/S) */
/*      ELSE */
/*        Q3 = 0.0D0 */
/*        Q4 = 0.0D0 */
/*      ENDIF */
/* Loop over all grid points. */
/*      DO 20 JZ = 1,MZ */
/*        ZDN = 30.0D0 + (JZ - 1.5D0)*DZ */
/*        ZUP = ZDN + DZ */
/*        CZDN = VDCO*EXP(0.2D0*ZDN) */
/*        CZUP = VDCO*EXP(0.2D0*ZUP) */
/*        IBLOK0 = (JZ-1)*MX */
/*        IDN = -MX */
/*        IF (JZ .EQ. 1) IDN = MX */
/*        IUP = MX */
/*        IF (JZ .EQ. MZ) IUP = -MX */
/*        DO 10 JX = 1,MX */
/*          IBLOK = IBLOK0 + JX */
/*          C1 = Y(1,IBLOK) */
/*          C2 = Y(2,IBLOK) */
/* Set kinetic rate terms. */
/*          QQ1 = Q1*C1*C3 */
/*          QQ2 = Q2*C1*C2 */
/*          QQ3 = Q3*C3 */
/*          QQ4 = Q4*C2 */
/*          RKIN1 = -QQ1 - QQ2 + 2.0D0*QQ3 + QQ4 */
/*          RKIN2 = QQ1 - QQ2 - QQ4 */
/* Set vertical diffusion terms. */
/*          C1DN = Y(1,IBLOK+IDN) */
/*          C2DN = Y(2,IBLOK+IDN) */
/*          C1UP = Y(1,IBLOK+IUP) */
/*          C2UP = Y(2,IBLOK+IUP) */
/*          VERTD1 = CZUP*(C1UP - C1) - CZDN*(C1 - C1DN) */
/*          VERTD2 = CZUP*(C2UP - C2) - CZDN*(C2 - C2DN) */
/* Set horizontal diffusion and advection terms. */
/*          ILEFT = -1 */
/*          IF (JX .EQ. 1) ILEFT = 1 */
/*          IRIGHT = 1 */
/*          IF (JX .EQ. MX) IRIGHT = -1 */
/*          C1LT = Y(1,IBLOK+ILEFT) */
/*          C2LT = Y(2,IBLOK+ILEFT) */
/*          C1RT = Y(1,IBLOK+IRIGHT) */
/*          C2RT = Y(2,IBLOK+IRIGHT) */
/*          HORD1 = HDCO*(C1RT - 2.0D0*C1 + C1LT) */
/*          HORD2 = HDCO*(C2RT - 2.0D0*C2 + C2LT) */
/*          HORAD1 = HACO*(C1RT - C1LT) */
/*          HORAD2 = HACO*(C2RT - C2LT) */
/* Load all terms into YDOT. */
/*          YDOT(1,IBLOK) = VERTD1 + HORD1 + HORAD1 + RKIN1 */
/*          YDOT(2,IBLOK) = VERTD2 + HORD2 + HORAD2 + RKIN2 */
/* 10       CONTINUE */
/* 20     CONTINUE */
/*      RETURN */
/*      END */

/*      SUBROUTINE JACBD (F, NEQ, T, Y, YSV, REWT, F0, F1, HRL1, */
/*     1                  BD, IPBD, IER, RPAR, IPAR) */
/*      EXTERNAL F */
/*      DOUBLE PRECISION T, Y, YSV, REWT, F0, F1, HRL1, BD, RPAR */
/*      DIMENSION Y(2, *), YSV(*), REWT(*), F0(*), F1(*), BD(2, 2, *), */
/*     1          IPBD(2, *) */
/*      DOUBLE PRECISION Q1,Q2,Q3,Q4, A3,A4, OM, C3, DZ, HDCO,VDCO,HACO */
/*      COMMON /PCOM/ Q1,Q2,Q3,Q4,A3,A4,OM,C3,DZ,HDCO,VDCO,HACO,MX,MZ,MM */
/*      DOUBLE PRECISION C1, C2, CZDN, CZUP, DIAG, ZDN, ZUP */

/* Compute diagonal Jacobian blocks, multiplied by -HRL1 */
/*   (using q3 and q4 values computed on last F call). */
/*      DO 20 JZ = 1,MZ */
/*        ZDN = 30.0D0 + (JZ - 1.5D0)*DZ */
/*        ZUP = ZDN + DZ */
/*        CZDN = VDCO*EXP(0.2D0*ZDN) */
/*        CZUP = VDCO*EXP(0.2D0*ZUP) */
/*        DIAG = -(CZDN + CZUP + 2.0D0*HDCO) */
/*        IBLOK0 = (JZ-1)*MX */
/*        DO 10 JX = 1,MX */
/*          IBLOK = IBLOK0 + JX */
/*          C1 = Y(1,IBLOK) */
/*          C2 = Y(2,IBLOK) */
/*          BD(1,1,IBLOK) = -HRL1*( (-Q1*C3 - Q2*C2) + DIAG ) */
/*          BD(1,2,IBLOK) = -HRL1*( -Q2*C1 + Q4 ) */
/*          BD(2,1,IBLOK) = -HRL1*( Q1*C3 - Q2*C2 ) */
/*          BD(2,2,IBLOK) = -HRL1*( (-Q2*C1 - Q4) + DIAG ) */
/* 10       CONTINUE */
/* 20     CONTINUE */
/* Add identity matrix and do LU decompositions on blocks. */
/*      DO 40 IBLOK = 1,MM */
/*        BD(1,1,IBLOK) = BD(1,1,IBLOK) + 1.0D0 */
/*        BD(2,2,IBLOK) = BD(2,2,IBLOK) + 1.0D0 */
/*        CALL DGEFA (BD(1,1,IBLOK), 2, 2, IPBD(1,IBLOK), IER) */
/*        IF (IER .NE. 0) RETURN */
/* 40     CONTINUE */
/*      RETURN */
/*      END */

/*      SUBROUTINE SOLBD (NEQ, T, Y, F0, WK, HRL1, BD, IPBD, V, LR, IER, */
/*     1                  RPAR, IPAR) */
/*      DOUBLE PRECISION T, Y, F0, WK, HRL1, BD, V, RPAR */
/*      DIMENSION BD(2,2,*), IPBD(2,*), V(2,*) */
/*      DOUBLE PRECISION Q1,Q2,Q3,Q4, A3,A4, OM, C3, DZ, HDCO,VDCO,HACO */
/*      COMMON /PCOM/ Q1,Q2,Q3,Q4,A3,A4,OM,C3,DZ,HDCO,VDCO,HACO,MX,MZ,MM */

/* Solve the block-diagonal system Px = v using LU factors stored in BD */
/* and pivot data in IPBD, and return the solution in V. */
/*      IER = 0 */
/*      DO 10 I = 1,MM */
/*        CALL DGESL (BD(1,1,I), 2, 2, IPBD(1,I), V(1,I), 0) */
/* 10     CONTINUE */
/*      RETURN */
/*      END */

/* The output of this program, on a Cray-1 in single precision, */
/* is as follows: */

/* t =  7.20e+03     no. steps =  194   order =  5   stepsize =  1.17e+02 */
/*  c1 (bot.left/middle/top rt.) =   1.047e+04   2.964e+04   1.119e+04 */
/*  c2 (bot.left/middle/top rt.) =   2.527e+11   7.154e+11   2.700e+11 */

/* t =  1.44e+04     no. steps =  227   order =  5   stepsize =  2.73e+02 */
/*  c1 (bot.left/middle/top rt.) =   6.659e+06   5.316e+06   7.301e+06 */
/*  c2 (bot.left/middle/top rt.) =   2.582e+11   2.057e+11   2.833e+11 */

/* t =  2.16e+04     no. steps =  252   order =  5   stepsize =  4.21e+02 */
/*  c1 (bot.left/middle/top rt.) =   2.665e+07   1.036e+07   2.931e+07 */
/*  c2 (bot.left/middle/top rt.) =   2.993e+11   1.028e+11   3.313e+11 */

/* t =  2.88e+04     no. steps =  291   order =  4   stepsize =  2.13e+02 */
/*  c1 (bot.left/middle/top rt.) =   8.702e+06   1.292e+07   9.650e+06 */
/*  c2 (bot.left/middle/top rt.) =   3.380e+11   5.029e+11   3.751e+11 */

/* t =  3.60e+04     no. steps =  321   order =  5   stepsize =  9.90e+01 */
/*  c1 (bot.left/middle/top rt.) =   1.404e+04   2.029e+04   1.561e+04 */
/*  c2 (bot.left/middle/top rt.) =   3.387e+11   4.894e+11   3.765e+11 */

/* t =  4.32e+04     no. steps =  374   order =  4   stepsize =  4.44e+02 */
/*  c1 (bot.left/middle/top rt.) =  -5.457e-09  -4.365e-09  -6.182e-09 */
/*  c2 (bot.left/middle/top rt.) =   3.382e+11   1.355e+11   3.804e+11 */

/* t =  5.04e+04     no. steps =  393   order =  5   stepsize =  5.22e+02 */
/*  c1 (bot.left/middle/top rt.) =   3.396e-12   2.798e-12   3.789e-12 */
/*  c2 (bot.left/middle/top rt.) =   3.358e+11   4.930e+11   3.864e+11 */

/* t =  5.76e+04     no. steps =  407   order =  5   stepsize =  3.54e+02 */
/*  c1 (bot.left/middle/top rt.) =   7.738e-12   6.455e-12   8.598e-12 */
/*  c2 (bot.left/middle/top rt.) =   3.320e+11   9.650e+11   3.909e+11 */

/* t =  6.48e+04     no. steps =  419   order =  5   stepsize =  5.90e+02 */
/*  c1 (bot.left/middle/top rt.) =  -2.018e-11  -1.680e-11  -2.243e-11 */
/*  c2 (bot.left/middle/top rt.) =   3.313e+11   8.922e+11   3.963e+11 */

/* t =  7.20e+04     no. steps =  432   order =  5   stepsize =  5.90e+02 */
/*  c1 (bot.left/middle/top rt.) =  -2.837e-11  -2.345e-11  -3.166e-11 */
/*  c2 (bot.left/middle/top rt.) =   3.330e+11   6.186e+11   4.039e+11 */

/* t =  7.92e+04     no. steps =  444   order =  5   stepsize =  5.90e+02 */
/*  c1 (bot.left/middle/top rt.) =  -4.861e-14  -4.433e-14  -5.162e-14 */
/*  c2 (bot.left/middle/top rt.) =   3.334e+11   6.669e+11   4.120e+11 */

/* t =  8.64e+04     no. steps =  456   order =  5   stepsize =  5.90e+02 */
/*  c1 (bot.left/middle/top rt.) =   2.511e-15   2.071e-15   2.802e-15 */
/*  c2 (bot.left/middle/top rt.) =   3.352e+11   9.107e+11   4.163e+11 */


/* Final statistics: */
/* RWORK size = 3861      IWORK size = 230 */
/* Number of steps        =  456     Number of f evals.     = 1317 */
/* Number of prec. evals. =   82     Number of prec. solves = 1226 */
/* Number of nonl. iters. =  571     Number of lin. iters.  =  743 */
/* Average Krylov subspace dimension (NLI/NNI)  =  1.3012 */
/* Number of conv. failures:  nonlinear =  0  linear =  0 */
/* ----------------------------------------------------------------------- */
/* Full description of user interface to DVODPK. */

/* The user interface to DVODPK consists of the following parts. */

/* i.   The call sequence to subroutine DVODPK, which is a driver */
/*      routine for the solver.  This includes descriptions of both */
/*      the call sequence arguments and of user-supplied routines. */
/*      Following these descriptions are */
/*        * a description of optional input available through the */
/*          call sequence, */
/*        * a description of optional output (in the work arrays), and */
/*        * instructions for interrupting and restarting a solution. */

/* ii.  Descriptions of other routines in the DVODPK package that may be */
/*      (optionally) called by the user.  These provide the ability to */
/*      alter error message handling, save and restore the internal */
/*      COMMON, and obtain specified derivatives of the solution y(t). */

/* iii. Descriptions of COMMON blocks to be declared in overlay */
/*      or similar environments. */

/* iv.  Description of two routines in the DVODPK package, either of */
/*      which the user may replace with the user's own version, if */
/*      desired.  These relate to the measurement of errors. */

/* ----------------------------------------------------------------------- */
/* Part i.  Call Sequence. */

/* The call sequence parameters used for input only are */
/*     F, NEQ, TOUT, ITOL, RTOL, ATOL, ITASK, IOPT, LRW, LIW, */
/*     JAC, PSOL, MF, */
/* and those used for both input and output are */
/*     Y, T, ISTATE. */
/* The work arrays RWORK and IWORK are also used for conditional and */
/* optional input and optional output.  (The term output here refers */
/* to the return from subroutine DVODPK to the user's calling program.) */

/* The legality of input parameters will be thoroughly checked on the */
/* initial call for the problem, but not checked thereafter unless a */
/* change in input parameters is flagged by ISTATE = 3 in the input. */

/* The descriptions of the call arguments are as follows. */

/* F      = The name of the user-supplied subroutine defining the */
/*          ODE system.  The system must be put in the first-order */
/*          form dy/dt = f(t,y), where f is a vector-valued function */
/*          of the scalar t and the vector y.  Subroutine F is to */
/*          compute the function f.  It is to have the form */
/*               SUBROUTINE F (NEQ, T, Y, YDOT, RPAR, IPAR) */
/*               DOUBLE PRECISION T, Y(NEQ), YDOT(NEQ), RPAR(*) */
/*               INTEGER IPAR(*) */
/*          where NEQ, T, and Y are input, and the array YDOT = f(t,y) */
/*          is output.  Y and YDOT are arrays of length NEQ. */
/*          (In the DIMENSION statement above, NEQ  can be replaced by */
/*          *  to make  Y  and  YDOT  assumed size arrays.) */
/*          Subroutine F should not alter Y or T. */
/*          F must be declared EXTERNAL in the calling program. */

/*          Subroutine F may access user-defined real and integer */
/*          work arrays RPAR and IPAR, which are to be dimensioned */
/*          in the user's calling (main) program. */

/*          If quantities computed in the F routine are needed */
/*          externally to DVODPK, an extra call to F should be made */
/*          for this purpose, for consistent and accurate results. */
/*          If only the derivative dy/dt is needed, use DVINDY instead. */

/* NEQ    = The size of the ODE system (number of first order */
/*          ordinary differential equations).  Used only for input. */
/*          NEQ may not be increased during the problem, but */
/*          can be decreased (with ISTATE = 3 in the input). */

/* Y      = A real array for the vector of dependent variables, of */
/*          length NEQ or more.  Used for both input and output on the */
/*          first call (ISTATE = 1), and only for output on other calls. */
/*          On the first call, Y must contain the vector of initial */
/*          values.  In the output, Y contains the computed solution */
/*          evaluated at T.  If desired, the Y array may be used */
/*          for other purposes between calls to the solver. */

/*          This array is passed as the Y argument in all calls to */
/*          F, JAC, and PSOL. */

/* T      = The independent variable.  In the input, T is used only on */
/*          the first call, as the initial point of the integration. */
/*          In the output, after each call, T is the value at which a */
/*          computed solution Y is evaluated (usually the same as TOUT). */
/*          On an error return, T is the farthest point reached. */

/* TOUT   = The next value of t at which a computed solution is desired. */
/*          Used only for input. */

/*          When starting the problem (ISTATE = 1), TOUT may be equal */
/*          to T for one call, then should .ne. T for the next call. */
/*          For the initial T, an input value of TOUT .ne. T is used */
/*          in order to determine the direction of the integration */
/*          (i.e. the algebraic sign of the step sizes) and the rough */
/*          scale of the problem.  Integration in either direction */
/*          (forward or backward in t) is permitted. */

/*          If ITASK = 2 or 5 (one-step modes), TOUT is ignored after */
/*          the first call (i.e. the first call with TOUT .ne. T). */
/*          Otherwise, TOUT is required on every call. */

/*          If ITASK = 1, 3, or 4, the values of TOUT need not be */
/*          monotone, but a value of TOUT which backs up is limited */
/*          to the current internal t interval, whose endpoints are */
/*          TCUR - HU and TCUR.  (See optional output, below, for */
/*          TCUR and HU.) */

/* ITOL   = An indicator for the type of error control.  See */
/*          description below under ATOL.  Used only for input. */

/* RTOL   = A relative error tolerance parameter, either a scalar or */
/*          an array of length NEQ.  See description below under ATOL. */
/*          Input only. */

/* ATOL   = An absolute error tolerance parameter, either a scalar or */
/*          an array of length NEQ.  Input only. */

/*          The input parameters ITOL, RTOL, and ATOL determine */
/*          the error control performed by the solver.  The solver will */
/*          control the vector e = (e(i)) of estimated local errors */
/*          in Y, according to an inequality of the form */
/*                      rms-norm of ( e(i)/EWT(i) )   .le.   1, */
/*          where       EWT(i) = RTOL(i)*abs(Y(i)) + ATOL(i), */
/*          and the rms-norm (root-mean-square norm) here is */
/*          rms-norm(v) = sqrt(sum v(i)**2 / NEQ).  Here EWT = (EWT(i)) */
/*          is a vector of weights which must always be positive, and */
/*          the values of RTOL and ATOL should all be non-negative. */
/*          The following table gives the types (scalar/array) of */
/*          RTOL and ATOL, and the corresponding form of EWT(i). */

/*             ITOL    RTOL       ATOL          EWT(i) */
/*              1     scalar     scalar     RTOL*ABS(Y(i)) + ATOL */
/*              2     scalar     array      RTOL*ABS(Y(i)) + ATOL(i) */
/*              3     array      scalar     RTOL(i)*ABS(Y(i)) + ATOL */
/*              4     array      array      RTOL(i)*ABS(Y(i)) + ATOL(i) */

/*          When either of these parameters is a scalar, it need not */
/*          be dimensioned in the user's calling program. */

/*          If none of the above choices (with ITOL, RTOL, and ATOL */
/*          fixed throughout the problem) is suitable, more general */
/*          error controls can be obtained by substituting */
/*          user-supplied routines for the setting of EWT and/or for */
/*          the norm calculation.  See Part iv below. */

/*          If global errors are to be estimated by making a repeated */
/*          run on the same problem with smaller tolerances, then all */
/*          components of RTOL and ATOL (i.e. of EWT) should be scaled */
/*          down uniformly. */

/* ITASK  = An index specifying the task to be performed. */
/*          Input only.  ITASK has the following values and meanings. */
/*          1  means normal computation of output values of y(t) at */
/*             t = TOUT (by overshooting and interpolating). */
/*          2  means take one step only and return. */
/*          3  means stop at the first internal mesh point at or */
/*             beyond t = TOUT and return. */
/*          4  means normal computation of output values of y(t) at */
/*             t = TOUT but without overshooting t = TCRIT. */
/*             TCRIT must be input as RWORK(1).  TCRIT may be equal to */
/*             or beyond TOUT, but not behind it in the direction of */
/*             integration.  This option is useful if the problem */
/*             has a singularity at or beyond t = TCRIT. */
/*          5  means take one step, without passing TCRIT, and return. */
/*             TCRIT must be input as RWORK(1). */

/*          Note:  If ITASK = 4 or 5 and the solver reaches TCRIT */
/*          (within roundoff), it will return T = TCRIT (exactly) to */
/*          indicate this (unless ITASK = 4 and TOUT comes before TCRIT, */
/*          in which case answers at T = TOUT are returned first). */

/* ISTATE = an index used for input and output to specify the */
/*          the state of the calculation. */

/*          In the input, the values of ISTATE are as follows. */
/*          1  means this is the first call for the problem */
/*             (initializations will be done).  See note below. */
/*          2  means this is not the first call, and the calculation */
/*             is to continue normally, with no change in any input */
/*             parameters except possibly TOUT and ITASK. */
/*             (If ITOL, RTOL, and/or ATOL are changed between calls */
/*             with ISTATE = 2, the new values will be used but not */
/*             tested for legality.) */
/*          3  means this is not the first call, and the */
/*             calculation is to continue normally, but with */
/*             a change in input parameters other than */
/*             TOUT and ITASK.  Changes are allowed in */
/*             NEQ, ITOL, RTOL, ATOL, IOPT, LRW, LIW, MF, */
/*             and any of the optional input except H0. */

/*          Note:  A preliminary call with TOUT = T is not counted */
/*          as a first call here, as no initialization or checking of */
/*          input is done.  (Such a call is sometimes useful to include */
/*          the initial conditions in the output.) */
/*          Thus the first call for which TOUT .ne. T requires */
/*          ISTATE = 1 in the input. */

/*          In the output, ISTATE has the following values and meanings. */
/*           1  means nothing was done, as TOUT was equal to T with */
/*              ISTATE = 1 in the input. */
/*           2  means the integration was performed successfully. */
/*          -1  means an excessive amount of work (more than MXSTEP */
/*              steps) was done on this call, before completing the */
/*              requested task, but the integration was otherwise */
/*              successful as far as T.  (MXSTEP is an optional input */
/*              and is normally 500.)  To continue, the user may */
/*              simply reset ISTATE to a value .gt. 1 and call again. */
/*              (The excess work step counter will be reset to 0.) */
/*              In addition, the user may increase MXSTEP to avoid */
/*              this error return.  (See optional input below.) */
/*          -2  means too much accuracy was requested for the precision */
/*              of the machine being used.  This was detected before */
/*              completing the requested task, but the integration */
/*              was successful as far as T.  To continue, the tolerance */
/*              parameters must be reset, and ISTATE must be set */
/*              to 3.  The optional output TOLSF may be used for this */
/*              purpose.  (Note: If this condition is detected before */
/*              taking any steps, then an illegal input return */
/*              (ISTATE = -3) occurs instead.) */
/*          -3  means illegal input was detected, before taking any */
/*              integration steps.  See written message for details. */
/*              Note:  If the solver detects an infinite loop of calls */
/*              to the solver with illegal input, it will cause */
/*              the run to stop. */
/*          -4  means there were repeated error test failures on */
/*              one attempted step, before completing the requested */
/*              task, but the integration was successful as far as T. */
/*              The problem may have a singularity, or the input */
/*              may be inappropriate. */
/*          -5  means there were repeated convergence test failures on */
/*              one attempted step, before completing the requested */
/*              task, but the integration was successful as far as T. */
/*              This may be caused by a poor preconditioner matrix. */
/*          -6  means EWT(i) became zero for some i during the */
/*              integration.  Pure relative error control (ATOL(i)=0.0) */
/*              was requested on a variable which has now vanished. */
/*              The integration was successful as far as T. */
/*          -7  means an unrecoverable error occurred in JAC or PSOL. */
/*              Either JAC returned IER .ne. 0, or PSOL returned */
/*              IER .lt. 0. */

/*          Note:  Since the normal output value of ISTATE is 2, */
/*          it does not need to be reset for normal continuation. */
/*          Also, since a negative input value of ISTATE will be */
/*          regarded as illegal, a negative output value requires the */
/*          user to change it, and possibly other input, before */
/*          calling the solver again. */

/* IOPT   = An integer flag to specify whether or not any optional */
/*          input is being used on this call.  Input only. */
/*          The optional input is listed separately below. */
/*          IOPT = 0 means no optional input is being used. */
/*                   Default values will be used in all cases. */
/*          IOPT = 1 means optional input is being used. */

/* RWORK  = A real working array (double precision). */
/*          The length of RWORK must be at least */
/*             20 + NYH*(MAXORD + 1) + 3*NEQ + LENK + LWP   where */
/*          NYH    = the initial value of NEQ, */
/*          MAXORD = 12 (if METH = 1) or 5 (if METH = 2) (unless a */
/*                   smaller value is given as an optional input), */
/*          LENK = length of work space for Krylov-related data: */
/*          LENK = 0                                 if MITER = 0, */
/*          LENK = NEQ*(MAXL+3+MIN(1,MAXL-KMP)) */
/*                  + (MAXL+3)*MAXL + 1              if MITER = 1, */
/*          LENK = 3*NEQ                             if MITER = 9. */
/*          LWP = length of real user work space for preconditioning. */
/*          (See JAC/PSOL.) */
/*          (See the MF description for METH and MITER.) */
/*          Thus if MAXORD etc. have default values and NEQ is constant, */
/*          this length is: */
/*             20 + 16*NEQ                    for MF = 10, */
/*             61 + 24*NEQ + LWP              for MF = 11, */
/*             20 + 19*NEQ + LWP              for MF = 19, */
/*             20 + 9*NEQ                     for MF = 20, */
/*             61 + 17*NEQ + LWP              for MF = 21, */
/*             20 + 12*NEQ + LWP              for MF = 29 */
/*          The first 20 words of RWORK are reserved for conditional */
/*          and optional input and optional output. */

/*          The following word in RWORK is a conditional input: */
/*            RWORK(1) = TCRIT = critical value of t which the solver */
/*                       is not to overshoot.  Required if ITASK is */
/*                       4 or 5, and ignored otherwise.  (See ITASK.) */

/* LRW    = The length of the array RWORK, as declared by the user. */
/*          (This will be checked by the solver.) */

/* IWORK  = An integer work array.  The length of IWORK must be at least */
/*             30        if MITER = 0  (MF = 10, 20), or */
/*             30 + LIWP  otherwise (MF = 11, 21, 19, 29). */
/*          LIWP = length of integer user work space for preconditioning. */
/*          (See conditional input list following). */

/*          The first 30 words of IWORK are reserved for conditional and */
/*          optional input and optional output. */

/*          The following 4 words in IWORK are conditional input, */
/*          required if MITER .ge. 1: */

/*          IWORK(1) = LWP  = length of real array WP for use in */
/*                     preconditioning (part of RWORK array). */
/*          IWORK(2) = LIWP = length of integer array IWP for use in */
/*                     preconditioning (part of IWORK array). */
/*                     The arrays WP and IWP are work arrays under the */
/*                     user's control, for use in the routines that */
/*                     perform preconditioning operations (JAC and PSOL). */
/*          IWORK(3) = JPRE = preconditioner type flag: */
/*                   = 0 for no preconditioning (P1 = P2 = I */
/*                   = 1 for left-only preconditioning (P2 = I) */
/*                   = 2 for right-only preconditioning (P1 = I) */
/*                   = 3 for two-sided preconditioning */
/*          IWORK(4) = JACFLG = flag for whether JAC is called. */
/*                   = 0 if JAC is not to be called, */
/*                   = 1 if JAC is to be called. */
/*                     Use JACFLG = 1 if JAC computes any nonconstant */
/*                     data needed in preconditioning operations, */
/*                     such as some of the Jacobian elements. */


/* LIW    = the length of the array IWORK, as declared by the user. */
/*          (This will be checked by the solver.) */

/* Note:  The work arrays must not be altered between calls to DVODPK */
/* for the same problem, except possibly for the conditional and */
/* optional input, and except for the last 3*NEQ words of RWORK. */
/* The latter space is used for internal scratch space, and so is */
/* available for use by the user outside DVODPK between calls, if */
/* desired (but not for use by F or JAC). */

/* JAC    = The name of the user-supplied routine (MITER = 1 or 9) to */
/*          compute the Jacobian matrix, df/dy, as a function of */
/*          the scalar t and the vector y.  It is to have the form */
/*             SUBROUTINE JAC (F, NEQ, T, Y, YSV, REWT, FTY, V, HRL1, */
/*            1                WP, IWP, IER, RPAR, IPAR) */
/*             EXTERNAL F */
/*             DOUBLE PRECISION T, Y(NEQ), YSV(NEQ), REWT(NEQ), FTY(NEQ), */
/*            1                 V(NEQ), HRL1, WP(*), RPAR(*) */
/*             INTEGER IWP(*), IPAR(*) */
/*          This routine must evaluate and preprocess any parts of the */
/*          Jacobian matrix df/dy used in the preconditioners P1, P2 . */
/*          The Y and FTY arrays contain the current values of y and */
/*          f(t,y), respectively, and YSV also contains the current */
/*          value of y.  The array V is work space of length */
/*          NEQ for use by JAC.  REWT is the array of reciprocal error */
/*          weights (1/ewt).  JAC must multiply all computed Jacobian */
/*          elements by the scalar -hrl1, add the identity matrix I and */
/*          do any factorization operations called for, in preparation */
/*          for solving linear systems with a coefficient matrix of */
/*          P1 or P2.  The matrix P1*P2 should be an approximation to */
/*          I - hrl1 * (df/dy), where hrl1 is stored in HRL1.  JAC should */
/*          return IER = 0 if successful, and IER .ne. 0 if not. */
/*          (If IER .ne. 0, a smaller time step will be tried.) */
/*          The arrays WP (of length LWP) and IWP (of length LIWP) */
/*          are for use by JAC and PSOL for work space and for storage */
/*          of data needed for the solution of the preconditioner */
/*          linear systems.  Their lengths and contents are under the */
/*          user's control. */
/*          The JAC routine may save relevant Jacobian elements (or */
/*          approximations) used in the preconditioners, along with the */
/*          value of hrl1, and use these to reconstruct preconditioner */
/*          matrices later without reevaluationg those elements. */
/*          This may be cost-effective if JAC is called with hrl1 */
/*          considerably different from its earlier value, indicating */
/*          that a corrector convergence failure has occurred because */
/*          of the change in hrl1, not because of changes in the */
/*          value of the Jacobian.  In doing this, use the saved and */
/*          current values of hrl1 to decide whether to use saved */
/*          or reevaluated elements. */
/*          JAC may alter V, but not Y, YSV, REWT, FTY, or HRL1. */
/*          JAC must be declared external in the calling program. */

/* PSOL   = the name of the user-supplied routine for the */
/*          solution of preconditioner linear systems. */
/*          It is to have the form */
/*             SUBROUTINE PSOL (NEQ, T, Y, FTY, WK, HRL1, WP, IWP, B, LR, */
/*            1                 IER, RPAR, IPAR) */
/*             DOUBLE PRECISION T, Y(NEQ), FTY(NEQ), WK(NEQ), HRL1, */
/*            1                 WP(*), B(NEQ), RPAR(*) */
/*             INTEGER  IWP(*), IPAR(*) */
/*          This routine must solve a linear system with b (stored in B) */
/*          as right-hand side and one of the preconditioning matrices, */
/*          P1 or P2, as coefficient matrix, and return the solution */
/*          vector in B.  LR is a flag concerning left vs. right */
/*          preconditioning, input to PSOL.  PSOL is to use P1 if LR = 1 */
/*          and P2 if LR = 2.  In the case MITER = 9 (no Krylov */
/*          iteration), LR will be 1 and then 2, according to JPRE, and */
/*          PSOL is to return in B the desired approximate solution to */
/*          A * x = b, where A = I - hrl1 * (df/dy).  (hrl1 is stored in */
/*          HRL1.)  PSOL can use data generated in the JAC routine and */
/*          stored in WP and IWP.  The Y and FTY arrays contain the */
/*          current values of y and f(t,y), respectively. */
/*          The array WK is work space of length NEQ for use by PSOL. */
/*          The argument HRL1 is the current value of the scalar appear- */
/*          ing in the linear system.  If the old value, as of the last */
/*          JAC call, is needed, it must have been saved by JAC in WP. */
/*          On return, PSOL should set the error flag IER as follows: */
/*            IER = 0 if PSOL was successful, */
/*            IER .gt. 0 on a recoverable error, meaning that the */
/*                   time step will be retried, */
/*            IER .lt. 0 on an unrecoverable error, meaning that the */
/*                   solver is to stop immediately. */
/*          PSOL may not alter Y, FTY, or HRL1. */
/*          PSOL must be declared external in the calling program. */

/* MF     = The method flag.  Used only for input.  The legal values of */
/*          MF are 10, 11, 19, 20, 21, 29 . */
/*          MF is a two-digit integer, MF = 10*METH + MITER . */
/*          METH indicates the basic linear multistep method: */
/*            METH = 1 means the implicit Adams method. */
/*            METH = 2 means the method based on backward */
/*                     differentiation formulas (BDF-s). */
/*          MITER indicates the corrector iteration method.  Currently, */
/*            the values of MITER have the following meanings: */

/*          0 means functional iteration is used (no Jacobian matrix */
/*            is involved). */

/*          1 means SPIGMR, a scaled, preconditioned, incomplete version */
/*            of GMRES, a generalized minimum residual method, is used. */
/*            This is the best choice in general. */

/*          9 means that only a user-supplied matrix P (approximating A) */
/*            will be used, with no Krylov iteration done internally to */
/*            DVODPK.  This option allows the user to provide the */
/*            complete linear system solution algorithm, if desired. */

/* The user can apply preconditioning to the linear system A*x = b, */
/* by means of arbitrary matrices (the preconditioners). */

/* RPAR     User-specified array used to communicate real parameters */
/*          to user-supplied subroutines.  If RPAR is a vector, then */
/*          it must be dimensioned in the user's main program.  If it */
/*          is unused or a scalar, then it need not be dimensioned. */

/* IPAR     User-specified array used to communicate integer parameter */
/*          to user-supplied subroutines.  The comments on dimensioning */
/*          RPAR apply to IPAR. */
/* ----------------------------------------------------------------------- */
/* Optional Inputs. */

/* The following is a list of the optional input provided for in the */
/* call sequence.  (See also Part ii.)  For each such input variable, */
/* this table lists its name as used in this documentation, its */
/* location in the call sequence, its meaning, and the default value. */
/* The use of any of this input requires IOPT = 1, and in that */
/* case all of this input is examined.  A value of zero for any */
/* of these optional input variables will cause the default value to be */
/* used.  Thus to use a subset of the optional input, simply preload */
/* locations 5 to 10 in RWORK and IWORK to 0.0 and 0 respectively, and */
/* then set those of interest to nonzero values. */

/* NAME    LOCATION      MEANING AND DEFAULT VALUE */

/* H0      RWORK(5)  The step size to be attempted on the first step. */
/*                   The default value is determined by the solver. */

/* HMAX    RWORK(6)  The maximum absolute step size allowed. */
/*                   The default value is infinite. */

/* HMIN    RWORK(7)  The minimum absolute step size allowed. */
/*                   The default value is 0.  (This lower bound is not */
/*                   enforced on the final step before reaching TCRIT */
/*                   when ITASK = 4 or 5.) */

/* DELT    RWORK(8)  Convergence test constant used in Krylov iteration */
/*                   algorithm.  The default value is 0.05. */

/* MAXORD  IWORK(5)  The maximum order to be allowed.  The default */
/*                   value is 12 if METH = 1, and 5 if METH = 2. */
/*                   If MAXORD exceeds the default value, it will */
/*                   be reduced to the default value. */
/*                   If MAXORD is changed during the problem, it may */
/*                   cause the current order to be reduced. */

/* MXSTEP  IWORK(6)  Maximum number of (internally defined) steps */
/*                   allowed during one call to the solver. */
/*                   The default value is 500. */

/* MXHNIL  IWORK(7)  Maximum number of messages printed (per problem) */
/*                   warning that T + H = T on a step (H = step size). */
/*                   This must be positive to result in a non-default */
/*                   value.  The default value is 10. */

/* MAXL    IWORK(8)  maximum number of iterations in the SPIGMR */
/*                   algorithm (.le. NEQ).  The default is */
/*                   MAXL = min(5, NEQ). */

/* KMP     IWORK(9)  number of vectors on which orthogonalization */
/*                   is done in the SPIGMR algorithm (.le. MAXL). */
/*                   The default is KMP = MAXL (complete GMRES method). */
/*                   See Ref. 2 for details on incomplete GMRES. */
/*                   Note:  When KMP .lt. MAXL and MITER = 1, the length */
/*                   of RWORK must be set accordingly.  See RWORK above. */
/* ----------------------------------------------------------------------- */
/* Optional Outputs. */

/* As optional additional output from DVODPK, the variables listed */
/* below are quantities related to the performance of DVODPK */
/* which are available to the user.  These are communicated by way of */
/* the work arrays, but also have internal mnemonic names as shown. */
/* Except where stated otherwise, all of this output is defined */
/* on any successful return from DVODPK, and on any return with */
/* ISTATE = -1, -2, -4, -5, -6, or -7.  On an illegal input return */
/* (ISTATE = -3), they will be unchanged from their existing values */
/* (if any), except possibly for TOLSF, LENRW, and LENIW. */
/* On any error return, output relevant to the error will be defined, */
/* as noted below. */

/* NAME    LOCATION      MEANING */

/* HU      RWORK(11) The step size in t last used (successfully). */

/* HCUR    RWORK(12) The step size to be attempted on the next step. */

/* TCUR    RWORK(13) The current value of the independent variable */
/*                   which the solver has actually reached, i.e. the */
/*                   current internal mesh point in t.  In the output, */
/*                   TCUR will always be at least as far from the */
/*                   initial value of t as the current argument T, */
/*                   but may be farther (if interpolation was done). */

/* TOLSF   RWORK(14) A tolerance scale factor, greater than 1.0, */
/*                   computed when a request for too much accuracy was */
/*                   detected (ISTATE = -3 if detected at the start of */
/*                   the problem, ISTATE = -2 otherwise).  If ITOL is */
/*                   left unaltered but RTOL and ATOL are uniformly */
/*                   scaled up by a factor of TOLSF for the next call, */
/*                   then the solver is deemed likely to succeed. */
/*                   (The user may also ignore TOLSF and alter the */
/*                   tolerance parameters in any other way appropriate.) */

/* NST     IWORK(11) The number of steps taken for the problem so far. */

/* NFE     IWORK(12) The number of f evaluations for the problem so far. */

/* NPE     IWORK(13) The number of preconditioner evaluations (JAC calls) */
/*                   so far. */

/* NQU     IWORK(14) The method order last used (successfully). */

/* NQCUR   IWORK(15) The order to be attempted on the next step. */

/* IMXER   IWORK(16) The index of the component of largest magnitude in */
/*                   the weighted local error vector ( e(i)/EWT(i) ), */
/*                   on an error return with ISTATE = -4 or -5. */

/* LENRW   IWORK(17) The length of RWORK actually required. */
/*                   This is defined on normal returns and on an illegal */
/*                   input return for insufficient storage. */

/* LENIW   IWORK(18) The length of IWORK actually required. */
/*                   This is defined on normal returns and on an illegal */
/*                   input return for insufficient storage. */

/* NNI     IWORK(20) The number of nonlinear iterations so far (each of */
/*                   which calls the Krylov iterative linear solver). */

/* NCFN    IWORK(21) The number of convergence failures of the nonlinear */
/*                   (Newton) iteration so far. */
/*                   Note: A measure of success is the overall rate of */
/*                   nonlinear convergence failures, NCFN/NST. */

/* NETF    IWORK(22) The number of error test failures of the integrator */
/*                   so far. */

/* NLI     IWORK(23) The number of linear iterations so far. */
/*                   Note: a measure of the success of SPIGMR algorithm */
/*                   is the average number of linear iterations per */
/*                   nonlinear iteration, given by NLI/NNI. */
/*                   If this is close to MAXL, MAXL may be too small. */

/* NPS     IWORK(24) The number of preconditioning solve operations */
/*                   (PSOL calls) so far. */

/* NCFL    IWORK(25) The number of convergence failures of the linear */
/*                   iteration so far. */
/*                   Note: A measure of success is the overall rate of */
/*                   linear convergence failures, NCFL/NNI. */

/* The following two arrays are segments of the RWORK array which */
/* may also be of interest to the user as optional output. */
/* For each array, the table below gives its internal name, */
/* its base address in RWORK, and its description. */

/* NAME    BASE ADDRESS      DESCRIPTION */

/* YH      21             The Nordsieck history array, of size NYH by */
/*                        (NQCUR + 1), where NYH is the initial value */
/*                        of NEQ.  For j = 0,1,...,NQCUR, column j+1 */
/*                        of YH contains HCUR**j/factorial(j) times */
/*                        the j-th derivative of the interpolating */
/*                        polynomial currently representing the */
/*                        solution, evaluated at t = TCUR. */

/* ACOR     LENRW-NEQ+1   Array of size NEQ used for the accumulated */
/*                        corrections on each step, scaled in the output */
/*                        to represent the estimated local error in Y on */
/*                        the last step.  This is the vector e in the */
/*                        description of the error control.  Defined */
/*                        only on a successful return from DVODPK. */

/* ----------------------------------------------------------------------- */
/* Interrupting and Restarting */

/* If the integration of a given problem by DVODPK is to be */
/* interrrupted and then later continued, such as when restarting */
/* an interrupted run or alternating between two or more ODE problems, */
/* the user should save, following the return from the last DVODPK call */
/* prior to the interruption, the contents of the call sequence */
/* variables and internal COMMON blocks, and later restore these */
/* values before the next DVODPK call for that problem.  To save */
/* and restore the COMMON blocks, use subroutine DVKSRC, as */
/* described below in Part ii. */

/* In addition, if non-default values for either LUN or MFLAG are */
/* desired, an extra call to XSETUN and/or XSETF should be made just */
/* before continuing the integration.  See Part ii below for details. */

/* ----------------------------------------------------------------------- */
/* Part ii.  Other Routines Callable. */

/* The following are optional calls which the user may make to */
/* gain additional capabilities in conjunction with DVODPK. */
/* (The routines XSETUN and XSETF are designed to conform to the */
/* SLATEC error handling package.) */

/*     FORM OF CALL                  FUNCTION */

/*  CALL XSETUN(LUN)           Set the logical unit number, LUN, for */
/*                             output of messages from DVODPK, if */
/*                             the default is not desired. */
/*                             The default value of LUN is 6. */

/*  CALL XSETF(MFLAG)          Set a flag to control the printing of */
/*                             messages by DVODPK. */
/*                             MFLAG = 0 means do not print. (Danger: */
/*                             This risks losing valuable information.) */
/*                             MFLAG = 1 means print (the default). */

/*                             Either of the above calls may be made at */
/*                             any time and will take effect immediately. */

/*  CALL DVKSRC(RSAV,ISAV,JOB) Saves and restores the contents of */
/*                             the internal COMMON blocks used by */
/*                             DVODPK. (See Part iii below.) */
/*                             RSAV must be a real array of length 52 */
/*                             or more, and ISAV must be an integer */
/*                             array of length 52 or more. */
/*                             JOB=1 means save COMMON into RSAV/ISAV. */
/*                             JOB=2 means restore COMMON from RSAV/ISAV. */

/*                                DVKSRC is useful if one is */
/*                             interrupting a run and restarting */
/*                             later, or alternating between two or */
/*                             more problems solved with DVODPK. */

/*  CALL DVINDY(,,,,,)         Provide derivatives of y, of various */
/*        (See below.)         orders, at a specified point T, if */
/*                             desired.  It may be called only after */
/*                             a successful return from DVODPK. */

/* The detailed instructions for using DVINDY are as follows. */
/* The form of the call is: */

/*  CALL DVINDY (T, K, RWORK(21), NYH, DKY, IFLAG) */

/* The input parameters are: */

/* T         = Value of independent variable where answers are desired */
/*             (normally the same as the T last returned by DVODPK). */
/*             For valid results, T must lie between TCUR - HU and TCUR. */
/*             (See optional output for TCUR and HU.) */
/* K         = Integer order of the derivative desired.  K must satisfy */
/*             0 .le. K .le. NQCUR, where NQCUR is the current order */
/*             (see optional output).  The capability corresponding */
/*             to K = 0, i.e. computing y(T), is already provided */
/*             by DVODPK directly.  Since NQCUR .ge. 1, the first */
/*             derivative dy/dt is always available with DVINDY. */
/* RWORK(21) = The base address of the history array YH. */
/* NYH       = Column length of YH, equal to the initial value of NEQ. */

/* The output parameters are: */

/* DKY       = A real array of length NEQ containing the computed value */
/*             of the K-th derivative of y(t). */
/* IFLAG     = Integer flag, returned as 0 if K and T were legal, */
/*             -1 if K was illegal, and -2 if T was illegal. */
/*             On an error return, a message is also written. */
/* ----------------------------------------------------------------------- */
/* Part iii.  COMMON Blocks. */
/* If DVODPK is to be used in an overlay situation, the user */
/* must declare, in the primary overlay, the variables in: */
/*   (1) the call sequence to DVODPK, */
/*   (2) the three internal COMMON blocks */
/*         /DVOD01/  of length  81  (48 double precision words */
/*                         followed by 33 integer words), */
/*         /DVOD02/  of length  9  (1 double precision word */
/*                         followed by 8 integer words), */
/*         /DVPK01/  of length 14 (3 double precision words */
/*                         followed by 11 integer words) */

/* If DVODPK is used on a system in which the contents of internal */
/* COMMON blocks are not preserved between calls, the user should */
/* declare the above three COMMON blocks in the calling (main) program */
/* to insure that their contents are preserved. */

/* ----------------------------------------------------------------------- */
/* Part iv.  Optionally Replaceable Solver Routines. */

/* Below are descriptions of two routines in the DVODPK package which */
/* relate to the measurement of errors.  Either routine can be */
/* replaced by a user-supplied version, if desired.  However, since such */
/* a replacement may have a major impact on performance, it should be */
/* done only when absolutely necessary, and only with great caution. */
/* (Note: The means by which the package version of a routine is */
/* superseded by the user's version may be system-dependent.) */

/* (a) DEWSET. */
/* The following subroutine is called just before each internal */
/* integration step, and sets the array of error weights, EWT, as */
/* described under ITOL/RTOL/ATOL above: */
/*     SUBROUTINE DEWSET (NEQ, ITOL, RTOL, ATOL, YCUR, EWT) */
/* where NEQ, ITOL, RTOL, and ATOL are as in the DVODPK call sequence, */
/* YCUR contains the current dependent variable vector, and */
/* EWT is the array of weights set by DEWSET. */

/* If the user supplies this subroutine, it must return in EWT(i) */
/* (i = 1,...,NEQ) a positive quantity suitable for comparison with */
/* errors in Y(i).  The EWT array returned by DEWSET is passed to the */
/* DVNORM routine (see below), and also used by DVODPK in the computation */
/* of the optional output IMXER, the diagonal Jacobian approximation, */
/* and the increments for difference quotient Jacobians. */

/* In the user-supplied version of DEWSET, it may be desirable to use */
/* the current values of derivatives of y.  Derivatives up to order NQ */
/* are available from the history array YH, described above under */
/* Optional Output.  In DEWSET, YH is identical to the YCUR array, */
/* extended to NQ + 1 columns with a column length of NYH and scale */
/* factors of h**j/factorial(j).  On the first call for the problem, */
/* given by NST = 0, NQ is 1 and H is temporarily set to 1.0. */
/* NYH is the initial value of NEQ.  The quantities NQ, H, and NST */
/* can be obtained by including in DEWSET the statements */
/*     COMMON /DVOD01/ RVOD(48), IVOD(33) */
/*     COMMON /DVOD02/ HU, NCFN, NETF, NFE, NPE, NLU, NNI, NQU, NST */
/*     NQ = IVOD(28) */
/*     H = RVOD(21) */
/* Thus, for example, the current value of dy/dt can be obtained as */
/* YCUR(NYH+i)/H  (i=1,...,NEQ)  (and the division by H is */
/* unnecessary when NST = 0). */

/* (b) DVNORM. */
/* The following is a real function routine which computes the weighted */
/* root-mean-square norm of a vector v: */
/*     D = DVNORM (N, V, W) */
/* where: */
/*   N = the length of the vector, */
/*   V = real array of length N containing the vector, */
/*   W = real array of length N containing weights, */
/*   D = sqrt( (1/N) * sum(V(i)*W(i))**2 ). */
/* DVNORM is called with N = NEQ and with W(i) = 1.0/EWT(i), where */
/* EWT is as set by subroutine DEWSET. */

/* If the user supplies this function, it should return a non-negative */
/* value of DVNORM suitable for use in the error control in DVODPK. */
/* None of the arguments should be altered by DVNORM. */
/* For example, a user-supplied DVNORM routine might: */
/*   -substitute a max-norm of (V(i)*W(i)) for the RMS-norm, or */
/*   -ignore some components of V in the norm, with the effect of */
/*    suppressing the error control on those components of Y. */
/* ----------------------------------------------------------------------- */

/* Revision History (YYYYMMDD) */
/* 19910315  DATE WRITTEN */
/* 19910415  Minor revisions to VODPK prologue. */
/* 19920715  In demo, corrected name R1MACH to D1MACH. */
/* 19921106  In VSTEP, added ETAQ and ETAQM1 to SAVE statement. */
/* 19930701  IN VNLSK, moved line setting HRL1 below statement 220. */
/* 19940502  Minor revisions to VODPK prologue and internal comments. */
/*           In VODPK, set JACFLG = 0 if 0 < MITER < 9 and JPRE = 0. */
/*           In VNLSK, add conditions on rescaling of correction vector. */
/* 19940504  In demo programs, fixed logic in SOLSBG involving LR. */
/* 19970515  Minor revisions to VODPK prologue and internal comments. */
/*           In VHIN, attached sign to H in second derivative estimation. */
/* 19981111  In VODPK, at end of Block B, when ISTATE = 3, jump to 200. */
/* 20020423  Major upgrade: Added *DECK lines.  Renamed all routines and */
/*           Common blocks for uniqueness across single/double prec. */
/*           versions and for sharing of routines with VODE and ODEPACK. */
/*           Changed names R1MACH/D1MACH to RUMACH/DUMACH. */
/*           Converted intrinsic names to generic form. */
/*           Numerous revisions to main prologue. */
/*           Revisions to demo program - formats, intrinsics, comments. */
/* 20020426  Converted upgraded single precision version to double prec. */

/* ----------------------------------------------------------------------- */
/* Other Routines in the DVODPK Package. */

/* In addition to subroutine DVODPK, the DVODPK package includes the */
/* following subroutines and function routines: */
/*  DVHIN    computes an approximate step size for the initial step. */
/*  DVINDY   computes an interpolated value of the y vector at t = TOUT. */
/*  DVSTEP   is the core integrator, which does one step of the */
/*           integration and the associated error control. */
/*  DVSET    sets all method coefficients and test constants. */
/*  DVJUST   adjusts the history array on a change of order. */
/*  DVNLSK   solves the underlying nonlinear system -- the corrector. */
/*  DVSLPK   manages solution of linear system in chord iteration. */
/*  DVSPIG   performs the SPIGMR algorithm. */
/*  DVATV    computes a scaled, preconditioned product (I-hrl1*J)*v. */
/*  DORTHOG  orthogonalizes a vector against previous basis vectors. */
/*  DHEQR    generates a QR factorization of a Hessenberg matrix. */
/*  DHELS    finds the least squares solution of a Hessenberg system. */
/*  DVUSOL   interfaces to the user's PSOL routine (MITER = 9). */
/*  DEWSET   sets the error weight vector EWT before each step. */
/*  DVNORM   computes the weighted r.m.s. norm of a vector. */
/*  DVKSRC   is a user-callable routine to save and restore */
/*           the contents of the internal COMMON blocks. */
/*  DAXPY, DCOPY, DDOT, DNRM2, and DSCAL are basic linear */
/*           algebra modules (BLAS) used by this package. */
/*  DUMACH   computes the unit roundoff in a machine-independent manner. */
/*  XERRWD, XSETUN, XSETF, IXSAV, and IUMACH handle the printing of all */
/*           error messages and warnings.  XERRWD is machine-dependent. */
/* Note:  DVNORM, DDOT, DNRM2, DUMACH, IXSAV, and IUMACH are function */
/* routines.  All the others are subroutines. */

/* ----------------------------------------------------------------------- */

/* Declarations for external routines and function subroutines called --- */

/* Declarations for local variables ------------------------------------- */

/* ----------------------------------------------------------------------- */
/* The following Fortran-77 declarations are to cause the values of the */
/* listed (local) variables to be saved between calls to DVODPK. */
/* ----------------------------------------------------------------------- */

/* Type declarations for labeled COMMON block DVOD01 -------------------- */


/* Type declarations for labeled COMMON block DVOD02 -------------------- */


/* Type declarations for labeled COMMON block DVPK01 -------------------- */


/* ----------------------------------------------------------------------- */
/* The following internal COMMON blocks contain variables which are */
/* communicated between subroutines in the DVODPK package, or which are */
/* to be saved between calls to DVODPK. */
/* In each block, real variables precede integers. */
/* The block /DVOD01/ appears in subroutines DVODPK, DVINDY, DVSTEP, */
/* DVSET, DVJUST, DVNLSK, DVSLPK, DVATV, and DVKSRC. */
/* The block /DVOD02/ appears in subroutines DVODPK, DVINDY, DVSTEP, */
/* DVNLSK, DVSLPK, DVATV, and DVKSRC. */
/* The block /DVPK01/ appears in subroutines DVODPK, DVNLSK, DVSLPK, */
/* and DVKSRC. */

/* The variables stored in the internal COMMON blocks are as follows: */

/* ACNRM  = Weighted r.m.s. norm of accumulated correction vectors. */
/* CCMXJ  = Threshhold on DRC for updating the Jacobian. (See DRC.) */
/* CONP   = The saved value of TQ(5). */
/* CRATE  = Estimated corrector convergence rate constant. */
/* DRC    = Relative change in H*RL1 since last VJAC call. */
/* EL     = Real array of integration coefficients.  See DVSET. */
/* ETA    = Saved tentative ratio of new to old H. */
/* ETAMAX = Saved maximum value of ETA to be allowed. */
/* H      = The step size. */
/* HMIN   = The minimum absolute value of the step size H to be used. */
/* HMXI   = Inverse of the maximum absolute value of H to be used. */
/*          HMXI = 0.0 is allowed and corresponds to an infinite HMAX. */
/* HNEW   = The step size to be attempted on the next step. */
/* HSCAL  = Stepsize in scaling of YH array. */
/* PRL1   = The saved value of RL1. */
/* RC     = Ratio of current H*RL1 to value on last VJAC call. */
/* RL1    = The reciprocal of the coefficient EL(1). */
/* TAU    = Real vector of past NQ step sizes, length 13. */
/* TQ     = A real vector of length 5 in which DVSET stores constants */
/*          used for the convergence test, the error test, and the */
/*          selection of H at a new order. */
/* TN     = The independent variable, updated on each step taken. */
/* UROUND = The machine unit roundoff.  The smallest positive real number */
/*          such that  1.0 + UROUND .ne. 1.0 */
/* ICF    = Integer flag for convergence failure in DVNLSK: */
/*            0 means no failures. */
/*            1 means convergence failure with out of date Jacobian */
/*                   (recoverable error). */
/*            2 means convergence failure with current Jacobian or */
/*                   singular matrix (unrecoverable error). */
/* INIT   = Saved integer flag indicating whether initialization of the */
/*          problem has been done (INIT = 1) or not. */
/* IPUP   = Saved flag to signal updating of Newton matrix. */
/* JCUR   = Output flag from VJAC showing Jacobian status: */
/*            JCUR = 0 means J is not current. */
/*            JCUR = 1 means J is current. */
/* JSTART = Integer flag used as input to DVSTEP: */
/*            0  means perform the first step. */
/*            1  means take a new step continuing from the last. */
/*            -1 means take the next step with a new value of MAXORD, */
/*                  HMIN, HMXI, N, METH, MITER, and/or matrix parameters. */
/*          On return, DVSTEP sets JSTART = 1. */
/* JSV    = Integer flag for Jacobian saving, = sign(MF). */
/* KFLAG  = A completion code from DVSTEP with the following meanings: */
/*               0      the step was succesful. */
/*              -1      the requested error could not be achieved. */
/*              -2      corrector convergence could not be achieved. */
/*              -3, -4  fatal error in VNLS. */
/* KUTH   = Input flag to DVSTEP showing whether H was reduced by the */
/*          driver.  KUTH = 1 if H was reduced, = 0 otherwise. */
/* L      = Integer variable, NQ + 1, current order plus one. */
/* LMAX   = MAXORD + 1 (used for dimensioning). */
/* LOCJS  = A pointer to the saved Jacobian, whose storage starts at */
/*          WM(LOCJS), if JSV = 1. */
/* LYH, LEWT, LACOR, LSAVF, LWM, LIWM = Saved integer pointers */
/*          to segments of RWORK and IWORK. */
/* MAXORD = The maximum order of integration method to be allowed. */
/* METH/MITER = The method flags.  See MF. */
/* MSBJ   = The maximum number of steps between J evaluations, = 50. */
/* MXHNIL = Saved value of optional input MXHNIL. */
/* MXSTEP = Saved value of optional input MXSTEP. */
/* N      = The number of first-order ODEs, = NEQ. */
/* NEWH   = Saved integer to flag change of H. */
/* NEWQ   = The method order to be used on the next step. */
/* NHNIL  = Saved counter for occurrences of T + H = T. */
/* NQ     = Integer variable, the current integration method order. */
/* NQNYH  = Saved value of NQ*NYH. */
/* NQWAIT = A counter controlling the frequency of order changes. */
/*          An order change is about to be considered if NQWAIT = 1. */
/* NSLJ   = The number of steps taken as of the last Jacobian update. */
/* NSLP   = Saved value of NST as of last Newton matrix update. */
/* NYH    = Saved value of the initial value of NEQ. */

/* HU     = The step size in t last used. */
/* NCFN   = Number of nonlinear convergence failures so far. */
/* NETF   = The number of error test failures of the integrator so far. */
/* NFE    = The number of f evaluations for the problem so far. */
/* NPE    = The number of preconditioner evaluations (JAC calls) so far. */
/* NLU    = The number of matrix LU decompositions so far. */
/* NNI    = Number of nonlinear iterations so far. */
/* NQU    = The method order last used. */
/* NST    = The number of steps taken for the problem so far. */

/* DELT   = Convergence test constant in Krylov iterations. */
/* SQRTN  = SQRT(NEQ), for use in weights in Krylov convergence tests. */
/* RSQRTN = 1.0/SQRTN, also for use in convergence weights. */
/* JPRE   = Preconditioner type flag. */
/* JACFLG = Indicator for presence of user-supplied JAC routine. */
/* LOCWP  = Location of start of user's WP array in WM work array. */
/* LOCIWP = Location of start of user's IWP array in IWM work array. */
/* LVSAV  = Saved pointer to VSAV array in RWORK. */
/* KMP    = Number of vectors on which orthogonalization is done in */
/*          Krylov iteration. */
/* MAXL   = Maximum dimension of Krylov subspace used. */
/* MNEWT  = Newton iteration index. */
/* NLI    = Number of linear (Krylov) iterations done. */
/* NPS    = Number of preconditioner solvers (PSOL calls) done. */
/* NCFL   = Number of convergence failures in Krylov iteration. */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --y;
    --rtol;
    --atol;
    --rwork;
    --iwork;
    --rpar;
    --ipar;

    /* Function Body */
/* ----------------------------------------------------------------------- */
/* Block A. */
/* This code block is executed on every call. */
/* It tests ISTATE and ITASK for legality and branches appropriately. */
/* If ISTATE .gt. 1 but the flag INIT shows that initialization has */
/* not yet been done, an error return occurs. */
/* If ISTATE = 1 and TOUT = T, jump to Block G and return immediately. */
/* ----------------------------------------------------------------------- */
    if (*istate < 1 || *istate > 3) {
	goto L601;
    }
    if (*itask < 1 || *itask > 5) {
	goto L602;
    }
    if (*istate == 1) {
	goto L10;
    }
    if (dvod01_1.init == 0) {
	goto L603;
    }
    if (*istate == 2) {
	goto L200;
    }
    goto L20;
L10:
    dvod01_1.init = 0;
    if (*tout == *t) {
	return 0;
    }
/* ----------------------------------------------------------------------- */
/* Block B. */
/* The next code block is executed for the initial call (ISTATE = 1), */
/* or for a continuation call with parameter changes (ISTATE = 3). */
/* It contains checking of all inputs and various initializations. */

/* First check legality of the non-optional inputs NEQ, ITOL, IOPT, MF. */
/* ----------------------------------------------------------------------- */
L20:
    if (*neq <= 0) {
	goto L604;
    }
    if (*istate == 1) {
	goto L25;
    }
    if (*neq > dvod01_1.n) {
	goto L605;
    }
L25:
    dvod01_1.n = *neq;
    if (*itol < 1 || *itol > 4) {
	goto L606;
    }
    if (*iopt < 0 || *iopt > 1) {
	goto L607;
    }
    dvod01_1.jsv = 0;
    dvod01_1.meth = *mf / 10;
    dvod01_1.miter = *mf - dvod01_1.meth * 10;
    if (dvod01_1.meth < 1 || dvod01_1.meth > 2) {
	goto L608;
    }
    if (dvod01_1.miter < 0) {
	goto L608;
    }
    if (dvod01_1.miter > 1 && dvod01_1.miter < 9) {
	goto L608;
    }
    if (dvod01_1.miter >= 1) {
	dvpk01_1.jpre = iwork[3];
    }
    dvpk01_1.jacflg = 0;
    if (dvod01_1.miter >= 1) {
	dvpk01_1.jacflg = iwork[4];
    }
    if (dvod01_1.miter >= 1 && dvod01_1.miter != 9 && dvpk01_1.jpre == 0) {
	dvpk01_1.jacflg = 0;
    }
/* Next process and check the optional inputs. -------------------------- */
    if (*iopt == 1) {
	goto L40;
    }
    dvod01_1.maxord = mord[dvod01_1.meth - 1];
    dvod01_1.mxstep = mxstp0;
    dvod01_1.mxhnil = mxhnl0;
    if (*istate == 1) {
	h0 = zero;
    }
    dvod01_1.hmxi = zero;
    dvod01_1.hmin = zero;
    dvpk01_1.maxl = min(5,dvod01_1.n);
    dvpk01_1.kmp = dvpk01_1.maxl;
    dvpk01_1.delt = pt05;
    goto L60;
L40:
    dvod01_1.maxord = iwork[5];
    if (dvod01_1.maxord < 0) {
	goto L611;
    }
    if (dvod01_1.maxord == 0) {
	dvod01_1.maxord = 100;
    }
/* Computing MIN */
    i__1 = dvod01_1.maxord, i__2 = mord[dvod01_1.meth - 1];
    dvod01_1.maxord = min(i__1,i__2);
    dvod01_1.mxstep = iwork[6];
    if (dvod01_1.mxstep < 0) {
	goto L612;
    }
    if (dvod01_1.mxstep == 0) {
	dvod01_1.mxstep = mxstp0;
    }
    dvod01_1.mxhnil = iwork[7];
    if (dvod01_1.mxhnil < 0) {
	goto L613;
    }
    if (dvod01_1.mxhnil == 0) {
	dvod01_1.mxhnil = mxhnl0;
    }
    if (*istate != 1) {
	goto L50;
    }
    h0 = rwork[5];
    if ((*tout - *t) * h0 < zero) {
	goto L614;
    }
L50:
    hmax = rwork[6];
    if (hmax < zero) {
	goto L615;
    }
    dvod01_1.hmxi = zero;
    if (hmax > zero) {
	dvod01_1.hmxi = one / hmax;
    }
    dvod01_1.hmin = rwork[7];
    if (dvod01_1.hmin < zero) {
	goto L616;
    }
    dvpk01_1.maxl = iwork[8];
    if (dvpk01_1.maxl == 0) {
	dvpk01_1.maxl = 5;
    }
    dvpk01_1.maxl = min(dvpk01_1.maxl,dvod01_1.n);
    dvpk01_1.kmp = iwork[9];
    if (dvpk01_1.kmp == 0 || dvpk01_1.kmp > dvpk01_1.maxl) {
	dvpk01_1.kmp = dvpk01_1.maxl;
    }
    dvpk01_1.delt = rwork[8];
    if (dvpk01_1.delt == 0.) {
	dvpk01_1.delt = pt05;
    }
/* ----------------------------------------------------------------------- */
/* Set work array pointers and check lengths lrw and liw. */
/* Pointers to segments of RWORK and iwork are named by prefixing l to */
/* the name of the segment.  e.g., the segment YH starts at RWORK(LYH). */
/* Segments of RWORK (in order) are  YH, WM, EWT, SAVF, VSAV, ACOR. */
/* Within WM, LOCWP is the location of the WP work array, */
/* and within IWM, LOCIWP is the location of the IWP work array. */
/* ----------------------------------------------------------------------- */
L60:
    dvod01_1.lyh = 21;
    if (*istate == 1) {
	dvod01_1.nyh = dvod01_1.n;
    }
    dvod01_1.lwm = dvod01_1.lyh + (dvod01_1.maxord + 1) * dvod01_1.nyh;
    if (dvod01_1.miter == 0) {
	lenwk = 0;
    }
    if (dvod01_1.miter == 1) {
/* Computing MIN */
	i__1 = 1, i__2 = dvpk01_1.maxl - dvpk01_1.kmp;
	lenwk = dvod01_1.n * (dvpk01_1.maxl + 2 + min(i__1,i__2)) + (
		dvpk01_1.maxl + 3) * dvpk01_1.maxl + 1;
    }
    if (dvod01_1.miter == 9) {
	lenwk = dvod01_1.n << 1;
    }
    lwp = 0;
    if (dvod01_1.miter >= 1) {
	lwp = iwork[1];
    }
    lenwm = lenwk + lwp;
    dvpk01_1.locwp = lenwk + 1;
    dvod01_1.lewt = dvod01_1.lwm + lenwm;
    dvod01_1.lsavf = dvod01_1.lewt + dvod01_1.n;
    dvpk01_1.lvsav = dvod01_1.lsavf + dvod01_1.n;
    dvod01_1.lacor = dvpk01_1.lvsav + dvod01_1.n;
    if (dvod01_1.miter == 0) {
	dvod01_1.lacor = dvpk01_1.lvsav;
    }
    lenrw = dvod01_1.lacor + dvod01_1.n - 1;
    iwork[17] = lenrw;
    dvod01_1.liwm = 31;
    leniwk = 0;
    liwp = 0;
    if (dvod01_1.miter >= 1) {
	liwp = iwork[2];
    }
    leniw = leniwk + 30 + liwp;
    dvpk01_1.lociwp = leniwk + 1;
    iwork[18] = leniw;
    if (lenrw > *lrw) {
	goto L617;
    }
    if (leniw > *liw) {
	goto L618;
    }
/* Check RTOL and ATOL for legality. ------------------------------------ */
    rtoli = rtol[1];
    atoli = atol[1];
    i__1 = dvod01_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (*itol >= 3) {
	    rtoli = rtol[i__];
	}
	if (*itol == 2 || *itol == 4) {
	    atoli = atol[i__];
	}
	if (rtoli < zero) {
	    goto L619;
	}
	if (atoli < zero) {
	    goto L620;
	}
/* L70: */
    }
/* Load SQRT(N) and its reciprocal in common. --------------------------- */
    dvpk01_1.sqrtn = sqrt((doublereal) dvod01_1.n);
    dvpk01_1.rsqrtn = one / dvpk01_1.sqrtn;
    if (*istate == 1) {
	goto L100;
    }
/* If ISTATE = 3, set flag to signal parameter changes to DVSTEP. ------- */
    dvod01_1.jstart = -1;
    if (dvod01_1.nq <= dvod01_1.maxord) {
	goto L200;
    }
/* MAXORD was reduced below NQ.  Copy YH(*,MAXORD+2) into SAVF. --------- */
    dcopy_(&dvod01_1.n, &rwork[dvod01_1.lwm], &c__1, &rwork[dvod01_1.lsavf], &
	    c__1);
    goto L200;
/* ----------------------------------------------------------------------- */
/* Block C. */
/* The next block is for the initial call only (ISTATE = 1). */
/* It contains all remaining initializations, the initial call to F, */
/* and the calculation of the initial step size. */
/* The error weights in EWT are inverted after being loaded. */
/* ----------------------------------------------------------------------- */
L100:
    dvod01_1.uround = dumach_();
    dvod01_1.tn = *t;
    if (*itask != 4 && *itask != 5) {
	goto L110;
    }
    tcrit = rwork[1];
    if ((tcrit - *tout) * (*tout - *t) < zero) {
	goto L625;
    }
    if (h0 != zero && (*t + h0 - tcrit) * h0 > zero) {
	h0 = tcrit - *t;
    }
L110:
    dvod01_1.jstart = 0;
    dvod01_1.ccmxj = pt2;
    dvod01_1.msbj = 50;
    dvod01_1.nhnil = 0;
    dvod02_1.nst = 0;
    nslast = 0;
    dvod02_1.hu = zero;
    dvod02_1.nqu = 0;
    dvod02_3.npe = 0;
    nli0 = 0;
    nni0 = 0;
    ncfn0 = 0;
    ncfl0 = 0;
    nwarn = 0;
    dvod02_1.nni = 0;
    dvpk01_1.nli = 0;
    dvpk01_1.nps = 0;
    dvod02_1.netf = 0;
    dvod02_1.ncfn = 0;
    dvpk01_1.ncfl = 0;
/* Initial call to F.  (LF0 points to YH(*,2).) ------------------------- */
    lf0 = dvod01_1.lyh + dvod01_1.nyh;
    (*f)(&dvod01_1.n, t, &y[1], &rwork[lf0], &rpar[1], &ipar[1]);
    dvod02_1.nfe = 1;
/* Load the initial value vector in YH. --------------------------------- */
    dcopy_(&dvod01_1.n, &y[1], &c__1, &rwork[dvod01_1.lyh], &c__1);
/* Load and invert the EWT array.  (H is temporarily set to 1.0.) ------- */
    dvod01_1.nq = 1;
    dvod01_1.h__ = one;
    dewset_(&dvod01_1.n, itol, &rtol[1], &atol[1], &rwork[dvod01_1.lyh], &
	    rwork[dvod01_1.lewt]);
    i__1 = dvod01_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (rwork[i__ + dvod01_1.lewt - 1] <= zero) {
	    goto L621;
	}
/* L120: */
	rwork[i__ + dvod01_1.lewt - 1] = one / rwork[i__ + dvod01_1.lewt - 1];
    }
    if (h0 != zero) {
	goto L180;
    }
/* Call DVHIN to set initial step size H0 to be attempted. -------------- */
    dvhin_(&dvod01_1.n, t, &rwork[dvod01_1.lyh], &rwork[lf0], (S_fp)f, &rpar[
	    1], &ipar[1], tout, &dvod01_1.uround, &rwork[dvod01_1.lewt], itol,
	     &atol[1], &y[1], &rwork[dvod01_1.lacor], &h0, &niter, &ier);
    dvod02_1.nfe += niter;
    if (ier != 0) {
	goto L622;
    }
/* Adjust H0 if necessary to meet HMAX bound. --------------------------- */
L180:
    rh = abs(h0) * dvod01_1.hmxi;
    if (rh > one) {
	h0 /= rh;
    }
/* Load H with H0 and scale YH(*,2) by H0. ------------------------------ */
    dvod01_1.h__ = h0;
    dscal_(&dvod01_1.n, &h0, &rwork[lf0], &c__1);
    goto L270;
/* ----------------------------------------------------------------------- */
/* Block D. */
/* The next code block is for continuation calls only (ISTATE = 2 or 3) */
/* and is to check stop conditions before taking a step. */
/* ----------------------------------------------------------------------- */
L200:
    nslast = dvod02_1.nst;
    dvod01_1.kuth = 0;
    nli0 = dvpk01_1.nli;
    nni0 = dvod02_1.nni;
    ncfn0 = dvod02_1.ncfn;
    ncfl0 = dvpk01_1.ncfl;
    nwarn = 0;
    switch (*itask) {
	case 1:  goto L210;
	case 2:  goto L250;
	case 3:  goto L220;
	case 4:  goto L230;
	case 5:  goto L240;
    }
L210:
    if ((dvod01_1.tn - *tout) * dvod01_1.h__ < zero) {
	goto L250;
    }
    dvindy_(tout, &c__0, &rwork[dvod01_1.lyh], &dvod01_1.nyh, &y[1], &iflag);
    if (iflag != 0) {
	goto L627;
    }
    *t = *tout;
    goto L420;
L220:
    tp = dvod01_1.tn - dvod02_1.hu * (one + hun * dvod01_1.uround);
    if ((tp - *tout) * dvod01_1.h__ > zero) {
	goto L623;
    }
    if ((dvod01_1.tn - *tout) * dvod01_1.h__ < zero) {
	goto L250;
    }
    goto L400;
L230:
    tcrit = rwork[1];
    if ((dvod01_1.tn - tcrit) * dvod01_1.h__ > zero) {
	goto L624;
    }
    if ((tcrit - *tout) * dvod01_1.h__ < zero) {
	goto L625;
    }
    if ((dvod01_1.tn - *tout) * dvod01_1.h__ < zero) {
	goto L245;
    }
    dvindy_(tout, &c__0, &rwork[dvod01_1.lyh], &dvod01_1.nyh, &y[1], &iflag);
    if (iflag != 0) {
	goto L627;
    }
    *t = *tout;
    goto L420;
L240:
    tcrit = rwork[1];
    if ((dvod01_1.tn - tcrit) * dvod01_1.h__ > zero) {
	goto L624;
    }
L245:
    hmx = abs(dvod01_1.tn) + abs(dvod01_1.h__);
    ihit = (d__1 = dvod01_1.tn - tcrit, abs(d__1)) <= hun * dvod01_1.uround * 
	    hmx;
    if (ihit) {
	goto L400;
    }
    tnext = dvod01_1.tn + dvod01_1.hnew * (one + four * dvod01_1.uround);
    if ((tnext - tcrit) * dvod01_1.h__ <= zero) {
	goto L250;
    }
    dvod01_1.h__ = (tcrit - dvod01_1.tn) * (one - four * dvod01_1.uround);
    dvod01_1.kuth = 1;
/* ----------------------------------------------------------------------- */
/* Block E. */
/* The next block is normally executed for all calls and contains */
/* the call to the one-step core integrator DVSTEP. */

/* This is a looping point for the integration steps. */

/* First check for too many steps being taken, */
/* check for poor Newton/Krylov performance, update EWT (if not at */
/* start of problem), check for too much accuracy being requested, and */
/* check for H below the roundoff level in T. */
/* ----------------------------------------------------------------------- */
L250:
    if (dvod02_1.nst - nslast >= dvod01_1.mxstep) {
	goto L500;
    }
    dewset_(&dvod01_1.n, itol, &rtol[1], &atol[1], &rwork[dvod01_1.lyh], &
	    rwork[dvod01_1.lewt]);
    nstd = dvod02_1.nst - nslast;
    nnid = dvod02_1.nni - nni0;
    if (nstd < 10 || nnid == 0) {
	goto L255;
    }
    avdim = (real) (dvpk01_1.nli - nli0) / (real) nnid;
    rcfn = (real) (dvod02_1.ncfn - ncfn0) / (real) nstd;
    rcfl = (real) (dvpk01_1.ncfl - ncfl0) / (real) nnid;
    lavd = avdim > dvpk01_1.maxl - pt05;
    lcfn = rcfn > pt9;
    lcfl = rcfl > pt9;
    lwarn = lavd || lcfn || lcfl;
    if (! lwarn) {
	goto L255;
    }
    ++nwarn;
    if (nwarn > 10) {
	goto L255;
    }
    if (lavd) {
	s_copy(msg, "DVODPK- Warning. Poor iterative algorithm performance   "
		, (ftnlen)80, (ftnlen)56);
	xerrwd_(msg, &c__56, &c__111, &c__0, &c__0, &c__0, &c__0, &c__0, &
		zero, &zero, (ftnlen)80);
	s_copy(msg, "      at T = R1. Average no. of linear iterations = R2  "
		, (ftnlen)80, (ftnlen)56);
	xerrwd_(msg, &c__56, &c__111, &c__0, &c__0, &c__0, &c__0, &c__2, &
		dvod01_1.tn, &avdim, (ftnlen)80);
    }
    if (lcfn) {
	s_copy(msg, "DVODPK- Warning. Poor iterative algorithm performance   "
		, (ftnlen)80, (ftnlen)56);
	xerrwd_(msg, &c__56, &c__112, &c__0, &c__0, &c__0, &c__0, &c__0, &
		zero, &zero, (ftnlen)80);
	s_copy(msg, "      at T = R1. Nonlinear convergence failure rate = R2"
		, (ftnlen)80, (ftnlen)56);
	xerrwd_(msg, &c__56, &c__112, &c__0, &c__0, &c__0, &c__0, &c__2, &
		dvod01_1.tn, &rcfn, (ftnlen)80);
    }
    if (lcfl) {
	s_copy(msg, "DVODPK- Warning. Poor iterative algorithm performance   "
		, (ftnlen)80, (ftnlen)56);
	xerrwd_(msg, &c__56, &c__113, &c__0, &c__0, &c__0, &c__0, &c__0, &
		zero, &zero, (ftnlen)80);
	s_copy(msg, "      at T = R1. Linear convergence failure rate = R2   "
		, (ftnlen)80, (ftnlen)56);
	xerrwd_(msg, &c__56, &c__113, &c__0, &c__0, &c__0, &c__0, &c__2, &
		dvod01_1.tn, &rcfl, (ftnlen)80);
    }
L255:
    i__1 = dvod01_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (rwork[i__ + dvod01_1.lewt - 1] <= zero) {
	    goto L510;
	}
/* L260: */
	rwork[i__ + dvod01_1.lewt - 1] = one / rwork[i__ + dvod01_1.lewt - 1];
    }
L270:
    tolsf = dvod01_1.uround * dvnorm_(&dvod01_1.n, &rwork[dvod01_1.lyh], &
	    rwork[dvod01_1.lewt]);
    if (tolsf <= one) {
	goto L280;
    }
    tolsf *= two;
    if (dvod02_1.nst == 0) {
	goto L626;
    }
    goto L520;
L280:
    if (dvod01_1.tn + dvod01_1.h__ != dvod01_1.tn) {
	goto L290;
    }
    ++dvod01_1.nhnil;
    if (dvod01_1.nhnil > dvod01_1.mxhnil) {
	goto L290;
    }
    s_copy(msg, "DVODPK-  Warning: internal T (=R1) and H (=R2) are", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__101, &c__1, &c__0, &c__0, &c__0, &c__0, &zero, &
	    zero, (ftnlen)80);
    s_copy(msg, "      such that in the machine, T + H = T on the next step  "
	    , (ftnlen)80, (ftnlen)60);
    xerrwd_(msg, &c__60, &c__101, &c__1, &c__0, &c__0, &c__0, &c__0, &zero, &
	    zero, (ftnlen)80);
    s_copy(msg, "      (H = step size). solver will continue anyway", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__101, &c__1, &c__0, &c__0, &c__0, &c__2, &
	    dvod01_1.tn, &dvod01_1.h__, (ftnlen)80);
    if (dvod01_1.nhnil < dvod01_1.mxhnil) {
	goto L290;
    }
    s_copy(msg, "DVODPK-  Above warning has been issued I1 times.  ", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__102, &c__1, &c__0, &c__0, &c__0, &c__0, &zero, &
	    zero, (ftnlen)80);
    s_copy(msg, "      it will not be issued again for this problem", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__102, &c__1, &c__1, &dvod01_1.mxhnil, &c__0, &
	    c__0, &zero, &zero, (ftnlen)80);
L290:
/* ----------------------------------------------------------------------- */
/*  CALL DVSTEP (Y, YH, NYH, YH, EWT, SAVF, VSAV, ACOR, WM, IWM, */
/*                                     F, JAC, PSOL, DVNLSK, RPAR, IPAR) */
/* ----------------------------------------------------------------------- */
    dvstep_(&y[1], &rwork[dvod01_1.lyh], &dvod01_1.nyh, &rwork[dvod01_1.lyh], 
	    &rwork[dvod01_1.lewt], &rwork[dvod01_1.lsavf], &rwork[
	    dvpk01_1.lvsav], &rwork[dvod01_1.lacor], &rwork[dvod01_1.lwm], &
	    iwork[dvod01_1.liwm], (S_fp)f, (U_fp)jac, (U_fp)psol, (U_fp)
	    dvnlsk_, &rpar[1], &ipar[1]);
    kgo = 1 - dvod01_1.kflag;
    switch (kgo) {
	case 1:  goto L300;
	case 2:  goto L530;
	case 3:  goto L540;
	case 4:  goto L550;
	case 5:  goto L555;
    }
/* ----------------------------------------------------------------------- */
/* Block F. */
/* The following block handles the case of a successful return from the */
/* core integrator (KFLAG = 0).  Test for stop conditions. */
/* ----------------------------------------------------------------------- */
L300:
    dvod01_1.init = 1;
    dvod01_1.kuth = 0;
    switch (*itask) {
	case 1:  goto L310;
	case 2:  goto L400;
	case 3:  goto L330;
	case 4:  goto L340;
	case 5:  goto L350;
    }
/* ITASK = 1.  if TOUT has been reached, interpolate. ------------------- */
L310:
    if ((dvod01_1.tn - *tout) * dvod01_1.h__ < zero) {
	goto L250;
    }
    dvindy_(tout, &c__0, &rwork[dvod01_1.lyh], &dvod01_1.nyh, &y[1], &iflag);
    *t = *tout;
    goto L420;
/* ITASK = 3.  Jump to exit if TOUT was reached. ------------------------ */
L330:
    if ((dvod01_1.tn - *tout) * dvod01_1.h__ >= zero) {
	goto L400;
    }
    goto L250;
/* ITASK = 4.  See if TOUT or TCRIT was reached.  Adjust H if necessary. */
L340:
    if ((dvod01_1.tn - *tout) * dvod01_1.h__ < zero) {
	goto L345;
    }
    dvindy_(tout, &c__0, &rwork[dvod01_1.lyh], &dvod01_1.nyh, &y[1], &iflag);
    *t = *tout;
    goto L420;
L345:
    hmx = abs(dvod01_1.tn) + abs(dvod01_1.h__);
    ihit = (d__1 = dvod01_1.tn - tcrit, abs(d__1)) <= hun * dvod01_1.uround * 
	    hmx;
    if (ihit) {
	goto L400;
    }
    tnext = dvod01_1.tn + dvod01_1.h__ * (one + four * dvod01_1.uround);
    if ((tnext - tcrit) * dvod01_1.h__ <= zero) {
	goto L250;
    }
    dvod01_1.h__ = (tcrit - dvod01_1.tn) * (one - four * dvod01_1.uround);
    dvod01_1.kuth = 1;
    goto L250;
/* ITASK = 5.  See if TCRIT was reached and jump to exit. --------------- */
L350:
    hmx = abs(dvod01_1.tn) + abs(dvod01_1.h__);
    ihit = (d__1 = dvod01_1.tn - tcrit, abs(d__1)) <= hun * dvod01_1.uround * 
	    hmx;
/* ----------------------------------------------------------------------- */
/* Block G. */
/* The following block handles all successful returns from DVODPK. */
/* If ITASK .ne. 1, Y is loaded from YH and T is set accordingly. */
/* ISTATE is set to 2, and the optional outputs are loaded into the */
/* work arrays before returning. */
/* ----------------------------------------------------------------------- */
L400:
    dcopy_(&dvod01_1.n, &rwork[dvod01_1.lyh], &c__1, &y[1], &c__1);
    *t = dvod01_1.tn;
    if (*itask != 4 && *itask != 5) {
	goto L420;
    }
    if (ihit) {
	*t = tcrit;
    }
L420:
    *istate = 2;
    rwork[11] = dvod02_1.hu;
    rwork[12] = dvod01_1.h__;
    rwork[13] = dvod01_1.tn;
    iwork[11] = dvod02_1.nst;
    iwork[12] = dvod02_1.nfe;
    iwork[13] = dvod02_3.npe;
    iwork[14] = dvod02_1.nqu;
    iwork[15] = dvod01_1.nq;
    iwork[20] = dvod02_1.nni;
    iwork[21] = dvod02_1.ncfn;
    iwork[22] = dvod02_1.netf;
    iwork[23] = dvpk01_1.nli;
    iwork[24] = dvpk01_1.nps;
    iwork[25] = dvpk01_1.ncfl;
    return 0;
/* ----------------------------------------------------------------------- */
/* Block H. */
/* The following block handles all unsuccessful returns other than */
/* those for illegal input.  First the error message routine is called. */
/* if there was an error test or convergence test failure, IMXER is set. */
/* Then Y is loaded from YH, and T is set to TN.  The optional outputs */
/* are loaded into the work arrays before returning. */
/* ----------------------------------------------------------------------- */
/* The maximum number of steps was taken before reaching TOUT. ---------- */
L500:
    s_copy(msg, "DVODPK-  At current T (=R1), MXSTEP (=I1) steps   ", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__201, &c__1, &c__0, &c__0, &c__0, &c__0, &zero, &
	    zero, (ftnlen)80);
    s_copy(msg, "      taken on this call before reaching TOUT     ", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__201, &c__1, &c__1, &dvod01_1.mxstep, &c__0, &
	    c__1, &dvod01_1.tn, &zero, (ftnlen)80);
    *istate = -1;
    goto L580;
/* EWT(i) .le. 0.0 for some i (not at start of problem). ---------------- */
L510:
    ewti = rwork[dvod01_1.lewt + i__ - 1];
    s_copy(msg, "DVODPK-  At T (=R1), EWT(I1) has become R2 .le. 0.", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__202, &c__1, &c__1, &i__, &c__0, &c__2, &
	    dvod01_1.tn, &ewti, (ftnlen)80);
    *istate = -6;
    goto L580;
/* Too much accuracy requested for machine precision. ------------------- */
L520:
    s_copy(msg, "DVODPK-  At T (=R1), too much accuracy requested  ", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__203, &c__1, &c__0, &c__0, &c__0, &c__0, &zero, &
	    zero, (ftnlen)80);
    s_copy(msg, "      for precision of machine:  See TOLSF (=R2)  ", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__203, &c__1, &c__0, &c__0, &c__0, &c__2, &
	    dvod01_1.tn, &tolsf, (ftnlen)80);
    rwork[14] = tolsf;
    *istate = -2;
    goto L580;
/* KFLAG = -1.  Error test failed repeatedly or with ABS(H) = HMIN. ----- */
L530:
    s_copy(msg, "DVODPK-  At T(=R1) and step size H(=R2), the error", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__204, &c__1, &c__0, &c__0, &c__0, &c__0, &zero, &
	    zero, (ftnlen)80);
    s_copy(msg, "      test failed repeatedly or with abs(H) = HMIN", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__204, &c__1, &c__0, &c__0, &c__0, &c__2, &
	    dvod01_1.tn, &dvod01_1.h__, (ftnlen)80);
    *istate = -4;
    goto L560;
/* KFLAG = -2.  Convergence failed repeatedly or with ABS(H) = HMIN. ---- */
L540:
    s_copy(msg, "DVODPK-  At T (=R1) and step size H (=R2), the    ", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__205, &c__1, &c__0, &c__0, &c__0, &c__0, &zero, &
	    zero, (ftnlen)80);
    s_copy(msg, "      corrector convergence failed repeatedly     ", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__205, &c__1, &c__0, &c__0, &c__0, &c__0, &zero, &
	    zero, (ftnlen)80);
    s_copy(msg, "      or with abs(H) = HMIN   ", (ftnlen)80, (ftnlen)30);
    xerrwd_(msg, &c__30, &c__205, &c__1, &c__0, &c__0, &c__0, &c__2, &
	    dvod01_1.tn, &dvod01_1.h__, (ftnlen)80);
    *istate = -5;
    goto L560;
/* KFLAG = -3.  Unrecoverable error from JAC. --------------------------- */
L550:
    s_copy(msg, "DVODPK-  at T (=R1) an unrecoverable error return ", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__206, &c__0, &c__0, &c__0, &c__0, &c__0, &zero, &
	    zero, (ftnlen)80);
    s_copy(msg, "      was made from subroutine JAC      ", (ftnlen)80, (
	    ftnlen)40);
    xerrwd_(msg, &c__40, &c__206, &c__0, &c__0, &c__0, &c__0, &c__1, &
	    dvod01_1.tn, &zero, (ftnlen)80);
    *istate = -7;
    goto L580;
/* KFLAG = -4.  Unrecoverable error from PSOL. -------------------------- */
L555:
    s_copy(msg, "DVODPK-  at T (=R1) an unrecoverable error return ", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__207, &c__0, &c__0, &c__0, &c__0, &c__0, &zero, &
	    zero, (ftnlen)80);
    s_copy(msg, "      was made from subroutine PSOL     ", (ftnlen)80, (
	    ftnlen)40);
    xerrwd_(msg, &c__40, &c__207, &c__0, &c__0, &c__0, &c__0, &c__1, &
	    dvod01_1.tn, &zero, (ftnlen)80);
    *istate = -7;
    goto L580;
/* Compute IMXER if relevant. ------------------------------------------- */
L560:
    big = zero;
    imxer = 1;
    i__1 = dvod01_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	size = (d__1 = rwork[i__ + dvod01_1.lacor - 1] * rwork[i__ + 
		dvod01_1.lewt - 1], abs(d__1));
	if (big >= size) {
	    goto L570;
	}
	big = size;
	imxer = i__;
L570:
	;
    }
    iwork[16] = imxer;
/* Set Y vector, T, and optional outputs. ------------------------------- */
L580:
    dcopy_(&dvod01_1.n, &rwork[dvod01_1.lyh], &c__1, &y[1], &c__1);
    *t = dvod01_1.tn;
    rwork[11] = dvod02_1.hu;
    rwork[12] = dvod01_1.h__;
    rwork[13] = dvod01_1.tn;
    iwork[11] = dvod02_1.nst;
    iwork[12] = dvod02_1.nfe;
    iwork[13] = dvod02_3.npe;
    iwork[14] = dvod02_1.nqu;
    iwork[15] = dvod01_1.nq;
    iwork[20] = dvod02_1.nni;
    iwork[21] = dvod02_1.ncfn;
    iwork[22] = dvod02_1.netf;
    iwork[23] = dvpk01_1.nli;
    iwork[24] = dvpk01_1.nps;
    iwork[25] = dvpk01_1.ncfl;
    return 0;
/* ----------------------------------------------------------------------- */
/* Block I. */
/* The following block handles all error returns due to illegal input */
/* (ISTATE = -3), as detected before calling the core integrator. */
/* Call the error message routine and then return. */
/* ----------------------------------------------------------------------- */
L601:
    s_copy(msg, "DVODPK-  ISTATE (=I1) illegal ", (ftnlen)80, (ftnlen)30);
    xerrwd_(msg, &c__30, &c__1, &c__1, &c__1, istate, &c__0, &c__0, &zero, &
	    zero, (ftnlen)80);
    if (*istate < 0) {
	goto L800;
    }
    goto L700;
L602:
    s_copy(msg, "DVODPK-  ITASK (=I1) illegal  ", (ftnlen)80, (ftnlen)30);
    xerrwd_(msg, &c__30, &c__2, &c__1, &c__1, itask, &c__0, &c__0, &zero, &
	    zero, (ftnlen)80);
    goto L700;
L603:
    s_copy(msg, "DVODPK-   ISTATE (=I1) .gt. 1 but DVODPK not initialized    "
	    , (ftnlen)80, (ftnlen)60);
    xerrwd_(msg, &c__60, &c__3, &c__1, &c__1, istate, &c__0, &c__0, &zero, &
	    zero, (ftnlen)80);
    goto L700;
L604:
    s_copy(msg, "DVODPK-  NEQ (=I1) .lt. 1     ", (ftnlen)80, (ftnlen)30);
    xerrwd_(msg, &c__30, &c__4, &c__1, &c__1, neq, &c__0, &c__0, &zero, &zero,
	     (ftnlen)80);
    goto L700;
L605:
    s_copy(msg, "DVODPK-  ISTATE = 3 and NEQ increased (I1 to I2)  ", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__5, &c__1, &c__2, &dvod01_1.n, neq, &c__0, &zero, 
	    &zero, (ftnlen)80);
    goto L700;
L606:
    s_copy(msg, "DVODPK-  ITOL (=I1) illegal   ", (ftnlen)80, (ftnlen)30);
    xerrwd_(msg, &c__30, &c__6, &c__1, &c__1, itol, &c__0, &c__0, &zero, &
	    zero, (ftnlen)80);
    goto L700;
L607:
    s_copy(msg, "DVODPK-  IOPT (=I1) illegal   ", (ftnlen)80, (ftnlen)30);
    xerrwd_(msg, &c__30, &c__7, &c__1, &c__1, iopt, &c__0, &c__0, &zero, &
	    zero, (ftnlen)80);
    goto L700;
L608:
    s_copy(msg, "DVODPK-  MF (=I1) illegal     ", (ftnlen)80, (ftnlen)30);
    xerrwd_(msg, &c__30, &c__8, &c__1, &c__1, mf, &c__0, &c__0, &zero, &zero, 
	    (ftnlen)80);
    goto L700;
L611:
    s_copy(msg, "DVODPK-  MAXORD (=I1) .lt. 0  ", (ftnlen)80, (ftnlen)30);
    xerrwd_(msg, &c__30, &c__11, &c__1, &c__1, &dvod01_1.maxord, &c__0, &c__0,
	     &zero, &zero, (ftnlen)80);
    goto L700;
L612:
    s_copy(msg, "DVODPK-  MXSTEP (=I1) .lt. 0  ", (ftnlen)80, (ftnlen)30);
    xerrwd_(msg, &c__30, &c__12, &c__1, &c__1, &dvod01_1.mxstep, &c__0, &c__0,
	     &zero, &zero, (ftnlen)80);
    goto L700;
L613:
    s_copy(msg, "DVODPK-  MXHNIL (=I1) .lt. 0  ", (ftnlen)80, (ftnlen)30);
    xerrwd_(msg, &c__30, &c__13, &c__1, &c__1, &dvod01_1.mxhnil, &c__0, &c__0,
	     &zero, &zero, (ftnlen)80);
    goto L700;
L614:
    s_copy(msg, "DVODPK-  TOUT (=R1) behind T (=R2)      ", (ftnlen)80, (
	    ftnlen)40);
    xerrwd_(msg, &c__40, &c__14, &c__1, &c__0, &c__0, &c__0, &c__2, tout, t, (
	    ftnlen)80);
    s_copy(msg, "      integration direction is given by H0 (=R1)  ", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__14, &c__1, &c__0, &c__0, &c__0, &c__1, &h0, &
	    zero, (ftnlen)80);
    goto L700;
L615:
    s_copy(msg, "DVODPK-  HMAX (=R1) .lt. 0.0  ", (ftnlen)80, (ftnlen)30);
    xerrwd_(msg, &c__30, &c__15, &c__1, &c__0, &c__0, &c__0, &c__1, &hmax, &
	    zero, (ftnlen)80);
    goto L700;
L616:
    s_copy(msg, "DVODPK-  HMIN (=R1) .lt. 0.0  ", (ftnlen)80, (ftnlen)30);
    xerrwd_(msg, &c__30, &c__16, &c__1, &c__0, &c__0, &c__0, &c__1, &
	    dvod01_1.hmin, &zero, (ftnlen)80);
    goto L700;
L617:
    s_copy(msg, "DVODPK-  RWORK length needed, LENRW (=I1), exceeds LRW (=I2)"
	    , (ftnlen)80, (ftnlen)60);
    xerrwd_(msg, &c__60, &c__17, &c__1, &c__2, &lenrw, lrw, &c__0, &zero, &
	    zero, (ftnlen)80);
    goto L700;
L618:
    s_copy(msg, "DVODPK-  IWORK length needed, LENIW (=I1), exceeds LIW (=I2)"
	    , (ftnlen)80, (ftnlen)60);
    xerrwd_(msg, &c__60, &c__18, &c__1, &c__2, &leniw, liw, &c__0, &zero, &
	    zero, (ftnlen)80);
    goto L700;
L619:
    s_copy(msg, "DVODPK-  RTOL(I1) is R1 .lt. 0.0        ", (ftnlen)80, (
	    ftnlen)40);
    xerrwd_(msg, &c__40, &c__19, &c__1, &c__1, &i__, &c__0, &c__1, &rtoli, &
	    zero, (ftnlen)80);
    goto L700;
L620:
    s_copy(msg, "DVODPK-  ATOL(I1) is R1 .lt. 0.0        ", (ftnlen)80, (
	    ftnlen)40);
    xerrwd_(msg, &c__40, &c__20, &c__1, &c__1, &i__, &c__0, &c__1, &atoli, &
	    zero, (ftnlen)80);
    goto L700;
L621:
    ewti = rwork[dvod01_1.lewt + i__ - 1];
    s_copy(msg, "DVODPK-  EWT(I1) is R1 .le. 0.0         ", (ftnlen)80, (
	    ftnlen)40);
    xerrwd_(msg, &c__40, &c__21, &c__1, &c__1, &i__, &c__0, &c__1, &ewti, &
	    zero, (ftnlen)80);
    goto L700;
L622:
    s_copy(msg, "DVODPK-  TOUT (=R1) too close to T(=R2) to start integration"
	    , (ftnlen)80, (ftnlen)60);
    xerrwd_(msg, &c__60, &c__22, &c__1, &c__0, &c__0, &c__0, &c__2, tout, t, (
	    ftnlen)80);
    goto L700;
L623:
    s_copy(msg, "DVODPK-  ITASK = I1 and TOUT (=R1) behind TCUR - HU (= R2)  "
	    , (ftnlen)80, (ftnlen)60);
    xerrwd_(msg, &c__60, &c__23, &c__1, &c__1, itask, &c__0, &c__2, tout, &tp,
	     (ftnlen)80);
    goto L700;
L624:
    s_copy(msg, "DVODPK-  ITASK = 4 or 5 and TCRIT (=R1) behind TCUR (=R2)   "
	    , (ftnlen)80, (ftnlen)60);
    xerrwd_(msg, &c__60, &c__24, &c__1, &c__0, &c__0, &c__0, &c__2, &tcrit, &
	    dvod01_1.tn, (ftnlen)80);
    goto L700;
L625:
    s_copy(msg, "DVODPK-  ITASK = 4 or 5 and TCRIT (=R1) behind TOUT (=R2)   "
	    , (ftnlen)80, (ftnlen)60);
    xerrwd_(msg, &c__60, &c__25, &c__1, &c__0, &c__0, &c__0, &c__2, &tcrit, 
	    tout, (ftnlen)80);
    goto L700;
L626:
    s_copy(msg, "DVODPK-  At start of problem, too much accuracy   ", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__26, &c__1, &c__0, &c__0, &c__0, &c__0, &zero, &
	    zero, (ftnlen)80);
    s_copy(msg, "      requested for precision of machine:  See TOLSF (=R1)  "
	    , (ftnlen)80, (ftnlen)60);
    xerrwd_(msg, &c__60, &c__26, &c__1, &c__0, &c__0, &c__0, &c__1, &tolsf, &
	    zero, (ftnlen)80);
    rwork[14] = tolsf;
    goto L700;
L627:
    s_copy(msg, "DVODPK-  Trouble from DVINDY. ITASK = I1, TOUT = R1         "
	    , (ftnlen)80, (ftnlen)60);
    xerrwd_(msg, &c__60, &c__27, &c__1, &c__1, itask, &c__0, &c__1, tout, &
	    zero, (ftnlen)80);

L700:
    *istate = -3;
    return 0;

L800:
    s_copy(msg, "DVODPK-  Run aborted: apparent infinite loop      ", (ftnlen)
	    80, (ftnlen)50);
    xerrwd_(msg, &c__50, &c__303, &c__2, &c__0, &c__0, &c__0, &c__0, &zero, &
	    zero, (ftnlen)80);
    return 0;
/* ----------------------- End of Subroutine DVODPK ---------------------- */
} /* dvodpk_ */

/* DECK DVNLSK */
/* Subroutine */ int dvnlsk_(doublereal *y, doublereal *yh, integer *ldyh, 
	doublereal *vsav, doublereal *savf, doublereal *ewt, doublereal *acor,
	 integer *iwm, doublereal *wm, S_fp f, S_fp jac, U_fp psol, integer *
	nflag, doublereal *rpar, integer *ipar)
{
    /* Initialized data */

    static doublereal ccmax = .3;
    static doublereal crdown = .3;
    static integer maxcor = 3;
    static integer msbp = 20;
    static doublereal rdiv = 2.;
    static doublereal one = 1.;
    static doublereal two = 2.;
    static doublereal zero = 0.;

    /* System generated locals */
    integer yh_dim1, yh_offset, i__1;
    doublereal d__1, d__2;

    /* Local variables */
    integer i__, m;
    doublereal del, hrl1, dcon, delp;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    integer ierpj, iersl;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    doublereal cscale;
    extern /* Subroutine */ int dvslpk_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, S_fp, U_fp, 
	    integer *, doublereal *, integer *);
    extern doublereal dvnorm_(integer *, doublereal *, doublereal *);


/* ----------------------------------------------------------------------- */
/* Call sequence input -- YH, LDYH, F, JAC, EWT, ACOR, IWM, WM, */
/*                        NFLAG, RPAR, IPAR */
/* Call sequence output -- Y, YH, VSAV, SAVF, ACOR, IWM, WM, NFLAG */
/* COMMON block variables accessed: */
/*        /DVOD01/  ACNRM, CRATE, DRC, H, ICF, IPUP, JCUR, JSTART, */
/*                  METH, MITER, N, NSLP, RC, RL1, TN, TQ */
/*        /DVOD02/  NFE, NNI, NPE, NST */
/*        /DVPK01/  JACFLG, LOCIWP, LOCWP, MNEWT */
/* Subroutines called: F, JAC, PSOL, DAXPY, DCOPY, DSCAL, DVSLPK */
/* Function subroutines called: DVNORM */
/* ----------------------------------------------------------------------- */
/* Subroutine DVNLSK is a nonlinear system solver, which uses either */
/* functional iteration (MITER = 0), or a combination of an inexact */
/* Newton method and preconditioned Krylov iteration (MITER .gt. 0) */
/* to solve the implicit system for the corrector y vector. */
/* It calls Subroutine JAC (user-supplied) for preprocessing the */
/* preconditioner, and Subroutine DVSLPK for the Krylov iteration. */

/* In addition to variables described elsewhere, communication with */
/* DVNLSK uses the following variables: */

/* Y          = The dependent variable, a vector of length N, input. */
/* YH         = The Nordsieck (Taylor) array, LDYH by LMAX, input */
/*              and output.  On input, it contains predicted values. */
/* LDYH       = A constant .ge. N, the first dimension of YH, input. */
/* VSAV       = A work array of length N. */
/* SAVF       = A work array of length N. */
/* EWT        = An error weight vector of length N, input. */
/* ACOR       = A work array of length N, used for the accumulated */
/*              corrections to the predicted y vector. */
/* WM,IWM     = Real and integer work arrays associated with matrix */
/*              operations in Newton iteration (MITER .ne. 0). */
/* F          = Dummy name for user-supplied routine for f. */
/* JAC        = Dummy name for user-supplied routine for Jacobian data */
/*              and associated preconditioner matrix. */
/* PSOL       = Dummy name for user-supplied subroutine to solve */
/*              preconditioner linear system. */
/* NFLAG      = Input/output flag, with values and meanings as follows: */
/*              INPUT */
/*                  0 first call for this time step. */
/*                 -1 convergence failure in previous call to DVNLSK. */
/*                 -2 error test failure in DVSTEP. */
/*              OUTPUT */
/*                  0 successful completion of nonlinear solver. */
/*                 -1 convergence failure or failure in JAC. */
/*                 -2 unrecoverable error in matrix preprocessing */
/*                    (cannot occur here). */
/*                 -3 unrecoverable error in PSOL. */
/* RPAR, IPAR = Dummy names for user's real and integer work arrays. */

/* IPUP       = Own variable flag with values and meanings as follows: */
/*              0,            do not update preconditioner. */
/*              MITER .ne. 0, update the preconditioner, because it is */
/*                            the initial step, user input changed, */
/*                            there was an error test failure, or an */
/*                            update is indicated by a change in the */
/*                            scalar RC or step counter NST. */

/* For more details, see comments in driver subroutine. */
/* ----------------------------------------------------------------------- */

/* Type declarations for labeled COMMON block DVOD01 -------------------- */


/* Type declarations for labeled COMMON block DVOD02 -------------------- */


/* Type declarations for labeled COMMON block DVPK01 -------------------- */



/* Type declarations for local variables -------------------------------- */


/* Type declaration for function subroutines called --------------------- */

/* ----------------------------------------------------------------------- */
/* The following Fortran-77 declarations are to cause the values of the */
/* listed (local) variables to be saved between calls to DVODPK. */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --y;
    yh_dim1 = *ldyh;
    yh_offset = 1 + yh_dim1;
    yh -= yh_offset;
    --vsav;
    --savf;
    --ewt;
    --acor;
    --iwm;
    --wm;
    --rpar;
    --ipar;

    /* Function Body */
/* ----------------------------------------------------------------------- */
/* Up to MAXCOR corrector iterations are taken.  A convergence test is */
/* made on the RMS norm of each correction, weighted by the error */
/* weight vector EWT.  The sum of the corrections is accumulated in the */
/* vector ACOR(*).  The YH array is not altered in the corrector loop. */
/* ----------------------------------------------------------------------- */
    if (dvod01_1.jstart == 0) {
	dvod01_1.nslp = 0;
    }
    if (*nflag == 0) {
	dvod01_1.icf = 0;
    }
    if (*nflag == -2) {
	dvod01_1.ipup = dvod01_1.miter;
    }
    if (dvod01_1.jstart == 0 || dvod01_1.jstart == -1) {
	dvod01_1.ipup = dvod01_1.miter;
    }
    if (dvpk01_1.jacflg == 0) {
	dvod01_1.ipup = 0;
	dvod01_1.crate = one;
	goto L220;
    }
    dvod01_1.drc = (d__1 = dvod01_1.rc - one, abs(d__1));
    if (dvod01_1.drc > ccmax || dvod02_1.nst >= dvod01_1.nslp + msbp) {
	dvod01_1.ipup = dvod01_1.miter;
    }
L220:
    m = 0;
    hrl1 = dvod01_1.h__ * dvod01_1.rl1;
    delp = zero;
    dvpk01_1.mnewt = 0;
    dcopy_(&dvod01_1.n, &yh[yh_dim1 + 1], &c__1, &y[1], &c__1);
    (*f)(&dvod01_1.n, &dvod01_1.tn, &y[1], &savf[1], &rpar[1], &ipar[1]);
    ++dvod02_1.nfe;
    if (dvod01_1.ipup <= 0) {
	goto L250;
    }
/* ----------------------------------------------------------------------- */
/* If indicated, the preconditioner matrix is reevaluated and */
/* preprocessed before starting the corrector iteration.  IPUP is set */
/* to 0 as an indicator that this has been done. */
/* ----------------------------------------------------------------------- */
    dvod01_1.jcur = 1;
    ierpj = 0;
    (*jac)((S_fp)f, &dvod01_1.n, &dvod01_1.tn, &y[1], &yh[yh_offset], &ewt[1],
	     &savf[1], &acor[1], &hrl1, &wm[dvpk01_1.locwp], &iwm[
	    dvpk01_1.lociwp], &ierpj, &rpar[1], &ipar[1]);
    ++dvod02_3.npe;
    dvod01_1.ipup = 0;
    dvod01_1.rc = one;
    dvod01_1.drc = zero;
    dvod01_1.crate = one;
    dvod01_1.nslp = dvod02_1.nst;
    if (ierpj != 0) {
	goto L420;
    }
L250:
    i__1 = dvod01_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L260: */
	acor[i__] = zero;
    }
L270:
    if (dvod01_1.miter != 0) {
	goto L350;
    }
/* ----------------------------------------------------------------------- */
/* In the case of functional iteration, update Y directly from */
/* the result of the last function evaluation. */
/* ----------------------------------------------------------------------- */
    i__1 = dvod01_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	savf[i__] = dvod01_1.rl1 * (dvod01_1.h__ * savf[i__] - yh[i__ + (
		yh_dim1 << 1)]);
/* L290: */
	y[i__] = savf[i__] - acor[i__];
    }
    del = dvnorm_(&dvod01_1.n, &y[1], &ewt[1]);
    i__1 = dvod01_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L300: */
	y[i__] = yh[i__ + yh_dim1] + savf[i__];
    }
    dcopy_(&dvod01_1.n, &savf[1], &c__1, &acor[1], &c__1);
    goto L400;
/* ----------------------------------------------------------------------- */
/* In the case of the Newton method, compute the corrector error, */
/* and solve the linear system with that as right-hand side and */
/* A as coefficient matrix.  In the case of Modified Newton iteration */
/* with BDF, the correction is scaled by the factor 2/(1+RC) to */
/* account for changes in H*RL1 since the last JAC call. */
/* ----------------------------------------------------------------------- */
L350:
    i__1 = dvod01_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L360: */
	vsav[i__] = hrl1 * savf[i__] - (dvod01_1.rl1 * yh[i__ + (yh_dim1 << 1)
		] + acor[i__]);
    }
    dvslpk_(&y[1], &savf[1], &vsav[1], &ewt[1], &wm[1], &iwm[1], (S_fp)f, (
	    U_fp)psol, &iersl, &rpar[1], &ipar[1]);
    ++dvod02_1.nni;
    if (dvod01_1.meth == 2 && dvpk01_1.jacflg == 1 && dvod01_1.miter == 9 && 
	    dvod01_1.rc != one) {
	cscale = two / (one + dvod01_1.rc);
	dscal_(&dvod01_1.n, &cscale, &vsav[1], &c__1);
    }
    if (iersl < 0) {
	goto L440;
    }
    if (iersl > 0) {
	goto L410;
    }
    del = dvnorm_(&dvod01_1.n, &vsav[1], &ewt[1]);
    daxpy_(&dvod01_1.n, &one, &vsav[1], &c__1, &acor[1], &c__1);
    i__1 = dvod01_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L380: */
	y[i__] = yh[i__ + yh_dim1] + acor[i__];
    }
/* ----------------------------------------------------------------------- */
/* Test for convergence.  If M.gt.0, an estimate of the convergence */
/* rate constant is stored in CRATE, and this is used in the test. */
/* ----------------------------------------------------------------------- */
L400:
    if (m != 0) {
/* Computing MAX */
	d__1 = crdown * dvod01_1.crate, d__2 = del / delp;
	dvod01_1.crate = max(d__1,d__2);
    }
    dcon = del * min(one,dvod01_1.crate) / dvod01_1.tq[3];
    if (dcon <= one) {
	goto L450;
    }
    ++m;
    if (m == maxcor) {
	goto L410;
    }
    if (m >= 2 && del > rdiv * delp) {
	goto L410;
    }
    dvpk01_1.mnewt = m;
    delp = del;
    (*f)(&dvod01_1.n, &dvod01_1.tn, &y[1], &savf[1], &rpar[1], &ipar[1]);
    ++dvod02_1.nfe;
    goto L270;

L410:
    if (dvod01_1.miter == 0 || dvod01_1.jcur == 1 || dvpk01_1.jacflg == 0) {
	goto L420;
    }
    dvod01_1.icf = 1;
    dvod01_1.ipup = dvod01_1.miter;
    goto L220;

L420:
    dvod01_1.icf = 2;
    *nflag = -1;
    return 0;
L440:
    *nflag = -3;
    return 0;
/* Return for successful step. ------------------------------------------ */
L450:
    *nflag = 0;
    dvod01_1.jcur = 0;
    dvod01_1.icf = 0;
    if (m == 0) {
	dvod01_1.acnrm = del;
    }
    if (m > 0) {
	dvod01_1.acnrm = dvnorm_(&dvod01_1.n, &acor[1], &ewt[1]);
    }
    return 0;
/* ----------------------- End of Subroutine DVNLSK ---------------------- */
} /* dvnlsk_ */

/* DECK DVSLPK */
/* Subroutine */ int dvslpk_(doublereal *y, doublereal *savf, doublereal *x, 
	doublereal *ewt, doublereal *wm, integer *iwm, S_fp f, U_fp psol, 
	integer *iersl, doublereal *rpar, integer *ipar)
{
    /* System generated locals */
    integer i__1, i__2;

    /* Local variables */
    integer lb, lq, lv, ldl, lwk;
    doublereal hrl1;
    integer lgmr, lhes, npsl, iflag;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    doublereal delta;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    integer maxlp1;
    extern /* Subroutine */ int dvspig_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *, 
	    integer *, integer *, doublereal *, doublereal *, integer *, 
	    integer *, S_fp, U_fp, integer *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *, 
	    doublereal *, doublereal *, doublereal *, integer *, integer *), 
	    dvusol_(integer *, doublereal *, doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, doublereal *, integer *,
	     integer *, U_fp, integer *, doublereal *, doublereal *, integer *
	    , doublereal *, doublereal *, integer *, integer *);

/* ----------------------------------------------------------------------- */
/* Call sequence input -- Y, SAVF, X, EWT, F, PSOL, RPAR, IPAR */
/* Call sequence output -- Y, SAVF, X, WM, IWM, IERSL */
/* COMMON block variables accessed: */
/*        /DVOD01/  H, RL1, TQ, TN, MITER, N */
/*        /DVPK01/  DELT, SQRTN, RSQRTN, JPRE, LOCIWP, LOCWP, */
/*                  KMP, MAXL, MNEWT, NLI, NPS, NCFL */
/* Subroutines called: F, PSOL, DCOPY, DSCAL, DVSPIG, DVUSOL */
/* ----------------------------------------------------------------------- */
/* This routine interfaces with  DVSPIG  or  DVUSOL  for the solution of */
/* the linear system arising from a Newton iteration (MITER .ne. 0). */

/* In addition to variables described elsewhere, communication with */
/* DVSLPK uses the following variables: */
/* WM    = real work space containing data for the algorithm */
/*         (Krylov basis vectors, Hessenberg matrix, etc.) */
/* IWM   = integer work space containing data for the algorithm */
/* X     = the right-hand side vector on input, and the solution vector */
/*         on output, of length N. */
/* IERSL = output flag (in COMMON): */
/*         IERSL =  0 means no trouble occurred. */
/*         IERSL =  1 means the iterative method failed to converge. */
/*                    If the preconditioner is out of date, the step */
/*                    is repeated with a new preconditioner.  Otherwise, */
/*                    the stepsize is reduced (forcing a new evalua- */
/*                    tion of the preconditioner) and the step is */
/*                    repeated. */
/*         IERSL = -1 means there was a nonrecoverable error in the */
/*                    iterative solver.  The stepsize is reduced in */
/*                    DVSTEP and the step is repeated. */
/* ----------------------------------------------------------------------- */

/* Type declarations for labeled COMMON block DVOD01 -------------------- */


/* Type declarations for labeled COMMON block DVPK01 -------------------- */


/* Type declarations for local variables -------------------------------- */

/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --ipar;
    --rpar;
    --iwm;
    --wm;
    --ewt;
    --x;
    --savf;
    --y;

    /* Function Body */
    *iersl = 0;
    hrl1 = dvod01_1.h__ * dvod01_1.rl1;
    delta = dvpk01_1.delt * dvod01_1.tq[3];
    if (dvod01_1.miter == 1) {
/* ----------------------------------------------------------------------- */
/* Use the SPIGMR algorithm to solve the linear system A*x = -f. */
/* ----------------------------------------------------------------------- */
	maxlp1 = dvpk01_1.maxl + 1;
	lv = 1;
	lb = lv + dvod01_1.n * dvpk01_1.maxl;
	lhes = lb + dvod01_1.n + 1;
	lq = lhes + dvpk01_1.maxl * maxlp1;
	lwk = lq + (dvpk01_1.maxl << 1);
/* Computing MIN */
	i__1 = 1, i__2 = dvpk01_1.maxl - dvpk01_1.kmp;
	ldl = lwk + min(i__1,i__2) * dvod01_1.n;
	dcopy_(&dvod01_1.n, &x[1], &c__1, &wm[lb], &c__1);
	dscal_(&dvod01_1.n, &dvpk01_1.rsqrtn, &ewt[1], &c__1);
	dvspig_(&dvod01_1.tn, &y[1], &savf[1], &wm[lb], &ewt[1], &dvod01_1.n, 
		&dvpk01_1.maxl, &maxlp1, &dvpk01_1.kmp, &delta, &hrl1, &
		dvpk01_1.jpre, &dvpk01_1.mnewt, (S_fp)f, (U_fp)psol, &npsl, &
		x[1], &wm[lv], &wm[lhes], &wm[lq], &lgmr, &wm[dvpk01_1.locwp],
		 &iwm[dvpk01_1.lociwp], &wm[lwk], &wm[ldl], &rpar[1], &ipar[1]
		, &iflag);
	dvpk01_1.nli += lgmr;
	dvpk01_1.nps += npsl;
	dscal_(&dvod01_1.n, &dvpk01_1.sqrtn, &ewt[1], &c__1);
	if (iflag != 0) {
	    ++dvpk01_1.ncfl;
	}
	if (iflag >= 2) {
	    *iersl = 1;
	}
	if (iflag < 0) {
	    *iersl = -1;
	}
	return 0;
    } else if (dvod01_1.miter == 9) {
/* ----------------------------------------------------------------------- */
/* Use DVUSOL, which interfaces to PSOL, to solve the linear system */
/* (No Krylov iteration). */
/* ----------------------------------------------------------------------- */
	lb = 1;
	lwk = lb + dvod01_1.n;
	dcopy_(&dvod01_1.n, &x[1], &c__1, &wm[lb], &c__1);
	dvusol_(&dvod01_1.n, &dvod01_1.tn, &y[1], &savf[1], &wm[lb], &ewt[1], 
		&delta, &hrl1, &dvpk01_1.jpre, &dvpk01_1.mnewt, (U_fp)psol, &
		npsl, &x[1], &wm[dvpk01_1.locwp], &iwm[dvpk01_1.lociwp], &wm[
		lwk], &rpar[1], &ipar[1], &iflag);
	dvpk01_1.nps += npsl;
	if (iflag != 0) {
	    ++dvpk01_1.ncfl;
	}
	if (iflag == 3) {
	    *iersl = 1;
	}
	if (iflag < 0) {
	    *iersl = -1;
	}
	return 0;
    }
/* ----------------------- End of Subroutine DVSLPK ---------------------- */
    return 0;
} /* dvslpk_ */

/* DECK DVSPIG */
/* Subroutine */ int dvspig_(doublereal *tn, doublereal *y, doublereal *savf, 
	doublereal *b, doublereal *wght, integer *n, integer *maxl, integer *
	maxlp1, integer *kmp, doublereal *delta, doublereal *hb0, integer *
	jpre, integer *mnewt, S_fp f, S_fp psol, integer *npsl, doublereal *x,
	 doublereal *v, doublereal *hes, doublereal *q, integer *lgmr, 
	doublereal *wp, integer *iwp, doublereal *wk, doublereal *dl, 
	doublereal *rpar, integer *ipar, integer *iflag)
{
    /* System generated locals */
    integer v_dim1, v_offset, hes_dim1, hes_offset, i__1, i__2, i__3;
    doublereal d__1;

    /* Local variables */
    doublereal c__;
    integer i__, j, k;
    doublereal s;
    integer i2, ll, ip1, ier;
    doublereal tem, rho;
    integer llp1, info;
    doublereal bnrm, prod, bnrm0;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *), dhels_(doublereal *, integer *, integer *, doublereal 
	    *, doublereal *), dheqr_(doublereal *, integer *, integer *, 
	    doublereal *, integer *, integer *);
    doublereal dlnrm;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *), dvatv_(doublereal *, doublereal *, 
	    doublereal *, doublereal *, doublereal *, S_fp, S_fp, doublereal *
	    , integer *, doublereal *, doublereal *, doublereal *, integer *, 
	    doublereal *, integer *, integer *, integer *), daxpy_(integer *, 
	    doublereal *, doublereal *, integer *, doublereal *, integer *);
    doublereal snormw;
    extern /* Subroutine */ int dorthog_(doublereal *, doublereal *, 
	    doublereal *, integer *, integer *, integer *, integer *, 
	    doublereal *);

/* ----------------------------------------------------------------------- */
/* Call sequence input --  TN, Y, SAVF, B, WGHT, N, MAXL, MAXLP1, DELTA, */
/*                         HB0, JPRE, MNEWT, F, PSOL, RPAR, IPAR */
/* Call sequence output -- B, KMP, DELTA, NPSL, X, V, HES, Q, LGMR, WP, */
/*                         IWP, WK, DL, RPAR, IPAR, IFLAG */
/* COMMON block variables accessed: None */
/* Subroutines called: F, DORTHOG, PSOL, DAXPY, DCOPY, DHELS, DHEQR, */
/*                      DSCAL, DVATV */
/* Function subroutines called: DNRM2 */
/* ----------------------------------------------------------------------- */
/* This routine solves the linear system  A * x = b using SPIGMR, */
/* a scaled preconditioned incomplete version of the generalized */
/* minimum residual method GMRES. */
/* An initial guess of x = 0 is assumed. */
/* ----------------------------------------------------------------------- */

/*      On entry */

/*           TN = current value of t. */

/*            Y = array containing current dependent variable vector. */

/*         SAVF = array containing current value of f(t,y). */

/*            B = the right hand side of the system A*x = b. */
/*                B is also used as work space when computing */
/*                the final approximation. */
/*                (B is the same as V(*,MAXL+1) in the call to DVSPIG.) */

/*         WGHT = the vector of length N containing the nonzero */
/*                elements of the diagonal scaling matrix. */

/*            N = the order of the matrix A, and the lengths */
/*                of the vectors WGHT, B and X. */

/*         MAXL = the maximum allowable order of the matrix HES. */

/*       MAXLP1 = MAXL + 1, used for dynamic dimensioning of HES. */

/*          KMP = the number of previous vectors the new vector VNEW */
/*                must be made orthogonal to.  KMP .le. MAXL. */

/*        DELTA = tolerance on residuals  b - A*x  in weighted RMS norm. */

/*          HB0 = current value of (step size h) * (coefficient beta0). */

/*         JPRE = preconditioner type flag. */

/*        MNEWT = Newton iteration counter (.ge. 0). */

/*           WK = real work array used by routine DVATV and PSOL. */

/*           DL = real work array used for calculation of the residual */
/*                norm rho when the method is incomplete (KMP.lt.MAXL). */

/*           WP = real work array used by preconditioner PSOL. */

/*          IWP = integer work array used by preconditioner PSOL. */

/*      On return */

/*         X    = the final computed approximation to the solution */
/*                of the system A*x = b. */

/*         LGMR = the number of iterations performed and the current */
/*                order of the upper Hessenberg matrix HES. */

/*         NPSL = the number of calls to PSOL. */

/*         V    = the N by (LGMR+1) array containing the LGMR */
/*                orthogonal vectors V(*,1) to V(*,LGMR). */

/*         HES  = the upper triangular factor of the QR decomposition */
/*                of the (LGMR+1) by LGMR upper Hessenberg matrix whose */
/*                entries are the scaled inner-products of A*V(*,i) */
/*                and V(*,k). */

/*         Q    = real array of length 2*MAXL containing the components */
/*                of the Givens rotations used in the QR decomposition */
/*                of HES.  It is loaded in DHEQR and used in DHELS. */

/*        IFLAG = integer error flag: */
/*                0 means convergence in LGMR iterations, LGMR.le.MAXL. */
/*                1 means the convergence test did not pass in MAXL */
/*                  iterations, but the residual norm is .lt. 1, */
/*                  or .lt. norm(b) if MNEWT = 0, and so x is computed. */
/*                2 means the convergence test did not pass in MAXL */
/*                  iterations, residual .gt. 1, and x is undefined. */
/*                3 means there was a recoverable error in PSOL */
/*                  caused by the preconditioner being out of date. */
/*               -1 means there was a nonrecoverable error in PSOL. */

/* ----------------------------------------------------------------------- */

/* Type declarations for local variables -------------------------------- */


/* Type declaration for function subroutines called --------------------- */


    /* Parameter adjustments */
    --y;
    --savf;
    --b;
    --wght;
    v_dim1 = *n;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    hes_dim1 = *maxlp1;
    hes_offset = 1 + hes_dim1;
    hes -= hes_offset;
    --x;
    --q;
    --wp;
    --iwp;
    --wk;
    --dl;
    --rpar;
    --ipar;

    /* Function Body */
    *iflag = 0;
    *lgmr = 0;
    *npsl = 0;
/* ----------------------------------------------------------------------- */
/* The initial residual is the vector b.  Apply scaling to b, and test */
/* for an immediate return with x = 0 or x = b. */
/* ----------------------------------------------------------------------- */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	v[i__ + v_dim1] = b[i__] * wght[i__];
    }
    bnrm0 = dnrm2_(n, &v[v_offset], &c__1);
    bnrm = bnrm0;
    if (bnrm0 > *delta) {
	goto L30;
    }
    if (*mnewt > 0) {
	goto L20;
    }
    dcopy_(n, &b[1], &c__1, &x[1], &c__1);
    return 0;
L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L25: */
	x[i__] = 0.;
    }
    return 0;
L30:
/* Apply inverse of left preconditioner to vector b. -------------------- */
    ier = 0;
    if (*jpre == 0 || *jpre == 2) {
	goto L55;
    }
    (*psol)(n, tn, &y[1], &savf[1], &wk[1], hb0, &wp[1], &iwp[1], &b[1], &
	    c__1, &ier, &rpar[1], &ipar[1]);
    *npsl = 1;
    if (ier != 0) {
	goto L300;
    }
/* Calculate norm of scaled vector V(*, 1) and normalize it. ------------ */
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L50: */
	v[i__ + v_dim1] = b[i__] * wght[i__];
    }
    bnrm = dnrm2_(n, &v[v_offset], &c__1);
    *delta *= bnrm / bnrm0;
L55:
    tem = 1. / bnrm;
    dscal_(n, &tem, &v[v_dim1 + 1], &c__1);
/* Zero out the HES array. ---------------------------------------------- */
    i__1 = *maxl;
    for (j = 1; j <= i__1; ++j) {
	i__2 = *maxlp1;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L60: */
	    hes[i__ + j * hes_dim1] = 0.;
	}
/* L65: */
    }
/* ----------------------------------------------------------------------- */
/* Main loop to compute the vectors V(*,2) to V(*,MAXL). */
/* The running product PROD is needed for the convergence test. */
/* ----------------------------------------------------------------------- */
    prod = 1.;
    i__1 = *maxl;
    for (ll = 1; ll <= i__1; ++ll) {
	*lgmr = ll;
/* ----------------------------------------------------------------------- */
/* Call routine DVATV to compute VNEW = Abar*v(ll), where Abar is */
/* the matrix A with scaling and inverse preconditioner factors applied. */
/* Call routine DORTHOG to orthogonalize the new vector VNEW = V(*,LL+1). */
/* Call routine DHEQR to update the factors of HES. */
/* ----------------------------------------------------------------------- */
	dvatv_(&y[1], &savf[1], &v[ll * v_dim1 + 1], &wght[1], &x[1], (S_fp)f,
		 (S_fp)psol, &rpar[1], &ipar[1], &v[(ll + 1) * v_dim1 + 1], &
		wk[1], &wp[1], &iwp[1], hb0, jpre, &ier, npsl);
	if (ier != 0) {
	    goto L300;
	}
	dorthog_(&v[(ll + 1) * v_dim1 + 1], &v[v_offset], &hes[hes_offset], n,
		 &ll, maxlp1, kmp, &snormw);
	hes[ll + 1 + ll * hes_dim1] = snormw;
	dheqr_(&hes[hes_offset], maxlp1, &ll, &q[1], &info, &ll);
	if (info == ll) {
	    goto L120;
	}
/* ----------------------------------------------------------------------- */
/* Update RHO, the estimate of the norm of the residual b - A*xl. */
/* If KMP .lt. MAXL, then the vectors V(*,1),...,V(*,LL+1) are not */
/* necessarily orthogonal for LL .gt. KMP.  The vector DL must then */
/* be computed, and its norm used in the calculation of RHO. */
/* ----------------------------------------------------------------------- */
	prod *= q[ll * 2];
	rho = (d__1 = prod * bnrm, abs(d__1));
	if (ll > *kmp && *kmp < *maxl) {
	    if (ll == *kmp + 1) {
		dcopy_(n, &v[v_dim1 + 1], &c__1, &dl[1], &c__1);
		i__2 = *kmp;
		for (i__ = 1; i__ <= i__2; ++i__) {
		    ip1 = i__ + 1;
		    i2 = i__ << 1;
		    s = q[i2];
		    c__ = q[i2 - 1];
		    i__3 = *n;
		    for (k = 1; k <= i__3; ++k) {
/* L70: */
			dl[k] = s * dl[k] + c__ * v[k + ip1 * v_dim1];
		    }
/* L75: */
		}
	    }
	    s = q[ll * 2];
	    c__ = q[(ll << 1) - 1] / snormw;
	    llp1 = ll + 1;
	    i__2 = *n;
	    for (k = 1; k <= i__2; ++k) {
/* L80: */
		dl[k] = s * dl[k] + c__ * v[k + llp1 * v_dim1];
	    }
	    dlnrm = dnrm2_(n, &dl[1], &c__1);
	    rho *= dlnrm;
	}
/* ----------------------------------------------------------------------- */
/* Test for convergence.  If passed, compute approximation xl. */
/* If failed and LL .lt. MAXL, then continue iterating. */
/* ----------------------------------------------------------------------- */
	if (rho <= *delta) {
	    goto L200;
	}
	if (ll == *maxl) {
	    goto L100;
	}
/* ----------------------------------------------------------------------- */
/* Rescale so that the norm of V(1,LL+1) is one. */
/* ----------------------------------------------------------------------- */
	tem = 1. / snormw;
	dscal_(n, &tem, &v[(ll + 1) * v_dim1 + 1], &c__1);
/* L90: */
    }
L100:
    if (rho <= 1.) {
	goto L150;
    }
    if (rho <= bnrm && *mnewt == 0) {
	goto L150;
    }
L120:
    *iflag = 2;
    return 0;
L150:
    *iflag = 1;
/* ----------------------------------------------------------------------- */
/* Compute the approximation xl to the solution. */
/* Since the vector X was used as work space, and the initial guess */
/* of the Newton correction is zero, X must be reset to zero. */
/* ----------------------------------------------------------------------- */
L200:
    ll = *lgmr;
    llp1 = ll + 1;
    i__1 = llp1;
    for (k = 1; k <= i__1; ++k) {
/* L210: */
	b[k] = 0.;
    }
    b[1] = bnrm;
    dhels_(&hes[hes_offset], maxlp1, &ll, &q[1], &b[1]);
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
/* L220: */
	x[k] = 0.;
    }
    i__1 = ll;
    for (i__ = 1; i__ <= i__1; ++i__) {
	daxpy_(n, &b[i__], &v[i__ * v_dim1 + 1], &c__1, &x[1], &c__1);
/* L230: */
    }
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L240: */
	x[i__] /= wght[i__];
    }
    if (*jpre <= 1) {
	return 0;
    }
    (*psol)(n, tn, &y[1], &savf[1], &wk[1], hb0, &wp[1], &iwp[1], &x[1], &
	    c__2, &ier, &rpar[1], &ipar[1]);
    ++(*npsl);
    if (ier != 0) {
	goto L300;
    }
    return 0;
/* ----------------------------------------------------------------------- */
/* This block handles error returns forced by routine PSOL. */
/* ----------------------------------------------------------------------- */
L300:
    if (ier < 0) {
	*iflag = -1;
    }
    if (ier > 0) {
	*iflag = 3;
    }

    return 0;
/* ----------------------- End of Subroutine DVSPIG ---------------------- */
} /* dvspig_ */

/* DECK DVATV */
/* Subroutine */ int dvatv_(doublereal *y, doublereal *savf, doublereal *v, 
	doublereal *wght, doublereal *ftem, S_fp f, S_fp psol, doublereal *
	rpar, integer *ipar, doublereal *z__, doublereal *vtem, doublereal *
	wp, integer *iwp, doublereal *hb0, integer *jpre, integer *ier, 
	integer *npsl)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__;
    doublereal fac;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    doublereal tempn, rnorm;

/* ----------------------------------------------------------------------- */
/* Call sequence input -- Y, SAVF, V, WGHT, F, PSOL, RPAR, IPAR, */
/*                        WP, IWP, HB0, NPSL */
/* Call sequence output --Z, IER, NPSL */
/* COMMON block variables accessed: */
/*        /DVOD01/  TN, N */
/*        /DVOD02/  NFE */
/* Subroutines called: F, PSOL, DCOPY */
/* Function subroutines called: DNRM2 */
/* ----------------------------------------------------------------------- */
/* This routine computes the product */

/*   (D-inverse)*(P1-inverse)*(I - hb0*df/dy)*(P2-inverse)*(D*v), */

/* where D is a diagonal scaling matrix, and P1 and P2 are the */
/* left and right preconditioning matrices, respectively. */
/* v is assumed to have L2 norm equal to 1. */
/* The product is stored in Z.  This is computed by a */
/* difference quotient, a call to F, and two calls to PSOL. */
/* ----------------------------------------------------------------------- */

/*      On entry */

/*            Y = array containing current dependent variable vector. */

/*         SAVF = array containing current value of f(t,y). */

/*            V = real array of length N (can be the same array as Z). */

/*         WGHT = array of length N containing scale factors. */
/*                1/WGHT(i) are the diagonal elements of the matrix D. */

/*         FTEM = work array of length N. */

/*         VTEM = work array of length N used to store the */
/*                unscaled version of v. */

/*           WP = real work array used by preconditioner PSOL. */

/*          IWP = integer work array used by preconditioner PSOL. */

/*          HB0 = current value of (step size h) * (coefficient beta0). */

/*         JPRE = preconditioner type flag. */


/*      On return */

/*            Z = array of length N containing desired scaled */
/*                matrix-vector product. */

/*          IER = error flag from PSOL. */

/*         NPSL = the number of calls to PSOL. */
/* ----------------------------------------------------------------------- */

/* Type declarations for labeled COMMON block DVOD01 -------------------- */


/* Type declarations for labeled COMMON block DVOD02 -------------------- */


/* Type declarations for local variables -------------------------------- */


/* Type declaration for function subroutines called --------------------- */


/* ----------------------------------------------------------------------- */
/* Set vtem = D * v. ---------------------------------------------------- */
    /* Parameter adjustments */
    --iwp;
    --wp;
    --vtem;
    --z__;
    --ipar;
    --rpar;
    --ftem;
    --wght;
    --v;
    --savf;
    --y;

    /* Function Body */
    i__1 = dvod01_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	vtem[i__] = v[i__] / wght[i__];
    }
    *ier = 0;
    if (*jpre >= 2) {
	goto L30;
    }

/* JPRE = 0 or 1.  Save y in Z and increment Y by VTEM. ----------------- */
    dcopy_(&dvod01_1.n, &y[1], &c__1, &z__[1], &c__1);
    i__1 = dvod01_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L20: */
	y[i__] = z__[i__] + vtem[i__];
    }
    fac = *hb0;
    goto L60;

/* JPRE = 2 or 3.  Apply inverse of right preconditioner to VTEM. ------- */
L30:
    (*psol)(&dvod01_1.n, &dvod01_1.tn, &y[1], &savf[1], &ftem[1], hb0, &wp[1],
	     &iwp[1], &vtem[1], &c__2, ier, &rpar[1], &ipar[1]);
    ++(*npsl);
    if (*ier != 0) {
	return 0;
    }
/* Calculate l-2 norm of (D-inverse) * VTEM. ---------------------------- */
    i__1 = dvod01_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L40: */
	z__[i__] = vtem[i__] * wght[i__];
    }
    tempn = dnrm2_(&dvod01_1.n, &z__[1], &c__1);
    rnorm = 1. / tempn;
/* Save y in Z and increment Y by VTEM/norm. ---------------------------- */
    dcopy_(&dvod01_1.n, &y[1], &c__1, &z__[1], &c__1);
    i__1 = dvod01_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L50: */
	y[i__] = z__[i__] + vtem[i__] * rnorm;
    }
    fac = *hb0 * tempn;

/* For all JPRE, call F with incremented Y argument, and restore Y. ----- */
L60:
    (*f)(&dvod01_1.n, &dvod01_1.tn, &y[1], &ftem[1], &rpar[1], &ipar[1]);
    ++dvod02_1.nfe;
    dcopy_(&dvod01_1.n, &z__[1], &c__1, &y[1], &c__1);
/* Set Z = (I - HB0*Jacobian) * VTEM, using difference quotient. -------- */
    i__1 = dvod01_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L70: */
	z__[i__] = ftem[i__] - savf[i__];
    }
    i__1 = dvod01_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L80: */
	z__[i__] = vtem[i__] - fac * z__[i__];
    }
/* Apply inverse of left preconditioner to Z, if nontrivial. ------------ */
    if (*jpre == 0 || *jpre == 2) {
	goto L85;
    }
    (*psol)(&dvod01_1.n, &dvod01_1.tn, &y[1], &savf[1], &ftem[1], hb0, &wp[1],
	     &iwp[1], &z__[1], &c__1, ier, &rpar[1], &ipar[1]);
    ++(*npsl);
    if (*ier != 0) {
	return 0;
    }
L85:
/* Apply D-inverse to Z and return. ------------------------------------- */
    i__1 = dvod01_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L90: */
	z__[i__] *= wght[i__];
    }
    return 0;
/* ----------------------- End of Subroutine DVATV ----------------------- */
} /* dvatv_ */

/* DECK DVUSOL */
/* Subroutine */ int dvusol_(integer *n, doublereal *tn, doublereal *y, 
	doublereal *savf, doublereal *b, doublereal *wght, doublereal *delta, 
	doublereal *hb0, integer *jpre, integer *mnewt, S_fp psol, integer *
	npsl, doublereal *x, doublereal *wp, integer *iwp, doublereal *wk, 
	doublereal *rpar, integer *ipar, integer *iflag)
{
    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, ier;
    doublereal bnrm;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    extern doublereal dvnorm_(integer *, doublereal *, doublereal *);

/* ----------------------------------------------------------------------- */
/* This routine solves the linear system A * x = b using only */
/* calls to the user-supplied routine PSOL (no Krylov iteration). */
/* If the norm of the right-hand side vector b is smaller than DELTA, */
/* the vector x returned is x = b (if MNEWT = 0) or x = 0 otherwise. */
/* PSOL is called with an LR argument of 1 (if JPRE = 1 or 3), */
/* then 2 (if JPRE = 2 or 3). */
/* ----------------------------------------------------------------------- */

/*      On entry */

/*          NEQ = problem size, passed to F and PSOL (NEQ(1) = N). */

/*           TN = current value of t. */

/*            Y = array containing current dependent variable vector. */

/*         SAVF = array containing current value of f(t,y). */

/*            B = the right hand side of the system A*x = b. */

/*         WGHT = the vector of length N containing the nonzero */
/*                elements of the diagonal scaling matrix. */

/*            N = the order of the matrix A, and the lengths */
/*                of the vectors WGHT, b and x. */

/*        DELTA = tolerance on residuals  b - A*x  in weighted RMS norm. */

/*          HB0 = current value of (step size h) * (coefficient beta0). */

/*         JPRE = preconditioner type flag. */

/*        MNEWT = Newton iteration counter (.ge. 0). */

/*           WK = real work array used by PSOL. */

/*           WP = real work array used by preconditioner PSOL. */

/*          IWP = integer work array used by preconditioner PSOL. */

/*      On return */

/*         X    = the final computed approximation to the solution */
/*                of the system A*x = b. */

/*         NPSL = the number of calls to PSOL. */

/*        IFLAG = integer error flag: */
/*                0 means no trouble occurred. */
/*                3 means there was a recoverable error in PSOL */
/*                  caused by the preconditioner being out of date. */
/*               -1 means there was a nonrecoverable error in PSOL. */

/* ----------------------------------------------------------------------- */

/* Type declarations for local variables -------------------------------- */


/* Type declaration for function subroutines called --------------------- */


    /* Parameter adjustments */
    --ipar;
    --rpar;
    --wk;
    --iwp;
    --wp;
    --x;
    --wght;
    --b;
    --savf;
    --y;

    /* Function Body */
    *iflag = 0;
    *npsl = 0;
/* ----------------------------------------------------------------------- */
/* Test for an immediate return with x = 0 or x = b. */
/* ----------------------------------------------------------------------- */
    bnrm = dvnorm_(n, &b[1], &wght[1]);
    if (bnrm > *delta) {
	goto L30;
    }
    if (*mnewt > 0) {
	goto L10;
    }
    dcopy_(n, &b[1], &c__1, &x[1], &c__1);
    return 0;
L10:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L20: */
	x[i__] = 0.;
    }
    return 0;
/* Apply inverse of left preconditioner to vector b. -------------------- */
L30:
    ier = 0;
    if (*jpre == 0 || *jpre == 2) {
	goto L40;
    }
    (*psol)(n, tn, &y[1], &savf[1], &wk[1], hb0, &wp[1], &iwp[1], &b[1], &
	    c__1, &ier, &rpar[1], &ipar[1]);
    *npsl = 1;
    if (ier != 0) {
	goto L100;
    }
/* Apply inverse of right preconditioner to result, and copy to X. ------ */
L40:
    if (*jpre <= 1) {
	goto L50;
    }
    (*psol)(n, tn, &y[1], &savf[1], &wk[1], hb0, &wp[1], &iwp[1], &b[1], &
	    c__2, &ier, &rpar[1], &ipar[1]);
    ++(*npsl);
    if (ier != 0) {
	goto L100;
    }
L50:
    dcopy_(n, &b[1], &c__1, &x[1], &c__1);
    return 0;
/* ----------------------------------------------------------------------- */
/* This block handles error returns forced by routine PSOL. */
/* ----------------------------------------------------------------------- */
L100:
    if (ier < 0) {
	*iflag = -1;
    }
    if (ier > 0) {
	*iflag = 3;
    }
    return 0;
/* ----------------------- End of Subroutine DVUSOL ---------------------- */
} /* dvusol_ */

/* DECK DVKSRC */
/* Subroutine */ int dvksrc_(doublereal *rsav, integer *isav, integer *job)
{
    /* Initialized data */

    static integer lenrv1 = 48;
    static integer leniv1 = 33;
    static integer lenrv2 = 1;
    static integer leniv2 = 8;
    static integer lrvk1 = 3;
    static integer livk1 = 11;

    /* System generated locals */
    integer i__1;

    /* Local variables */
    integer i__, ioff;

/* ----------------------------------------------------------------------- */
/* Call sequence input -- RSAV, ISAV, JOB */
/* Call sequence output -- RSAV, ISAV */
/* COMMON block variables accessed: all of /DVOD01/, /DVOD02/, /DVPK01/ */

/* Subroutines/functions called by DVKSRC: None */
/* ----------------------------------------------------------------------- */
/* This routine saves or restores (depending on JOB) the contents of the */
/* COMMON blocks DVOD01, DVOD02, DVPK01, used internally by DVODPK. */

/* RSAV = real array of length 52 or more. */
/* ISAV = integer array of length 52 or more. */
/* JOB  = flag indicating to save or restore the COMMON blocks: */
/*        JOB  = 1 if COMMON is to be saved (written to RSAV/ISAV). */
/*        JOB  = 2 if COMMON is to be restored (read from RSAV/ISAV). */
/*        A call with JOB = 2 presumes a prior call with JOB = 1. */
/* ----------------------------------------------------------------------- */

/* Type declarations for labeled COMMON block DVOD01 -------------------- */


/* Type declarations for labeled COMMON block DVOD02 -------------------- */


/* Type declarations for labeled COMMON block DVPK01 -------------------- */


/* Type declarations for local variables -------------------------------- */

/* ----------------------------------------------------------------------- */
/* The following Fortran-77 declaration is to cause the values of the */
/* listed (local) variables to be saved between calls to this integrator. */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --isav;
    --rsav;

    /* Function Body */

    if (*job == 2) {
	goto L100;
    }
    i__1 = lenrv1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
	rsav[i__] = dvod01_2.rvod1[i__ - 1];
    }
    i__1 = lenrv2;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L12: */
	rsav[lenrv1 + i__] = dvod02_2.rvod2[i__ - 1];
    }
    ioff = lenrv1 + lenrv2;
    i__1 = lrvk1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L14: */
	rsav[ioff + i__] = dvpk01_2.rvpk1[i__ - 1];
    }

    i__1 = leniv1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L20: */
	isav[i__] = dvod01_2.ivod1[i__ - 1];
    }
    i__1 = leniv2;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L22: */
	isav[leniv1 + i__] = dvod02_2.ivod2[i__ - 1];
    }
    ioff = leniv1 + leniv2;
    i__1 = livk1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L24: */
	isav[ioff + i__] = dvpk01_2.ivpk1[i__ - 1];
    }

    return 0;

L100:
    i__1 = lenrv1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L110: */
	dvod01_2.rvod1[i__ - 1] = rsav[i__];
    }
    i__1 = lenrv2;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L112: */
	dvod02_2.rvod2[i__ - 1] = rsav[lenrv1 + i__];
    }
    ioff = lenrv1 + lenrv2;
    i__1 = lrvk1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L114: */
	dvpk01_2.rvpk1[i__ - 1] = rsav[ioff + i__];
    }

    i__1 = leniv1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L120: */
	dvod01_2.ivod1[i__ - 1] = isav[i__];
    }
    i__1 = leniv2;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L122: */
	dvod02_2.ivod2[i__ - 1] = isav[leniv1 + i__];
    }
    ioff = leniv1 + leniv2;
    i__1 = livk1;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L124: */
	dvpk01_2.ivpk1[i__ - 1] = isav[ioff + i__];
    }

    return 0;
/* ----------------------- End of Subroutine DVKSRC ---------------------- */
} /* dvksrc_ */

#ifndef VODE_COMMON_DEFINED
/* The following functions are shared between vode.c, vodpk.c, and zvode.c.
   When VODE_COMMON_DEFINED is set, these are provided by vode_common.c */

/* DECK DVHIN */
/* Subroutine */ int dvhin_(integer *n, doublereal *t0, doublereal *y0, 
	doublereal *ydot, S_fp f, doublereal *rpar, integer *ipar, doublereal 
	*tout, doublereal *uround, doublereal *ewt, integer *itol, doublereal 
	*atol, doublereal *y, doublereal *temp, doublereal *h0, integer *
	niter, integer *ier)
{
    /* Initialized data */

    static doublereal half = .5;
    static doublereal hun = 100.;
    static doublereal pt1 = .1;
    static doublereal two = 2.;

    /* System generated locals */
    integer i__1;
    doublereal d__1, d__2;

    /* Builtin functions */
    double sqrt(doublereal), d_sign(doublereal *, doublereal *);

    /* Local variables */
    doublereal h__;
    integer i__;
    doublereal t1, hg, afi, hlb, hub, hrat, hnew;
    integer iter;
    doublereal delyi, atoli, tdist, yddnrm;
    extern doublereal dvnorm_(integer *, doublereal *, doublereal *);
    doublereal tround;

/* ----------------------------------------------------------------------- */
/* Call sequence input -- N, T0, Y0, YDOT, F, RPAR, IPAR, TOUT, UROUND, */
/*                        EWT, ITOL, ATOL, Y, TEMP */
/* Call sequence output -- H0, NITER, IER */
/* COMMON block variables accessed -- None */

/* Subroutines called by DVHIN:  F */
/* Function routines called by DVHI: DVNORM */
/* ----------------------------------------------------------------------- */
/* This routine computes the step size, H0, to be attempted on the */
/* first step, when the user has not supplied a value for this. */

/* First we check that TOUT - T0 differs significantly from zero.  Then */
/* an iteration is done to approximate the initial second derivative */
/* and this is used to define h from w.r.m.s.norm(h**2 * yddot / 2) = 1. */
/* A bias factor of 1/2 is applied to the resulting h. */
/* The sign of H0 is inferred from the initial values of TOUT and T0. */

/* Communication with DVHIN is done with the following variables: */

/* N      = Size of ODE system, input. */
/* T0     = Initial value of independent variable, input. */
/* Y0     = Vector of initial conditions, input. */
/* YDOT   = Vector of initial first derivatives, input. */
/* F      = Name of subroutine for right-hand side f(t,y), input. */
/* RPAR, IPAR = Dummy names for user's real and integer work arrays. */
/* TOUT   = First output value of independent variable */
/* UROUND = Machine unit roundoff */
/* EWT, ITOL, ATOL = Error weights and tolerance parameters */
/*                   as described in the driver routine, input. */
/* Y, TEMP = Work arrays of length N. */
/* H0     = Step size to be attempted, output. */
/* NITER  = Number of iterations (and of f evaluations) to compute H0, */
/*          output. */
/* IER    = The error flag, returned with the value */
/*          IER = 0  if no trouble occurred, or */
/*          IER = -1 if TOUT and T0 are considered too close to proceed. */
/* ----------------------------------------------------------------------- */

/* Type declarations for local variables -------------------------------- */


/* Type declaration for function subroutines called --------------------- */

/* ----------------------------------------------------------------------- */
/* The following Fortran-77 declaration is to cause the values of the */
/* listed (local) variables to be saved between calls to this integrator. */
/* ----------------------------------------------------------------------- */
    /* Parameter adjustments */
    --temp;
    --y;
    --atol;
    --ewt;
    --ipar;
    --rpar;
    --ydot;
    --y0;

    /* Function Body */

    *niter = 0;
    tdist = (d__1 = *tout - *t0, abs(d__1));
/* Computing MAX */
    d__1 = abs(*t0), d__2 = abs(*tout);
    tround = *uround * max(d__1,d__2);
    if (tdist < two * tround) {
	goto L100;
    }

/* Set a lower bound on h based on the roundoff level in T0 and TOUT. --- */
    hlb = hun * tround;
/* Set an upper bound on h based on TOUT-T0 and the initial Y and YDOT. - */
    hub = pt1 * tdist;
    atoli = atol[1];
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
	if (*itol == 2 || *itol == 4) {
	    atoli = atol[i__];
	}
	delyi = pt1 * (d__1 = y0[i__], abs(d__1)) + atoli;
	afi = (d__1 = ydot[i__], abs(d__1));
	if (afi * hub > delyi) {
	    hub = delyi / afi;
	}
/* L10: */
    }

/* Set initial guess for h as geometric mean of upper and lower bounds. - */
    iter = 0;
    hg = sqrt(hlb * hub);
/* If the bounds have crossed, exit with the mean value. ---------------- */
    if (hub < hlb) {
	*h0 = hg;
	goto L90;
    }

/* Looping point for iteration. ----------------------------------------- */
L50:
/* Estimate the second derivative as a difference quotient in f. -------- */
    d__1 = *tout - *t0;
    h__ = d_sign(&hg, &d__1);
    t1 = *t0 + h__;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L60: */
	y[i__] = y0[i__] + h__ * ydot[i__];
    }
    (*f)(n, &t1, &y[1], &temp[1], &rpar[1], &ipar[1]);
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L70: */
	temp[i__] = (temp[i__] - ydot[i__]) / h__;
    }
    yddnrm = dvnorm_(n, &temp[1], &ewt[1]);
/* Get the corresponding new value of h. -------------------------------- */
    if (yddnrm * hub * hub > two) {
	hnew = sqrt(two / yddnrm);
    } else {
	hnew = sqrt(hg * hub);
    }
    ++iter;
/* ----------------------------------------------------------------------- */
/* Test the stopping conditions. */
/* Stop if the new and previous h values differ by a factor of .lt. 2. */
/* Stop if four iterations have been done.  Also, stop with previous h */
/* if HNEW/HG .gt. 2 after first iteration, as this probably means that */
/* the second derivative value is bad because of cancellation error. */
/* ----------------------------------------------------------------------- */
    if (iter >= 4) {
	goto L80;
    }
    hrat = hnew / hg;
    if (hrat > half && hrat < two) {
	goto L80;
    }
    if (iter >= 2 && hnew > two * hg) {
	hnew = hg;
	goto L80;
    }
    hg = hnew;
    goto L50;

/* Iteration done.  Apply bounds, bias factor, and sign.  Then exit. ---- */
L80:
    *h0 = hnew * half;
    if (*h0 < hlb) {
	*h0 = hlb;
    }
    if (*h0 > hub) {
	*h0 = hub;
    }
L90:
    d__1 = *tout - *t0;
    *h0 = d_sign(h0, &d__1);
    *niter = iter;
    *ier = 0;
    return 0;
/* Error return for TOUT - T0 too small. -------------------------------- */
L100:
    *ier = -1;
    return 0;
/* ----------------------- End of Subroutine DVHIN ----------------------- */
} /* dvhin_ */

/* DECK DVINDY */
/* Subroutine */ int dvindy_(doublereal *t, integer *k, doublereal *yh, 
	integer *ldyh, doublereal *dky, integer *iflag)
{
    /* Initialized data */

    static doublereal hun = 100.;
    static doublereal zero = 0.;

    /* System generated locals */
    integer yh_dim1, yh_offset, i__1, i__2;
    doublereal d__1;

    /* Builtin functions */
    double d_sign(doublereal *, doublereal *), pow_di(doublereal *, integer *)
	    ;
    /* Subroutine */ int s_copy(char *, char *, ftnlen, ftnlen);

    /* Local variables */
    doublereal c__;
    integer i__, j;
    doublereal r__, s;
    integer ic, jb, jj;
    doublereal tp;
    integer jb2, jj1, jp1;
    doublereal tn1;
    char msg[80];
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    doublereal tfuzz;
    extern /* Subroutine */ int xerrwd_(char *, integer *, integer *, integer 
	    *, integer *, integer *, integer *, integer *, doublereal *, 
	    doublereal *, ftnlen);

/* ----------------------------------------------------------------------- */
/* Call sequence input -- T, K, YH, LDYH */
/* Call sequence output -- DKY, IFLAG */
/* COMMON block variables accessed: */
/*     /DVOD01/ --  H, TN, UROUND, L, N, NQ */
/*     /DVOD02/ --  HU */

/* Subroutines called by DVINDY: DSCAL, XERRWD */
/* Function routines called by DVINDY: None */
/* ----------------------------------------------------------------------- */
/* DVINDY computes interpolated values of the K-th derivative of the */
/* dependent variable vector y, and stores it in DKY.  This routine */
/* is called within the package with K = 0 and T = TOUT, but may */
/* also be called by the user for any K up to the current order. */
/* (See detailed instructions in the usage documentation.) */
/* ----------------------------------------------------------------------- */
/* The computed values in DKY are gotten by interpolation using the */
/* Nordsieck history array YH.  This array corresponds uniquely to a */
/* vector-valued polynomial of degree NQCUR or less, and DKY is set */
/* to the K-th derivative of this polynomial at T. */
/* The formula for DKY is: */
/*              q */
/*  DKY(i)  =  sum  c(j,K) * (T - TN)**(j-K) * H**(-j) * YH(i,j+1) */
/*             j=K */
/* where  c(j,K) = j*(j-1)*...*(j-K+1), q = NQCUR, TN = TCUR, H = HCUR. */
/* The quantities  NQ = NQCUR, L = NQ+1, N, TN, and H are */
/* communicated by COMMON.  The above sum is done in reverse order. */
/* IFLAG is returned negative if either K or T is out of bounds. */

/* Discussion above and comments in driver explain all variables. */
/* ----------------------------------------------------------------------- */

/* Type declarations for labeled COMMON block DVOD01 -------------------- */


/* Type declarations for labeled COMMON block DVOD02 -------------------- */


/* Type declarations for local variables -------------------------------- */

/* ----------------------------------------------------------------------- */
/* The following Fortran-77 declaration is to cause the values of the */
/* listed (local) variables to be saved between calls to this integrator. */
/* ----------------------------------------------------------------------- */


    /* Parameter adjustments */
    yh_dim1 = *ldyh;
    yh_offset = 1 + yh_dim1;
    yh -= yh_offset;
    --dky;

    /* Function Body */

    *iflag = 0;
    if (*k < 0 || *k > dvod01_1.nq) {
	goto L80;
    }
    d__1 = abs(dvod01_1.tn) + abs(dvod02_3.hu);
    tfuzz = hun * dvod01_1.uround * d_sign(&d__1, &dvod02_3.hu);
    tp = dvod01_1.tn - dvod02_3.hu - tfuzz;
    tn1 = dvod01_1.tn + tfuzz;
    if ((*t - tp) * (*t - tn1) > zero) {
	goto L90;
    }

    s = (*t - dvod01_1.tn) / dvod01_1.h__;
    ic = 1;
    if (*k == 0) {
	goto L15;
    }
    jj1 = dvod01_1.l - *k;
    i__1 = dvod01_1.nq;
    for (jj = jj1; jj <= i__1; ++jj) {
/* L10: */
	ic *= jj;
    }
L15:
    c__ = (real) ic;
    i__1 = dvod01_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L20: */
	dky[i__] = c__ * yh[i__ + dvod01_1.l * yh_dim1];
    }
    if (*k == dvod01_1.nq) {
	goto L55;
    }
    jb2 = dvod01_1.nq - *k;
    i__1 = jb2;
    for (jb = 1; jb <= i__1; ++jb) {
	j = dvod01_1.nq - jb;
	jp1 = j + 1;
	ic = 1;
	if (*k == 0) {
	    goto L35;
	}
	jj1 = jp1 - *k;
	i__2 = j;
	for (jj = jj1; jj <= i__2; ++jj) {
/* L30: */
	    ic *= jj;
	}
L35:
	c__ = (real) ic;
	i__2 = dvod01_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L40: */
	    dky[i__] = c__ * yh[i__ + jp1 * yh_dim1] + s * dky[i__];
	}
/* L50: */
    }
    if (*k == 0) {
	return 0;
    }
L55:
    i__1 = -(*k);
    r__ = pow_di(&dvod01_1.h__, &i__1);
    dscal_(&dvod01_1.n, &r__, &dky[1], &c__1);
    return 0;

L80:
    s_copy(msg, "DVINDY-- K (=I1) illegal      ", (ftnlen)80, (ftnlen)30);
    xerrwd_(msg, &c__30, &c__51, &c__1, &c__1, k, &c__0, &c__0, &zero, &zero, 
	    (ftnlen)80);
    *iflag = -1;
    return 0;
L90:
    s_copy(msg, "DVINDY-- T (=R1) illegal      ", (ftnlen)80, (ftnlen)30);
    xerrwd_(msg, &c__30, &c__52, &c__1, &c__0, &c__0, &c__0, &c__1, t, &zero, 
	    (ftnlen)80);
    s_copy(msg, "      T not in interval TCUR - HU (= R1) to TCUR (=R2)      "
	    , (ftnlen)80, (ftnlen)60);
    xerrwd_(msg, &c__60, &c__52, &c__1, &c__0, &c__0, &c__0, &c__2, &tp, &
	    dvod01_1.tn, (ftnlen)80);
    *iflag = -2;
    return 0;
/* ----------------------- End of Subroutine DVINDY ---------------------- */
} /* dvindy_ */

#endif /* VODE_COMMON_DEFINED - end guard for dvhin_ and dvindy_ */

#ifndef VODE_COMMON_DEFINED
/* dvstep_, dvset_, and dvjust_ are also shared with vode.c.
   When VODE_COMMON_DEFINED is set, use the versions from vode.c */

/* DECK DVSTEP */
/* Subroutine */ int dvstep_(doublereal *y, doublereal *yh, integer *ldyh, 
	doublereal *yh1, doublereal *ewt, doublereal *savf, doublereal *vsav, 
	doublereal *acor, doublereal *wm, integer *iwm, S_fp f, U_fp jac, 
	U_fp psol, S_fp vnls, doublereal *rpar, integer *ipar)
{
    /* Initialized data */

    static integer kfc = -3;
    static integer kfh = -7;
    static integer mxncf = 10;
    static doublereal addon = 1e-6;
    static doublereal bias1 = 6.;
    static doublereal bias2 = 6.;
    static doublereal bias3 = 10.;
    static doublereal etacf = .25;
    static doublereal etamin = .1;
    static doublereal etamxf = .2;
    static doublereal etamx1 = 1e4;
    static doublereal etamx2 = 10.;
    static doublereal etamx3 = 10.;
    static doublereal onepsm = 1.00001;
    static doublereal thresh = 1.5;
    static doublereal one = 1.;
    static doublereal zero = 0.;

    /* System generated locals */
    integer yh_dim1, yh_offset, i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double pow_dd(doublereal *, doublereal *), pow_di(doublereal *, integer *)
	    ;

    /* Local variables */
    integer i__, j;
    doublereal r__;
    integer i1, i2, jb;
    doublereal ddn;
    integer ncf;
    doublereal dsm, dup;
    static doublereal etaq;
    doublereal told;
    integer iback;
    extern /* Subroutine */ int dscal_(integer *, doublereal *, doublereal *, 
	    integer *);
    integer nflag;
    extern /* Subroutine */ int dcopy_(integer *, doublereal *, integer *, 
	    doublereal *, integer *);
    doublereal flotl;
    extern /* Subroutine */ int dvset_(void), daxpy_(integer *, doublereal *, 
	    doublereal *, integer *, doublereal *, integer *);
    static doublereal etaqm1;
    doublereal etaqp1;
    extern doublereal dvnorm_(integer *, doublereal *, doublereal *);
    doublereal cnquot;
    extern /* Subroutine */ int dvjust_(doublereal *, integer *, integer *);

/* ----------------------------------------------------------------------- */
/* Call sequence input -- Y, YH, LDYH, YH1, EWT, SAVF, VSAV, */
/*                        ACOR, WM, IWM, F, JAC, PSOL, VNLS, RPAR, IPAR */
/* Call sequence output -- YH, ACOR, WM, IWM */
/* COMMON block variables accessed: */
/*     /DVOD01/  ACNRM, EL(13), H, HMIN, HMXI, HNEW, HSCAL, RC, TAU(13), */
/*               TQ(5), TN, JCUR, JSTART, KFLAG, KUTH, */
/*               L, LMAX, MAXORD, N, NEWQ, NQ, NQWAIT */
/*     /DVOD02/  HU, NCFN, NETF, NFE, NQU, NST */

/* Subroutines called by DVSTEP: F, DAXPY, DCOPY, DSCAL, */
/*                               DVJUST, VNLS, DVSET */
/* Function routines called by DVSTEP: DVNORM */
/* ----------------------------------------------------------------------- */
/* DVSTEP performs one step of the integration of an initial value */
/* problem for a system of ordinary differential equations. */
/* DVSTEP calls subroutine VNLS for the solution of the nonlinear system */
/* arising in the time step.  Thus it is independent of the problem */
/* Jacobian structure and the type of nonlinear system solution method. */
/* DVSTEP returns a completion flag KFLAG (in COMMON). */
/* A return with KFLAG = -1 or -2 means either ABS(H) = HMIN or 10 */
/* consecutive failures occurred.  On a return with KFLAG negative, */
/* the values of TN and the YH array are as of the beginning of the last */
/* step, and H is the last step size attempted. */

/* Communication with DVSTEP is done with the following variables: */

/* Y      = An array of length N used for the dependent variable vector. */
/* YH     = An LDYH by LMAX array containing the dependent variables */
/*          and their approximate scaled derivatives, where */
/*          LMAX = MAXORD + 1.  YH(i,j+1) contains the approximate */
/*          j-th derivative of y(i), scaled by H**j/factorial(j) */
/*          (j = 0,1,...,NQ).  On entry for the first step, the first */
/*          two columns of YH must be set from the initial values. */
/* LDYH   = A constant integer .ge. N, the first dimension of YH. */
/*          N is the number of ODEs in the system. */
/* YH1    = A one-dimensional array occupying the same space as YH. */
/* EWT    = An array of length N containing multiplicative weights */
/*          for local error measurements.  Local errors in y(i) are */
/*          compared to 1.0/EWT(i) in various error tests. */
/* SAVF   = An array of working storage, of length N. */
/*          also used for input of YH(*,MAXORD+2) when JSTART = -1 */
/*          and MAXORD .lt. the current order NQ. */
/* VSAV   = A work array of length N passed to subroutine VNLS. */
/* ACOR   = A work array of length N, used for the accumulated */
/*          corrections.  On a successful return, ACOR(i) contains */
/*          the estimated one-step local error in y(i). */
/* WM,IWM = Real and integer work arrays associated with matrix */
/*          operations in VNLS. */
/* F      = Dummy name for the user supplied subroutine for f. */
/* JAC    = Dummy name for the user supplied Jacobian subroutine. */
/* PSOL   = Dummy name for the subroutine passed to VNLS, for */
/*          possible use there. */
/* VNLS   = Dummy name for the nonlinear system solving subroutine, */
/*          whose real name is dependent on the method used. */
/* RPAR, IPAR = Dummy names for user's real and integer work arrays. */
/* ----------------------------------------------------------------------- */

/* Type declarations for labeled COMMON block DVOD01 -------------------- */


/* Type declarations for labeled COMMON block DVOD02 -------------------- */


/* Type declarations for local variables -------------------------------- */


/* Type declaration for function subroutines called --------------------- */

/* ----------------------------------------------------------------------- */
/* The following Fortran-77 declaration is to cause the values of the */
/* listed (local) variables to be saved between calls to this integrator. */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    --y;
    yh_dim1 = *ldyh;
    yh_offset = 1 + yh_dim1;
    yh -= yh_offset;
    --yh1;
    --ewt;
    --savf;
    --vsav;
    --acor;
    --wm;
    --iwm;
    --rpar;
    --ipar;

    /* Function Body */

    dvod01_1.kflag = 0;
    told = dvod01_1.tn;
    ncf = 0;
    dvod01_1.jcur = 0;
    nflag = 0;
    if (dvod01_1.jstart > 0) {
	goto L20;
    }
    if (dvod01_1.jstart == -1) {
	goto L100;
    }
/* ----------------------------------------------------------------------- */
/* On the first call, the order is set to 1, and other variables are */
/* initialized.  ETAMAX is the maximum ratio by which H can be increased */
/* in a single step.  It is normally 10, but is larger during the */
/* first step to compensate for the small initial H.  If a failure */
/* occurs (in corrector convergence or error test), ETAMAX is set to 1 */
/* for the next increase. */
/* ----------------------------------------------------------------------- */
    dvod01_1.lmax = dvod01_1.maxord + 1;
    dvod01_1.nq = 1;
    dvod01_1.l = 2;
    dvod01_1.nqnyh = dvod01_1.nq * *ldyh;
    dvod01_1.tau[0] = dvod01_1.h__;
    dvod01_1.prl1 = one;
    dvod01_1.rc = zero;
    dvod01_1.etamax = etamx1;
    dvod01_1.nqwait = 2;
    dvod01_1.hscal = dvod01_1.h__;
    goto L200;
/* ----------------------------------------------------------------------- */
/* Take preliminary actions on a normal continuation step (JSTART.GT.0). */
/* If the driver changed H, then ETA must be reset and NEWH set to 1. */
/* If a change of order was dictated on the previous step, then */
/* it is done here and appropriate adjustments in the history are made. */
/* On an order decrease, the history array is adjusted by DVJUST. */
/* On an order increase, the history array is augmented by a column. */
/* On a change of step size H, the history array YH is rescaled. */
/* ----------------------------------------------------------------------- */
L20:
    if (dvod01_1.kuth == 1) {
/* Computing MIN */
	d__1 = dvod01_1.eta, d__2 = dvod01_1.h__ / dvod01_1.hscal;
	dvod01_1.eta = min(d__1,d__2);
	dvod01_1.newh = 1;
    }
L50:
    if (dvod01_1.newh == 0) {
	goto L200;
    }
    if (dvod01_1.newq == dvod01_1.nq) {
	goto L150;
    }
    if (dvod01_1.newq < dvod01_1.nq) {
	dvjust_(&yh[yh_offset], ldyh, &c_n1);
	dvod01_1.nq = dvod01_1.newq;
	dvod01_1.l = dvod01_1.nq + 1;
	dvod01_1.nqwait = dvod01_1.l;
	goto L150;
    }
    if (dvod01_1.newq > dvod01_1.nq) {
	dvjust_(&yh[yh_offset], ldyh, &c__1);
	dvod01_1.nq = dvod01_1.newq;
	dvod01_1.l = dvod01_1.nq + 1;
	dvod01_1.nqwait = dvod01_1.l;
	goto L150;
    }
/* ----------------------------------------------------------------------- */
/* The following block handles preliminaries needed when JSTART = -1. */
/* If N was reduced, zero out part of YH to avoid undefined references. */
/* If MAXORD was reduced to a value less than the tentative order NEWQ, */
/* then NQ is set to MAXORD, and a new H ratio ETA is chosen. */
/* Otherwise, we take the same preliminary actions as for JSTART .gt. 0. */
/* In any case, NQWAIT is reset to L = NQ + 1 to prevent further */
/* changes in order for that many steps. */
/* The new H ratio ETA is limited by the input H if KUTH = 1, */
/* by HMIN if KUTH = 0, and by HMXI in any case. */
/* Finally, the history array YH is rescaled. */
/* ----------------------------------------------------------------------- */
L100:
    dvod01_1.lmax = dvod01_1.maxord + 1;
    if (dvod01_1.n == *ldyh) {
	goto L120;
    }
    i1 = (dvod01_1.newq + 1) * *ldyh + 1;
    i2 = (dvod01_1.maxord + 1) * *ldyh;
    if (i1 > i2) {
	goto L120;
    }
    i__1 = i2;
    for (i__ = i1; i__ <= i__1; ++i__) {
/* L110: */
	yh1[i__] = zero;
    }
L120:
    if (dvod01_1.newq <= dvod01_1.maxord) {
	goto L140;
    }
    flotl = (real) dvod01_1.lmax;
    if (dvod01_1.maxord < dvod01_1.nq - 1) {
	ddn = dvnorm_(&dvod01_1.n, &savf[1], &ewt[1]) / dvod01_1.tq[0];
	d__1 = bias1 * ddn;
	d__2 = one / flotl;
	dvod01_1.eta = one / (pow_dd(&d__1, &d__2) + addon);
    }
    if (dvod01_1.maxord == dvod01_1.nq && dvod01_1.newq == dvod01_1.nq + 1) {
	dvod01_1.eta = etaq;
    }
    if (dvod01_1.maxord == dvod01_1.nq - 1 && dvod01_1.newq == dvod01_1.nq + 
	    1) {
	dvod01_1.eta = etaqm1;
	dvjust_(&yh[yh_offset], ldyh, &c_n1);
    }
    if (dvod01_1.maxord == dvod01_1.nq - 1 && dvod01_1.newq == dvod01_1.nq) {
	ddn = dvnorm_(&dvod01_1.n, &savf[1], &ewt[1]) / dvod01_1.tq[0];
	d__1 = bias1 * ddn;
	d__2 = one / flotl;
	dvod01_1.eta = one / (pow_dd(&d__1, &d__2) + addon);
	dvjust_(&yh[yh_offset], ldyh, &c_n1);
    }
    dvod01_1.eta = min(dvod01_1.eta,one);
    dvod01_1.nq = dvod01_1.maxord;
    dvod01_1.l = dvod01_1.lmax;
L140:
    if (dvod01_1.kuth == 1) {
/* Computing MIN */
	d__2 = dvod01_1.eta, d__3 = (d__1 = dvod01_1.h__ / dvod01_1.hscal, 
		abs(d__1));
	dvod01_1.eta = min(d__2,d__3);
    }
    if (dvod01_1.kuth == 0) {
/* Computing MAX */
	d__1 = dvod01_1.eta, d__2 = dvod01_1.hmin / abs(dvod01_1.hscal);
	dvod01_1.eta = max(d__1,d__2);
    }
/* Computing MAX */
    d__1 = one, d__2 = abs(dvod01_1.hscal) * dvod01_1.hmxi * dvod01_1.eta;
    dvod01_1.eta /= max(d__1,d__2);
    dvod01_1.newh = 1;
    dvod01_1.nqwait = dvod01_1.l;
    if (dvod01_1.newq <= dvod01_1.maxord) {
	goto L50;
    }
/* Rescale the history array for a change in H by a factor of ETA. ------ */
L150:
    r__ = one;
    i__1 = dvod01_1.l;
    for (j = 2; j <= i__1; ++j) {
	r__ *= dvod01_1.eta;
	dscal_(&dvod01_1.n, &r__, &yh[j * yh_dim1 + 1], &c__1);
/* L180: */
    }
    dvod01_1.h__ = dvod01_1.hscal * dvod01_1.eta;
    dvod01_1.hscal = dvod01_1.h__;
    dvod01_1.rc *= dvod01_1.eta;
    dvod01_1.nqnyh = dvod01_1.nq * *ldyh;
/* ----------------------------------------------------------------------- */
/* This section computes the predicted values by effectively */
/* multiplying the YH array by the Pascal triangle matrix. */
/* DVSET is called to calculate all integration coefficients. */
/* RC is the ratio of new to old values of the coefficient H/EL(2)=h/l1. */
/* ----------------------------------------------------------------------- */
L200:
    dvod01_1.tn += dvod01_1.h__;
    i1 = dvod01_1.nqnyh + 1;
    i__1 = dvod01_1.nq;
    for (jb = 1; jb <= i__1; ++jb) {
	i1 -= *ldyh;
	i__2 = dvod01_1.nqnyh;
	for (i__ = i1; i__ <= i__2; ++i__) {
/* L210: */
	    yh1[i__] += yh1[i__ + *ldyh];
	}
/* L220: */
    }
    dvset_();
    dvod01_1.rl1 = one / dvod01_1.el[1];
    dvod01_1.rc *= dvod01_1.rl1 / dvod01_1.prl1;
    dvod01_1.prl1 = dvod01_1.rl1;

/* Call the nonlinear system solver. ------------------------------------ */

    (*vnls)(&y[1], &yh[yh_offset], ldyh, &vsav[1], &savf[1], &ewt[1], &acor[1]
	    , &iwm[1], &wm[1], (S_fp)f, (U_fp)jac, (U_fp)psol, &nflag, &rpar[
	    1], &ipar[1]);

    if (nflag == 0) {
	goto L450;
    }
/* ----------------------------------------------------------------------- */
/* The VNLS routine failed to achieve convergence (NFLAG .NE. 0). */
/* The YH array is retracted to its values before prediction. */
/* The step size H is reduced and the step is retried, if possible. */
/* Otherwise, an error exit is taken. */
/* ----------------------------------------------------------------------- */
    ++ncf;
    ++dvod02_3.ncfn;
    dvod01_1.etamax = one;
    dvod01_1.tn = told;
    i1 = dvod01_1.nqnyh + 1;
    i__1 = dvod01_1.nq;
    for (jb = 1; jb <= i__1; ++jb) {
	i1 -= *ldyh;
	i__2 = dvod01_1.nqnyh;
	for (i__ = i1; i__ <= i__2; ++i__) {
/* L420: */
	    yh1[i__] -= yh1[i__ + *ldyh];
	}
/* L430: */
    }
    if (nflag < -1) {
	goto L680;
    }
    if (abs(dvod01_1.h__) <= dvod01_1.hmin * onepsm) {
	goto L670;
    }
    if (ncf == mxncf) {
	goto L670;
    }
    dvod01_1.eta = etacf;
/* Computing MAX */
    d__1 = dvod01_1.eta, d__2 = dvod01_1.hmin / abs(dvod01_1.h__);
    dvod01_1.eta = max(d__1,d__2);
    nflag = -1;
    goto L150;
/* ----------------------------------------------------------------------- */
/* The corrector has converged (NFLAG = 0).  The local error test is */
/* made and control passes to statement 500 if it fails. */
/* ----------------------------------------------------------------------- */
L450:
    dsm = dvod01_1.acnrm / dvod01_1.tq[1];
    if (dsm > one) {
	goto L500;
    }
/* ----------------------------------------------------------------------- */
/* After a successful step, update the YH and TAU arrays and decrement */
/* NQWAIT.  If NQWAIT is then 1 and NQ .lt. MAXORD, then ACOR is saved */
/* for use in a possible order increase on the next step. */
/* If ETAMAX = 1 (a failure occurred this step), keep NQWAIT .ge. 2. */
/* ----------------------------------------------------------------------- */
    dvod01_1.kflag = 0;
    ++dvod02_3.nst;
    dvod02_3.hu = dvod01_1.h__;
    dvod02_3.nqu = dvod01_1.nq;
    i__1 = dvod01_1.nq;
    for (iback = 1; iback <= i__1; ++iback) {
	i__ = dvod01_1.l - iback;
/* L470: */
	dvod01_1.tau[i__] = dvod01_1.tau[i__ - 1];
    }
    dvod01_1.tau[0] = dvod01_1.h__;
    i__1 = dvod01_1.l;
    for (j = 1; j <= i__1; ++j) {
	daxpy_(&dvod01_1.n, &dvod01_1.el[j - 1], &acor[1], &c__1, &yh[j * 
		yh_dim1 + 1], &c__1);
/* L480: */
    }
    --dvod01_1.nqwait;
    if (dvod01_1.l == dvod01_1.lmax || dvod01_1.nqwait != 1) {
	goto L490;
    }
    dcopy_(&dvod01_1.n, &acor[1], &c__1, &yh[dvod01_1.lmax * yh_dim1 + 1], &
	    c__1);
    dvod01_1.conp = dvod01_1.tq[4];
L490:
    if (dvod01_1.etamax != one) {
	goto L560;
    }
    if (dvod01_1.nqwait < 2) {
	dvod01_1.nqwait = 2;
    }
    dvod01_1.newq = dvod01_1.nq;
    dvod01_1.newh = 0;
    dvod01_1.eta = one;
    dvod01_1.hnew = dvod01_1.h__;
    goto L690;
/* ----------------------------------------------------------------------- */
/* The error test failed.  KFLAG keeps track of multiple failures. */
/* Restore TN and the YH array to their previous values, and prepare */
/* to try the step again.  Compute the optimum step size for the */
/* same order.  After repeated failures, H is forced to decrease */
/* more rapidly. */
/* ----------------------------------------------------------------------- */
L500:
    --dvod01_1.kflag;
    ++dvod02_3.netf;
    nflag = -2;
    dvod01_1.tn = told;
    i1 = dvod01_1.nqnyh + 1;
    i__1 = dvod01_1.nq;
    for (jb = 1; jb <= i__1; ++jb) {
	i1 -= *ldyh;
	i__2 = dvod01_1.nqnyh;
	for (i__ = i1; i__ <= i__2; ++i__) {
/* L510: */
	    yh1[i__] -= yh1[i__ + *ldyh];
	}
/* L520: */
    }
    if (abs(dvod01_1.h__) <= dvod01_1.hmin * onepsm) {
	goto L660;
    }
    dvod01_1.etamax = one;
    if (dvod01_1.kflag <= kfc) {
	goto L530;
    }
/* Compute ratio of new H to current H at the current order. ------------ */
    flotl = (real) dvod01_1.l;
    d__1 = bias2 * dsm;
    d__2 = one / flotl;
    dvod01_1.eta = one / (pow_dd(&d__1, &d__2) + addon);
/* Computing MAX */
    d__1 = dvod01_1.eta, d__2 = dvod01_1.hmin / abs(dvod01_1.h__), d__1 = max(
	    d__1,d__2);
    dvod01_1.eta = max(d__1,etamin);
    if (dvod01_1.kflag <= -2 && dvod01_1.eta > etamxf) {
	dvod01_1.eta = etamxf;
    }
    goto L150;
/* ----------------------------------------------------------------------- */
/* Control reaches this section if 3 or more consecutive failures */
/* have occurred.  It is assumed that the elements of the YH array */
/* have accumulated errors of the wrong order.  The order is reduced */
/* by one, if possible.  Then H is reduced by a factor of 0.1 and */
/* the step is retried.  After a total of 7 consecutive failures, */
/* an exit is taken with KFLAG = -1. */
/* ----------------------------------------------------------------------- */
L530:
    if (dvod01_1.kflag == kfh) {
	goto L660;
    }
    if (dvod01_1.nq == 1) {
	goto L540;
    }
/* Computing MAX */
    d__1 = etamin, d__2 = dvod01_1.hmin / abs(dvod01_1.h__);
    dvod01_1.eta = max(d__1,d__2);
    dvjust_(&yh[yh_offset], ldyh, &c_n1);
    dvod01_1.l = dvod01_1.nq;
    --dvod01_1.nq;
    dvod01_1.nqwait = dvod01_1.l;
    goto L150;
L540:
/* Computing MAX */
    d__1 = etamin, d__2 = dvod01_1.hmin / abs(dvod01_1.h__);
    dvod01_1.eta = max(d__1,d__2);
    dvod01_1.h__ *= dvod01_1.eta;
    dvod01_1.hscal = dvod01_1.h__;
    dvod01_1.tau[0] = dvod01_1.h__;
    (*f)(&dvod01_1.n, &dvod01_1.tn, &y[1], &savf[1], &rpar[1], &ipar[1]);
    ++dvod02_3.nfe;
    i__1 = dvod01_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L550: */
	yh[i__ + (yh_dim1 << 1)] = dvod01_1.h__ * savf[i__];
    }
    dvod01_1.nqwait = 10;
    goto L200;
/* ----------------------------------------------------------------------- */
/* If NQWAIT = 0, an increase or decrease in order by one is considered. */
/* Factors ETAQ, ETAQM1, ETAQP1 are computed by which H could */
/* be multiplied at order q, q-1, or q+1, respectively. */
/* The largest of these is determined, and the new order and */
/* step size set accordingly. */
/* A change of H or NQ is made only if H increases by at least a */
/* factor of THRESH.  If an order change is considered and rejected, */
/* then NQWAIT is set to 2 (reconsider it after 2 steps). */
/* ----------------------------------------------------------------------- */
/* Compute ratio of new H to current H at the current order. ------------ */
L560:
    flotl = (real) dvod01_1.l;
    d__1 = bias2 * dsm;
    d__2 = one / flotl;
    etaq = one / (pow_dd(&d__1, &d__2) + addon);
    if (dvod01_1.nqwait != 0) {
	goto L600;
    }
    dvod01_1.nqwait = 2;
    etaqm1 = zero;
    if (dvod01_1.nq == 1) {
	goto L570;
    }
/* Compute ratio of new H to current H at the current order less one. --- */
    ddn = dvnorm_(&dvod01_1.n, &yh[dvod01_1.l * yh_dim1 + 1], &ewt[1]) / 
	    dvod01_1.tq[0];
    d__1 = bias1 * ddn;
    d__2 = one / (flotl - one);
    etaqm1 = one / (pow_dd(&d__1, &d__2) + addon);
L570:
    etaqp1 = zero;
    if (dvod01_1.l == dvod01_1.lmax) {
	goto L580;
    }
/* Compute ratio of new H to current H at current order plus one. ------- */
    d__1 = dvod01_1.h__ / dvod01_1.tau[1];
    cnquot = dvod01_1.tq[4] / dvod01_1.conp * pow_di(&d__1, &dvod01_1.l);
    i__1 = dvod01_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L575: */
	savf[i__] = acor[i__] - cnquot * yh[i__ + dvod01_1.lmax * yh_dim1];
    }
    dup = dvnorm_(&dvod01_1.n, &savf[1], &ewt[1]) / dvod01_1.tq[2];
    d__1 = bias3 * dup;
    d__2 = one / (flotl + one);
    etaqp1 = one / (pow_dd(&d__1, &d__2) + addon);
L580:
    if (etaq >= etaqp1) {
	goto L590;
    }
    if (etaqp1 > etaqm1) {
	goto L620;
    }
    goto L610;
L590:
    if (etaq < etaqm1) {
	goto L610;
    }
L600:
    dvod01_1.eta = etaq;
    dvod01_1.newq = dvod01_1.nq;
    goto L630;
L610:
    dvod01_1.eta = etaqm1;
    dvod01_1.newq = dvod01_1.nq - 1;
    goto L630;
L620:
    dvod01_1.eta = etaqp1;
    dvod01_1.newq = dvod01_1.nq + 1;
    dcopy_(&dvod01_1.n, &acor[1], &c__1, &yh[dvod01_1.lmax * yh_dim1 + 1], &
	    c__1);
/* Test tentative new H against THRESH, ETAMAX, and HMXI, then exit. ---- */
L630:
    if (dvod01_1.eta < thresh || dvod01_1.etamax == one) {
	goto L640;
    }
    dvod01_1.eta = min(dvod01_1.eta,dvod01_1.etamax);
/* Computing MAX */
    d__1 = one, d__2 = abs(dvod01_1.h__) * dvod01_1.hmxi * dvod01_1.eta;
    dvod01_1.eta /= max(d__1,d__2);
    dvod01_1.newh = 1;
    dvod01_1.hnew = dvod01_1.h__ * dvod01_1.eta;
    goto L690;
L640:
    dvod01_1.newq = dvod01_1.nq;
    dvod01_1.newh = 0;
    dvod01_1.eta = one;
    dvod01_1.hnew = dvod01_1.h__;
    goto L690;
/* ----------------------------------------------------------------------- */
/* All returns are made through this section. */
/* On a successful return, ETAMAX is reset and ACOR is scaled. */
/* ----------------------------------------------------------------------- */
L660:
    dvod01_1.kflag = -1;
    goto L720;
L670:
    dvod01_1.kflag = -2;
    goto L720;
L680:
    if (nflag == -2) {
	dvod01_1.kflag = -3;
    }
    if (nflag == -3) {
	dvod01_1.kflag = -4;
    }
    goto L720;
L690:
    dvod01_1.etamax = etamx3;
    if (dvod02_3.nst <= 10) {
	dvod01_1.etamax = etamx2;
    }
/* L700: */
    r__ = one / dvod01_1.tq[1];
    dscal_(&dvod01_1.n, &r__, &acor[1], &c__1);
L720:
    dvod01_1.jstart = 1;
    return 0;
/* ----------------------- End of Subroutine DVSTEP ---------------------- */
} /* dvstep_ */

/* DECK DVSET */
/* Subroutine */ int dvset_(void)
{
    /* Initialized data */

    static doublereal cortes = .1;
    static doublereal one = 1.;
    static doublereal six = 6.;
    static doublereal two = 2.;
    static doublereal zero = 0.;

    /* System generated locals */
    integer i__1, i__2;
    doublereal d__1;

    /* Local variables */
    integer i__, j;
    doublereal s, t1, t2, t3, t4, t5, t6, em[13], xi, em0;
    integer jp1;
    doublereal elp, rxi;
    integer nqm1, nqm2;
    doublereal csum, hsum, rxis, alph0, cnqm1;
    integer iback;
    doublereal floti, flotl, ahatn0, flotnq;

/* ----------------------------------------------------------------------- */
/* Call sequence communication: None */
/* COMMON block variables accessed: */
/*     /DVOD01/ -- EL(13), H, TAU(13), TQ(5), L(= NQ + 1), */
/*                 METH, NQ, NQWAIT */

/* Subroutines called by DVSET: None */
/* Function routines called by DVSET: None */
/* ----------------------------------------------------------------------- */
/* DVSET is called by DVSTEP and sets coefficients for use there. */

/* For each order NQ, the coefficients in EL are calculated by use of */
/*  the generating polynomial lambda(x), with coefficients EL(i). */
/*      lambda(x) = EL(1) + EL(2)*x + ... + EL(NQ+1)*(x**NQ). */
/* For the backward differentiation formulas, */
/*                                     NQ-1 */
/*      lambda(x) = (1 + x/xi*(NQ)) * product (1 + x/xi(i) ) . */
/*                                     i = 1 */
/* For the Adams formulas, */
/*                              NQ-1 */
/*      (d/dx) lambda(x) = c * product (1 + x/xi(i) ) , */
/*                              i = 1 */
/*      lambda(-1) = 0,    lambda(0) = 1, */
/* where c is a normalization constant. */
/* In both cases, xi(i) is defined by */
/*      H*xi(i) = t sub n  -  t sub (n-i) */
/*              = H + TAU(1) + TAU(2) + ... TAU(i-1). */


/* In addition to variables described previously, communication */
/* with DVSET uses the following: */
/*   TAU    = A vector of length 13 containing the past NQ values */
/*            of H. */
/*   EL     = A vector of length 13 in which vset stores the */
/*            coefficients for the corrector formula. */
/*   TQ     = A vector of length 5 in which vset stores constants */
/*            used for the convergence test, the error test, and the */
/*            selection of H at a new order. */
/*   METH   = The basic method indicator. */
/*   NQ     = The current order. */
/*   L      = NQ + 1, the length of the vector stored in EL, and */
/*            the number of columns of the YH array being used. */
/*   NQWAIT = A counter controlling the frequency of order changes. */
/*            An order change is about to be considered if NQWAIT = 1. */
/* ----------------------------------------------------------------------- */

/* Type declarations for labeled COMMON block DVOD01 -------------------- */


/* Type declarations for local variables -------------------------------- */


/* ----------------------------------------------------------------------- */
/* The following Fortran-77 declaration is to cause the values of the */
/* listed (local) variables to be saved between calls to this integrator. */
/* ----------------------------------------------------------------------- */



    flotl = (real) dvod01_1.l;
    nqm1 = dvod01_1.nq - 1;
    nqm2 = dvod01_1.nq - 2;
    switch (dvod01_1.meth) {
	case 1:  goto L100;
	case 2:  goto L200;
    }

/* Set coefficients for Adams methods. ---------------------------------- */
L100:
    if (dvod01_1.nq != 1) {
	goto L110;
    }
    dvod01_1.el[0] = one;
    dvod01_1.el[1] = one;
    dvod01_1.tq[0] = one;
    dvod01_1.tq[1] = two;
    dvod01_1.tq[2] = six * dvod01_1.tq[1];
    dvod01_1.tq[4] = one;
    goto L300;
L110:
    hsum = dvod01_1.h__;
    em[0] = one;
    flotnq = flotl - one;
    i__1 = dvod01_1.l;
    for (i__ = 2; i__ <= i__1; ++i__) {
/* L115: */
	em[i__ - 1] = zero;
    }
    i__1 = nqm1;
    for (j = 1; j <= i__1; ++j) {
	if (j != nqm1 || dvod01_1.nqwait != 1) {
	    goto L130;
	}
	s = one;
	csum = zero;
	i__2 = nqm1;
	for (i__ = 1; i__ <= i__2; ++i__) {
	    csum += s * em[i__ - 1] / (real) (i__ + 1);
/* L120: */
	    s = -s;
	}
	dvod01_1.tq[0] = em[nqm1 - 1] / (flotnq * csum);
L130:
	rxi = dvod01_1.h__ / hsum;
	i__2 = j;
	for (iback = 1; iback <= i__2; ++iback) {
	    i__ = j + 2 - iback;
/* L140: */
	    em[i__ - 1] += em[i__ - 2] * rxi;
	}
	hsum += dvod01_1.tau[j - 1];
/* L150: */
    }
/* Compute integral from -1 to 0 of polynomial and of x times it. ------- */
    s = one;
    em0 = zero;
    csum = zero;
    i__1 = dvod01_1.nq;
    for (i__ = 1; i__ <= i__1; ++i__) {
	floti = (real) i__;
	em0 += s * em[i__ - 1] / floti;
	csum += s * em[i__ - 1] / (floti + one);
/* L160: */
	s = -s;
    }
/* In EL, form coefficients of normalized integrated polynomial. -------- */
    s = one / em0;
    dvod01_1.el[0] = one;
    i__1 = dvod01_1.nq;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L170: */
	dvod01_1.el[i__] = s * em[i__ - 1] / (real) i__;
    }
    xi = hsum / dvod01_1.h__;
    dvod01_1.tq[1] = xi * em0 / csum;
    dvod01_1.tq[4] = xi / dvod01_1.el[dvod01_1.l - 1];
    if (dvod01_1.nqwait != 1) {
	goto L300;
    }
/* For higher order control constant, multiply polynomial by 1+x/xi(q). - */
    rxi = one / xi;
    i__1 = dvod01_1.nq;
    for (iback = 1; iback <= i__1; ++iback) {
	i__ = dvod01_1.l + 1 - iback;
/* L180: */
	em[i__ - 1] += em[i__ - 2] * rxi;
    }
/* Compute integral of polynomial. -------------------------------------- */
    s = one;
    csum = zero;
    i__1 = dvod01_1.l;
    for (i__ = 1; i__ <= i__1; ++i__) {
	csum += s * em[i__ - 1] / (real) (i__ + 1);
/* L190: */
	s = -s;
    }
    dvod01_1.tq[2] = flotl * em0 / csum;
    goto L300;

/* Set coefficients for BDF methods. ------------------------------------ */
L200:
    i__1 = dvod01_1.l;
    for (i__ = 3; i__ <= i__1; ++i__) {
/* L210: */
	dvod01_1.el[i__ - 1] = zero;
    }
    dvod01_1.el[0] = one;
    dvod01_1.el[1] = one;
    alph0 = -one;
    ahatn0 = -one;
    hsum = dvod01_1.h__;
    rxi = one;
    rxis = one;
    if (dvod01_1.nq == 1) {
	goto L240;
    }
    i__1 = nqm2;
    for (j = 1; j <= i__1; ++j) {
/* In EL, construct coefficients of (1+x/xi(1))*...*(1+x/xi(j+1)). ------ */
	hsum += dvod01_1.tau[j - 1];
	rxi = dvod01_1.h__ / hsum;
	jp1 = j + 1;
	alph0 -= one / (real) jp1;
	i__2 = jp1;
	for (iback = 1; iback <= i__2; ++iback) {
	    i__ = j + 3 - iback;
/* L220: */
	    dvod01_1.el[i__ - 1] += dvod01_1.el[i__ - 2] * rxi;
	}
/* L230: */
    }
    alph0 -= one / (real) dvod01_1.nq;
    rxis = -dvod01_1.el[1] - alph0;
    hsum += dvod01_1.tau[nqm1 - 1];
    rxi = dvod01_1.h__ / hsum;
    ahatn0 = -dvod01_1.el[1] - rxi;
    i__1 = dvod01_1.nq;
    for (iback = 1; iback <= i__1; ++iback) {
	i__ = dvod01_1.nq + 2 - iback;
/* L235: */
	dvod01_1.el[i__ - 1] += dvod01_1.el[i__ - 2] * rxis;
    }
L240:
    t1 = one - ahatn0 + alph0;
    t2 = one + (real) dvod01_1.nq * t1;
    dvod01_1.tq[1] = (d__1 = alph0 * t2 / t1, abs(d__1));
    dvod01_1.tq[4] = (d__1 = t2 / (dvod01_1.el[dvod01_1.l - 1] * rxi / rxis), 
	    abs(d__1));
    if (dvod01_1.nqwait != 1) {
	goto L300;
    }
    cnqm1 = rxis / dvod01_1.el[dvod01_1.l - 1];
    t3 = alph0 + one / (real) dvod01_1.nq;
    t4 = ahatn0 + rxi;
    elp = t3 / (one - t4 + t3);
    dvod01_1.tq[0] = (d__1 = elp / cnqm1, abs(d__1));
    hsum += dvod01_1.tau[dvod01_1.nq - 1];
    rxi = dvod01_1.h__ / hsum;
    t5 = alph0 - one / (real) (dvod01_1.nq + 1);
    t6 = ahatn0 - rxi;
    elp = t2 / (one - t6 + t5);
    dvod01_1.tq[2] = (d__1 = elp * rxi * (flotl + one) * t5, abs(d__1));
L300:
    dvod01_1.tq[3] = cortes * dvod01_1.tq[1];
    return 0;
/* ----------------------- End of Subroutine DVSET ----------------------- */
} /* dvset_ */

/* DECK DVJUST */
/* Subroutine */ int dvjust_(doublereal *yh, integer *ldyh, integer *iord)
{
    /* Initialized data */

    static doublereal one = 1.;
    static doublereal zero = 0.;

    /* System generated locals */
    integer yh_dim1, yh_offset, i__1, i__2;

    /* Local variables */
    integer i__, j;
    doublereal t1, xi;
    integer jp1, lp1, nqm1, nqm2, nqp1;
    doublereal prod, hsum, alph0, alph1;
    integer iback;
    doublereal xiold;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);

/* ----------------------------------------------------------------------- */
/* Call sequence input -- YH, LDYH, IORD */
/* Call sequence output -- YH */
/* COMMON block input -- NQ, METH, LMAX, HSCAL, TAU(13), N */
/* COMMON block variables accessed: */
/*     /DVOD01/ -- HSCAL, TAU(13), LMAX, METH, N, NQ, */

/* Subroutines called by DVJUST: DAXPY */
/* Function routines called by DVJUST: None */
/* ----------------------------------------------------------------------- */
/* This subroutine adjusts the YH array on reduction of order, */
/* and also when the order is increased for the stiff option (METH = 2). */
/* Communication with DVJUST uses the following: */
/* IORD  = An integer flag used when METH = 2 to indicate an order */
/*         increase (IORD = +1) or an order decrease (IORD = -1). */
/* HSCAL = Step size H used in scaling of Nordsieck array YH. */
/*         (If IORD = +1, DVJUST assumes that HSCAL = TAU(1).) */
/* See References 1 and 2 for details. */
/* ----------------------------------------------------------------------- */

/* Type declarations for labeled COMMON block DVOD01 -------------------- */


/* Type declarations for local variables -------------------------------- */

/* ----------------------------------------------------------------------- */
/* The following Fortran-77 declaration is to cause the values of the */
/* listed (local) variables to be saved between calls to this integrator. */
/* ----------------------------------------------------------------------- */


    /* Parameter adjustments */
    yh_dim1 = *ldyh;
    yh_offset = 1 + yh_dim1;
    yh -= yh_offset;

    /* Function Body */

    if (dvod01_1.nq == 2 && *iord != 1) {
	return 0;
    }
    nqm1 = dvod01_1.nq - 1;
    nqm2 = dvod01_1.nq - 2;
    switch (dvod01_1.meth) {
	case 1:  goto L100;
	case 2:  goto L200;
    }
/* ----------------------------------------------------------------------- */
/* Nonstiff option... */
/* Check to see if the order is being increased or decreased. */
/* ----------------------------------------------------------------------- */
L100:
    if (*iord == 1) {
	goto L180;
    }
/* Order decrease. ------------------------------------------------------ */
    i__1 = dvod01_1.lmax;
    for (j = 1; j <= i__1; ++j) {
/* L110: */
	dvod01_1.el[j - 1] = zero;
    }
    dvod01_1.el[1] = one;
    hsum = zero;
    i__1 = nqm2;
    for (j = 1; j <= i__1; ++j) {
/* Construct coefficients of x*(x+xi(1))*...*(x+xi(j)). ----------------- */
	hsum += dvod01_1.tau[j - 1];
	xi = hsum / dvod01_1.hscal;
	jp1 = j + 1;
	i__2 = jp1;
	for (iback = 1; iback <= i__2; ++iback) {
	    i__ = j + 3 - iback;
/* L120: */
	    dvod01_1.el[i__ - 1] = dvod01_1.el[i__ - 1] * xi + dvod01_1.el[
		    i__ - 2];
	}
/* L130: */
    }
/* Construct coefficients of integrated polynomial. --------------------- */
    i__1 = nqm1;
    for (j = 2; j <= i__1; ++j) {
/* L140: */
	dvod01_1.el[j] = (real) dvod01_1.nq * dvod01_1.el[j - 1] / (real) j;
    }
/* Subtract correction terms from YH array. ----------------------------- */
    i__1 = dvod01_1.nq;
    for (j = 3; j <= i__1; ++j) {
	i__2 = dvod01_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L160: */
	    yh[i__ + j * yh_dim1] -= yh[i__ + dvod01_1.l * yh_dim1] * 
		    dvod01_1.el[j - 1];
	}
/* L170: */
    }
    return 0;
/* Order increase. ------------------------------------------------------ */
/* Zero out next column in YH array. ------------------------------------ */
L180:
    lp1 = dvod01_1.l + 1;
    i__1 = dvod01_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L190: */
	yh[i__ + lp1 * yh_dim1] = zero;
    }
    return 0;
/* ----------------------------------------------------------------------- */
/* Stiff option... */
/* Check to see if the order is being increased or decreased. */
/* ----------------------------------------------------------------------- */
L200:
    if (*iord == 1) {
	goto L300;
    }
/* Order decrease. ------------------------------------------------------ */
    i__1 = dvod01_1.lmax;
    for (j = 1; j <= i__1; ++j) {
/* L210: */
	dvod01_1.el[j - 1] = zero;
    }
    dvod01_1.el[2] = one;
    hsum = zero;
    i__1 = nqm2;
    for (j = 1; j <= i__1; ++j) {
/* Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). --------------- */
	hsum += dvod01_1.tau[j - 1];
	xi = hsum / dvod01_1.hscal;
	jp1 = j + 1;
	i__2 = jp1;
	for (iback = 1; iback <= i__2; ++iback) {
	    i__ = j + 4 - iback;
/* L220: */
	    dvod01_1.el[i__ - 1] = dvod01_1.el[i__ - 1] * xi + dvod01_1.el[
		    i__ - 2];
	}
/* L230: */
    }
/* Subtract correction terms from YH array. ----------------------------- */
    i__1 = dvod01_1.nq;
    for (j = 3; j <= i__1; ++j) {
	i__2 = dvod01_1.n;
	for (i__ = 1; i__ <= i__2; ++i__) {
/* L240: */
	    yh[i__ + j * yh_dim1] -= yh[i__ + dvod01_1.l * yh_dim1] * 
		    dvod01_1.el[j - 1];
	}
/* L250: */
    }
    return 0;
/* Order increase. ------------------------------------------------------ */
L300:
    i__1 = dvod01_1.lmax;
    for (j = 1; j <= i__1; ++j) {
/* L310: */
	dvod01_1.el[j - 1] = zero;
    }
    dvod01_1.el[2] = one;
    alph0 = -one;
    alph1 = one;
    prod = one;
    xiold = one;
    hsum = dvod01_1.hscal;
    if (dvod01_1.nq == 1) {
	goto L340;
    }
    i__1 = nqm1;
    for (j = 1; j <= i__1; ++j) {
/* Construct coefficients of x*x*(x+xi(1))*...*(x+xi(j)). --------------- */
	jp1 = j + 1;
	hsum += dvod01_1.tau[jp1 - 1];
	xi = hsum / dvod01_1.hscal;
	prod *= xi;
	alph0 -= one / (real) jp1;
	alph1 += one / xi;
	i__2 = jp1;
	for (iback = 1; iback <= i__2; ++iback) {
	    i__ = j + 4 - iback;
/* L320: */
	    dvod01_1.el[i__ - 1] = dvod01_1.el[i__ - 1] * xiold + dvod01_1.el[
		    i__ - 2];
	}
	xiold = xi;
/* L330: */
    }
L340:
    t1 = (-alph0 - alph1) / prod;
/* Load column L + 1 in YH array. --------------------------------------- */
    lp1 = dvod01_1.l + 1;
    i__1 = dvod01_1.n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L350: */
	yh[i__ + lp1 * yh_dim1] = t1 * yh[i__ + dvod01_1.lmax * yh_dim1];
    }
/* Add correction terms to YH array. ------------------------------------ */
    nqp1 = dvod01_1.nq + 1;
    i__1 = nqp1;
    for (j = 3; j <= i__1; ++j) {
	daxpy_(&dvod01_1.n, &dvod01_1.el[j - 1], &yh[lp1 * yh_dim1 + 1], &
		c__1, &yh[j * yh_dim1 + 1], &c__1);
/* L370: */
    }
    return 0;
/* ----------------------- End of Subroutine DVJUST ---------------------- */
} /* dvjust_ */

#endif /* VODE_COMMON_DEFINED - end guard for dvstep_, dvset_, dvjust_ */

#ifndef VODE_COMMON_DEFINED
/* The following utility functions are provided by vode_common.c */

/* DECK DEWSET */
/* Subroutine */ int dewset_(integer *n, integer *itol, doublereal *rtol, 
	doublereal *atol, doublereal *ycur, doublereal *ewt)
{
    /* System generated locals */
    integer i__1;
    doublereal d__1;

    /* Local variables */
    integer i__;

/* ***BEGIN PROLOGUE  DEWSET */
/* ***SUBSIDIARY */
/* ***PURPOSE  Set error weight vector. */
/* ***TYPE      DOUBLE PRECISION (SEWSET-S, DEWSET-D) */
/* ***AUTHOR  Hindmarsh, Alan C., (LLNL) */
/* ***DESCRIPTION */

/*  This subroutine sets the error weight vector EWT according to */
/*      EWT(i) = RTOL(i)*ABS(YCUR(i)) + ATOL(i),  i = 1,...,N, */
/*  with the subscript on RTOL and/or ATOL possibly replaced by 1 above, */
/*  depending on the value of ITOL. */

/* ***SEE ALSO  DLSODE */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791129  DATE WRITTEN */
/*   890501  Modified prologue to SLATEC/LDOC format.  (FNF) */
/*   890503  Minor cosmetic changes.  (FNF) */
/*   930809  Renamed to allow single/double precision versions. (ACH) */
/* ***END PROLOGUE  DEWSET */
/* **End */

/* ***FIRST EXECUTABLE STATEMENT  DEWSET */
    /* Parameter adjustments */
    --ewt;
    --ycur;
    --rtol;
    --atol;

    /* Function Body */
    switch (*itol) {
	case 1:  goto L10;
	case 2:  goto L20;
	case 3:  goto L30;
	case 4:  goto L40;
    }
L10:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L15: */
	ewt[i__] = rtol[1] * (d__1 = ycur[i__], abs(d__1)) + atol[1];
    }
    return 0;
L20:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L25: */
	ewt[i__] = rtol[1] * (d__1 = ycur[i__], abs(d__1)) + atol[i__];
    }
    return 0;
L30:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L35: */
	ewt[i__] = rtol[i__] * (d__1 = ycur[i__], abs(d__1)) + atol[1];
    }
    return 0;
L40:
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L45: */
	ewt[i__] = rtol[i__] * (d__1 = ycur[i__], abs(d__1)) + atol[i__];
    }
    return 0;
/* ----------------------- END OF SUBROUTINE DEWSET ---------------------- */
} /* dewset_ */

/* DECK DVNORM */
doublereal dvnorm_(integer *n, doublereal *v, doublereal *w)
{
    /* System generated locals */
    integer i__1;
    doublereal ret_val, d__1;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer i__;
    doublereal sum;

/* ***BEGIN PROLOGUE  DVNORM */
/* ***SUBSIDIARY */
/* ***PURPOSE  Weighted root-mean-square vector norm. */
/* ***TYPE      DOUBLE PRECISION (SVNORM-S, DVNORM-D) */
/* ***AUTHOR  Hindmarsh, Alan C., (LLNL) */
/* ***DESCRIPTION */

/*  This function routine computes the weighted root-mean-square norm */
/*  of the vector of length N contained in the array V, with weights */
/*  contained in the array W of length N: */
/*    DVNORM = SQRT( (1/N) * SUM( V(i)*W(i) )**2 ) */

/* ***SEE ALSO  DLSODE */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   791129  DATE WRITTEN */
/*   890501  Modified prologue to SLATEC/LDOC format.  (FNF) */
/*   890503  Minor cosmetic changes.  (FNF) */
/*   930809  Renamed to allow single/double precision versions. (ACH) */
/* ***END PROLOGUE  DVNORM */
/* **End */

/* ***FIRST EXECUTABLE STATEMENT  DVNORM */
    /* Parameter adjustments */
    --w;
    --v;

    /* Function Body */
    sum = 0.;
    i__1 = *n;
    for (i__ = 1; i__ <= i__1; ++i__) {
/* L10: */
/* Computing 2nd power */
	d__1 = v[i__] * w[i__];
	sum += d__1 * d__1;
    }
    ret_val = sqrt(sum / *n);
    return ret_val;
/* ----------------------- END OF FUNCTION DVNORM ------------------------ */
} /* dvnorm_ */

#endif /* VODE_COMMON_DEFINED - end guard for dewset_ and dvnorm_ */

/* The following GMRES helper routines are VODPK-specific and not shared */

/* DECK DORTHOG */
/* Subroutine */ int dorthog_(doublereal *vnew, doublereal *v, doublereal *
	hes, integer *n, integer *ll, integer *ldhes, integer *kmp, 
	doublereal *snormw)
{
    /* System generated locals */
    integer v_dim1, v_offset, hes_dim1, hes_offset, i__1, i__2;
    doublereal d__1, d__2, d__3;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    integer i__, i0;
    doublereal arg, tem;
    extern doublereal ddot_(integer *, doublereal *, integer *, doublereal *, 
	    integer *);
    doublereal vnrm;
    extern doublereal dnrm2_(integer *, doublereal *, integer *);
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);
    doublereal sumdsq;

/* ----------------------------------------------------------------------- */
/* This routine orthogonalizes the vector VNEW against the previous */
/* KMP vectors in the V array.  It uses a modified Gram-Schmidt */
/* orthogonalization procedure with conditional reorthogonalization. */
/* This is the version of 28 may 1986. */
/* ----------------------------------------------------------------------- */

/*      On entry */

/*         VNEW = the vector of length N containing a scaled product */
/*                of the Jacobian and the vector V(*,LL). */

/*         V    = the N x l array containing the previous LL */
/*                orthogonal vectors v(*,1) to v(*,LL). */

/*         HES  = an LL x LL upper Hessenberg matrix containing, */
/*                in HES(i,k), k.lt.LL, scaled inner products of */
/*                A*V(*,k) and V(*,i). */

/*        LDHES = the leading dimension of the HES array. */

/*         N    = the order of the matrix A, and the length of VNEW. */

/*         LL   = the current order of the matrix HES. */

/*          KMP = the number of previous vectors the new vector VNEW */
/*                must be made orthogonal to (KMP .le. MAXL). */


/*      On return */

/*         VNEW = the new vector orthogonal to V(*,i0) to V(*,LL), */
/*                where i0 = MAX(1, LL-KMP+1). */

/*         HES  = upper Hessenberg matrix with column LL filled in with */
/*                scaled inner products of A*V(*,LL) and V(*,i). */

/*       SNORMW = L-2 norm of VNEW. */

/* ----------------------------------------------------------------------- */

/* Get norm of unaltered VNEW for later use. ---------------------------- */
    /* Parameter adjustments */
    --vnew;
    v_dim1 = *n;
    v_offset = 1 + v_dim1;
    v -= v_offset;
    hes_dim1 = *ldhes;
    hes_offset = 1 + hes_dim1;
    hes -= hes_offset;

    /* Function Body */
    vnrm = dnrm2_(n, &vnew[1], &c__1);
/* ----------------------------------------------------------------------- */
/* Do modified Gram-Schmidt on VNEW = A*v(LL). */
/* Scaled inner products give new column of HES. */
/* Projections of earlier vectors are subtracted from VNEW. */
/* ----------------------------------------------------------------------- */
/* Computing MAX */
    i__1 = 1, i__2 = *ll - *kmp + 1;
    i0 = max(i__1,i__2);
    i__1 = *ll;
    for (i__ = i0; i__ <= i__1; ++i__) {
	hes[i__ + *ll * hes_dim1] = ddot_(n, &v[i__ * v_dim1 + 1], &c__1, &
		vnew[1], &c__1);
	tem = -hes[i__ + *ll * hes_dim1];
	daxpy_(n, &tem, &v[i__ * v_dim1 + 1], &c__1, &vnew[1], &c__1);
/* L10: */
    }
/* ----------------------------------------------------------------------- */
/* Compute SNORMW = norm of VNEW. */
/* If VNEW is small compared to its input value (in norm), then */
/* reorthogonalize VNEW to V(*,1) through V(*,LL). */
/* Correct if relative correction exceeds 1000*(unit roundoff). */
/* finally, correct SNORMW using the dot products involved. */
/* ----------------------------------------------------------------------- */
    *snormw = dnrm2_(n, &vnew[1], &c__1);
    if (vnrm + *snormw * .001 != vnrm) {
	return 0;
    }
    sumdsq = 0.;
    i__1 = *ll;
    for (i__ = i0; i__ <= i__1; ++i__) {
	tem = -ddot_(n, &v[i__ * v_dim1 + 1], &c__1, &vnew[1], &c__1);
	if (hes[i__ + *ll * hes_dim1] + tem * .001 == hes[i__ + *ll * 
		hes_dim1]) {
	    goto L30;
	}
	hes[i__ + *ll * hes_dim1] -= tem;
	daxpy_(n, &tem, &v[i__ * v_dim1 + 1], &c__1, &vnew[1], &c__1);
/* Computing 2nd power */
	d__1 = tem;
	sumdsq += d__1 * d__1;
L30:
	;
    }
    if (sumdsq == 0.) {
	return 0;
    }
/* Computing MAX */
/* Computing 2nd power */
    d__3 = *snormw;
    d__1 = 0., d__2 = d__3 * d__3 - sumdsq;
    arg = max(d__1,d__2);
    *snormw = sqrt(arg);

    return 0;
/* ----------------------- End of Subroutine DORTHOG --------------------- */
} /* dorthog_ */

/* DECK DHEQR */
/* Subroutine */ int dheqr_(doublereal *a, integer *lda, integer *n, 
	doublereal *q, integer *info, integer *ijob)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Builtin functions */
    double sqrt(doublereal);

    /* Local variables */
    doublereal c__;
    integer i__, j, k;
    doublereal s, t, t1, t2;
    integer iq, km1, kp1, nm1;

/* ----------------------------------------------------------------------- */
/*     This routine performs a QR decomposition of an upper */
/*     Hessenberg matrix A.  There are two options available: */

/*          (1)  performing a fresh decomposition */
/*          (2)  updating the QR factors by adding a row and a */
/*               column to the matrix A. */
/* ----------------------------------------------------------------------- */
/*     DHEQR decomposes an upper Hessenberg matrix by using Givens */
/*     rotations. */

/*     On entry */

/*        A       DOUBLE PRECISION(LDA, N) */
/*                the matrix to be decomposed. */

/*        LDA     INTEGER */
/*                the leading dimension of the array  A . */

/*        N       INTEGER */
/*                A is an (N+1) by N Hessenberg matrix. */

/*        IJOB    INTEGER */
/*                = 1     means that a fresh decomposition of the */
/*                        matrix A is desired. */
/*                .ge. 2  means that the current decomposition of A */
/*                        will be updated by the addition of a row */
/*                        and a column. */
/*     On return */

/*        A       the upper triangular matrix R. */
/*                The factorization can be written Q*A = R, where */
/*                Q is a product of Givens rotations and R is upper */
/*                triangular. */

/*        Q       DOUBLE PRECISION(2*N) */
/*                the factors c and s of each Givens rotation used */
/*                in decomposing A. */

/*        INFO    INTEGER */
/*                = 0  normal value. */
/*                = k  if  A(k,k) .eq. 0.0 .  This is not an error */
/*                     condition for this subroutine, but it does */
/*                     indicate that DHELS will divide by zero */
/*                     if called. */

/*     Modification of LINPACK, by Peter Brown, LLNL. */
/*     Written 1/13/86.  This version dated 6/20/01. */
/* ----------------------------------------------------------------------- */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --q;

    /* Function Body */
    if (*ijob > 1) {
	goto L70;
    }

/* A new facorization is desired. */

/*     QR decomposition without pivoting */

    *info = 0;
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	km1 = k - 1;
	kp1 = k + 1;

/*           Compute kth column of R. */
/*           First, multiply the kth column of A by the previous */
/*           k-1 Givens rotations. */

	if (km1 < 1) {
	    goto L20;
	}
	i__2 = km1;
	for (j = 1; j <= i__2; ++j) {
	    i__ = (j - 1 << 1) + 1;
	    t1 = a[j + k * a_dim1];
	    t2 = a[j + 1 + k * a_dim1];
	    c__ = q[i__];
	    s = q[i__ + 1];
	    a[j + k * a_dim1] = c__ * t1 - s * t2;
	    a[j + 1 + k * a_dim1] = s * t1 + c__ * t2;
/* L10: */
	}

/*           Compute Givens components c and s */

L20:
	iq = (km1 << 1) + 1;
	t1 = a[k + k * a_dim1];
	t2 = a[kp1 + k * a_dim1];
	if (t2 != 0.) {
	    goto L30;
	}
	c__ = 1.;
	s = 0.;
	goto L50;
L30:
	if (abs(t2) < abs(t1)) {
	    goto L40;
	}
	t = t1 / t2;
	s = -1. / sqrt(t * t + 1.);
	c__ = -s * t;
	goto L50;
L40:
	t = t2 / t1;
	c__ = 1. / sqrt(t * t + 1.);
	s = -c__ * t;
L50:
	q[iq] = c__;
	q[iq + 1] = s;
	a[k + k * a_dim1] = c__ * t1 - s * t2;
	if (a[k + k * a_dim1] == 0.) {
	    *info = k;
	}
/* L60: */
    }
    return 0;

/* The old factorization of A will be updated.  A row and a column */
/* has been added to the matrix A. */
/* N by N-1 is now the old size of the matrix. */

L70:
    nm1 = *n - 1;

/* Multiply the new column by the N previous Givens rotations. */

    i__1 = nm1;
    for (k = 1; k <= i__1; ++k) {
	i__ = (k - 1 << 1) + 1;
	t1 = a[k + *n * a_dim1];
	t2 = a[k + 1 + *n * a_dim1];
	c__ = q[i__];
	s = q[i__ + 1];
	a[k + *n * a_dim1] = c__ * t1 - s * t2;
	a[k + 1 + *n * a_dim1] = s * t1 + c__ * t2;
/* L100: */
    }

/* Complete update of decomposition by forming last Givens rotation, */
/* and multiplying it times the column vector (A(N,N), A(N+1,N)). */

    *info = 0;
    t1 = a[*n + *n * a_dim1];
    t2 = a[*n + 1 + *n * a_dim1];
    if (t2 != 0.) {
	goto L110;
    }
    c__ = 1.;
    s = 0.;
    goto L130;
L110:
    if (abs(t2) < abs(t1)) {
	goto L120;
    }
    t = t1 / t2;
    s = -1. / sqrt(t * t + 1.);
    c__ = -s * t;
    goto L130;
L120:
    t = t2 / t1;
    c__ = 1. / sqrt(t * t + 1.);
    s = -c__ * t;
L130:
    iq = (*n << 1) - 1;
    q[iq] = c__;
    q[iq + 1] = s;
    a[*n + *n * a_dim1] = c__ * t1 - s * t2;
    if (a[*n + *n * a_dim1] == 0.) {
	*info = *n;
    }
    return 0;
/* ----------------------- End of Subroutine DHEQR ----------------------- */
} /* dheqr_ */

/* DECK DHELS */
/* Subroutine */ int dhels_(doublereal *a, integer *lda, integer *n, 
	doublereal *q, doublereal *b)
{
    /* System generated locals */
    integer a_dim1, a_offset, i__1, i__2;

    /* Local variables */
    doublereal c__;
    integer k;
    doublereal s, t, t1, t2;
    integer kb, iq, kp1;
    extern /* Subroutine */ int daxpy_(integer *, doublereal *, doublereal *, 
	    integer *, doublereal *, integer *);

/* ----------------------------------------------------------------------- */
/* This is part of the LINPACK routine DGESL with changes */
/* due to the fact that A is an upper Hessenberg matrix. */
/* ----------------------------------------------------------------------- */
/*     DHELS solves the least squares problem */

/*           min (b-A*x, b-A*x) */

/*     using the factors computed by DHEQR. */

/*     On entry */

/*        A       DOUBLE PRECISION(LDA, N) */
/*                the output from DHEQR which contains the upper */
/*                triangular factor R in the QR decomposition of A. */

/*        LDA     INTEGER */
/*                the leading dimension of the array  A . */

/*        N       INTEGER */
/*                A is originally an (N+1) by N matrix. */

/*        Q       DOUBLE PRECISION(2*N) */
/*                The coefficients of the N givens rotations */
/*                used in the QR factorization of A. */

/*        B       DOUBLE PRECISION(N+1) */
/*                the right hand side vector. */

/*     On return */

/*        B       the solution vector  x . */

/*     Modification of LINPACK, by Peter Brown, LLNL. */
/*     Written 1/13/86.  This version dated 6/20/01. */

/*     BLAS called: DAXPY */
/* ----------------------------------------------------------------------- */

/*        Minimize (b-A*x, b-A*x) */
/*        First form Q*b. */

    /* Parameter adjustments */
    a_dim1 = *lda;
    a_offset = 1 + a_dim1;
    a -= a_offset;
    --q;
    --b;

    /* Function Body */
    i__1 = *n;
    for (k = 1; k <= i__1; ++k) {
	kp1 = k + 1;
	iq = (k - 1 << 1) + 1;
	c__ = q[iq];
	s = q[iq + 1];
	t1 = b[k];
	t2 = b[kp1];
	b[k] = c__ * t1 - s * t2;
	b[kp1] = s * t1 + c__ * t2;
/* L20: */
    }

/*        Now solve  R*x = Q*b. */

    i__1 = *n;
    for (kb = 1; kb <= i__1; ++kb) {
	k = *n + 1 - kb;
	b[k] /= a[k + k * a_dim1];
	t = -b[k];
	i__2 = k - 1;
	daxpy_(&i__2, &t, &a[k * a_dim1 + 1], &c__1, &b[1], &c__1);
/* L40: */
    }
    return 0;
/* ----------------------- End of Subroutine DHELS ----------------------- */
} /* dhels_ */

#ifndef VODE_COMMON_DEFINED
/* The following utility functions are provided by vode_common.c */

/* DECK XERRWD */
/* Subroutine */ int xerrwd_(char *msg, integer *nmes, integer *nerr, integer 
	*level, integer *ni, integer *i1, integer *i2, integer *nr, 
	doublereal *r1, doublereal *r2, ftnlen msg_len)
{
    /* Format strings */
    static char fmt_10[] = "(1x,a)";
    static char fmt_20[] = "(6x,\002In above message,  I1 =\002,i10)";
    static char fmt_30[] = "(6x,\002In above message,  I1 =\002,i10,3x,\002I"
	    "2 =\002,i10)";
    static char fmt_40[] = "(6x,\002In above message,  R1 =\002,d21.13)";
    static char fmt_50[] = "(6x,\002In above,  R1 =\002,d21.13,3x,\002R2 "
	    "=\002,d21.13)";

    /* Builtin functions */
    integer s_wsfe(cilist *), do_fio(integer *, char *, ftnlen), e_wsfe(void);
    /* Subroutine */ int s_stop(char *, ftnlen);

    /* Local variables */
    extern integer ixsav_(integer *, integer *, logical *);
    integer lunit, mesflg;

    /* Fortran I/O blocks */
    static cilist io___270 = { 0, 0, 0, fmt_10, 0 };
    static cilist io___271 = { 0, 0, 0, fmt_20, 0 };
    static cilist io___272 = { 0, 0, 0, fmt_30, 0 };
    static cilist io___273 = { 0, 0, 0, fmt_40, 0 };
    static cilist io___274 = { 0, 0, 0, fmt_50, 0 };


/* ***BEGIN PROLOGUE  XERRWD */
/* ***SUBSIDIARY */
/* ***PURPOSE  Write error message with values. */
/* ***CATEGORY  R3C */
/* ***TYPE      DOUBLE PRECISION (XERRWV-S, XERRWD-D) */
/* ***AUTHOR  Hindmarsh, Alan C., (LLNL) */
/* ***DESCRIPTION */

/*  Subroutines XERRWD, XSETF, XSETUN, and the function routine IXSAV, */
/*  as given here, constitute a simplified version of the SLATEC error */
/*  handling package. */

/*  All arguments are input arguments. */

/*  MSG    = The message (character array). */
/*  NMES   = The length of MSG (number of characters). */
/*  NERR   = The error number (not used). */
/*  LEVEL  = The error level.. */
/*           0 or 1 means recoverable (control returns to caller). */
/*           2 means fatal (run is aborted--see note below). */
/*  NI     = Number of integers (0, 1, or 2) to be printed with message. */
/*  I1,I2  = Integers to be printed, depending on NI. */
/*  NR     = Number of reals (0, 1, or 2) to be printed with message. */
/*  R1,R2  = Reals to be printed, depending on NR. */

/*  Note..  this routine is machine-dependent and specialized for use */
/*  in limited context, in the following ways.. */
/*  1. The argument MSG is assumed to be of type CHARACTER, and */
/*     the message is printed with a format of (1X,A). */
/*  2. The message is assumed to take only one line. */
/*     Multi-line messages are generated by repeated calls. */
/*  3. If LEVEL = 2, control passes to the statement   STOP */
/*     to abort the run.  This statement may be machine-dependent. */
/*  4. R1 and R2 are assumed to be in double precision and are printed */
/*     in D21.13 format. */

/* ***ROUTINES CALLED  IXSAV */
/* ***REVISION HISTORY  (YYMMDD) */
/*   920831  DATE WRITTEN */
/*   921118  Replaced MFLGSV/LUNSAV by IXSAV. (ACH) */
/*   930329  Modified prologue to SLATEC format. (FNF) */
/*   930407  Changed MSG from CHARACTER*1 array to variable. (FNF) */
/*   930922  Minor cosmetic change. (FNF) */
/* ***END PROLOGUE  XERRWD */

/* *Internal Notes: */

/* For a different default logical unit number, IXSAV (or a subsidiary */
/* routine that it calls) will need to be modified. */
/* For a different run-abort command, change the statement following */
/* statement 100 at the end. */
/* ----------------------------------------------------------------------- */
/* Subroutines called by XERRWD.. None */
/* Function routine called by XERRWD.. IXSAV */
/* ----------------------------------------------------------------------- */
/* **End */

/*  Declare arguments. */


/*  Declare local variables. */


/*  Get logical unit number and message print flag. */

/* ***FIRST EXECUTABLE STATEMENT  XERRWD */
    lunit = ixsav_(&c__1, &c__0, &c_false);
    mesflg = ixsav_(&c__2, &c__0, &c_false);
    if (mesflg == 0) {
	goto L100;
    }

/*  Write the message. */

    io___270.ciunit = lunit;
    s_wsfe(&io___270);
    do_fio(&c__1, msg, msg_len);
    e_wsfe();
    if (*ni == 1) {
	io___271.ciunit = lunit;
	s_wsfe(&io___271);
	do_fio(&c__1, (char *)&(*i1), (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (*ni == 2) {
	io___272.ciunit = lunit;
	s_wsfe(&io___272);
	do_fio(&c__1, (char *)&(*i1), (ftnlen)sizeof(integer));
	do_fio(&c__1, (char *)&(*i2), (ftnlen)sizeof(integer));
	e_wsfe();
    }
    if (*nr == 1) {
	io___273.ciunit = lunit;
	s_wsfe(&io___273);
	do_fio(&c__1, (char *)&(*r1), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }
    if (*nr == 2) {
	io___274.ciunit = lunit;
	s_wsfe(&io___274);
	do_fio(&c__1, (char *)&(*r1), (ftnlen)sizeof(doublereal));
	do_fio(&c__1, (char *)&(*r2), (ftnlen)sizeof(doublereal));
	e_wsfe();
    }

/*  Abort the run if LEVEL = 2. */

L100:
    if (*level != 2) {
	return 0;
    }
    s_stop("", (ftnlen)0);
/* ----------------------- End of Subroutine XERRWD ---------------------- */
    return 0;
} /* xerrwd_ */

/* DECK XSETF */
/* Subroutine */ int xsetf_(integer *mflag)
{
    integer junk;
    extern integer ixsav_(integer *, integer *, logical *);

/* ***BEGIN PROLOGUE  XSETF */
/* ***PURPOSE  Reset the error print control flag. */
/* ***CATEGORY  R3A */
/* ***TYPE      ALL (XSETF-A) */
/* ***KEYWORDS  ERROR CONTROL */
/* ***AUTHOR  Hindmarsh, Alan C., (LLNL) */
/* ***DESCRIPTION */

/*   XSETF sets the error print control flag to MFLAG: */
/*      MFLAG=1 means print all messages (the default). */
/*      MFLAG=0 means no printing. */

/* ***SEE ALSO  XERRWD, XERRWV */
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  IXSAV */
/* ***REVISION HISTORY  (YYMMDD) */
/*   921118  DATE WRITTEN */
/*   930329  Added SLATEC format prologue. (FNF) */
/*   930407  Corrected SEE ALSO section. (FNF) */
/*   930922  Made user-callable, and other cosmetic changes. (FNF) */
/* ***END PROLOGUE  XSETF */

/* Subroutines called by XSETF.. None */
/* Function routine called by XSETF.. IXSAV */
/* ----------------------------------------------------------------------- */
/* **End */

/* ***FIRST EXECUTABLE STATEMENT  XSETF */
    if (*mflag == 0 || *mflag == 1) {
	junk = ixsav_(&c__2, mflag, &c_true);
    }
    return 0;
/* ----------------------- End of Subroutine XSETF ----------------------- */
} /* xsetf_ */

/* DECK XSETUN */
/* Subroutine */ int xsetun_(integer *lun)
{
    integer junk;
    extern integer ixsav_(integer *, integer *, logical *);

/* ***BEGIN PROLOGUE  XSETUN */
/* ***PURPOSE  Reset the logical unit number for error messages. */
/* ***CATEGORY  R3B */
/* ***TYPE      ALL (XSETUN-A) */
/* ***KEYWORDS  ERROR CONTROL */
/* ***DESCRIPTION */

/*   XSETUN sets the logical unit number for error messages to LUN. */

/* ***AUTHOR  Hindmarsh, Alan C., (LLNL) */
/* ***SEE ALSO  XERRWD, XERRWV */
/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  IXSAV */
/* ***REVISION HISTORY  (YYMMDD) */
/*   921118  DATE WRITTEN */
/*   930329  Added SLATEC format prologue. (FNF) */
/*   930407  Corrected SEE ALSO section. (FNF) */
/*   930922  Made user-callable, and other cosmetic changes. (FNF) */
/* ***END PROLOGUE  XSETUN */

/* Subroutines called by XSETUN.. None */
/* Function routine called by XSETUN.. IXSAV */
/* ----------------------------------------------------------------------- */
/* **End */

/* ***FIRST EXECUTABLE STATEMENT  XSETUN */
    if (*lun > 0) {
	junk = ixsav_(&c__1, lun, &c_true);
    }
    return 0;
/* ----------------------- End of Subroutine XSETUN ---------------------- */
} /* xsetun_ */

/* DECK IXSAV */
integer ixsav_(integer *ipar, integer *ivalue, logical *iset)
{
    /* Initialized data */

    static integer lunit = -1;
    static integer mesflg = 1;

    /* System generated locals */
    integer ret_val;

    /* Local variables */
    extern integer iumach_(void);

/* ***BEGIN PROLOGUE  IXSAV */
/* ***SUBSIDIARY */
/* ***PURPOSE  Save and recall error message control parameters. */
/* ***CATEGORY  R3C */
/* ***TYPE      ALL (IXSAV-A) */
/* ***AUTHOR  Hindmarsh, Alan C., (LLNL) */
/* ***DESCRIPTION */

/*  IXSAV saves and recalls one of two error message parameters: */
/*    LUNIT, the logical unit number to which messages are printed, and */
/*    MESFLG, the message print flag. */
/*  This is a modification of the SLATEC library routine J4SAVE. */

/*  Saved local variables.. */
/*   LUNIT  = Logical unit number for messages.  The default is obtained */
/*            by a call to IUMACH (may be machine-dependent). */
/*   MESFLG = Print control flag.. */
/*            1 means print all messages (the default). */
/*            0 means no printing. */

/*  On input.. */
/*    IPAR   = Parameter indicator (1 for LUNIT, 2 for MESFLG). */
/*    IVALUE = The value to be set for the parameter, if ISET = .TRUE. */
/*    ISET   = Logical flag to indicate whether to read or write. */
/*             If ISET = .TRUE., the parameter will be given */
/*             the value IVALUE.  If ISET = .FALSE., the parameter */
/*             will be unchanged, and IVALUE is a dummy argument. */

/*  On return.. */
/*    IXSAV = The (old) value of the parameter. */

/* ***SEE ALSO  XERRWD, XERRWV */
/* ***ROUTINES CALLED  IUMACH */
/* ***REVISION HISTORY  (YYMMDD) */
/*   921118  DATE WRITTEN */
/*   930329  Modified prologue to SLATEC format. (FNF) */
/*   930915  Added IUMACH call to get default output unit.  (ACH) */
/*   930922  Minor cosmetic changes. (FNF) */
/*   010425  Type declaration for IUMACH added. (ACH) */
/* ***END PROLOGUE  IXSAV */

/* Subroutines called by IXSAV.. None */
/* Function routine called by IXSAV.. IUMACH */
/* ----------------------------------------------------------------------- */
/* **End */
/* ----------------------------------------------------------------------- */
/* ----------------------------------------------------------------------- */
/* The following Fortran-77 declaration is to cause the values of the */
/* listed (local) variables to be saved between calls to this routine. */
/* ----------------------------------------------------------------------- */

/* ***FIRST EXECUTABLE STATEMENT  IXSAV */
    if (*ipar == 1) {
	if (lunit == -1) {
	    lunit = iumach_();
	}
	ret_val = lunit;
	if (*iset) {
	    lunit = *ivalue;
	}
    }

    if (*ipar == 2) {
	ret_val = mesflg;
	if (*iset) {
	    mesflg = *ivalue;
	}
    }

    return ret_val;
/* ----------------------- End of Function IXSAV ------------------------- */
} /* ixsav_ */

/* DECK IUMACH */
integer iumach_(void)
{
    /* System generated locals */
    integer ret_val;

/* ***BEGIN PROLOGUE  IUMACH */
/* ***PURPOSE  Provide standard output unit number. */
/* ***CATEGORY  R1 */
/* ***TYPE      INTEGER (IUMACH-I) */
/* ***KEYWORDS  MACHINE CONSTANTS */
/* ***AUTHOR  Hindmarsh, Alan C., (LLNL) */
/* ***DESCRIPTION */
/* *Usage: */
/*        INTEGER  LOUT, IUMACH */
/*        LOUT = IUMACH() */

/* *Function Return Values: */
/*     LOUT : the standard logical unit for Fortran output. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  (NONE) */
/* ***REVISION HISTORY  (YYMMDD) */
/*   930915  DATE WRITTEN */
/*   930922  Made user-callable, and other cosmetic changes. (FNF) */
/* ***END PROLOGUE  IUMACH */

/* *Internal Notes: */
/*  The built-in value of 6 is standard on a wide range of Fortran */
/*  systems.  This may be machine-dependent. */
/* **End */
/* ***FIRST EXECUTABLE STATEMENT  IUMACH */
    ret_val = 6;

    return ret_val;
/* ----------------------- End of Function IUMACH ------------------------ */
} /* iumach_ */

/* DECK DUMACH */
doublereal dumach_(void)
{
    /* System generated locals */
    doublereal ret_val;

    /* Local variables */
    doublereal u, comp;
    extern /* Subroutine */ int dumsum_(doublereal *, doublereal *, 
	    doublereal *);

/* ***BEGIN PROLOGUE  DUMACH */
/* ***PURPOSE  Compute the unit roundoff of the machine. */
/* ***CATEGORY  R1 */
/* ***TYPE      DOUBLE PRECISION (RUMACH-S, DUMACH-D) */
/* ***KEYWORDS  MACHINE CONSTANTS */
/* ***AUTHOR  Hindmarsh, Alan C., (LLNL) */
/* ***DESCRIPTION */
/* *Usage: */
/*        DOUBLE PRECISION  A, DUMACH */
/*        A = DUMACH() */

/* *Function Return Values: */
/*     A : the unit roundoff of the machine. */

/* *Description: */
/*     The unit roundoff is defined as the smallest positive machine */
/*     number u such that  1.0 + u .ne. 1.0.  This is computed by DUMACH */
/*     in a machine-independent manner. */

/* ***REFERENCES  (NONE) */
/* ***ROUTINES CALLED  DUMSUM */
/* ***REVISION HISTORY  (YYYYMMDD) */
/*   19930216  DATE WRITTEN */
/*   19930818  Added SLATEC-format prologue.  (FNF) */
/*   20030707  Added DUMSUM to force normal storage of COMP.  (ACH) */
/* ***END PROLOGUE  DUMACH */

/* ***FIRST EXECUTABLE STATEMENT  DUMACH */
    u = 1.;
L10:
    u *= .5;
    dumsum_(&c_b811, &u, &comp);
    if (comp != 1.) {
	goto L10;
    }
    ret_val = u * 2.;
    return ret_val;
/* ----------------------- End of Function DUMACH ------------------------ */
} /* dumach_ */

/* Subroutine */ int dumsum_(doublereal *a, doublereal *b, doublereal *c__)
{
/*     Routine to force normal storing of A + B, for DUMACH. */
    *c__ = *a + *b;
    return 0;
} /* dumsum_ */

#endif /* VODE_COMMON_DEFINED */

